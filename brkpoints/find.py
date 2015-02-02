"""
Find locations of translocation breakpoints within the sequence of
translocated genes by aligning junction sequences curated in the
TICdb database (http://www.unav.es/genetica/TICdb).

Inputs
------
1) translocations_file
Tab delimited file of translocation data from TICdb format described
below

2) sequence_file
FASTA format file of exons and introns derived from ensembl modelled
transcripts for translocated genes of interest

3) details_file
YAML format file holding transcript structure and coordinated
details cross-referenced to the sequence_file

4) workdir
Work area for blast

5) blast_bin_dir
Path to the directory containing blast+ executables in particular
makeblastdb and blastn required by this script

Output
------
1) breakpoints_file
YAML dump of the translocations data, with alignments and breakpoints
included
"""

import sys
import os
import shutil
import yaml
import subprocess
import optparse


def usage():
    return (
        "Usage: {prog} " +
        "<translocations_file> " +
        "<sequence_file> " +
        "<details_file> " +
        "<workdir> " +
        "<blast_bin_dir> " +
        "<breakpoints_file>"
        )


def parse_translocations_file(translocations_file):
    """
    Parse translocations file
    -------------------------
    The file from TICdb is in tab separated values format. Each line
    represents a translocation, with four items:
    1. Identifier for the 5' partner gene (HGNC symbol)
    2. Identifier for the 3' partner gene (HGNC symbol)
    3. A reference
    4. A junction sequence

    The file is parsed into a list of translocation objects, which is
    returned to the calling routine.

    Each translocation object is a dict with keys

        'partners' : []            array length 2 of gene objects
        'junction_sequence' : str  the junction sequence
        'xref' : str               the reference

    Each gene object is a dict with keys

        'hgnc_symbol' : str        identifier
        'disp' : str               '5prime' or '3prime'
        'aln' : {}                 alignment + breakpoint, to be populated
    """

    input_fh = open(translocations_file)
    lines = input_fh.readlines()
    input_fh.close()

    translocations = []

    for line in lines:
        if line.startswith('#'):
            continue
        stripped = line.rstrip()
        if not stripped:
            continue
        parts = stripped.split("\t")
        if len(parts) != 4:
            sys.stderr.write("Error - unexpected format on line " + line)
            sys.exit(1)
        translocation = {
            'partners': [{'hgnc_symbol': parts[0], 'disp': '5prime'},
                         {'hgnc_symbol': parts[1], 'disp': '3prime'}],
            'junction_sequence': parts[3],
            'xref': parts[2]
            }
        translocations.append(translocation)

    print "Found %d translocations in %s" % \
        (len(translocations), translocations_file)
    return translocations


def locate_breakpoints(translocation, transcripts, workdir,
                       database_file, blastn_exec):
    """
    Locate translocation breakpoints
    ----------------------------------
    Aligns the junction sequence for a given translocation to the
    genomic sequence in the locations of the translocated genes.

    translocation
    -------------
    Object containing the translocation data

    transcripts
    -----------
    A dict of keyed by HGNC gene symbol; trascripts['XX'] is a list
    of transcript objects, corresponding to the gene 'XX'.

    Each transcript object is a dict with keys

        'name' : str               the ENSEMBL transcript id
        'components' : []          list of exon/intron features
        'chromosome_name': str     the ENSEMBL chromosome_name
        'start_position' : int     chromosome start position
        'end_position' : int       chromosome end position
        'strand' : int             direction of transcript

    Each component (exon/intron feature) is a dict with keys

        'uid' : str                unique sequence id
        'name' : str               exon id (ENSEMBL) or constructed
        'start_position' : int     chromosome start position
        'end_position' : int       chromosome end position
        'type' : str               'exon' or 'intron'

    workdir
    -------
    Work area for blasting

    database_file
    -------------
    BLAST database file, assumed to have been formatted by
    makeblastdb for searches with blastn.

    blastn_exec
    -----------
    BLAST executable
    """

    # Write the junction sequence to a query file
    query_file = workdir + "/%s.fa" % (translocation['xref'])
    output_fh = open(query_file, 'w')
    output_fh.write(">query\n")
    for pos in range(0, len(translocation['junction_sequence']), 60):
        output_fh.write(translocation['junction_sequence'][pos:pos+60]
                        + "\n")
    output_fh.close()

    # Run blastn
    results_file = workdir + "/%s.alns" % (translocation['xref'])
    subprocess.check_call([
        blastn_exec,
        '-query', query_file,
        '-db', database_file,
        '-dust', 'no',                 # No filtering
        '-word_size', '7',             # Capture short matches
        '-ungapped',                   # Ungapped alignment only
        '-perc_identity', '100',       # Identical within aligned region
        '-evalue', '0.001',
        '-outfmt', '7 std sstrand',    # Simple output format
        '-out', results_file           # Output file
        ])

    # Load the results
    input_fh = open(results_file)
    lines = input_fh.readlines()
    input_fh.close()

    # Cleanup
    os.unlink(query_file)
    os.unlink(results_file)

    # Loop through results twice to pull out best (first) alignment
    # of sequences relevant to each partner; note that the database
    # file contains exon/intron sequences for other genes and these
    # may appear in the blast results; relevant sequences are
    # identified via the transcripts data
    #
    for partner in translocation['partners']:

        # Locate relevant sequences
        uids = {}
        aln = None
        transcript_lookup = {}
        if partner['hgnc_symbol'] in transcripts:
            # Transcripts exist for this gene
            for transcript in transcripts[partner['hgnc_symbol']]:
                for component in transcript['components']:
                    transcript_lookup[component['uid']] = transcript
                    uids[component['uid']] = component

        # Examine the BLAST output
        # Break at the first alignment referring to a relevant sequence
        for line in lines:
            if line.startswith('#'):
                continue
            stripped = line.rstrip()
            if not stripped:
                continue
            parts = stripped.split("\t")
            if len(parts) != 13:
                sys.stderr.write("Error - unexpected number of " +
                                 "alignment data\n")
                sys.exit(1)
            if parts[1] in uids:
                aln = {
                    'subject_id': parts[1],
                    'query_start': int(parts[6]),
                    'query_end': int(parts[7]),
                    'subject_start': int(parts[8]),
                    'subject_end': int(parts[9]),
                    'sstrand': parts[12]
                    }
                break
        if not aln:
            continue

        # Establish genomic coordinates of the alignment
        # Start with the aligned coordinates of the exone or intron
        start = aln['subject_start']
        end = aln['subject_end']
        strand = 1
        if aln['sstrand'] == 'minus':
            strand = -1
        # BLAST returns inverted coordinates if the match is on the
        # minus strand - so swap them
        if aln['sstrand'] == 'minus':
            (start, end) = (end, start)
        # Reflect the region if the subject sequence is reverse
        # complemented
        component = uids[aln['subject_id']]
        if transcript_lookup[aln['subject_id']]['strand'] == -1:
            strand = -strand
            length = component['end_position'] - \
                component['start_position'] + 1
            (start, end) = (length - end + 1, length - start + 1)
        # Get genomic coordinates for the alignment by offsetting
        # from the genomic start position of the subject sequence
        partner['aln'] = {
            'query_start': aln['query_start'],
            'query_end': aln['query_end'],
            'chromosome_name':
                transcript_lookup[aln['subject_id']]['chromosome_name'],
            'start_position': component['start_position'] + start - 1,
            'end_position': component['start_position'] + end - 1,
            'strand': strand
            }
        # Locate the break point
        partner['aln']['breakpoint'] = partner['aln']['start_position']
        if (partner['disp'] == '5prime' and strand == 1 or
                partner['disp'] == '3prime' and strand == -1):
            partner['aln']['breakpoint'] = partner['aln']['end_position']

    return


def main(args=None):
    """ Obtain from sys.argv if no args passed """

    if not args and len(sys.argv) > 1:
        args = sys.argv[1:]

    if len(args) != 6:
        print usage()
        sys.exit(1)

    paths = []
    for arg in args:
        paths.append(os.path.abspath(arg))

    (translocations_file,
     sequence_file,
     details_file,
     workdir,
     blast_bin_dir,
     breakpoints_file) = paths

    for file_name in [translocations_file,
                      sequence_file,
                      details_file]:
        if not os.path.isfile(file_name):
            sys.stderr.write("Error - file " + file_name + "\n")
            sys.exit(1)
    for dir_name in [workdir,
                     blast_bin_dir]:
        if not os.path.isdir(dir_name):
            sys.stderr.write("Error - directory " + dir_name + "\n")
            sys.exit(1)

    # Parse the translocations file
    translocations = parse_translocations_file(translocations_file)

    # Load the transcript details_file (YAML)
    stream = open(details_file, 'r')
    transcripts = yaml.load(stream)

    # Locate the makeblastdb executable
    makeblastdb_exec = blast_bin_dir + '/makeblastdb'
    blastn_exec = blast_bin_dir + '/blastn'
    for executable in [makeblastdb_exec, blastn_exec]:
        if not os.path.isfile(executable):
            sys.stderr.write("Error - could not find " + executable + "\n")
            sys.exit(1)

    # Make a copy of the sequence_file in the working area and format
    # it for blastn
    database_file = workdir + '/target.fa'
    log_file = workdir + '/makeblastdb.log'
    shutil.copy(sequence_file, database_file)
    subprocess.check_call([
        makeblastdb_exec,
        '-dbtype', 'nucl',     # Nucleotide sequences
        '-in', database_file,
        '-logfile', log_file
        ])

    # Loop over translocations and attempt to align the junction sequence
    # with the genomic sequence in the loci of the 5' and 3' genes
    mapped_translocations = []
    for translocation in translocations:
        mapped = 0
        for partner in translocation['partners']:
            if partner['hgnc_symbol'] in transcripts:
                mapped += 1
        if mapped < 2:
            continue

        locate_breakpoints(translocation, transcripts,
                           workdir, database_file, blastn_exec)

        mapped_translocations.append(translocation)

        print "%s" % (translocation['xref'])
        print "\tgene\tdisp\tquery\tstrand\tbreakpoint\tsubject"
        for partner in translocation['partners']:
            if 'aln' in partner:
                print "\t%s\t%s\t%d-%d\t%d\t%d\t%s:%d-%d" % (
                    partner['hgnc_symbol'],
                    partner['disp'],
                    partner['aln']['query_start'],
                    partner['aln']['query_end'],
                    partner['aln']['strand'],
                    partner['aln']['breakpoint'],
                    partner['aln']['chromosome_name'],
                    partner['aln']['start_position'],
                    partner['aln']['end_position']
                    )

    # Write mapped translocation results
    output_fh = open(breakpoints_file, 'w')
    output_fh.write(yaml.dump(mapped_translocations))
    output_fh.close()


if __name__ == "__main__":
    main()
