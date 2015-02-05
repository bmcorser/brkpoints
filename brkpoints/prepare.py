"""
Split an input file of Ensembl modelled transcripts sequences into
separate exon and intron sequences, keeping track of structure and
coordinates.

Inputs
------
1) transcripts_file
A FASTA format file containing ensembl modelled transcripts and the
unspliced transcript sequence (see function split_transcript below
for expected format of FASTA header).

Outputs
-------
1) sequences_file
A FASTA format file with all the modelled intron and exon sequences
split out independently

2) details_file
A YAML format file recording the transcript structure.

Sequences in the output sequences_file have unique ids cross references
to data in the output details_file.
"""


import sys
import os
import yaml


def usage():
    """ Return usage """
    return (
        "Usage: {prog} " +
        "<transcripts_file> " +
        "<sequence_file> " +
        "<details_file> "
        )


def start_position(component):
    """
    Support for sorting array of objects on key 'start_position'
    """
    return component['start_position']


def split_transcript(fasta_seq, output_fh, transcripts, seq_counter):
    """
    Split a single transcript
    -------------------------
    This routine excises the constituent intron and exon sequences
    and writes them independently to the output file handle.
    Coordinates and other data are collected and stored within the
    transcripts object.

    fasta_seq
    ---------
    The input fasta_seq object (keys 'header', 'body') holds
    the structure and (unspliced) sequence of an ensembl modelled
    transcript.

    fasta_seq header
    ----------------
    Contains transcript specific info followed by matched lists of
    exon id, exon start and exon end e.g.
    ENST00000393293|ABL1|9|130713980|130854122|1|ENSE00001751015;\
    ENSE00001514748|130854067;130713980|130854122;130714455

    fasta_seq sequence
    ------------------
    Is the (unspliced) transcript sequence, reverse complemented if
    the transcript runs on the -1 strand.
    """

    # Parse the header:
    # Create a transcript object to hold coordinates and other
    # transcript specific data; include within it an (initially empty)
    # array to hold the gene model components (exons and introns)
    # that we will excise from the input sequence.
    parts = fasta_seq['header'].split('|')
    hgnc_symbol = parts[1]
    transcript = {
        'name': parts[0],
        'components': [],
        'chromosome_name': parts[2],
        'start_position': int(parts[3]),
        'end_position': int(parts[4]),
        'strand': int(parts[5]),
        }

    # Check sequence length equals coordinate span
    transcript_length = transcript['end_position'] - \
        transcript['start_position'] + 1
    if len(fasta_seq['body']) != transcript_length:
        sys.stderr.write("Error - length mismatch in " +
                         fasta_seq['header'] + "\n")
        sys.exit(1)

    # Deal with the exons, make sure they match up and check
    # for anomalies.
    exon_ids = parts[6].split(';')
    exon_starts = map(int, parts[7].split(';'))
    exon_ends = map(int, parts[8].split(';'))
    exon_count = len(exon_ids)
    if exon_count != len(exon_starts) or exon_count != len(exon_ends):
        sys.stderr.write("Error - exon mismatches in " +
                         fasta_seq['header'] + "\n")
        sys.exit(1)
    for i in range(0, exon_count):
        if exon_starts[i] >= exon_ends[i]:
            sys.stderr.write("Error - exon coordinate inversion\n")
            sys.exit(1)
        if (exon_starts[i] < transcript['start_position'] or
                exon_ends[i] > transcript['end_position']):
            sys.stderr.write("Error - exon coordinates outsite " +
                             "transcript coordinates\n")
            sys.exit(1)

    # Create an exon object for each exon, append to components
    # array within the transcript object
    for i in range(0, exon_count):
        exon = {
            'name': exon_ids[i],
            'start_position': exon_starts[i],
            'end_position': exon_ends[i],
            'type': 'exon'
            }
        transcript['components'].append(exon)

    # Sort those exons on start position
    transcript['components'].sort(key=start_position)
    for i in range(1, exon_count):
        if (transcript['components'][i]['start_position'] <=
                transcript['components'][i - 1]['end_position']):
            sys.stderr.write("Error - exon overlap\n")
            sys.exit(1)

    # Create an intron object between each of those exons and
    # append to the components array within the transcript object
    for i in range(1, exon_count):
        intron = {
            'name': "%s:%d" % (transcript['name'], i),
            'start_position':
                transcript['components'][i-1]['end_position'] + 1,
            'end_position':
                transcript['components'][i]['start_position'] - 1,
            'type': 'intron'
            }
        transcript['components'].append(intron)

    # Write the exons and intron components to the database
    # file - but sort again on start position so the final
    # ordering is sensible
    transcript['components'].sort(key=start_position)
    for component in transcript['components']:
        # Assign a unique id for this component
        seq_counter += 1
        component['uid'] = '%d' % (seq_counter)
        # Write to FASTA header; include other info but only for
        # readibility in command line work
        output_fh.write(">%s %s|%s|%s|%d|%d|%d|%s\n" % (
            component['uid'],
            component['name'],
            hgnc_symbol,
            transcript['chromosome_name'],
            component['start_position'],
            component['end_position'],
            transcript['strand'],
            component['type']
            ))
        offset = component['start_position'] - transcript['start_position']
        length = component['end_position'] - component['start_position'] + 1
        if length < 0:
            sys.stderr.write("Error - inverted coordinates\n")
            sys.exit(1)
        if transcript['strand'] == -1:
            end = component['end_position'] - \
                transcript['start_position'] + 1
            offset = transcript_length - end
        component_sequence = fasta_seq['body'][offset:offset+length]
        if len(component_sequence) != length:
            sys.stderr.write("Error - invalid length\n")
            sys.exit(1)
        for pos in range(0, length, 60):
            output_fh.write(component_sequence[pos:pos+60] + "\n")

    if hgnc_symbol not in transcripts:
        transcripts[hgnc_symbol] = []
    transcripts[hgnc_symbol].append(transcript)

    return seq_counter


def process_transcripts(transcripts_file, sequences_file, details_file):
    """
    Process the ensembl transcripts file
    ------------------------------------
    The input transcripts_file contains the (unspliced) sequence and
    exon/intron structure information for a collection of ensembl
    modelled transcripts, in FASTA format.

    To avoid holding the sequence data in memory, the subroutine
    iterates through the input file line by line and calls a second
    processing routine (split_transcript) once it has finished
    reading each sequence.

    The split_transcript routine inserts coordinates and other
    information extracted from the FASTA headers into a transcripts
    object written to the details_file at the end of the subroutine.
    """

    # Dict to collect the transcript information from split_transcript
    # Keyed on hgnc_symbol
    transcripts = {}

    # Initialise the fasta_seq object
    fasta_seq = None

    # Open input and output files
    input_fh = open(transcripts_file)
    output_fh = open(sequences_file, 'w')

    # Initialise counter (used to assign output sequence ids)
    seq_counter = 0

    # Read input file
    while True:
        line = input_fh.readline()
        if not line:
            break
        stripped = line.rstrip()
        if not stripped:
            continue
        if stripped.startswith('>'):
            if fasta_seq:
                seq_counter = split_transcript(
                    fasta_seq, output_fh, transcripts, seq_counter)
            fasta_seq = {'header': stripped[1:], 'body': ''}
        else:
            if not fasta_seq:
                sys.stderr.write("Error - file format\n")
                sys.exit(1)
            fasta_seq['body'] += stripped
    if fasta_seq:
        seq_counter = split_transcript(
            fasta_seq, output_fh, transcripts, seq_counter)

    output_fh.close()
    print "Wrote %d exons/introns to %s" % (seq_counter, sequences_file)

    # Write transcript details to details_file as YAML
    output_fh = open(details_file, 'w')
    output_fh.write(yaml.dump(transcripts))
    output_fh.close()

    return


def main(args=None):
    """ Obtain from sys.argv if no args passed """

    if not args and len(sys.argv) > 1:
        args = sys.argv[1:]

    if not args or len(args) != 3:
        print usage()
        sys.exit(1)

    paths = []
    for arg in args:
        paths.append(os.path.abspath(arg))

    transcripts_file, sequences_file, details_file = paths

    if not os.path.isfile(transcripts_file):
        sys.stderr.write("Error - file " + transcripts_file + "\n")
        sys.exit(1)

    process_transcripts(transcripts_file, sequences_file, details_file)


if __name__ == "__main__":
    main()
