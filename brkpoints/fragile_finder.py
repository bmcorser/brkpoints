"""Name: fragile_finder.py
By: Neil Humphryes 18-Jan-2015
Purpose: To identify genomic coordinates of fragile Site origins of
fused sequences resulting from translocations
Input: Tab-delimited TICdb output
Output: de-duplicated csv file of most likely fragile site locations"""

####################################################
# Modules
import csv
import optparse
import os
import subprocess

from Bio.Blast import NCBIXML

####################################################
# where the local biodata lives on the host system
BIODATA = os.environ.get('BIODATA', '/opt/biodata')
DEFAULT_DB = '/usr/local/share/Blast/db/human_genomic'


# Functions
def parse_ticdb_file(ticdb_file):
    """Parse tab-delimited TICdb file and convert into fasta format
    because fasta is convenient for BLAST input"""
    with open(ticdb_file) as file:
        contents = file.readlines()
    tic_entries = [line.strip().split('\t') for line in contents]
    fasta = []
    for entry in tic_entries:
        fasta.append(">" + '_'.join(entry[0:3]) + '\n' + entry[3])
    f = open('TIC.fasta', 'w')
    f.write('\n'.join(i for i in fasta))
    f.close()


def blast_seqs(input_fastafile, outfile_xml, db_path=DEFAULT_DB):
    """Run nucleotide Blast on fasta entries and output as XML, used blast
    to get best alignments rather than relying on perfect alignments or
    allowing for mismatches using BWA. Search returns 100 alignments, which
    should be sufficient to cover both sides of translocation junction"""
    # if running on a local server:
    subprocess.call(['blastn', '-query', input_fastafile, '-db',
                     db_path, '-outfmt', '5',
                     '-out', outfile_xml, '-num_alignments', '100'])


def parse_blast(blast_xml_file):
    """Parse BLAST XML output and compile relevant information. Query and
    subject positions are obtained from BLAST output, and these alignments
    are sorted into which side of the translocation junction they align to.
    The first (best) alignments from these lists are compiled into a list,
    one for each fragile site (two per translocation) """
    blast_xml = open(blast_xml_file)
    blast_results = {}
    n = 1
    # extract relevant data from XML
    for result in NCBIXML.parse(blast_xml):
        name = str(n) + '_' + str(result.query)
        n = n+1
        length = result.query_length
        alns = []
        for aln in result.alignments:
            if aln.hsps[0].expect < 5:
                if aln.hit_def.find('Homo sapiens chromosome ') >= 0:
                    alns.append([name, length, aln.hit_id, aln.hit_def,
                     aln.hsps[0].query_start, aln.hsps[0].query_end,
                     aln.hsps[0].sbjct_start, aln.hsps[0].sbjct_end,
                     aln.hsps[0].expect])
                else:
                    continue
        # make dictionary, one entry per translocation
        blast_results[name] = alns
    # compile best alignments into lists
    fragiles = []
    for hits in blast_results:
        best = {}
        # need to obtain best alignment for each side of fusion junction
        for fragment in blast_results[hits]:
            # Left-side alignment
            if fragment[4] <= 2:
                # only get best alignment, don't duplicate
                if best.get('left'):
                    continue
                else:
                    best['left'] = fragment
                    # get Chromosome number
                    chr1 = str(fragment[3]).split('Homo sapiens chromosome ')
                    chr2 = chr1[1].split(' ')
                    chr3 = chr2[0].split(',')
                    # determine strand and compile data
                    if fragment[7] - fragment[6] > 0:
                        fragiles.append([fragment[0]+'_a', chr3[0],
                                           fragment[6], fragment[7], 1])
                    else:
                        fragiles.append([fragment[0]+'_a', chr3[0],
                                           fragment[7], fragment[6], -1])
            # Right-side alignment
            elif fragment[5] >= fragment[1]-5:
                # only get best alignment, don't duplicate
                if best.get('right'):
                    continue
                else:
                    best['right'] = fragment
                    # get Chromosome number
                    chr1 = str(fragment[3]).split('Homo sapiens chromosome ')
                    chr2 = chr1[1].split(' ')
                    chr3 = chr2[0].split(',')
                    # determine strand and compile data:
                    # name (+_a for left, _b for right), chomosome number,
                    # abs start position,  abs end position, strand
                    if fragment[7] - fragment[6] > 0:
                        fragiles.append([fragment[0]+'_b', chr3[0],
                                           fragment[6], fragment[7], 1])
                    else:
                        fragiles.append([fragment[0]+'_b', chr3[0],
                                           fragment[7], fragment[6], -1])
    return fragiles


def remove_duplicates(fragile_list, distance):
    """Remove duplicated sites from list, which are defined by sites
    that are located less than the specified distance in kb apart,
    default = 100kb"""
    X = int(distance)*1000
    by_chr = {}
    # organise fragile sites by chromosome
    # one entry of by_chr dict per chromosome
    for site in fragile_list[1:]:
        name = site[0].split('_')
        if name[4] == 'a':
            gene = name[1]
        elif name[4] == 'b':
            gene = name[2]
        if by_chr.get(gene + '_' + site[1]):
            by_chr[gene + '_' + site[1]].append(site)
        else:
            by_chr[gene + '_' + site[1]] = [site]
    # find duplicates, one chromosome at a time
    dups = []
    for chr in by_chr:
        positions = []
        for i in by_chr[chr]:
            # make a list of sites within range of each site
            # i = each site to go through one by one
            # j = same sites as i to compare each site to i
            # add to positions list if midpoint of j is less than Xbp from i
            # if duplicates are found they will make identical lists for each
            # site in i
            positions = [j for j in by_chr[chr] if
                     (int(i[2])+int(i[3]))/2 >= (int(j[2])+int(j[3]))/2 - X and
                      (int(i[2])+int(i[3]))/2 <= (int(j[2])+int(j[3]))/2 + X]
            # form one entry from each list covering full range of site
            # k = each site in the duplicated positions list.
            # obtain minimum and maximum genomic coordinates across all
            # duplicated entries.
            dups.append([chr,
             '|'.join([positions[k][0] for k in range(0, len(positions))]),
             positions[0][1],
             min([min([int(positions[k][2]) for k in range(0, len(positions))]),
                  min([int(positions[k][3]) for k in range(0, len(positions))])
                  ]), max([
              max([int(positions[k][2]) for k in range(0, len(positions))]),
              max([int(positions[k][3]) for k in range(0, len(positions))])]),
             max([int(positions[k][4]) for k in range(0, len(positions))])])
    # remove duplicates
    # duplicated sites will form identical entries from previous loop
    # so can easily be removed by making a dictionary of unique entries
    no_dup = {}
    for site in dups:
        if no_dup.get(site[1]):
            continue
        else:
            no_dup[site[1]] = site
    return no_dup


def write_csv(fragile_dict, out_file):
    """Output fragile sites dictionary to csv file with header. Each unique
    fragile site is given a name, details of original sites that localise to
    this feature, chromosome number, absolute start, absolute end, strand"""
    with open(out_file, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['Fragile_Site_Identifier', 'Fragile_site_entries',
                         'Chromosome', 'Start', 'End', 'Strand'])
        writer.writerows([[fragile_dict[i][0], fragile_dict[i][1],
                            fragile_dict[i][2], fragile_dict[i][3],
                            fragile_dict[i][4], fragile_dict[i][5]]
                            for i in fragile_dict])


def main(ticdb_file, out_file='Fragile_Sites.csv',
         distance=100, blast_db=DEFAULT_DB):
    """Parse a file of translocations, split the translocations into putative
    fragile sites and map the fragile sites to the genome. Output the resulting
    de-duplicated list to a csv file.
    """
    # Parse input data
    parse_ticdb_file(ticdb_file)
    # Run BLAST
    blast_seqs('TIC.fasta', 'TICdb_BLASTout.xml', blast_db)
    # Parse BLAST output
    fragiles = parse_blast('TICdb_BLASTout.xml')
    # remove duplicates
    no_dup = remove_duplicates(fragiles, distance)
    # write to .csv
    write_csv(no_dup, out_file)


####################################################
# Main
if __name__ == '__main__':
    # parse object for managing input options
    parser = optparse.OptionParser()

    # essential data, defines commandline options
    parser.add_option('-t', dest='ticdb_file', default='',
                      help="This input is the tab-delimited TICdb output file.")

    parser.add_option('-d', dest='distance', default='100',
                       help="This input is the threshold distance in kb \
                       to determine if two fragile sites are from the same \
                       fragile region.")

    parser.add_option('--genome', dest='blastdb', default=DEFAULT_DB,
                      help="Path to the blast genome database to search")

    parser.add_option('-o', dest='out_file', default='Fragile_Sites.csv',
                      help="This input is the name of the .csv output file.")
    (options, args) = parser.parse_args()

    main(options.ticdb_file,
         options.out_file,
         options.distance,
         options.blastdb)
