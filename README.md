# Breakage Point Finder

##Overview
This is a "real" Readme file.

## Installation
1. Make a virtual environment (checkout virtualenvwrapper)
——not sure how to do this——

2. ```python setup.py develop```

### Dependencies
The following dependencies are required:
1. Biopython - Developed using version 1.64

2. Blastn
> Download and install BLAST commandline software using the NCBI FTP service,
suitable for the operating platform:
(this was developed using BLAST version 2.2.30)
http://www.ncbi.nlm.nih.gov/books/NBK1763/
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/

> by default, blast will install in:
/usr/local/ncbi/blast/

> create a folder for databases using linux commandline:
sudo mkdir /usr/local/ncbi/blast/db

3. A local copy of the human_genomic BLAST database
> download all the human_genomic database files (human_genomic.??.tar.gz)
using the NCBI FTP service,
ftp://ftp.ncbi.nlm.nih.gov/blast/db/
> Ensure these files are in the /usr/local/ncbi/blast/db folder and unzip all files:
This needs to be done individually for all files:

tar zxvpf human_genomic.00.tar.gz

and may require root access:

sudo tar zxvpf human_genomic.00.tar.gz

## Usage
fragile_finder can be run from a linux commandline using the following
command and options:

python path/fragile_finder.py -t ticdb_file -d distance -o out_file

> ticdb_file = tab-delimited output file from http://www.unav.es/genetica/TICdb/
> distance = maximum distance in kb between fragile sites to identify them as duplicates
	default = 100kb
> out_file = name of .csv output file in the form filename.csv
> path = specifies location of python script

Run the test suite with:
py.test -q --ticdb_file=fragile_finder_test_site.txt --output_file=ff_out.csv --distance=100

This will run 3 tests: test_main, test_parse_blast and test_remove_duplicates
