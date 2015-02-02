
Breakpoints assignment
======================

Where is the 'tab separated file' ...? 

Ah - there is a reference to some database ... of translocations? 

Googling suggests TICdb ... and yes there is a tab separated file
http://www.unav.es/genetica/allseqs_TICdb.txt

This seems to correspond with the description in the assignment,
one translocation per line, etc.

A few experiments:
BLATing a few of the junction sequences at UCSC
http://genome.ucsc.edu/cgi-bin/hgBlat
seems to give the expected thing; hits on different chromosomes, 
one in the vicinity of what looks like the 5' gene, the other 
in the vicinity of the 3' matching the rest, matching the
5' and 3' ends of the junction sequence.

But that's not a great approach - could hit lots of other stuff
that way, pseudogenes; need a map between the gene id and the
genomic sequence.

Had a look at the paper:
Hmm... maybe not a bad idea to download genomic regions for
each partner gene from ensembl and search the junction 
sequence in just those areas; 
this would depend on the gene ids being OK, but they look 
alright ... a few non HGNC symbols such as 'Ig' in there, but
the IDs must have been sorted out reasonably well around the 
time of the original paper ...

Approach 1
----------
Extracted some sequences from ensembl using Biomart 
http://www.ensembl.org/biomart 
... obtained genomic sequence for each locus (test set of 
50 genes)
... presumably at deskgen there would be a local data 
repository for this? 
Or would they use a web API? 
Whatever, this will do for now.

Ah...
looking again at the paper...
it seems some of the junction sequences are fusion mRNAs and
others are genomic fusion sequences
... not at all sure how to tell
doesn't seem to be any way, looking at the translocations
file ... length?

Looks like the paper authors had a quite a faff for the
fusion mRNAs - breakpoints can only be isolated to the
flanking intron ...

!!!!!!!
So the 'end of the alignment' doesn't necessarily correspond
to the fragile site!! Probably doesn't for fusion mRNAs

How to detect that? Perhaps it's easy, the breakpoint is 
near an exon end ... just ignore those? 

Approach 2
----------
The original authors of the paper looked at the exon/intron
boundaries. We'll possibly need similar to make sense of 
the alignments...
So back to Biomart, download the transcripts with the
exon/intron structure in the header.

So a plan
=========
Two programs:
1) to prepare a database of exon and intron sequences,
based on the ensembl downloaded transcripts file
2) to blast the junction sequence for each translocation
agains the exons and introns of the partner genes

Then: if a junction sequence hits introns, but no exons
the alignment break is possibly a reliable indication of
a fragile site.

Question in assignment: why blastn?
----------------------------------
1) Lots of controls - it can be a bit fiddly getting the 
alignments for very short sequences
2) Nice simple output for progamming - the tabular one, 
no real need for a wrapper
3) Defensible - good to implement something well
established that can be explained easily 
BUT maybe nice alternatives would be
1) exonerate suite (Birney et al.)
SW mode may be useful for the short junction sequeces
nice structured output
2) fasta
Also SW mode, quick to set up - no formatted database
3) Don't know about bwa-mem

Programs
========
1) prepare.py
mogrify the ensembl transcripts file into introns and exons, 
output structures in YAML
2) find.py
scan each junction sequence, identify the alignment end points 
(and perhaps do some exon/intron logic to refine)

Procedure
=========
Selected genes for around 30 translocations from the TICdb 
file and extracted corresponding transcript sequences from 
ensembl Biomart, downloaded to:
transcripts_file: ensembl_GRCh38.fa
147 transcript sequences (lots of alternative structures)

Ran:
./prepare.py ensembl_GRCh38.fa sequence_file.fa details_file.yaml

sequences.fa: 2183 introns/exons
details.yaml: corresponding coordiante info etc.

Ran:
./find.py allseqs_TICdb_3.3.txt sequences.fa details.yaml tmp /usr/bin/ find.yaml > find.out

From find.out:
1) 
18594527
	gene	disp	query	strand	breakpoint	subject
	ACSL3	5prime	1-36	1	222900780	2:222900745-222900780
	ETV1	3prime	35-71	-1	13935898	7:13935862-13935898
Tested with BLAT:
junction: ctgtgtcacaccaccttagcctcttgatcgaggaagTGCCTATGATCAGAAGCCACAAGTGGGAATGAGGC
Exact agreement for genomic regions
The query regions overlap by 2bps...
-> Consistent

2) 
17268511
	gene	disp	query	strand	breakpoint	subject
	ABL1	3prime	36-57	1	130854064	9:130854064-130854085
Didn't align the 5'?
junction aacaaggacaggggtcttcggagtCctgaactcagaagcccttcagcggccagtag
BLAT aligns to 33-56 9:130854061-130854084 so quite similar
BLAT also aligns the 3' to 39-56 15:45857800-45857819
but the 5' is in chromosome 22
-> Consistent

3)
AJ131466
	gene	disp	query	strand	breakpoint	subject
	BCR	5prime	17-93	1	23290413	22:23290337-23290413
	ABL1	3prime	94-210	1	130854064	9:130854064-130854180
nice junction ...
junction sequence TGACCATCAATAAGGAAGATGATGAGTCTCCGGGGCTCTATGGGTTTCTGAATGTCATCGTCCACTCAGCCACTGGATTTAAGCAGAGTTCAAAAGCCCTTCAGCGGCCAGTAGCATCTGACTTTGAGCCTCAGGGTCTGAGTGAAGCCGCTCGTTGGAACTCCAAGGAAAACCTTCTCGCTGGACCCAGTGAAAATGACCCCAACCTTT
BLAT agrees exactly at 5'
It's consistent at 3'' but has a large gap in the alignment
-> Consistent

4)
AY034078
	gene	disp	query	strand	breakpoint	subject
	ASPSCR1	5prime	1-87	1	81996846	17:81996760-81996846
	TFE3	3prime	193-300	-1	49034251	X:49034144-49034251
junction AAGAAGTCCAAGTCGGGCCAGGATCCCCAGCAGGAGCAGGAGCAGGAGCGGGAGCGGGATCCCCAGCAGGAGCAGGAGCGGGAGCGGATTGATGATGTCATTGATGAGATCATCAGCCTGGAGTCCAGTTACAATGATGAAATGCTCAGCTATCTGCCCGGAGGCACCACAGGACTGCAGCTCCCCAGCACGCTGCCTGTGTCAGGGAATCTGCTTGATGTGTACAGTAGTCAAGGCGTGGCCACACCAGCCATCACTGTCAGCAACTCCTGCCCAGCTGAGCTGCCCAACATCAAACGG
BLAT agrees on chr17
On chrX it gets  87-300 X:49034144-49038115 
which looks like the rest of it, but there is a BIG gap 
in the alignment.
I wonder if there was another suboptimal alignment for the
rest of the 3' end.
Perhaps the -ungapped flag for blastn is a bit strong.
-> Not inconsistent, but suggests changes to blastn parameters

5) 
AJ131466
	gene	disp	query	strand	breakpoint	subject
	BCR	5prime	17-93	1	23290413	22:23290337-23290413
	ABL1	3prime	94-210	1	130854064	9:130854064-130854180
AJ131467
	gene	disp	query	strand	breakpoint	subject
	BCR	5prime	1-88	1	23289621	22:23289534-23289621
	ABL1	3prime	89-210	1	130854064	9:130854064-130854185
These show the same breakpoint in the 3' gene, but different
break points in the 5' gene.
Does 'fragile' correspond to an extended region?

CONCLUSION
==========

Well looks like the basic program is working OK.

Not sure what to do about the fusion mRNAs. The simple 
association of a break in the alignment with the location 
of a fragile site doesn't work for these.



