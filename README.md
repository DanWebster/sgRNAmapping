We first generated a multi-fasta file of the hg19 genome using Bedtools getFasta. These regions and their reverse complements were parsed for spCas9 PAM sites (NGG) and then filtered based on two main criteria: no TTTTT allowed (this is a polI terminator), and no off-target effects for the identified 23-mer gRNA. Off-target determination was established withBowtie2 using the parameters first described in Kearns et al.:
bowtie2 -f -x HG19_GENOME --local -f -k 10 --very-sensitive-local -L 9 -N 1 -U GRNA_23MERS -S GRNA_HITS.sam

# sgRNAmapping
This script will take in a tab-delimited text file of sgRNA data in a defined format
and returns the genomic location (assuming no double-hits in the genome and bowtie2 
mapping was possible) together with the 30mer genomic context.

Usage: perl quilt20merGRNAto30mer.pl -i INPUT_FILENAME -t LIBRARY_TYPE -o OUTPUT_FILENAME

Parameters:

**Note** Bowtie2 and Bedtools must be installed and executable from here

-i INPUT_FILENAME Input is list of gRNA sequences in Quilt format with at least unique name and 20mer

Format below:
targetGene	gRNA.name	gRNA.20merSeq	oligo.plate.F	oligo.plate.R	oligo.library
A1BG	A1BG_zhangGeckoV2_0	GTCGCTGAGCTCCGATTCGA	ACCGGTCGCTGAGCTCCGATTCGAAACTCGAATCGGAGCTCAGCGAC	GGAAAGGACGAAACACCGGTCGCTGAGCTCCGATTCGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC
A1BG	A1BG_zhangGeckoV2_1	ACCTGTAGTTGCCGGCGTGC	ACCGACCTGTAGTTGCCGGCGTGAAACGCACGCCGGCAACTACAGGT	GGAAAGGACGAAACACCGACCTGTAGTTGCCGGCGTGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC
A1BG	A1BG_zhangGeckoV2_2	CGTCAGCGTCACATTGGCCA	ACCGCGTCAGCGTCACATTGGCCAAACTGGCCAATGTGACGCTGACG	GGAAAGGACGAAACACCGCGTCAGCGTCACATTGGCCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC

-t type of library (will write to a separate field with this string, i.e. sp_Cas9_cut)

-o OUTPUT_FILENAME

-x INDEX_PATH_FOR_BOWTIE2
(optional, default = /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome assuming Biowulf2 server)

-g GENOMEWIDE_FASTA_PATH
(optional, default = /fdb/genome/hg19/chr_all.fa)
