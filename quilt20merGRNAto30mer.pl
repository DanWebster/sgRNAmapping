use strict;
use warnings;
use List::Util;
use Data::Dumper;
my $start_time = time();

=head1 
This script will take in a tab-delimited text file of sgRNA data in a defined format
and returns the genomic location (assuming no double-hits in the genome and bowtie2 
mapping was possible) together with the 30mer genomic context 

Written by: Dan Webster
Last Updated/Checked: 7/27/2017
=cut

my $VERSION_NUMBER = "1.0";
my $LAST_UPDATED = "7/27/2017";
my $PROGRAM_CALL = "@ARGV\n\n";

my $USAGE_DESCRIPTION = "
Usage: perl quilt20merGRNAto30mer.pl -i INPUT_FILENAME -t LIBRARY_TYPE -o OUTPUT_FILENAME

Parameters:

**Note** Bowtie2 and Bedtools must be installed and executable from here

-i INPUT_FILENAME Input is list of gRNA sequences in Quilt format with at least unique name and 20mer

Format below:
targetGene	gRNA.name	gRNA.20merSeq	oligo.plate.F	oligo.plate.R	oligo.library
A1BG	A1BG_zhangGeckoV2_0	GTCGCTGAGCTCCGATTCGA	ACCGGTCGCTGAGCTCCGATTCGAAACTCGAATCGGAGCTCAGCGAC	GGAAAGGACGAAACACCGGTCGCTGAGCTCCGATTCGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC
A1BG	A1BG_zhangGeckoV2_1	ACCTGTAGTTGCCGGCGTGC	ACCGACCTGTAGTTGCCGGCGTGAAACGCACGCCGGCAACTACAGGT	GGAAAGGACGAAACACCGACCTGTAGTTGCCGGCGTGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC
A1BG	A1BG_zhangGeckoV2_2	CGTCAGCGTCACATTGGCCA	ACCGCGTCAGCGTCACATTGGCCAAACTGGCCAATGTGACGCTGACG	GGAAAGGACGAAACACCGCGTCAGCGTCACATTGGCCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC

-t type of library (will write to a separate field with this string)

-o OUTPUT_FILENAME

-x INDEX_PATH_FOR_BOWTIE2
(optional, default = /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome assuming Biowulf2 server)

-g GENOMEWIDE_FASTA_PATH
(optional, default = /fdb/genome/hg19/chr_all.fa)

";

############# PARAMETER VALUES ##############

my $INPUT_FILENAME = "";
my $OUTPUT_FILENAME = "";
my $INDEX_PATH_FOR_BOWTIE2 = "/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome";
my $PATH_TO_GENOMEWIDE_FASTA = "/fdb/genome/hg19/chr_all.fa";

my $fastaForBowtieTemp = "fastaForBowtie.fa";
my $bowtieOutput = "gRNAhits.sam";
my $bedtoolsInput = "bedtoolsInput.bed";
my $bedtoolsOutput = "bedtoolsOutput.fa";
my $libraryType = "";

############# GLOBAL VARIABLES ##############

my %gRNAnameToMetadata = ();

############ MAIN EXECUTION ##############

print STDERR "Your parameters: $PROGRAM_CALL\n";
print STDERR "Running version $VERSION_NUMBER last updated $LAST_UPDATED\n";
&parseArguments();

&readInInputFile($INPUT_FILENAME,1000,"INPUT1");

&writeOutFastaForBowtie();

&callBowtie2();

&readInInputFile($bowtieOutput,1000,"BOWTIE_OUTPUT");

#cleanup
`rm $fastaForBowtieTemp`;
`rm $bowtieOutput`;

&getFastaForEachRegion();

&printOutputForSuccessfulGRNAs();

#cleanup
`rm bedtoolsOutput.fa`;
`rm bedtoolsInput.bed`;

&printRunTime();
############## SUBROUTINES ###############

#parse arguments for the program
sub parseArguments
{
	if (scalar(@ARGV) == 0) {die "No arguments provided: See usage description below: \n$USAGE_DESCRIPTION";}
	
	for(my $i = 0; $i < scalar(@ARGV); $i++)
	{
		if ($ARGV[$i] =~/^\-i$/) {$INPUT_FILENAME = $ARGV[$i+1];}
		elsif ($ARGV[$i] =~/^\-o$/) {$OUTPUT_FILENAME = $ARGV[$i+1];}
		elsif ($ARGV[$i] =~/^\-x$/) {$INDEX_PATH_FOR_BOWTIE2 = $ARGV[$i+1];}
		elsif ($ARGV[$i] =~/^\-g$/) {$INDEX_PATH_FOR_BOWTIE2 = $ARGV[$i+1];}
		elsif ($ARGV[$i] =~/^\-t$/) {$libraryType = $ARGV[$i+1];}
	}
	unless($INPUT_FILENAME && $OUTPUT_FILENAME && $libraryType) {die "Didn't have all parameters defined, see usage description below:\n$USAGE_DESCRIPTION";}
	&replaceCarriageReturnsWithNewlines($INPUT_FILENAME);
}

=head1
target.gene	gRNA.name	gRNA.20merSeq	oligo.plate.F	oligo.plate.R	oligo.library
A1BG	A1BG_zhangGeckoV2_0	GTCGCTGAGCTCCGATTCGA	ACCGGTCGCTGAGCTCCGATTCGAAACTCGAATCGGAGCTCAGCGAC	GGAAAGGACGAAACACCGGTCGCTGAGCTCCGATTCGAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC
A1BG	A1BG_zhangGeckoV2_1	ACCTGTAGTTGCCGGCGTGC	ACCGACCTGTAGTTGCCGGCGTGAAACGCACGCCGGCAACTACAGGT	GGAAAGGACGAAACACCGACCTGTAGTTGCCGGCGTGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC
A1BG	A1BG_zhangGeckoV2_2	CGTCAGCGTCACATTGGCCA	ACCGCGTCAGCGTCACATTGGCCAAACTGGCCAATGTGACGCTGACG	GGAAAGGACGAAACACCGCGTCAGCGTCACATTGGCCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC
=cut
sub parseInputLineINPUT1
{
	my $input = shift(@_);
	return if $input =~ /^targetGene/;
	my @fields = split(/\t/,$input);
	die "Expecting 6 fields from this input line: $input" unless scalar(@fields) == 6;
	my $gRNAname = $fields[1];
	my $gRNA20mer = $fields[2];
	unless ($gRNA20mer)
	{
		print STDERR "No valid 20mer for this line: $gRNAname\t$gRNA20mer\n";
		return;
	}
	$gRNAnameToMetadata{$gRNAname}->{"gRNA20mer"} = $gRNA20mer;
}

sub writeOutFastaForBowtie
{
	open(FASTA,">",$fastaForBowtieTemp) or die "Can't open $fastaForBowtieTemp to write to!\n";
	foreach my $gRNAname (keys %gRNAnameToMetadata)
	{
		my $gRNA20mer = $gRNAnameToMetadata{$gRNAname}->{"gRNA20mer"};
		print FASTA ">$gRNAname\n$gRNA20mer\n";
	}
	close FASTA;
}

sub callBowtie2
{
	print STDERR "Executing Bowtie2 to find potential off-target effects...\n";
	`bowtie2 -f -p 16 -x $INDEX_PATH_FOR_BOWTIE2 --local -f -k 10 --very-sensitive-local -L 9 -N 1 -U $fastaForBowtieTemp -S $bowtieOutput`;
}

=head1
@SQ	SN:chr21	LN:48129895
@SQ	SN:chr22	LN:51304566
@SQ	SN:chrX	LN:155270560
@SQ	SN:chrY	LN:59373566
@PG	ID:bowtie2	PN:bowtie2	VN:2.2.6	CL:"/usr/local/apps/bowtie/2-2.2.6/bowtie2-align-s --wrapper 
basic-0 -f -p 16 -x /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome --local -f -k 10 --very-sensiti
ve-local -L 9 -N 1 -S gRNAhits.sam -U gRNA23mers.fa"
TCONS_00000901_12	0	chr1	2476803	255	23M	*	0	0	ATGGCTCAGGAGGACCCTGCAGG	IIIIIIIIIIIIIIIIIIIIIII	AS:i:46	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:23	YT:Z:UU
TCONS_00000901_13	16	chr1	2476818	32	23M	*	0	0	CCTGCAGGCACCCACTACCATGC	IIIIIIIIIIIIIIIIIIIIIII	AS:i:46	XS:i:42	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:23	YT:Z:UU
TCONS_00000901_13	272	chr6	108700353	32	2S21M	*	0	0	CCTGCAGGCACCCACTACCATGC	IIIIIIIIIIIIIIIIIIIIIII	AS:i:42	XS:i:42	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:21	YT:Z:UU
Note that the only values that *should be observed for the flag field are: 0,16,256,272
These correspond to:
0: forward primary read
16: reverse primary read
256: forward secondary read
272: reverse secondary read

In short, don't consider any gRNAs that have any other alignment in the genome with up to 1 mismatch.
=cut
sub parseInputLineSAM
{
	my $input = shift(@_);
	return if $input =~ /^\@/;
	my @fields = split(/\t/,$input);
	#die "Expecting 12 fields from this input line: $input" unless scalar(@fields) == 12;
	
	my $strand = ".";
	my $gRNAname = $fields[0];
	my $flag = $fields[1];
	
	my $chr = $fields[2];
	print STDERR "Expecting chromosome field in 3rd column of SAM line below, but got $chr\n$input\n" unless $chr =~ /^chr/;
	
	my $position = $fields[3];
	
	#only consider alignments with 19 or greater perfect matches
	my $cigar = $fields[5];
	$cigar =~ /(\d\d)M/;
	my $numMatchingBases = $1;
	return unless $numMatchingBases >= 19;
	
	if ($flag == 0) {$strand = "+";}
	elsif ($flag == 16) {$strand = "-";}
	
	#found a non-primary read, meaning potential for off-target effects, get rid of it!
	elsif ($flag == 256 || $flag == 272) 
	{
		delete($gRNAnameToMetadata{$gRNAname});
		return;
	
	}
	
	else {print STDERR "Unaccounted for flag: $flag was observed. Only looking for 0,16,256,272 currently\n"; return;}
	
	unless ($strand eq "+" || $strand eq "-") {print STDERR "For some reason, didn't have strand recognized after passing all filters for: $input\n";return;}
	
	my $start = $position - 1; #because of 0-indexing to go and get this from the UCSC genome browser later
	my $stop = "";
	
	#determine the bed-formatted start and stop, extending in the correct direction to get the Root Score-formatted 30mer: #NNNN-20mer-NGG-NNN
	if ($strand eq "+")
	{
		$stop = $start + 26;
		$start = $start - 4;
	}
	elsif($strand eq "-")
	{
		$stop = $start + 24;
		$start = $start - 6;
	}

	
	#add the strand at the end of the name for use later in determining the correct sequence
	my $bedData = "$chr\t$start\t$stop\t$gRNAname\_$strand";
	
	#You've made it through this stage of filters, you can store metadata and move on!
	if (exists($gRNAnameToMetadata{$gRNAname}))
	{
		$gRNAnameToMetadata{$gRNAname}->{"BED_DATA"} = $bedData;
		$gRNAnameToMetadata{$gRNAname}->{"STRAND"} = $strand;
	}
	
}

sub getFastaForEachRegion
{
	open(BED_OUT,">",$bedtoolsInput) or die "Can't open $bedtoolsInput for writing!\n";
	foreach my $gRNAname (sort keys %gRNAnameToMetadata)
	{
		my $bed = $gRNAnameToMetadata{$gRNAname}->{"BED_DATA"};
		print BED_OUT "$bed\n";
	}
	close BED_OUT;
	
	print STDERR "Calling BedTools, this could take a minute...\n";
	`bedtools getfasta -tab -name -fi $PATH_TO_GENOMEWIDE_FASTA -bed $bedtoolsInput -fo $bedtoolsOutput`;
	&readInInputFile($bedtoolsOutput,100,"GETFASTA");
	
}

#returns the fasta sequence and also its reverseComplement sequence
sub parseInputLineGETFASTA
{
	my $input = shift(@_);
	my @fields = split(/\t/,$input);
	die "Expecting 2 fields from this input line: $input" unless scalar(@fields) == 2;
	my ($gRNAname,$seq) = @fields;
	my $rcSeq = &reverseComplement($seq);
	
	my $gRNA30mer = "";
	if ($gRNAname =~ /_\+$/)
	{
		#remove the strand and then store
		chop($gRNAname);
		chop($gRNAname);
		$gRNAnameToMetadata{$gRNAname}->{"gRNA30mer"} = $seq;
	}
	elsif($gRNAname =~ /_-$/)
	{
		chop($gRNAname);
		chop($gRNAname);
		#This may need to be changed later, but I think getFasta may account for this?
		$gRNAnameToMetadata{$gRNAname}->{"gRNA30mer"} = $rcSeq;
	}
	else 
	{
		print STDERR "Didn't match the gRNA with appended strand\n";
	}
	
	#valid 30mer has format: #NNNN-20mer-NGG-NNN
	#if ($seq =~ /([ATGCatgc]{25}GG[ATGCatgc]{3})/gi) {$valid30mer = $seq;}
	
}

sub printOutputForSuccessfulGRNAs
{
	open(OUT,">",$OUTPUT_FILENAME) or die "Can't open $OUTPUT_FILENAME!\n";
	print OUT "sgRNAname\tsgRNA20merSeq\tsgRNAgenomicContext\tchr_hg19\tstart_hg19\tstop_hg19\tstrand\tcripsrType\n";
	foreach my $gRNAname (sort keys %gRNAnameToMetadata)
	{
		next unless exists($gRNAnameToMetadata{$gRNAname}->{"gRNA30mer"});
		
		my $gRNA20mer = $gRNAnameToMetadata{$gRNAname}->{"gRNA20mer"};
		my $gRNA30mer = $gRNAnameToMetadata{$gRNAname}->{"gRNA30mer"};
		my $bedData = $gRNAnameToMetadata{$gRNAname}->{"BED_DATA"};
		my @beds = split("\t",$bedData);
		my $chr = $beds[0];
		my $start = $beds[1];
		my $stop = $beds[2];
		my $strand = $gRNAnameToMetadata{$gRNAname}->{"STRAND"};
		print OUT "$gRNAname\t$gRNA20mer\t$gRNA30mer\t$chr\t$start\t$stop\t$strand\t$libraryType\n";
	}
	close OUT;
}

#Read in input file 
sub readInInputFile
{
	my $FILE_TO_READ = shift(@_);
	my $COUNTER_DENOMINATOR = shift(@_); #every this many lines, print output to user (a dot)
	my $FLAG = shift(@_); #signal for which parsing subroutine to use
	
	die "Please supply filename and counter denominator to readInInputFile. This was provided: $FILE_TO_READ and $COUNTER_DENOMINATOR\n" unless ($FILE_TO_READ && $COUNTER_DENOMINATOR);
	
	print STDERR "Reading in $FILE_TO_READ..";
	my $counter = 0;
	
	open(IN,"<",$FILE_TO_READ) or die "Can't open $FILE_TO_READ\n";
	while(<IN>)
	{
		chomp;
		if ($counter % $COUNTER_DENOMINATOR == 0) {print STDERR ".";}
		if ($FLAG eq "INPUT1")
		{
			&parseInputLineINPUT1($_);	
		}
		if ($FLAG eq "BOWTIE_OUTPUT")
		{
			&parseInputLineSAM($_);	
		}
		if ($FLAG eq "GETFASTA")
		{
			&parseInputLineGETFASTA($_);	
		}
		
		$counter++;	
	}
	print STDERR "done!\n";
	close IN;
}

sub reverseComplement
{
	my $dna = shift(@_);
	# reverse the DNA sequence
	my $revcomp = reverse($dna);
	# complement the reversed DNA sequence
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}

sub parseInputLineINPUT2
{
	my $input = shift(@_);
	my @fields = split(/\t/,$input);
	die "Expecting [THIS_MANY] fields from this input line: $input" unless scalar(@fields) == 1;
}

sub replaceCarriageReturnsWithNewlines
{
	my $fileToChange = shift(@_);
	`perl -ne 's/\r/\n/g;print;' $fileToChange > $fileToChange\1`;
	`mv $fileToChange\1 $fileToChange`;
}


sub printRunTime
{
	my $end_time = time();
	my $run_time = $end_time - $start_time;
	print STDERR "Script Done!\nThis job took $run_time seconds.\n";
}