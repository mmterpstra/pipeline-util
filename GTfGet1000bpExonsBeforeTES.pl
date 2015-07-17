#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my $use = <<"END";
	simple extractor of the exon annotations at 500 bp distance to the Transcription end site (TES) extracts / converts so that the regions match 500bp or less
	input
	perl $0 ~/Downloads/Homo_sapiens.GRCh37.75.gtf > ~/Downloads/Homo_sapiens.GRCh37.75.TES500.gff
	no validation of input....
END

my $TESoffset=1000;

my $exampleGTF = <<"END";
	#examplefile

chr1	HAVANA	transcript	12010	13670	.	+	.	transcript_id "ENST00000450305.2"; locus_id "chr1:11869-14412W"; gene_id "ENSG00000223972.4"; reads 0.000000; length 632; RPKM 0.000000
chr1	HAVANA	transcript	11869	14409	.	+	.	transcript_id "ENST00000456328.2"; locus_id "chr1:11869-14412W"; gene_id "ENSG00000223972.4"; reads 232.000000; length 1657; RPKM 6.344208
chr1	ENSEMBL	transcript	11872	14412	.	+	.	transcript_id "ENST00000515242.2"; locus_id "chr1:11869-14412W"; gene_id "ENSG00000223972.4"; reads 0.000000; length 1653; RPKM 0.000000
chr1	ENSEMBL	transcript	11874	14409	.	+	.	transcript_id "ENST00000518655.2"; locus_id "chr1:11869-14412W"; gene_id "ENSG00000223972.4"; reads 0.000000; length 1483; RPKM 0.000000
chr1	ENSEMBL	transcript	30366	30503	.	+	.	transcript_id "ENST00000408384.1"; locus_id "chr1:29554-31109W"; gene_id "ENSG00000221311.1"; reads 0.000000; length 138; RPKM 0.000000
chr1	HAVANA	transcript	30267	31109	.	+	.	transcript_id "ENST00000469289.1"; locus_id "chr1:29554-31109W"; gene_id "ENSG00000243485.1"; reads 0.000000; length 535; RPKM 0.000000
chr1	HAVANA	transcript	29554	31097	.	+	.	transcript_id "ENST00000473358.1"; locus_id "chr1:29554-31109W"; gene_id "ENSG00000243485.1"; reads 56.000000; length 712; RPKM 3.563855
chr1	HAVANA	transcript	62948	63887	.	+	.	transcript_id "ENST00000492842.1"; locus_id "chr1:62948-63887W"; gene_id "ENSG00000240361.1"; reads 0.000000; length 940; RPKM 0.000000
chr1	HAVANA	transcript	69091	70008	.	+	.	transcript_id "ENST00000335137.3"; locus_id "chr1:69091-70008W"; gene_id "ENSG00000186092.4"; reads 36.000000; length 918; RPKM 1.776936
chr1	HAVANA	transcript	131104	133923	.	+	.	transcript_id "ENST00000442987.2"; locus_id "chr1:131104-133923W"; gene_id "ENSG00000233750.2"; reads 522.000000; length 2820; RPKM 8.387516
END

my $script;
################################################################################
#no arg cathching
if (scalar(@ARGV)==0){
print $use;
exit(0);
}
#file input#####################################################
my $file = shift(@ARGV);
#my $basename = $file
#$basename =~ s/.fluxCapacitor.txt//;
#check max amount of fields#####################################################
warn " ## info: interation for exon->transcript linking.\n";

open(my $in,'<', $file)or die "$0: failed to open script:'$file'!!!!\n$!\n";
my $transcriptExonLengths;
while(<$in>){
	next if (m/^#/);
	my %gtfline = &readFlux($_);
	$gtfline{'exon_number'}=~s/"//g if($gtfline{'ANNOTATIONTYPE'} eq 'exon');
	#die Dumper(\%gtfline);
	$transcriptExonLengths->{$gtfline{'transcript_id'}}->{$gtfline{'exon_number'}}=$gtfline{'END'}-$gtfline{'START'}+1 if($gtfline{'ANNOTATIONTYPE'} eq 'exon');
	#warn Dumper($transcriptExonMax) if $. > 5000;
	#warn Dumper(\%gtfline);
}
close $in or die "write error cannot close $file\n";


warn " ## info: looking for offsets/exons.\n";

my $transcriptExonMax;
#$transcriptExonMax->{$gtfline{'transcript_id'}}=$gtfline{'exon_number'} if($gtfline{'ANNOTATIONTYPE'} eq 'exon' &&(not(defined($transcriptExonMax->{$gtfline{'transcript_id'}}))||$transcriptExonMax->{$gtfline{'transcript_id'}} < $gtfline{'exon_number'}));

for my $transcript (keys(%{$transcriptExonLengths})){
	my $transcriptTesOffset = 0;
	my @exons=sort { $b <=> $a }(keys(%{$transcriptExonLengths->{$transcript}}));
	my $transcriptTesExon=$exons[0];
	for my $exon (@exons){
		if($transcriptTesOffset < $TESoffset){
			$transcriptTesOffset=$transcriptTesOffset + $transcriptExonLengths->{$transcript}->{$exon};
			$transcriptTesExon=$exon; 
		}
	}
	$transcriptExonMax->{$transcript}->{'exon'}=$transcriptTesExon;
	$transcriptExonMax->{$transcript}->{'offset'}=$transcriptTesOffset;
}

warn " ## info: printing TES offset annotations.\n";


open($in,'<', $file)or die "$0: failed to open script:'$file'!!!!\n$!\n";
while(<$in>){
	next if (m/^#/);
	my %gtfline = &readFlux($_);
	if($gtfline{'ANNOTATIONTYPE'} eq 'exon'){
		my $exon_nr=$gtfline{'exon_number'};
		$exon_nr=~s/"//g;
		print writeFlux(%gtfline) if($transcriptExonMax->{$gtfline{'transcript_id'}}->{'exon'} < $exon_nr);
		if($transcriptExonMax->{$gtfline{'transcript_id'}}->{'exon'} == $exon_nr){
			if($gtfline{"STRAND"} eq '+' && $transcriptExonMax->{$gtfline{'transcript_id'}}->{'offset'} > $TESoffset){
				$gtfline{"START"}=$gtfline{"START"}+($transcriptExonMax->{$gtfline{'transcript_id'}}->{'offset'}-$TESoffset);
			}elsif($gtfline{"STRAND"} eq '-' && $transcriptExonMax->{$gtfline{'transcript_id'}}->{'offset'} > $TESoffset){
				$gtfline{"END"}=$gtfline{"END"}-($transcriptExonMax->{$gtfline{'transcript_id'}}->{'offset'}-$TESoffset);
			}
			print writeFlux(%gtfline);
		}
		
	}
	#$gtfline{'exon_number'}=~s/"//g if($gtfline{'ANNOTATIONTYPE'} eq 'exon');
}
close $in or die "write error cannot close $file\n";

sub readFlux {
	my $inputline = shift @_;
	chomp $inputline;
	
	my @tabdelim;
	@tabdelim = split("\t",$inputline);
	
	my %gtfParsed;
	$gtfParsed{"CHROM"} = $tabdelim[0];
	$gtfParsed{"ANNOTATIONSCOURCE"} = $tabdelim[1];
	$gtfParsed{"ANNOTATIONTYPE"} = $tabdelim[2];
	$gtfParsed{"START"} = $tabdelim[3];
	$gtfParsed{"END"} = $tabdelim[4];
	$gtfParsed{"COL6"} = $tabdelim[5];
	$gtfParsed{"STRAND"} = $tabdelim[6];
	$gtfParsed{"COL8"} = $tabdelim[7];
	$gtfParsed{"INFO"} = $tabdelim[8];
	
	my @semicolomndelimINFO;
	@semicolomndelimINFO = split(';',$gtfParsed{"INFO"});
	
	for my $element (@semicolomndelimINFO){
		my @spacedelim = split(' ',$element);
		if($spacedelim[0] ne "" && $spacedelim[1] ne "" && not($spacedelim[2]) ){
			$gtfParsed{$spacedelim[0]} = $spacedelim[1];
		}elsif($spacedelim[0] eq "" && $spacedelim[1] ne "" && $spacedelim[2] ne "" ){
			$gtfParsed{$spacedelim[1]} = $spacedelim[2];
		}elsif($spacedelim[0] eq "" && $spacedelim[1] ne "" && $spacedelim[2] ne "" ){
			$gtfParsed{$spacedelim[1]} = $spacedelim[2];
		}else{
			die "incomplete datapairs in the INFO column!!!!on line: '$.'\non infotag:'".$gtfParsed{"INFO"}."'\n$!";
		}	
	}
	return %gtfParsed;
}
sub writeFlux {
	my %gtfParsed = @_;
	
	my $info="";
	
	for my $key (keys(%gtfParsed)){
		
		next if($key eq "CHROM" || $key eq "ANNOTATIONSCOURCE" || $key eq "ANNOTATIONTYPE" || $key eq "START" || $key eq "END" || $key eq "COL6" || $key eq "STRAND" || $key eq "COL8" || $key eq "INFO");
		if($info ne ""){
			$info=$info."; ";
		}
		$info=$info.$key.' '.$gtfParsed{$key};
	}
	
	#warn $info;
	return join("\t",($gtfParsed{"CHROM"},$gtfParsed{"ANNOTATIONSCOURCE"},$gtfParsed{"ANNOTATIONTYPE"},$gtfParsed{"START"},$gtfParsed{"END"},$gtfParsed{"COL6"},$gtfParsed{"STRAND"},$gtfParsed{"COL8"},$info))."\n";
}