#!/usr/bin/env perl
use warnings;
use strict;
use Vcf;

use Data::Dumper;

my $use =<<"END";
use: perl $0 dbNSFP.txt.gz in.vcf > out.vcf
uses the vcftools Vcf.pm library so please add it to PERL5LIB enviroment variable
Notes:
USE FOR SNV ONLY! At the moment dbNSFP only contains Single Nucleotide Variants, no indels.
This Adds fields from dbNSFP. dbNSFP should be the tabix indexed one from SNPsift or you can index it yourselfs(see the snpSift webpage how to prep this file).
Has support for multiallelic sites with the | field separator
END
#
$ARGV[0]='/home/terpstramm/Documents/gccclusterdata/data/resources/dbNSFPv2.3/dbNSFP2.3.txt.gz';
$ARGV[1]='/home/terpstramm/Documents/gccclusterdata/data/projects/exomeseq/rccprivsmeta/FFPE_exome_ferronika_2.7.2/vcf/FFPE_exome_ferronika.snpSift.vcf';

warn "## ".localtime(time())." ## started $0 with @ARGV='".join('\'',@ARGV)."'.\n";
die "ERROR:Please check the dbNSFP/vcf file specified!! '$ARGV[0]'\n $use\n$!" if(not(-e $ARGV[0]) && not(-e $ARGV[1]));

my $dbNSFPFile = $ARGV[0];
my $dbNSFPVer = getDbNsfpVersion($dbNSFPFile);#regex finds the dbNSFP version

my $vcfFile = $ARGV[1];

my $vcf = Vcf->new(file=>$vcfFile) or die "cannot open vcf file $vcfFile\n";
$vcf->parse_header();
my $cmd = "gzip -dc $dbNSFPFile";
open( my $headerHandle,'-|', $cmd) or die "Cannot open '$cmd'";
$_ = <$headerHandle>;
chomp($_);
s/\r//g;
my @dbNsfpHeader = split("\t",$_);
close($headerHandle);
#chr	pos(1-coor)	ref	alt	aaref	aaalt	hg18_pos(1-coor)	genename	Uniprot_acc	Uniprot_id	Uniprot_aapos	Interpro_domain	cds_strand	refcodon	SLR_test_statistic 	codonpos	fold-degenerate	Ancestral_allele	Ensembl_geneid	Ensembl_transcriptid	aapos	aapos_SIFT	aapos_FATHMM	SIFT_score	SIFT_score_converted	SIFT_pred	Polyphen2_HDIV_score	Polyphen2_HDIV_pred	Polyphen2_HVAR_score	Polyphen2_HVAR_pred	LRT_score	LRT_score_converted	LRT_pred	MutationTaster_score	MutationTaster_score_converted	MutationTaster_pred	MutationAssessor_score	MutationAssessor_score_convertedMutationAssessor_pred	FATHMM_score	FATHMM_score_converted	FATHMM_pred	RadialSVM_score	RadialSVM_score_converted	RadialSVM_pred	LR_score	LR_pred	Reliability_index	GERP++_NR	GERP++_RS	phyloP	29way_pi	29way_logOdds	LRT_Omega	UniSNP_ids	1000Gp1_AC	1000Gp1_AF	1000Gp1_AFR_AC	1000Gp1_AFR_AF	1000Gp1_EUR_AC	1000Gp1_EUR_AF	1000Gp1_AMR_AC	1000Gp1_AMR_AF	1000Gp1_ASN_AC	1000Gp1_ASN_AF	ESP6500_AA_AF	ESP6500_EA_AF
for my $i (4..(scalar(@dbNsfpHeader)-1)){
				$dbNsfpHeader[$i]=~s/\)|\(| |,|;|:/_/g;
				$vcf->add_header_line({key=>'INFO', ID=>"dbNSFP${dbNSFPVer}_$dbNsfpHeader[$i]",Number=>-1,Type=>'String',Description=>"$dbNsfpHeader[$i] from dbNSFP$dbNSFPVer file $ARGV[0]"});
}
print $vcf->format_header();
while (my $x=$vcf->next_data_hash()){
		warn "## ".localtime(time())." ## running $0 with @ARGV='".join('\'',@ARGV)."' at line $. .\n" if($. =~ m/0000$/);
		
		my $dbNSFP;
		for my $alt (@{$x->{'ALT'}}){
			#even  better get tabixGetMatch to support an array of alt
			my @tabixRet=tabixGetMatch('dbNSFP'=> $dbNSFPFile,'chrom'=> $x->{'CHROM'},'pos'=> $x->{'POS'},'ref'=> $x->{'REF'},'alt'=> $alt);
			if(defined($tabixRet[0])){
				for my $i (4..(scalar(@dbNsfpHeader)-1)){
					push(@{$dbNSFP->{"dbNSFP${dbNSFPVer}_$dbNsfpHeader[$i]"}}, $tabixRet[$i]);
					#$x->{'INFO'}{"dbNSFP__$dbNsfpHeader[$i]"}="$tabixRet[$i]"; 
				}
			}
		}
		#wierd elegance but still elegant.....
		for my $i (4..(scalar(@dbNsfpHeader)-1)){
			#warn "#error at variant with vcf line '$.' (is the data all snps??):".$vcf->format_line($x) if(not(defined($dbNSFP->{"dbNSFP${dbNSFPVer}_$dbNsfpHeader[$i]"})) && defined($x->{'INFO'}{'SNPEFF_IMPACT'}) &&  not( $x->{'INFO'}{'SNPEFF_IMPACT'} eq 'MODIFIER'|| $x->{'INFO'}{'SNPEFF_IMPACT'} eq 'LOW'));
			$x->{'INFO'}{"dbNSFP2.3_$dbNsfpHeader[$i]"}=join('|',@{$dbNSFP->{"dbNSFP${dbNSFPVer}_$dbNsfpHeader[$i]"}}) if(defined($dbNSFP->{"dbNSFP${dbNSFPVer}_$dbNsfpHeader[$i]"})); 
		}
		print $vcf->format_line($x);
		exit if ($. == 10);
}
warn "## ".localtime(time())." ## ended $0 with @ARGV='".join('\'',@ARGV)."'.\n";

sub tabixGetMatch {
	my %input = @_;
	my $cmd = "/home/terpstramm/software/tabix-0.2.6/tabix $input{'dbNSFP'} $input{'chrom'}:$input{'pos'}-$input{'pos'}";
	open(my $in,'-|',$cmd) or die "invalid read from command $cmd";
	my @tabixRet;
	while (<$in>){
		chomp;
		s/\r//g;
		@tabixRet=split("\t");
		die 'Error ref:vcf does not match ref:tabix!!'."line: $. ;cmd: $cmd; ref vcf: $input{'ref'}; alt vcf: $input{'alt'};tabix ref: $tabixRet[2];" if(length($input{'ref'}) == 1 && uc($input{'ref'}) ne uc($tabixRet[2]));
		#warn 'matches'.join(',',@tabixRet) if($input{'alt'} eq $tabixRet[2]);	
		for my $index (0..(scalar(@tabixRet)-1)){
			$tabixRet[$index]=~s/,| /_/g;
			#strip trailing '_' or ';'
			$tabixRet[$index]=~s/^_*|^;*|_*$|;*$//g;
			#also remove unneccery _
			$tabixRet[$index]=~s/;_|_;/;/g;
			#convert ';' to ','
			$tabixRet[$index]=~s/;/,/g;
		}
		#warn 'matches'.join(',',@tabixRet) if($input{'ref'} eq $tabixRet[2]);	
		return @tabixRet if($input{'alt'} eq $tabixRet[3]);	
	}
	return ();
}
sub getDbNsfpVersion {
	my $tmp = $_[0];
	$tmp =~ s/^.*\///g;
	$tmp =~ s/\.txt\.gz$|^dbNSFP//g;
	return $tmp;
}
