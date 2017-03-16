#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Getopt::Std;

main();

sub main {
	my $opts;
	warn "[INFO] cmdline:".join(' ',($0,@ARGV))."\n";
	getopts('r:b:i:o:s:n:', \%{$opts});
	#warn Dumper($opts);
	if(not( defined($opts -> {'b'}) && -e $opts -> {'b'} && defined($opts -> {'o'}) && 
	 ((defined($opts -> {'r'}) && defined($opts -> {'i'}) && -e $opts -> {'i'} ) ||  (defined($opts-> {'s'}) && -e $opts -> {'s'} )) )){
		die "[FATAL] invalid arguments \n".  usemsg();
	}
	
	#gues some more params
	$opts -> {'bin'} = $0;
	$opts -> {'bin'} =~ s!^(.*\/).*!$1!g;
	
	#tmpdir
	if(not(defined($opts -> {'t'}))){
		$opts -> {'t'} = $opts -> {'o'};
		$opts -> {'t'} =~ s!^(.*\/.*)(\.fastq\.gz|\.fq\.gz|\.fqgz)!$1!g;
		
		#warn $opts -> {'t'};
	}
	#warn $opts -> {'bin'};
	
	wrapper($opts);
}
sub wrapper {
	my $opts = shift @_;
	my $createSamFile;
	if(defined($opts -> {'r'}) && defined($opts -> {'i'})){
		$createSamFile = "bwa mem ".$opts -> {'r'}." ".$opts -> {'i'}." > " . $opts -> {'t'} . ".tmp.sam";
	}
	if($opts -> {'s'} && -e $opts -> {'s'}){
		$createSamFile = "cp -v " . $opts -> {'s'} . " " . $opts -> {'t'} . ".tmp.sam";		
	}elsif($opts -> {'s'} && ! -e $opts -> {'s'}){
		die "invalid -s samfile option set this should be a valid/existing samfile.";
	}
	
	if(-e $opts -> {'t'}.".samfifo"){
		warn "[INFO] Removing samfifo '".$opts -> {'t'}.".samfifo'\n";
		unlink($opts -> {'t'}.".samfifo");
	}
	#good luck at debug
	my $cmd="set -x;set -e && " . $createSamFile . "&&  mkfifo ".$opts -> {'t'}.".samfifo;".
	" samtools view -Sb ".$opts -> {'t'}.".tmp.sam |".
	" bedtools intersect -wao -bed -a -  -b ".$opts -> {'b'}."  | ".
	"perl ".$opts -> {'bin'}."tickerRefine.pl - | ".
	"paste - ".$opts -> {'t'}.".samfifo |".
	"perl ".$opts -> {'bin'}."tickertape.pl -1 " . $opts -> {'o'} . "_R1.fq.gz  -2 ".$opts -> {'o'}."_R2.fq.gz -U ".$opts -> {'o'}.".fq.gz -s ".$opts -> {'n'}." &".
	"perl ".$opts -> {'bin'}."refineSam.pl  ".$opts -> {'t'}.".tmp.sam >  ".$opts -> {'t'}.".samfifo; wait && rm -v "  . $opts -> {'t'} . ".tmp.sam "  . $opts -> {'t'} . ".samfifo";
	
	warn "[INFO] system call:". $cmd."\n";
	my $ret;
	@{$ret} = `($cmd )2>&1`;
	
	warn "[INFO] system call results:\n";
	for (@{$ret}){
		$_ =~ s!^!   !g;
		print STDERR $_;
	}
}
sub usemsg {
	my $use =<<"END";
use: $0 [-s SAMALIGNMENT| -r BWAINDEXBASE -i INFASTQGZ] -b TRIMBED -o OUTPREFIX
        BWAINDEXBASE    Index generated with bwa on reference
        TRIMBED         Regions to trim in bed6 format (bed file with the CHR, START, END, ID, SCORE and STRAND fields).
        INFASTQGZ       Input gzipped Fastq file.
        OUTPREFIX      Output gzipped Fastq prefix. Suffix is .fq.gz or _R1.fq.gz/_R2.fq.gz for single or paired end respectively.
	SAMALIGNMENT	Alignment in SAM Format. Specifying this option skips the BWA alignment step and the need for a BWA index / program.
notes
Requires bwa, bedtools, samtools and a unix environment. Takes a Fastq file and trims it based on alignment position respective to specified intervals.
Due to no validity checking also check that these are for the same reference file!

END
	return $use;
}
