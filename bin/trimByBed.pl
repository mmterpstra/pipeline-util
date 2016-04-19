#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Getopt::Std;

main();

sub main {
	my $opts;
	warn "[INFO] cmdline:".join(' ',($0,@ARGV))."\n";
	getopts('r:b:i:o:s:', \%{$opts});
	#warn Dumper($opts);
	if(not( defined($opts -> {'b'}) && -e $opts -> {'b'} && defined($opts -> {'o'}) && defined($opts -> {'i'}) && -e $opts -> {'i'} &&
	 ((defined($opts -> {'r'}) ) ||  (defined($opts-> {'s'}) && -e $opts -> {'s'} )) )){
		die "[FATAL] invalid arguments \n".  usemsg();
	}
	
	#gues some more params
	$opts -> {'bin'} = $0;
	$opts -> {'bin'} =~ s!^(.*\/).*!$1!g;
	
	#tmpdir
	if(not(defined($opts -> {'t'}))){
		$opts -> {'t'} = $opts -> {'o'};
		$opts -> {'t'} =~ s!^(.*\/.*)(\.fastq\.gz|\.fq\.gz|\.fqgz)!$1!g;
		
		warn $opts -> {'t'};
	}
	#warn $opts -> {'bin'};
	
	wrapper($opts);
}
sub wrapper {
	my $opts = shift @_;
	
	my $createSamFile = "bwa mem ".$opts -> {'r'}." ".$opts -> {'i'}." > " . $opts -> {'t'} . ".tmp.sam";
	
	if($opts -> {'s'} && -e $opts -> {'s'}){
		$createSamFile = "cp -v " . $opts -> {'s'} . " " . $opts -> {'t'} . ".tmp.sam";		
	}elsif($opts -> {'s'} && ! -e $opts -> {'s'}){
		die "invalid -s samfile option set this should be a valid/existing samfile.";
	}
	
	#good luck at debug
	my $cmd="set -x;set -e && " . $createSamFile . "&&  mkfifo ".$opts -> {'t'}.".samfifo;".
	" samtools view -Sb ".$opts -> {'t'}.".tmp.sam |".
	" bedtools intersect -wao -bed -a -  -b ".$opts -> {'b'}." -S | ".
	"perl ".$opts -> {'bin'}."tickerRefine.pl - | ".
	"paste - ".$opts -> {'t'}.".samfifo |".
	"perl ".$opts -> {'bin'}."tickertape.pl |".
	"gzip -c > " . $opts -> {'o'} ."&".
	"grep -vP '^\@'  ".$opts -> {'t'}.".tmp.sam >  ".$opts -> {'t'}.".samfifo; wait && rm -v "  . $opts -> {'t'} . ".tmp.sam "  . $opts -> {'t'} . ".samfifo";
	
	warn "[INFO] system call:". $cmd."\n";
	my $ret;
	@{$ret} =`($cmd )2>&1`;
	
	warn "[INFO] system call results:\n";
	for (@{$ret}){
		$_ =~ s!^!   !g;
		print STDERR $_;
	}
}
sub usemsg {
	my $use =<<"END";
use: $0 [-s SAMALIGNMENT| -r BWAINDEXBASE] -b TRIMBED -i INFASTQGZ -o OUTFASTQGZ
        BWAINDEXBASE    Index generated with bwa on reference
        TRIMBED         Regions to trim in bed6 format (bed file with the CHR, START, END, ID, SCORE and STRAND fields).
        INFASTQGZ       Input gzipped Fastq file.
        OUTFASTQGZ      Output gzipped Fastq file.
	SAMALIGNMENT	Alignment in SAM Format. Specifying this option skips the BWA alignment step and the need for a BWA index.
notes
Requires bwa, bedtools, samtools and a unix environment. Takes a Fastq file and trims it based on alignment position respective to specified intervals.
Due to no validity checking also check that these are for the same reference file!

END
	return $use;
}
