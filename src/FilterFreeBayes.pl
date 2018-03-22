#!/usr/bin/perl
use warnings;
use strict;
use Vcf;
use Data::Dumper;

my $use = <<"END";
	use $0 ratio in.vcf out.vcf
		Filters for variants only having an maximum alternate depth <3 in one of the samples and configurable ratio in one of the samples.
	Recommended usage:
		$0 0.1 freebayes.vcf > out.vcf
END

my $ratio = shift @ARGV or die " No ratio option supplied. $use";

die "no valid 'in.vcf' specified on command line\n$use" if(not(defined($ARGV[0])) ||not -e $ARGV[0]);
#main
main($ARGV[0],$ARGV[1],$ratio);

sub main{
	my $vcfin = $_[0];#'/home/terpstramm/Downloads/s12-193.annotated.vcf';
	my $out;
	if($_[1] && not(-e $_[1])){
		my $vcfout = $_[1];
		open($out,'>',$vcfout) or die 'Cannot write vcf file';
	}else{
		warn localtime(time())." [WARN] $0: no output vcf file specified printing to STDOUT"."\n";
		$out =*STDOUT;
	}
	my $ratio = $_[2];
	
	
	warn localtime(time())." [INFO] $0: Parsing header.";
	my $vcf = Vcf->new(file=>$vcfin);
	#warn "parse header";
	$vcf->parse_header();
	
	#interpret fields from cmdline
	die "invalid ratio" if($ratio < 0 || $ratio > 1); 
	
	#Update headers for FORMAT info and also add FT annotation and print header as vcf
	print {$out} $vcf->format_header();
	#die "header check";
	
	#iterate file
	warn localtime(time())." [INFO] $0: Iterating records."; my $records=0;my $filtered = 0; my $filteredbin; my $recordsbin;
	while (my $x=$vcf->next_data_hash()){
		if($. % 100000 ==0){
			warn localtime(time())." [INFO]  $0: at line $. stats for last $recordsbin records and filtered $filteredbin records aka ".sprintf('%.3f',1-$filtered/$records).".";
			$recordsbin = 0;
			$filteredbin = 0;
		}
		#die Dumper $x;
		my $filter = FilterFreeBayesRecord('record' => $x,'mingtdepth'=> 3, 'mingtratio' => $ratio);
		$filtered += $filter;
		$filteredbin += $filter;
		print {$out} $vcf->format_line($x) if(not($filter));
		#die Dumper $x if(FilterFreeBayesRecord('record' => $x,'mingtdepth'=> 3, 'mingtratio' => $ratio));
		$records++;
		$recordsbin++;
	}
	$vcf->close();
	if($records> 0 && $filtered >0){
		 warn localtime(time())." [INFO] $0: Done. Processed $records and filtered $filtered records aka kept ".sprintf('%.3f',1-$filtered/$records).".";
	}else{
		warn localtime(time())." [INFO] $0: Done. Processed $records and filtered $filtered records.";
	}
}

sub SelfRequire {
	my $self;
        %{$self}= @_;
	#die Dumper($self)." ";
	my $reqs = $self -> {'req'} or die "no requirement array as input";
	for my $tag (@{$reqs}){
		defined($self -> {$tag}) or die "no requirement '$tag' as input";
	}
}
sub FilterFreeBayesRecord {
	my $self;
	%{$self}= @_;
	SelfRequire(%{$self}, 'req'=> ['record','mingtdepth','mingtratio']);
	#$filter=1;
	#iterate genotypes
	for my $sample (keys(%{$self -> {'record'} -> {'gtypes'}})){
		#filter
		my @adValues = split(',',$self -> {'record'} -> {'gtypes'} -> {$sample} -> {'AD'}) if(defined($self -> {'record'} -> {'gtypes'} -> {$sample} -> {'AD'}));
		my $adsum = Sum(@adValues);
		next if($adsum == 0);
		#remove ref and test count/ratio
		shift @adValues;
		
		for my $count (@adValues){
			next if($count == 0);
			my $ratio=$adsum/$count;
			# = './.' if($field eq "GT" && $record -> {'gtypes'} -> {$sample} -> {$field} eq '.');
			#do not filter
			
			return 0 if($count/$adsum > $self -> {'mingtratio'} && $count > $self -> {'mingtdepth'});
			#$record -> {'gtypes'} -> {$sample} -> {$caller.$field} = $record -> {'gtypes'} -> {$sample} -> {$field};
		}
	}
	#filter?
	return 1;
}
sub Sum {
	my $sum = 0;
	for my $val ( @_ ) {
		$sum += $val if $val ne '.';#else +0
	}
	return $sum;
}

sub GetAlt{
	my $r = shift @_;
	if(defined($r->{'ALT'})){
		return $r->{'ALT'};
	}
	die "Vcf record doesn't have ALT Field or is strangely annotated!".Dumper($r);
}
