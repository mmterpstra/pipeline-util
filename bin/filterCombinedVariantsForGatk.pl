#!/usr/bin/perl
use warnings;
use strict;
my $last = '1'; 
my $last4;
my $_last;
my $_store;
while(<>){
	my @F = split("\t",$_);
	#if pos are equal to last then remove one of the pos keeping/preferring the set from GATK toolkit
	if(defined($F[1]) && $F[1] eq $last ){
		if($_ =~ m/set=GATK/){
			$_last =$_;
			$_ = '';
		}else{
			$_ = '';
		}
	};
	print $_last if(defined($_last));
	$_store=$_;
	$last=$F[1]; $_last=$_;$last4=$F[4];
}
print $_store;
