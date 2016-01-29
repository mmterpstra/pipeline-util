#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
main();

sub main{
	#modify for getting the part of the sequence not overhanging probe:
	#get the overhang.
	#scan the cigar to see how mutch that is in sequence space
	#trim fastq
	#print the fastq
	#while
	#get some coffee x 5423456987634567890987654345678 ^ 234567890987654323456789 ^ 2345678909876545678 ^ 3456789097654
	
	#use:
	#paste \
	# <(bedtools intersect -wao -bed -a <(\
	#  samtools view -Sb 3294553_IN-14_1545_A1_1_simple.sam) \
	#  -b /groups/umcg-oncogenetics/tmp04/res/probe_ET1262F_2_182.ensembl.bed -S |\
	# perl tickerRefine.pl - ) <(grep -vP '^\@'  3294553_IN-14_1545_A1_1_simple.sam  ) |\
	# tickertape.pl
	#in one line:
	# paste <(bedtools intersect -wao -bed -a <(samtools view -Sb 3294553_IN-14_1545_A1_1_simple.sam) -b /groups/umcg-oncogenetics/tmp04/res/probe_ET1262F_2_182.ensembl.bed -S | perl tickerRefine.pl - ) <(grep -vP '^\@'  3294553_IN-14_1545_A1_1_simple.sam  ) | perl tickertape.pl 
	
	
	warn "Commandline: $0 ".join(' ',@ARGV)."\n";
	
	my $recordLast;
	my $buffer;
	my $itstimetodie=5;
	my $record;
	
	while(<>){
		$record = DumbReader($_);
		
		if(my $fq = TrimReadByProbe($record)){
			print WriteFastq($fq);
		}
		#die Dumper($record) if($. > 100);
	}
		
}

sub EqualReadRecords{
	my $r=shift(@_);
	my $r2=shift(@_);
	if((GetChrRead($r) eq GetChrRead($r2)) 
		&& (GetStartRead($r) == GetStartRead($r2)) 
		&& (GetEndRead($r) == GetEndRead($r2)) 
		&& (GetNameRead($r) eq GetNameRead($r2)) 
		&& (GetStrandRead($r) eq GetStrandRead($r2))){
		return 1;
	}
	return 0;
}
sub GetChrRead{
	my $r=shift(@_);
	return($r -> [0]);# or die 'Record does not contain this many fields!'.Dumper($r);
}
sub GetStartRead{
	my $r=shift(@_);
	return($r -> [1]);# or die 'Record does not contain this many fields!'.Dumper($r);
}
sub GetEndRead{
	my $r=shift(@_);
	return($r -> [2]);# or die 'Record does not contain this many fields!'.Dumper($r);
}
sub GetNameRead{
	my $r=shift(@_);
	return($r -> [3]);# or die 'Record does not contain this many fields!'.Dumper($r);	
}
sub GetStrandRead{
	my $r=shift(@_);
	return($r -> [5]);# or die 'Record does not contain this many fields!'.Dumper($r);	
}
sub GetChrProbe{
	my $r=shift(@_);
	return($r -> [12]);# or die 'Record does not contain this many fields!'.Dumper($r);
}
sub GetStartProbe{
	my $r=shift(@_);
	return($r -> [13]);# or die 'Record does not contain this many fields!'.Dumper($r);
}
sub GetEndProbe{
	my $r=shift(@_);
	return($r -> [14]);# or die 'Record does not contain this many fields!'.Dumper($r);
}
sub GetNameProbe{
	my $r=shift(@_);
	return($r -> [15]);# or die 'Record does not contain this many fields!'.Dumper($r);	
}
sub GetStrandProbe{
	my $r=shift(@_);
	return($r -> [17]);# or die 'Record does not contain this many fields!'.Dumper($r);	
}
sub GetQualRead{
	my $r=shift(@_);
	my $ret;
	if($r -> [0] ne "."){
		$ret = $r -> [29];# or die 'Record does not contain this many fields!';
	}else{
		#warn $record -> [20];
		$ret = $r -> [30];# or die 'Record does not contain this many fields!';
	}
	return($ret);	
}
sub GetSeqRead{
	my $r=shift(@_);
	
	my $ret;
	if($r -> [0] ne "."){
		$ret = $r -> [28];# or die 'Record does not contain this many fields!';
	}else{

		$ret = $r -> [29];# or die 'Record does not contain this many fields!';
	}
	#
	die $ret.Dumper($r)."$." if(! ($ret =~ /^[ATCGNatcgn]*$/));
	return($ret);
	
	#return($r -> [28]);# or die 'Record does not contain this many fields!'.Dumper($r);	
}
sub GetFlagRead{
	my $r=shift(@_);
	# $r -> [20]."\n".Dumper($r)."\nqual=".GetQualRead($r)."\nseq=".GetSeqRead($r) if(! looks_like_number($r -> [20]));
	my $ret;
	if($r -> [0] ne "."){
		$ret = $r -> [20];# or die 'Record does not contain this many fields!';
	}else{

		$ret = $r -> [21];# or die 'Record does not contain this many fields!';
	}
	#
	return($ret);
	#return($r -> [20]);# or die 'Record does not contain this many fields!'.Dumper($r);	
}
sub GetCigarRead{
	my $r=shift(@_);
	my $ret;
	if($r -> [0] ne "."){
		$ret = $r -> [24];# or die 'Record does not contain this many fields!';
	}else{

		$ret = $r -> [25];# or die 'Record does not contain this many fields!';
	}
	#
	return($ret);
}

sub IsPrimaryAlignment {
	my $r=shift(@_);
	#warn "test".(GetFlagRead($r) & 256);
	if((GetFlagRead($r) & 256)){
		
		return 0;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return 1;
	}
}

sub DumbReader{
	$_= shift(@_);
	chomp;
	my $record;
	@{$record}=split("\t");
	return $record;
}
sub Writer{
	$_= shift(@_);
	join("\t",@{$_})."\n";
}

sub GetHeaderRead{
	my $record = shift(@_);
	#warn $record -> [3];
	return($record -> [3]);# or die 'Record does not contain this many fields!'.Dumper($record);
}
sub GetH2{
	my $record = shift(@_);
	#warn $record -> [19];
	if($record -> [0] ne "."){
		return($record -> [19]);# or die 'Record does not contain this many fields!';
	}else{
		#warn $record -> [20];
		return($record -> [20]);# or die 'Record does not contain this many fields!';
	}
}

sub RefineOverlapsBuffer {
	my $b = shift(@_);
	
	my $r1= shift(@{$b});
	if(scalar(@{$b}) > 0){
		while(my $r2= shift(@{$b})){
			$r1=GetBestOverlap($r1,$r2);
		}
	}
	push(@{$b},$r1);
	return $r1;
}
sub GetBestOverlap{
	my $r1= shift(@_);
	my $r2= shift(@_);
	
	#has two record or more:
	if(Get3PrimeOverlap($r1)<Get3PrimeOverlap($r2)){
		$r1=$r2;
	}
	#die Dumper($r1,$r2)."#";
	return $r1;
}

#notes for sub below
#'12',
#          '56493411',
#          '56493512',
#          'HISEQ-MFG:688:C7Y7WACXX_TCCAAT:1:1101:12442:2813',
#          '60',
#          '-',


#'12',
#          '56493381',
#          '56493431',
#          'ERBB3|NM_001982_exon_23_0_chr12_56493432_f_0_1937047',
#          '1',
#          '+',

#'12',
#          '56493507',
#          '56493557',
#          'ERBB3|NM_001982_exon_24_0_chr12_56493622_f_0_2022598',


#56493381	56493411			56493431	56493507	56493512	56493557
#|			|					|			|			|			|
#>>>>>>>>>>>>>>p1>>>>>>>>>>>>>>>>
#			<<<<<<<<<<<<<<<<<<<<<<read<<<<<<<<<<<<<<<<<<<
#											>>>>>>>>>>>>>>p2>>>>>>>>>
#
##has any 3p overlap -
#pend >= readstart + 1 
#pstart - wiggle <= readstart
#
#aka
#GetEndProbe($r) >= GetStartRead($r) + 1;
#GetStartProbe($r) - $wiggle <= GetStartRead($r);
#$overlap=GetEndProbe($r) - GetStartRead($r) + 1;
#56493381	56493411			56493431	56493507	56493512	56493557
#|			|					|			|			|			|
#<<<<<<<<<<<<<<p1<<<<<<<<<<<<<<<<
#			>>>>>>>>>>>>>>>>>>>>>>read>>>>>>>>>>>>>>>>>>>
#											<<<<<<<<<<<<<<p2<<<<<<<<<
#
##has any 3p overlap +
#pend + wiggle >= readend
#pstart + 1<= readend
#
#aka
#GetEndProbe($r) + $wiggle >= GetEndRead($r);
#GetStartProbe($r) +1 <= GetEndRead($r);
#$overlap=GetEndRead($r)-GetStartProbe($r) + 1;

sub Get3PrimeOverlap{
	my $r= shift(@_);
	my $overlap=0;
	my $wiggle = 0;
	
	if(	GetStrandRead($r) ne GetStrandProbe($r) 
		&& GetStrandRead($r) ne '.'
		&& GetStrandProbe($r) ne '.'){
		
		if(GetStrandRead($r) eq '+'
			&& GetEndProbe($r) + $wiggle >= GetEndRead($r)
			&& GetStartProbe($r) +1 <= GetEndRead($r)){
				
			#--->read>
			#-----<probe<
			
			$overlap = GetEndRead($r) - GetStartProbe($r);
			#die Dumper(\$overlap,\$overlap,$r)."#";
		}elsif(GetStrandRead($r) eq '-'
			&& GetEndProbe($r) >= GetStartRead($r) + 1 
			&& GetStartProbe($r) - $wiggle <= GetStartRead($r)){
				
			#---<read<
			#->probe>
			
			$overlap = GetEndProbe($r) - GetStartRead($r);
			#die Dumper(\$overlap,\$overlap,$r);
		}
	}else{
		#warn Dumper($r, \$overlap);
		#die Dumper($r, \$overlap)."#" if($overlap < 0);
	}
	#warn Dumper($r, \$overlap)."#";
	#die Dumper($r, \$overlap)if($overlap > 0);
	
	return $overlap;
}

sub TrimReadByProbe{
	my $r= shift(@_);
	my $fq;
	my $overlap = Get3PrimeOverlap($r);
	if(IsPrimaryAlignment($r)> 0){
		$fq->[0] = GetHeaderRead($r);
		$fq->[1] = GetSeqRead($r);
		$fq->[2] = GetQualRead($r);
		if($overlap){
			my $trimOffset = CalcTrim($overlap,$r);
			TrimFq($trimOffset,$fq);
			#die Dumper(\$overlap,$r,$fq);
		}
		return $fq;
	}else{
		return undef;
	}
}
sub TrimFq{
	my $trim = shift(@_);
	my $fq = shift(@_);
	
	$fq->[1]=substr($fq->[1],0,length($fq->[1])-$trim);
	
	$fq->[2]=substr($fq->[2],0,length($fq->[2])-$trim);
	
	return $fq;
}
sub CalcTrim{
	my $overlap = shift(@_);
	my $r= shift(@_);
	my $cigar;
	$cigar= CigarReader($r);
	my $trimbasesright=0;
	my $ref;
	#warn "####ref#overl".$ref.'#'.$overlap;
	
	#die id=zfsmsljhrbhfxkjzdhsbrvfkjhzkjhfxdzkjhbf
	#fix this should
	
	while((my $ref = pop(@{$cigar})) && $overlap >= 0){
		my $operation;
		my $amount;
		#warn "ref#overl".$ref.'#'.$overlap;
		($operation,$amount) = %$ref;
		#warn "TATATA".join("\t",($operation,$amount));
		if($operation =~ /^[M\=X]$/){
			if($overlap-$amount>=0){
				$trimbasesright +=$amount;
				$overlap-=$amount;
			}else{
				$trimbasesright +=$overlap;
				$overlap-=$amount;
			}
			
		}elsif($operation =~ /^[IS]$/){
			$trimbasesright += $amount;
		}
	}
	
	#die Dumper(\$overlap,\$trimbasesright,$r,$cigar)if(GetCigarRead($r) =~ /S|I/ && $. > 38);
	return $trimbasesright;
}
sub CigarReader {
	my $r= shift(@_);
	my $c;
	my $cigar;
	@{$c} = split(/([A-Z=])/,GetCigarRead($r));
	while(my $tag = shift(@{$c})){
		if(looks_like_number($tag)){
			my %h = (shift(@{$c}) => $tag);
			push(@{$cigar}, \%h);
		}else{
			my %h = ($tag => 1);
			push(@{$cigar}, \%h);
		}
	}
	return $cigar;
}


sub WriteFastq {
	#my $fqHandle = STDOUT;
	my $fq = shift @_;
	
	return "\@".$fq->[0]."\n".$fq->[1]."\n"."\+\n".$fq->[2]."\n";
	
	#print $fqHandle "\@".$fq->[0]."\n";
	#print $fqHandle $fq->[1]."\n";
	#print $fqHandle "\+\n";
	#print $fqHandle $fq->[2]."\n";
	
	#warn "\@".$fq->[0]."\n";
	#print $fqHandle $fq->[1]."\n";
	#print $fqHandle "\+\n";
	#print $fqHandle $fq->[2]."\n";
	#$fastq->[1]=$seq;
	#$fastq->[0]=$seqHeader;
	#$fastq->[2]=$qual;
}

