#!/usr/bin/env perl -w
use warnings;
use strict;
use Data::Dumper;
use Carp;
main();

sub main{
	#multihits on probes something needs to care of this. This is it.
	#use:
	#paste \
	# <(bedtools intersect -wao -bed -a <(\
	#  samtools view -Sb 3294553_IN-14_1545_A1_1_simple.sam) \
	#  -b /groups/umcg-oncogenetics/tmp04/res/probe_ET1262F_2_182.ensembl.bed -S |\
	# perl tickerRefine.pl - ) <(grep -vP '^\@'  3294553_IN-14_1545_A1_1_simple.sam  ) |\
	#tail -n 40
	
	#tail -n 40 || tickertape.pl
	
	warn "Commandline: $0 ".join(' ',@ARGV)."\n";
	
	my $recordsLast;
	my $buffer;
	my $itstimetodie=5;
	
	while(<>){
		my $record;
		$record = DumbReader($_);

		if(defined($recordsLast) && EqualReadRecords($recordsLast,$record)){

			#buffer some
			push(@{$buffer}, @{$recordsLast});

		}elsif(defined($buffer) && not(EqualReadRecords($recordsLast,$record))){

			push(@{$buffer}, @{$recordsLast});
			
			#warn Dumper($buffer)."#";
			
			#die Dumper($buffer).$. ;#if(scalar(@{$buffer}) > 2);
			$recordsLast = RefineOverlapsBuffer($buffer);
			
			#die Dumper($buffer,$record,$recordLast)if($itstimetodie<0);
			#$itstimetodie--;
			
			print Writer($recordsLast);
			$recordsLast=undef;
			
			$buffer=undef;
						
		
		#also read not(EqualReadRecords($recordLast,$record)) below, but not inserted because efficiency
		}elsif(not(defined($buffer)) && defined($recordsLast)){

			#push(@{$buffer}, $recordLast);
			print Writer($recordsLast);
			$recordsLast=undef;
		}elsif(not(defined($buffer)) && not(defined($recordsLast))){
			#just exit and assign $recordLast = $record;
		}else{

			die 'Invalid something';

		}
		
		$recordsLast -> [0] = $record;
	}
	#purge buffer
	push(@{$buffer}, @{$recordsLast});
	$recordsLast = RefineOverlapsBuffer($buffer);
	print Writer($recordsLast);
}

sub EqualReadRecords{
	my $r=shift(@_);
	my $r2=shift(@_);
	#pe is modified#
	if(scalar(@{$r->[0]})){
		#warn "scalar ref reassign of because of nesting ".Dumper($r,$r2);
		$r = $r -> [0];
		#warn "after:". Dumper($r,$r2);
		#die;
	};
	if((GetNameReadNoEnd($r) eq GetNameReadNoEnd($r2))){
		return 1;
	}
	return 0;
}
#sub EqualReadRecords{
#	my $r=shift(@_);
#	my $r2=shift(@_);#
#	if((GetChrRead($r) eq GetChrRead($r2)) 
#		&& (GetStartRead($r) == GetStartRead($r2)) 
#		&& (GetEndRead($r) == GetEndRead($r2)) 
#		&& (GetNameRead($r) eq GetNameRead($r2)) 
#		&& (GetStrandRead($r) eq GetStrandRead($r2))){
#		return 1;
#	}
#	return 0;
#}
sub GetChrRead{
	my $r=shift(@_);
	return($r -> [0]);
}
sub GetStartRead{
	my $r=shift(@_);
	return($r -> [1]);
}
sub GetEndRead{
	my $r=shift(@_);
	return($r -> [2]);
}
sub GetNameRead{
	my $r=shift(@_);
	#warn $r;
	#confess "dem stacktrace".$r if($r eq 12);
	return($r -> [3]);	
}
sub GetNameReadNoEnd{
	my $r=shift(@_);
	my $name;
	($name,undef)= split('/',GetNameRead($r));
	return($name);	
}

sub GetStrandRead{
	my $r=shift(@_);
	return($r -> [5]);	
}

sub GetChrProbe{
	my $r=shift(@_);
	return($r -> [12]);
}
sub GetStartProbe{
	my $r=shift(@_);
	return($r -> [13]);
}
sub GetEndProbe{
	my $r=shift(@_);
	if($r -> [0] eq "." || $r -> [1] == -1){
		return($r -> [15]);# or die 'Record does not contain this many fields!';
	}else{
		#warn $record -> [20];
		return($r -> [14]);# or die 'Record does not contain this many fields!';
	}
}
sub GetNameProbe{
	my $r=shift(@_);
	if($r -> [0] eq "." || $r -> [1] == -1){
		#warn $r -> [16];
		
		return($r -> [16]);# or die 'Record does not contain this many fields!';
	}else{
		return($r -> [15]);# or die 'Record does not contain this many fields!';
	}
}
sub GetStrandProbe{
	my $r=shift(@_);
	if($r -> [0] eq "." || $r -> [1] == -1){
		return($r -> [18]);# or die 'Record does not contain this many fields!';
	}else{
		#warn $record -> [20];
		return($r -> [17]);# or die 'Record does not contain this many fields!';
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
	#warn Dumper($_);
	for my $r (@{$_}){
		print join("\t",@{$r})."\n";
	}
}

sub GetHeader{
	my $record = shift(@_);
	#warn $record -> [3];
	return($record -> [3]);# or die 'Record does not contain this many fields!'.Dumper($record);
}
sub GetH2{
	my $record = shift(@_);
	#warn $record -> [19];
	if($record -> [0] eq "." || $record -> [1] == -1){
		return($record -> [20]);# or die 'Record does not contain this many fields!';
	}else{
		#warn $record -> [20];
		return($record -> [19]);# or die 'Record does not contain this many fields!';
	}
}

sub RefineOverlapsBuffer {
	my $b = shift(@_);
	#warn $b->[0];
	if(HasPairedEndTag($b -> [0])){
		#run for /2 then /1
		my $r1first;
		my $r2first;
		if($r2first=ShiftBufferReadEnd2($b)){
			if(scalar(@{$b}) > 0){
				while(my $r2second = ShiftBufferReadEnd2($b)){
					$r2first=GetBestOverlap($r2first,$r2second);
				}
			}	
		}
		
		#also limit by $r2first probename if found
		#refactor to GetmatchingProbeR1 or something
		my $r1firstNext;
		if($r1firstNext=ShiftBufferReadEnd1($b)){
			#
			
			$r1first=$r1firstNext;
			if(GetNameProbe($r2first) ne '.' && GetNameProbe($r2first) ne GetNameProbe($r1first)){
				while(defined($r1firstNext) && not((GetNameProbe($r2first) eq GetNameProbe($r1first) || '\.' eq GetNameProbe($r1first)))){
					$r1first=$r1firstNext if($r1firstNext = ShiftBufferReadEnd1($b));					
				}
			}
			
			
			
			confess Dumper($r1first,$b,$r2first)." is undef l$." if(not(defined($r1first)));
			
			if(scalar(@{$b}) > 0){
				while(my $r1second = ShiftBufferReadEnd1($b)){
					next if(GetNameProbe($r2first) ne '.' && (GetNameProbe($r1first) ne GetNameProbe($r2first) && '.' ne GetNameProbe($r1first)));
					$r1first=GetBestOverlap($r1first,$r1second);
				}
			}
			confess "$r1first is undef now.$. $!" if(not(defined($r1first)));
		}
		
		#confess Dumper($b).Dumper($r2first).Dumper($r1first).$. if(GetR1Overlap($r1first) > 0||GetR2Overlap($r2first) > 0);
		#push(@{$b},$r1first);
		#push(@{$b},$r2first);
		my $bNew;
		@{$bNew} = ($r1first,$r2first);
		
		return($bNew);
		
	}else{
		my $rfirst= shift(@{$b});
		if(scalar(@{$b}) > 0){
			while(my $rsecond = shift(@{$b})){
				$rfirst=GetBestOverlap($rfirst,$rsecond);
			}
		}
		push(@{$b},@{$rfirst});
		
		my $bNew;
		@{$bNew} = ($rfirst);
		
		return($bNew);
	}

	
}
sub ShiftBufferReadEnd {
	my $b = shift(@_);
	#
	my $readEnd = shift(@_);#1 or 2
	#die "Funky data yes\n".Dumper($b) if($readEnd != 1 || $readEnd != 2 );
	#
	#warn "buffer contains:".Dumper $b;
	for(my $bIndex = 0; $bIndex < scalar(@{$b}); $bIndex++){
	#	#die "Funky data yes\n".Dumper($b) if(GetReadEnd($b -> [$bIndex]) != 1 || GetReadEnd($b -> [$bIndex]) != 2 );
	#	#return(splice(@{$b},$bIndex,1))if(GetReadEnd($b -> [$bIndex]) eq $readEnd);
		#my @tmpBuffer = shift(@{$b});
		if(GetReadEnd($b -> [$bIndex]) eq $readEnd){
	#		warn "buffer return contains: \n".Dumper($b).Dumper($b -> [$bIndex],$.);
			#too experimental... goes with return			
			return(splice(@{$b},$bIndex,1));
		}
	}
	
	return(undef);
	
}
sub ShiftBufferReadEnd1 {
	return(ShiftBufferReadEnd($_[0],1));
}
sub ShiftBufferReadEnd2 {
	return(ShiftBufferReadEnd($_[0],2));
}
sub IsPairedEnd {
	return(HasPairedEndTag(@_));
}
sub HasPairedEndTag {
	my $r= shift(@_);
	#test
	my @split = split('/',GetNameRead($r));
	if(scalar(@split) > 1){;
		return(1);
	}else{
		return(0);
		
	}
}
sub GetReadEnd {
	my $r= shift(@_);
	my $readend;
	(undef,$readend) = split('/',GetNameRead($r));
	return($readend);
}
sub GetBestOverlap{
	my $r1= shift(@_);
	my $r2= shift(@_);
	my $wiggle = 6;
	
	#has two record or more:
	#GetReadEnd($b -> [$bIndex]) eq $readEnd)
	if(not(HasPairedEndTag($r1))){
		if(GetR1Overlap($r1)<GetR1Overlap($r2)){
			$r1=$r2;
		}
	}else{
		if(GetReadEnd($r1) eq 1){
			#confess "logic not implemented";
			#best part is < abs(overlap - 40 - wiggle/pseudo homology )
			my $overlapR2=GetR1Overlap($r2);
			my $overlapR1=GetR1Overlap($r1);
			#warn "Check this out M01785:319:000000000-APBEA:1:2106:17925:3408 $.:".abs($overlapR1 - 40 - $wiggle ) ."R1/2R". abs($overlapR2 - 40 - $wiggle ).Dumper(\$overlapR1,\$overlapR2,$r1,$r2)if(GetNameRead($r1) =~  m/M01785:319:000000000\-APBEA:1:2106:24472:4693/);
			if(abs($overlapR1 - 40 - $wiggle ) > abs($overlapR2 - 40 - $wiggle )){#wants 40 or less
				$r1=$r2;
			}
		}elsif(GetReadEnd($r1) eq 2){
			#confess "logic not implemented";
			
			my $overlapR2=GetR2Overlap($r2);
			my $overlapR1=GetR2Overlap($r1);
			#keep smallest of pair =5 46
			if(abs($overlapR1 - 40 - $wiggle ) > abs($overlapR2 - 40 - $wiggle )){#wants 40
				#warn "#####################################Check this out M01785:319:000000000-APBEA:1:2106:17925:3408 $.:".Dumper(\$overlapR2,\$overlapR1,$r1,$r2)if(GetNameRead($r1) =~  m/M01785:319:000000000\-APBEA:1:2106:24472:4693/);
				
				$r1=$r2;
			}
		}
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

sub GetR1Overlap{
	my $r = shift(@_);
	my $R2probe;
	$R2probe = shift @_ if(scalar(@_));
	my $overlap=0;
	my $wiggle = 6;
	
	if($R2probe && $R2probe eq GetNameProbe($r)){
		return 0;
	}
	
	
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
		#there should not be overlap here
		#warn Dumper($r, \$overlap);
		#die Dumper($r, \$overlap)."#" if($overlap > 41);
	}
	#warn "Is this an error? ".Dumper($r, \$overlap)."#";
	#die "plz check for errors:".Dumper($r, \$overlap)if($overlap > 46);
	
	return $overlap;
}

#
#has R2?
#find probeoverlap as in:
#	
sub GetR2Overlap{
	my $r= shift(@_);
	my $overlap=0;
	my $wiggleR2 = 6;
	
	
	#strands should be eq
	
	if(	GetStrandRead($r) eq GetStrandProbe($r) 
		&& GetStrandRead($r) ne '.'
		&& GetStrandProbe($r) ne '.'){
		
		if(GetStrandRead($r) eq '+' &&
			GetStartRead($r) - $wiggleR2  <= GetStartProbe($r) && 
			GetEndRead($r) >= GetStartProbe($r) ){
				
			#--->--read-->
			#--- >probe>
			
			$overlap = GetEndProbe($r) - GetStartRead($r);
			#die Dumper(\$overlap,\$overlap,$r)."#";
		}elsif(GetStrandRead($r) eq '-' &&
			GetStartRead($r) <= GetEndProbe($r) && 
			GetEndRead($r) + $wiggleR2  >= GetEndProbe($r)){
				
			#---<read<---
			#--<probe<---
			
			$overlap = GetEndRead($r) - GetStartProbe($r);
			#die Dumper(\$overlap,\$overlap,$r) if($overlap > 41);
		}
	}else{
		#warn Dumper($r, \$overlap);
		#die Dumper($r, \$overlap)."#" if($overlap > 41);
	}
	#warn Dumper($r, \$overlap)."#";
	#die Dumper($r, \$overlap)if($overlap > 41);
	
	return $overlap;
}
#put GetR2Overlap into script
#find a way to trim/select most concardant R1 read whit best R2 read or just Trim R2
#if no R2 return R1
#create nice fastq/sam intermediates
