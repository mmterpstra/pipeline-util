#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
use Carp;
#use Getopt::ArgParse;

 
#my $ap = Getopt::ArgParse->new_parser(
#	prog        => 'CollectNugeneLandingProbeMetrics.pl ',
#	description => 'Collects metrics from bam file when trimmed by trimbybed.pl and optionall duplicate marked',
#	epilog      => 'Examples: samtools view sample.bam | CollectNugeneLandingProbeMetrics.pl > sample.metrics.log',
#);
#$ap->add_arg('--bool', '-b', type => 'Bool', dest => 'boo');
 
# Parse an option: '--foo value' or '-f value'
#$ap->add_arg('--duplicates', '-d', 'required' => 0,'dest' => 'duplicates','help' => "Toggle off filtering for duplicates");
#my $opts = $parser->parse_args();


use Getopt::Std;
my $opts;
getopts('d', $opts);


main();

sub main{
	warn "Commandline: $0 ".join(' ',@ARGV)."\n";
	my $landingProbeMetrics;
	while(<>){
		next if(m/^\@/);
		my $record;
		$record = DumbReader($_);
		next if(not($opts->{'d'}) && IsDuplicate($record));
		$landingProbeMetrics = CollectLandingProbeInfo($landingProbeMetrics,$record);
	}
	#purge buffer

	print LandingProbeMetricsAsString($landingProbeMetrics);
}

sub CollectLandingProbeInfo {
	my $landingProbeMetrics = shift @_;
	my $r = shift @_;
	if(IsFirstInPair($r)){
		$landingProbeMetrics -> {'FIRST_IN_PAIR'} -> { GetNugeneProbeName($r) } -> { 'count' }++;
	}elsif(IsSecondInPair($r)){
		$landingProbeMetrics -> {'SECOND_IN_PAIR'} -> { GetNugeneProbeName($r) } -> { 'count' }++;
	}elsif(not(HasPairedEndTag($r))){
		$landingProbeMetrics -> {'UNPAIRED'} -> { GetNugeneProbeName($r) } -> { 'count' }++;
	}else{
		die "poorly formatted sam record???\n".Dumper($r)." ";
	}
	if(not(IsSupplementalAlignment($r)) && HasProperReadPair($r)){
		$landingProbeMetrics -> {'PAIRED'} -> { GetNugeneProbeName($r) } -> { 'count' }++;
	}
	return $landingProbeMetrics;

}

sub LandingProbeMetricsAsString {
	my $landingProbeMetrics = shift @_;
	my $string = "LandingProbe\tFIRST_IN_PAIR_Count\tSECOND_IN_PAIR_Count\tUNPAIRED_Count\tPAIRED_Count\n";
	
	my @lps; 
	for my $orientation (sort(keys(%{$landingProbeMetrics}))){
		for my $lp (sort(keys(%{$landingProbeMetrics -> { $orientation }}))){
			#$string.="$lp\t$orientation\t".$landingProbeMetrics -> { $orientation } -> { $lp } -> {'count'}."\n";
			my $found = 0;
			#warn $lp;
			for my $lpseen (@lps){
				if($lpseen eq $lp){
					$found = 1;
					push(@lps,$lp);
					last;
				}
				
				last if($found == 1);
			}
			if($found == 0){
				push(@lps,$lp);
			}
		}
	}
	my @oris = ('FIRST_IN_PAIR','SECOND_IN_PAIR','UNPAIRED','PAIRED');
	#warn "here is the list of probes ".Dumper(\@lps,\@oris);
	for my $lp (@lps){
		$string .= "$lp"; 
		for my $orientation (@oris){
			if(defined($landingProbeMetrics -> { $orientation } -> { $lp } -> {'count'})){
				$string .= "\t".$landingProbeMetrics -> { $orientation } -> { $lp } -> {'count'};
			}else{
				$string	.= "\t0";
			}
		}
		$string .= "\n";
	}
	return $string;
}

sub GetNugeneProbeName {
	my $r=shift(@_);
	
	my $nugeneprobetag='NU:Z:'; 
	for my $samtag (@{$r}[11..($#{$r})]){
		if(substr($samtag,0,length($nugeneprobetag)) eq $nugeneprobetag ){
			return substr($samtag,length($nugeneprobetag));
		}
	}
	return "NU_TAG_NOT_SET";
}

sub EqualReadsAndMappings{
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
	return($r -> [2]);
}
sub GetStartRead{
	my $r=shift(@_);
	return($r -> [3]);
}
sub GetEndRead{
	my $r=shift(@_);
	my $rEnd;
	if(GetCigarRead($r) ne '*'){
		$rEnd = GetStartRead($r) + GetCigarLength($r) - 1;
	}else{
		$rEnd = 0;
	}
	#die $rEnd;
	return($rEnd);
}
sub GetNameRead{
	my $r=shift(@_);
	return($r -> [0]);	
}

sub GetStrandRead{
	my $r=shift(@_);
	my $ret;
	if(IsReverseAlignment($r)){
		$ret='-';
	}else{
		$ret='+';
	}
	return($ret);	
}
sub IsReverseAlignment {
	my $r=shift(@_);
	
	if((GetFlagRead($r) & 16)){
		
		return 1;
	}else{
		#"die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return 0;
	}
}

sub HasPairedEndTag {
	my $r=shift(@_);
	
	if((GetFlagRead($r) & 1)){
		
		return 1;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return 0;
	}

}
sub HasProperReadPair {
	my $r=shift(@_);
	
	if((GetFlagRead($r) & 1) && (GetFlagRead($r) & 2)){
		return 1;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);
		return 0;
	}

}

sub IsFirstInPair {
	my $r=shift(@_);

		if((GetFlagRead($r) & 1 && GetFlagRead($r) & 64)){
		return 1;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);
		return 0;
	}
}

sub IsSecondInPair {
	my $r=shift(@_);

	if((GetFlagRead($r) & 1 && GetFlagRead($r) & 128)){
		return 1;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);
		return 0;
	}
}

sub IsDuplicate {
	my $r=shift(@_);
	if((GetFlagRead($r) & 1024)){
		return 1;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);
		return 0;
	}
}


sub GetFlagRead{
	my $r=shift(@_);
	my $ret;
	$ret = $r -> [1];# or die 'Record does not contain this many fields!';
	#
	defined($ret) or confess "Invalid record at line ". $. .": ".Dumper($r);
	return($ret);
	#return($r -> [20]);# or die 'Record does not contain this many fields!'.Dumper($r);	
}
sub GetCigarRead{
	my $r=shift(@_);
	my $ret;
	$ret = $r -> [5] or die "Invalid record at line ". $. .": ".Dumper($r);
	#
	return($ret);
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
	confess Dumper($_) if(not(defined($_ -> [0])));
	for my $r (@{$_}){
		print join("\t",@{$r})."\n";
	}
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

sub IsSupplementalAlignment {
	my $r=shift(@_);
	#warn "test".(GetFlagRead($r) & 2048);
	if((GetFlagRead($r) & 2048)){

		return 1;
	}else{
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);
		return 0;
	}
}


sub GetBestOverlap{
	my $r1= shift(@_);
	my $r2= shift(@_);
	
	#has two records or more:
	if(GetFlagRead($r1)>GetFlagRead($r2)){
		$r1=$r2;
	}
	#die Dumper($r1,$r2)."#";
	return $r1;
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
sub GetCigarLength {
	my $r= shift(@_);
	my $clen = 0;
	my $cigar= CigarReader($r);
	for my $tag (@{$cigar}){
		#die Dumper($tag);
		my ($operation)= keys(%{$tag});
		my ($amount)= values(%{$tag});
		#warn "op:'$operation';amnt:'$amount';";
		if($operation =~ m/^[MN\=X]$/){
			$clen += $amount;
			
		}elsif($operation =~ m/^[IS]$/){
			#nothing
		}elsif($operation =~ m/^[D]$/){
			#only reduce overlap not increment trimbases left
			$clen += $amount;
		}else{
			die "invalid cigar operator '$operation'" . Dumper($r).Dumper($tag);
		}
	}
	#warn $clen ."#".Dumper($cigar)."here";
	
	return $clen;	
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

sub GetReadEnd {
	my $r= shift(@_);
	my $readend;
	
	if((GetFlagRead($r) & 64)){
		
		return 1;
	}elsif((GetFlagRead($r) & 128)){
		#die Dumper($r) or die 'Record does not contain this many fields!'.Dumper($r);	
		return 2;
	}else{
		confess "Invalid Pe Tag in flag! Flag=".GetFlagRead($r). Dumper($r)."line".$.;
	}
}

