package pipeline::util;

use 5.020002;
use strict;
use warnings;
use Vcf;
use Data::Dumper;
use Carp qw(verbose confess carp croak);
require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use pipeline-util ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
#our %EXPORT_TAGS = ( 'all' => [ qw(
#	'TargetVcfReAnnotator'
#) ] );

#our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	LinePosLastHandleLogger
	TargetVcfReAnnotator
);

our $VERSION = "0.8.8";


# Preloaded methods go here.
# Below is stub documentation for your module. You'd better edit it!

sub TargetVcfReAnnotator{
	my $self ;%{$self}= @_ ;
	SelfRequire(%{$self},'req'=> ['targetvcf', 'vcflist']);
	my $targetvcf = Vcf->new(file=> $self -> {'targetvcf'});
	$targetvcf->parse_header();
	print $targetvcf->format_header();
	
	
	my $vcfs;
	for my $vcf (@{$self -> {'vcflist'}}){
		push(@{$vcfs},{'handle'=>Vcf->new(file=> $vcf)});
	        $vcfs -> [-1] -> {handle}->parse_header();
		$vcfs -> [-1] -> {buffer} = {'current' => [],'next' => []};
	}
	
	#iterate target file and collect data for one position across all files
	my $targetposbuffer = {'current' => [],'next' => []};
	
	my $vcfbuffers;
	#warn Dumper($targetvcf->next_data_hash()). ' ';
	while ($targetposbuffer=GetBufferedPosNext('vcf' =>\$targetvcf,'buffer' => $targetposbuffer)){
		#my $loc = $targetposbuffer -> {'current'} -> [0];
		
		for my $vcf (@{$vcfs}){
			warn '60606060';
			$vcf -> {buffer} = GetBufferedPosByLoc('vcf' => \$vcf -> {handle},'buffer' => $vcf -> {buffer}, 'loc'=>  $targetposbuffer -> {'current'} -> [0]);
			#do some annotations
			warn 'target='.
				$targetposbuffer -> {'current'} -> [0] -> {'CHROM'}.':'.
				$targetposbuffer -> {'current'} -> [0] -> {'POS'}.
				" vcf=".$vcf -> {'buffer'} -> {'current'} -> [0] -> {'CHROM'}.':'.
				$vcf -> {'buffer'} -> {'current'} -> [0] -> {'POS'}
		}
		 
	}
	#dump rest of buffer
	#insert code

}

sub SelfRequire {
        #calls like SelfRequire(%{$self}, 'req'=> ['vcf', 'fields']);
        my $self;
        %{$self}=@_ or  confess "input not parseable as hash:".Dumper(@_). ' ';
        #die Dumper($self)." ";
        my $reqs = $self -> {'req'} or die "no requirement array as input";
        for my $tag (@{$reqs}){
                defined($self -> {$tag}) or confess "no requirement '$tag' as input";
        }
}

sub GetBufferedPosNext {
	my $self; %{$self}=@_;
	SelfRequire(%{$self}, 'req'=> ['vcf', 'buffer']);
	my $targetvcf = ${$self -> {'vcf'}};
	         #warn Dumper($targetvcf->next_data_hash()). ' ';

	if(scalar(@{$self -> {'buffer'} -> {'next'}}) > 0){
		$self -> {'buffer'} -> {'current'}=$self -> {'buffer'} -> {'next'};
		@{$self -> {'buffer'} -> {'next'}}= ();
	}
	#what works : defined($self -> {'buffer'} -> {'current'}) or defined($self -> {'buffer'} -> {'current'} -> [0] ? )
	#warn '## Buffer ##'. Dumper($self -> {'buffer'})." ";
	my $x;
	warn "next amount ".scalar(@{$self -> {'buffer'} -> {'next'}})." ";
	while( scalar(@{$self -> {'buffer'} -> {'next'}}) == 0 and $x=$targetvcf -> next_data_hash()){
		#warn '## VC ## '. Dumper($x). ' ';
		if(scalar(@{$self -> {'buffer'} -> {'current'}}) == 0){
			#init and add first val of array
			$self -> {'buffer'} -> {'current'} = [$x] or confess Dumper($x);
			#push(@{$self -> {'buffer'} -> {'current'}},$x);
			warn 'Init';
		}elsif(scalar(@{$self -> {'buffer'} -> {'current'}} ) && 
			$x -> {'POS'} eq $self -> {'buffer'} -> {'current'} -> [0] ->  {'POS'} && 
			$x -> {'CHROM'} eq $self -> {'buffer'} -> {'current'} -> [0] ->  {'CHROM'}){
			warn 'Grow';
			#grow array
			push(@{$self -> {'buffer'} -> {'current'}},$x);
		}else{
			#chrom / pos not the same so dump as warning
			#croak(Dumper($self -> {'buffer'}));
			warn 'Next';
			#assign to next
			$self -> {'buffer'} -> {'next'} = [$x];

			#process next
                }

	}


		#croak(Dumper($self -> {'buffer'})) if (@{$self -> {'buffer'} -> {'current'}} >1);
		#carp "x=".Dumper($x)."vcf=".Dumper($targetvcf).".result## ".scalar(@{$self -> {'buffer'} -> {'current'}})." ";
	if(scalar(@{$self -> {'buffer'} -> {'current'}}) > 0){
		return $self -> {'buffer'} ;
	}else{
		return undef;
	}
}

sub GetBufferedPosByLoc {
	my $self; %{$self}=@_;
	SelfRequire(%{$self}, 'req'=> ['vcf', 'buffer','loc']);
	
	#$self->{'buffer'}=GetBufferedPosNext('vcf'=>$self->{'vcf'},'buffer'=>$self->{'buffer'});
	warn LocGetChromPosAsString('loc'=>  $self -> {'loc'})." ";
	#scalar(@{$self -> {'buffer'} -> {'next'}}) == 0)
	my $count = 0;
	while( (ChromPosIsNext('loc1'=>$self->{'loc'},'loc2'=>$self->{'buffer'} -> {current} -> [0],'vcf'=>$self->{'vcf'}) || 
		not(ChromPosIsEq('loc1'=>$self->{'loc'},'loc2'=>$self->{'buffer'} -> {current} -> [0]))) && 
		scalar($self->{'buffer'}->{'current'}) > 0){
		$count++;
		#croak "test";
		$self->{'buffer'}=GetBufferedPosNext('vcf'=>$self->{'vcf'},'buffer'=>$self->{'buffer'});
		last if(not(ref($self->{'buffer'} -> {current}) eq 'ARRAY' ) );		
		warn join (",\n",($count, LocGetChromPosAsString('loc'=>$self -> {'buffer'} -> {'current'} -> [0]) , LocGetChromPosAsString('loc'=>$self -> {'loc'})))." ";		
		die "problems here" if($count > 10);
	}
	
	if(not(defined(@{$self->{'buffer'} -> {current}})) ){
		return undef;
	}elsif(ChromPosIsEq('loc1'=>$self->{'loc'},'loc2'=>$self->{'buffer'} -> {current} -> [0])){
		#return if equal
		return $self -> {'buffer'};
	}elsif(ChromPosIsNext('loc1'=>$self->{'loc'},'loc2'=>$self->{'buffer'} -> {current} -> [0],'vcf'=>$self->{'vcf'})){
		#return if higher then pos #note the same behavior as above, im not sure what to do here
		return $self -> {'buffer'};
	}else{
		return undef;
	}
}

sub LinePosLastHandleLogger {
	my $self; %{$self}=@_;
	SelfRequire(%{$self}, 'req'=> ['lineinterval', 'timeinterval']);
	if($. % $self -> {'modulo'} == 0){
		warn localtime(time())."\tINFO\t$0:\tIterating records currently at line $.\n" if(time % $self -> {timeinterval} == 0 );
	}
}

sub GetContigsAsArray {
	my $self; %{$self} = @_;
	my @contigs;
	my $vcf = ${$self -> {'vcf'}};
	#warn Dumper($vcf -> {'header_lines'})." ";
	for my $hdata (@{$vcf -> {'header_lines'}}){
		if($hdata -> {'key'} eq 'contig'){
			push(@contigs, $hdata ->{'ID'});
		}
	}
	#die Dumper(@contigs);
	#die Dumper(${$self -> {vcf}} -> get_header_line('key'=>'contig'))." ";
	return @contigs;
}

sub ContigIsNext {
	my $self ;%{$self} = @_;
	
	my $contigcur=0;
	my $contigidx1=-1;my $contigidx2=-0.5;
	#warn Dumper($self -> { 'vcf'}). ' ' if($. > 300);
	for my $contig (GetContigsAsArray('vcf' => $self -> { 'vcf'})){
		$contigidx1=$contigcur if(defined $self -> { 'loc1'} -> {'CHROM'} && $self -> { 'loc1'} -> {'CHROM'} eq $contig);
		$contigidx2=$contigcur if(defined $self -> { 'loc2'} -> {'CHROM'} && $self -> { 'loc2'} -> {'CHROM'} eq $contig);
		$contigcur++;
	}
	return 1 if($contigidx1 <= $contigidx2 );

	return 0;
}

sub ChromPosIsEq {
	my $self;%{$self} = @_;
	if(defined($self -> { 'loc2'})&& $self -> { 'loc1'} -> {'CHROM'} eq $self -> { 'loc2'} -> {'CHROM'} && $self -> { 'loc1'} -> {'POS'} eq $self -> { 'loc2'} -> {'POS'}){
		return 1;
	}else{
		return 0;
	}
}

sub ChromPosIsNext  {
	my $self; %{$self} = @_;
	#carp "loc1".Dumper($self -> { 'loc1'})."loc2".Dumper($self -> { 'loc2'});
	if(defined( $self -> { 'loc2'})&& defined( $self -> { 'loc2'} -> {'CHROM'}) && ($self -> { 'loc1'} -> {'CHROM'} eq $self -> { 'loc2'} -> {'CHROM'} && 
			$self -> { 'loc1'} -> {'POS'} >= $self -> { 'loc2'} -> {'POS'}) || 
		ContigIsNext('loc1'=> $self -> { 'loc1'},'loc2' => $self -> { 'loc2'}, 'vcf' => $self -> {'vcf'}) ){
		warn 1;
		return 1;
	}else{
		return 0;
	}
}
sub LocGetChromPosAsString{
	my $self; %{$self} = @_;
	return $self -> {'loc'} -> {'CHROM'}.":".$self -> {'loc'} -> {'POS'};
}
1;
__END__

=head1 NAME

pipeline::util - Perl cmdline tools for hts pipeline

=head1 SYNOPSIS

execute scripts from bash

=head1 DESCRIPTION

read readme on github

=head2 EXPORT

None by default.

=head1 SEE ALSO

	https://github.com/mmterpstra/pipeline-util

=head1 AUTHOR

m.m.terpstra

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2018 by m.m.terpstra

LGPL licence

=cut

=head2 TargetVcfReAnnotator
    About   : function which lifts over annotations which are lost by merges
    Usage   : TargetVcfReAnnotator('targetvcf'=>'test.vcf','vcflist'=>['a.vcf,b.vcf']);               # (from a script)
    Args    : File names of A the target vcf to annotate and B the list of vcfs to resque annotations from.
=cut


=head2 SelfRequire
    About   : minimal validation tool of subroutine input
    Usage   : SelfRequire(%{$self}, 'req'=> ['vcf', 'fields']); assumes input of self like a hash and checks it the fields are present in this hash
    Args    : 'req' = Fields to present in self hash
=cut

=head2 GetBufferedPosNext
    About   : Buffer all vcf records by position
    Usage   : GetBufferedPosNext(vcf=$vcf;buffer=$buffer); assumes input of self like a hash and checks it the fields are present in this hash
    Args    :   my $targetvcf = Vcf->new(file=> $self -> {'targetvcf'})
	        $targetvcf->parse_header();
	        print $targetvcf->format_header();
	
		#iterate target file and collect data for one position across all files
	        my $posbuffer;
	
	        while ($posbuffer=GetBufferedPosNext('vcf' => $targetvcf,'buffer' => $posbuffer)){
			for my $x in @{$posbuffer-> {'current'}} {
				#do something
			}
		}

=cut

=head2 GetBufferedPosByLoc 
    About   : Get the next parsed vcf records in the given location. 
    Usage   : GetPositionBufferedPos(vcf=>$vcf,buffer=>$buffer,loc=[$chr,$pos]); assumes input of self like a hash and checks it the fields are present in this hash
    Args    : my $targetvcf = Vcf->new(file=> $self -> {'targetvcf'})
	        $targetvcf->parse_header();
	        print $targetvcf->format_header();
	
		#iterate target file and collect data for one position across all files
	        my $posbuffer;
		my $loc = {'CHROM' => 'Y','POS'=> 1000};
	        while ($posbuffer=GetBufferedPosByLoc('vcf' => $targetvcf,'buffer' => $posbuffer,$loc)){
			for my $x in @{$posbuffer-> {'current'}} {
				#do something
			}
		}

=cut

=head2 
    About   : LinePosLastHandleLogger
    Usage   :  LinePosLastHandleLogger();
    Args    :  LinePosLastHandleLogger(lineinterval => 1000,timeinterval=20). 
		Just uses time, localtime and $. to get what it needs. Checks every modulo lineinterval if time matches timeinterval and gives a warning.
=cut
