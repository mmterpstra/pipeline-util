#!/usr/bin/env perl
use warnings;
use strict;
use Vcf;
use Getopt::Std;
use Data::Dumper;
&main();
sub main {
	#use List::Util qw/max/; my $in;
	
	my $use="$0 -v vcf sample1,fgbio.ReviewConsensusVariants.txt [sample,fgbio.ReviewConsensusVariants.txt] > annotated.vcf
	quick n dirtay fragmentcounts into your variants
";
	my $opts;
	getopts('v:o:', \%{$opts});
	
	warn "## ". localtime(time()). " ## INFO ## Reading into mem the frag counts.\n";
	
	my $sampleVariantToFragCount;#contains variant and frag count
	for my $sampleandtsv (@ARGV){
		my $VariantToFragCount;
		my ($samplename,$tsv) = split(',',$sampleandtsv);

		open(my $tsvHandle,"<",$tsv) 
			or die "[FATAL] Cannot read open tsv file ".$tsv;
		
		my $headerline = <$tsvHandle>;
		warn "## ". localtime(time()). " ## header: ".$headerline;
		my $originalheader="chrom	pos	ref	genotype	filters	A	C	G	T	N	consensus_read	consensus_insert	consensus_call	consensus_qual	a	c	g	t	n\n";
		die "## ". localtime(time()). " ## Not matching headers" if($headerline ne $originalheader);
		#my @header < maybe improve by adding header
		while(my $line = <$tsvHandle>){
			chomp $line;
			my @tsv = split("\t",$line);
			my ($fragLoc,$fragOrientation) = split(' \| ',$tsv[11]);
			$VariantToFragCount -> {$tsv[0].':'.$tsv[1].'-'.$tsv[2]} -> {$tsv[12]} -> {$fragLoc} -> {$fragOrientation}++;
			warn "## ". localtime(time()). " ## This was all that you can see".Dumper($VariantToFragCount)." " if($. == 10);
		}
		warn "## ". localtime(time()). " ## This was all that you can see".Dumper($VariantToFragCount)." ";
		
		for my $loc (keys(%{$VariantToFragCount})){
			for my $alt (keys(%{$VariantToFragCount -> {$loc}})){
				my $fragcount;
				my $fragByOriCount;
				for my $fragment (keys(%{$VariantToFragCount -> {$loc} -> {$alt}})){
					$sampleVariantToFragCount -> {$loc} -> {$samplename} -> {'frags'} ->  {$alt} ++;
					for my $ori (keys(%{$VariantToFragCount -> {$loc} -> {$alt} -> {$fragment}})){
						$sampleVariantToFragCount -> {$loc} -> {$samplename} -> {'fragsbyori'} -> {$alt} -> {$ori} ++
					}
				}
			}
			
		}
		
		
	}
	warn "## ". localtime(time()). " ## This was all that you can see".Dumper($sampleVariantToFragCount)." ";
	
	#warn localtime(time())." [INFO] $0: Parsing header.";
	my $vcf = Vcf->new(file=>$opts -> {'v'}) 
		or die "## ". localtime(time()). " ## ERROR cannot open input vcf. ".$opts -> {'v'}." ";

	my $out;
	if(defined($opts -> {'o'})){
		my $vcfout = $opts -> {'o'};
		open($out,'>',$vcfout) or die 'Cannot write vcf file';
	}else{
		warn localtime(time())." [WARN] $0: no output vcf file specified printing to STDOUT"."\n";
		$out =*STDOUT;
	}

	#warn "parse header";
	$vcf->parse_header();
	$vcf -> add_header_line({'key' => 'FORMAT', 'ID' => "FDP", 'Number' => 'A','Type' => 'String','Description' => "Fragment DP. Using the consensus read workflow these are the unique by start and end position consensus fragment depth supporting the alternate alleles."});
	$vcf -> add_header_line({'key' => 'FORMAT', 'ID' => "FDPF1R2", 'Number' => 'A','Type' => 'String','Description' => "FDP for read pair orientation F1R2."});
	$vcf -> add_header_line({'key' => 'FORMAT', 'ID' => "FDPF2R1", 'Number' => 'A','Type' => 'String','Description' => "FDP for read pair orientation F2R1."});
		
	print {$out} $vcf->format_header();
	warn localtime(time())." [INFO] $0: Iterating records."; my $records=0;
	while (my $x = $vcf->next_data_hash()){
		#if(defined($annheaderline) && scalar(@{$annheaderline})){
		#	FillSNPEFFANNFields('record' => $x,'ann' => $annotations);
		#}
		#warn "## ". localtime(time()). " ## Test record before. ".Dumper($x)." ";
		if(defined($sampleVariantToFragCount -> {$x -> {'CHROM'}.':'.$x -> {'POS'}.'-'.$x -> {'REF'}})){	

			for my $gt (keys(%{$x -> {'gtypes'}})){
				if(defined($sampleVariantToFragCount -> {$x -> {'CHROM'}.':'.$x -> {'POS'}.'-'.$x -> {'REF'}} -> {$gt})){
					my @fdp = ();
					my @fdpF1R2 = ();	my @fdpF2R1 = ();					
					for my $alt (@{$x -> {'ALT'}}){
						if(defined($sampleVariantToFragCount -> {$x -> {'CHROM'}.':'.$x -> {'POS'}.'-'.$x -> {'REF'}} -> {$gt} -> {'frags'} -> {$alt})){
							push(@fdp,$sampleVariantToFragCount -> {$x -> {'CHROM'}.':'.$x -> {'POS'}.'-'.$x -> {'REF'}} -> {$gt} -> {'frags'} -> {$alt});
						}else{ 
							push(@fdp,'.');

						}
						

						if(defined($sampleVariantToFragCount -> {$x -> {'CHROM'}.':'.$x -> {'POS'}.'-'.$x -> {'REF'}} -> {$gt} -> {'fragsbyori'} -> {$alt} -> {'F1R2'})){
							push(@fdpF1R2,$sampleVariantToFragCount -> {$x -> {'CHROM'}.':'.$x -> {'POS'}.'-'.$x -> {'REF'}} -> {$gt} -> {'fragsbyori'} -> {$alt} -> {'F1R2'});
						}else{ 
							push(@fdpF1R2,'.');
						}
						
						if(defined($sampleVariantToFragCount -> {$x -> {'CHROM'}.':'.$x -> {'POS'}.'-'.$x -> {'REF'}} -> {$gt} -> {'fragsbyori'} -> {$alt} -> {'F2R1'})){
							push(@fdpF2R1,$sampleVariantToFragCount -> {$x -> {'CHROM'}.':'.$x -> {'POS'}.'-'.$x -> {'REF'}} -> {$gt} -> {'fragsbyori'} -> {$alt} -> {'F2R1'});
						}else{ 
							push(@fdpF2R1,'.');
						}
					}
					$vcf->add_format_field($x,'FDP'); $x -> {gtypes} -> { $gt } -> {'FDP'}=join(',',@fdp);
					$vcf->add_format_field($x,'FDPF1R2'); $x -> {gtypes} -> { $gt } -> {'FDPF1R2'}=join(',',@fdpF1R2);
					$vcf->add_format_field($x,'FDPF2R1'); $x -> {gtypes} -> { $gt } -> {'FDPF2R1'}=join(',',@fdpF2R1);
				}else{
					#. < default to dot or do nothing due to sparseness
				}
			}
		}
		#die "## ". localtime(time()). " ## Test record after. ".Dumper($x)." ";
		print {$out} $vcf->format_line($x);
		$records++;
#               if(scalar(@{GetAlt($x)})>1){
#				die Dumper($x)."\n".GetAlt($x)."\n".$. ;
#		}

	}
	$vcf->close();
	warn localtime(time())." [INFO] $0: Done. Processed $records";
}


