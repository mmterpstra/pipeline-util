#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use Data::Dumper;

my %opts;
my $use =<<"END";
use:
perl $0 <vcf.table> > <vcf.oneSampleVariantPerLine.table>
prints out two lines if a single Variant is in two samples in the <vcf.table> and generates the Sample table field.

END

for my $file (@ARGV){
	open my $in, '<', $file or die "Read error on file: $file";
	#
	##header
	#
	$_ = <$in>;
	chomp;
	my @header = split("\t");
	
	my %headerToIndex;
	my $i = 0;
	map{my $tmp = $_;$headerToIndex{$tmp} = $i; $i++;}(@header);
	
	my %samples;
	my %GtFields;
	my %generalFields;
	my $generalFieldsIndex = 0;
	for my $field (reverse(@header)){
		#warn $field;
		#exac 
		my $lastdotindexfield = GetLastIndex($field,'.');
		if($lastdotindexfield > 0 && not(((index($field,'exac') >= 0)|| (index($field,'1000gPhase1') >= 0)|| (index($field,'cosmic') >= 0 )||(index($field,'dbNSFP') >= 0 )|| (index($field,'dbsnp') >= 0 ) || (index($field,'DBSNP') >= 0 ) || (index($field,'dbSNP') >= 0 ) || (index($field,'dbSnp') >= 0 ) ||  (index($field,'ID') >= 0 ) || (index($field,'dbSNPBuildID') >= 0 ) || (index($field,'1000g') >= 0 ) || (index($field,'VariantType') >= 0 ) || (index($field,'Samples') >= 0)))){
			my $sample = substr($field,0,$lastdotindexfield);
			$samples{$sample}++;
			my $GtField = substr($field,$lastdotindexfield+1);
			$GtFields{$GtField}++;
		}else{
			$generalFields{$generalFieldsIndex} = $field;
			$generalFieldsIndex++;
		}
	}
	
	#warn Dumper(\%generalFields,\%GtFields);
	#warn "##########################################################\n";
	warn "## $0 WARN:Found Genotype Fields (as a dump):\n".Dumper(\%GtFields)."\n";
	warn "## $0 WARN:Found General Fields (as a dump):\n".Dumper(\%generalFields)."\n";
	
	my @sampleheader;
	map{$_ =~ s/^.*\./SAMPLE./g;push @sampleheader,$_}(@header[($generalFieldsIndex)..($generalFieldsIndex+scalar(values(%GtFields))-1)]);
	print join("\t",(@header[0..($generalFieldsIndex-1)],'Sample',sort(@sampleheader)))."\n";
	#
	##body
	#
	while (<$in>){
		chomp;
		my @line = split("\t");
		for my $sample (sort(keys(%samples))){
			my @lineout;
			# for printing print otherfields:
			for my $field (sort { $a <=> $b } (keys(%generalFields))){
				push(@lineout, $line[$field]);
			}
			if(defined($headerToIndex{'Samples'})){
				#if($line[$headerToIndex{'Samples'}] =~ m/(^$sample$|^$sample,|,$sample$|,$sample,)/){
				
					push @lineout,$sample;
					for my $gtField (sort(keys(%GtFields))){
						#print gt fields
						if($headerToIndex{$sample.'.'.$gtField}){
							push @lineout,$line[$headerToIndex{$sample.'.'.$gtField}]; 
						}
					}
					#warn join("\t",(@line))."\n";
					print join("\t",(@lineout))."\n";
				#}
			}else{
				die "'Samples' field not found check for renaming of a gatk original samples field and/or other errors"
			}
		}
	}
	close $in;
}


sub GetLastIndex {
	my $field = shift @_ or die "Needs  field input";
	my $search = shift @_ or die "Needs input fo search the field";
	
	
	my $fieldrev= $field;
	$fieldrev=reverse($fieldrev);
	
	my $lastindex=length( $field ) - index( $fieldrev , $search ) -1;

	#warn "ERROR: in field='$field';search='$search' forward idx " . index($field,$search) . " reverse idx " . $lastindex  if( index( $field, $search ) > 0);
	
	return $lastindex if(index( $fieldrev , $search ) > 0);
	#replicat the alternative codes
	return index( $field , $search );
}
