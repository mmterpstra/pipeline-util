# Pipeline-util

![alt master](https://travis-ci.org/mmterpstra/pipeline-util.svg?branch=master)

Misc tools for pipelines written in molgenis-compute-5.

These are too a degree expermental tools, before using them on your own data for analysis please check/validate the output. In the design the preference was made for 'more flexable' > 'less safe'.

---

## Perl dependencies

perl scripts require the following libraries:

- **Data::Dumper**
- **Getopt::Std**
- **List::MoreUtils** qw(uniq)
- **List::Util** qw(sum min max)
- **Scalar::Util** qw /looks_like_number/ 
- **Vcf** (the vcftools library for reading vcf files)

## INSTALLATION

To install this module type the following:

   perl Makefile.PL # optionally PREFIX=~
   make
   make test
   make install

Testing /commiting: 

   #edit stuff
   git add file
   git commit -m "description"
   perl Makefile.PL
   make test
   git add Makefile.PL
   git commit -m "sync PM version to git"
   git push
   git push --tags

## DEPENDENCIES

This module requires these other modules and libraries:

software (executables should be available in $PATH):
- R
- DNAcopy (R library)
- samtools
- htslib
- bwa
- vcftools (This contains the Vcf.pm the tools complain about)

perl:  
- **Data::Dumper**
- **Getopt::Std**
- **List::MoreUtils** qw(uniq)
- **List::Util** qw(sum min max)
- **Scalar::Util** qw /looks_like_number/
- not on cpan: **Vcf** (the vcftools library for reading vcf files)


COPYRIGHT AND LICENCE

Put the correct copyright and licence information here.

Copyright (C) 2018 by m.m.terpstra.cluster@gmail.com

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.20.2 or,
at your option, any later version of Perl 5 you may have available.

----

## Script descriptions

What more can I say. Read the descriptions. Also the usemessages from the scripts should be helpful.

#### AddAdBasedAnnotations.pl

- Adds Z scores based on the AD annotation fields. Needs cleanup.

#### dbNSFPAnnotator.pl

- DIY dbnsfp annotator: experimental: works similarly like SNPSifts dbnsfp annotator needs tabix and bgzip.

#### ForGerardTeMeermanVcfAnnotation.pl

- Pipeline annotator script should work out of the box on the dataset

#### GTfGet1000bpExonsBeforeTES.pl

- GTF extraction tool: extracts the 1000 bp in exon space from the genes/transcripts and dumps them to a file: this is 
made especially for 3 prime mrna-seq (lexogen) and derived protocols

#### PlotFloatsOnInterVals.pl

- A beast of legency code used for production: Does some reformatting on Varscan2 output files to perform a basic aCGH 
 like analysis within R using the package 'DNAcopy', making nice plots and making a SEG file for downstream analysis

#### tickerRefine.pl

- Utility script used by 'trimByBed.pl'. See 'trimByBed.pl' for description

#### trimByBed.pl

- Trims sequences based on mapping location and specified bed files. Neatly trims off 'landing probes' used in the SPET
(Nugene) protocol. Wrapper script using 'tickerRefine.pl' and 'tickertape.pl' they need to be in the same folder for it
 to work, also 'bwa mem' , samtools and bedtools should be available in path 

 - Now has an experimental option to specify samfile instead of letting 'bwa mem' doing the dirty work, not everything
 ending with '.sam' should work though, also see version/commit messages on tested aligners.

#### AddAlleleFrequenciesToSeg.pl

- Adds F scores based on the AD annotation fields. Needs cleanup.

#### DumpFieldsForVariantsToTable.pl

- For easy conversion of your vcf to table format using the GATK VariantsToTable, also read the use message.

#### GenerateTableDescriptionByVcfHeader.pl

- It walks trought your VCF file and prints the annotations and the filtering results in a tsv like format.

#### NugeneDigitalSplitter.pl

- This crudly demultipexes data based on insertion of Sequencing info in the header of your file

#### tickertape.pl

- Utility script used by 'trimByBed.pl'. See 'trimByBed.pl' for description

#### VcfaddAdBasedZScores.pl

- Adds Z scores based on the AD annotation fields, needs cleanup. 

#### filterCombinedVariantsForGatk.pl 
 - This filters after CombineVariants giving prio to GATK if two records are found on the same position

#### RecoverSampleAnnotationsAfterCombineVariantsByPosWalk.pl
 - Tries to annotate sample level annotations from other vcfs with the same samples. Used after CombineVariants to get a more exhaustive vcf so that the multiple alt could be annotated.
 - pseudo workflow (code does not work for briefness): 
    - gatkCombinevariants gatkHCcaller.vcf Freebayes.vcf -o combined.vcf
    - RecoverSampleAnnotationsAfterCombineVariants complex.out.vcf combined.vcf gatkHCcaller.vcf Freebayes.vcf > merge.no_complex_anno.vcf
    - gatkHCcaller --genotype_given_alles --alles complex.out.vcf -o regeno.vcf
    - RecoverSampleAnnotationsAfterCombineVariants reallycomplex.out.vcf merge.no_complex_anno.vcf regeno.vcf merge.final.vcf

#### RecoverSampleAnnotationsAfterCombineVariants.pl

 - similar tool as RecoverSampleAnnotationsAfterCombineVariantsByPosWalk.pl but depreciated

#### multiIntersectSeg.pl

 - an .seg adaptation of the `bedtools multiinter` tool from the bedtools suite also generates plots and all possible combinations of plots from the seg files

### RenameChromosomes.pl

 - sloppy liftover tool if you are switching between 'chr' and non 'chr' prefixed chromosome versions of the same build this is your friend (mainly grch38 derived builds)

### AdFilter.pl
 - filters on AD by count and frequency cutoffs the filters are applied directly on the sample or on top of other samples their AD values(eg. if one of the control samples has 8 reads and a frequency of 0.5 then with count 4 and freq 0.1 the effective filter is 12 reads and 0.6 freq minimum or else it is filtered). Wip considering to add a either a filter based on SD in control samples or adding a freequency multiplier...



-----
