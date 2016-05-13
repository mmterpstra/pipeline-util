#Pipeline-util

Misc tools for pipelines written in molgenis-compute-5.

These are too a degree expermental tools, before using them on your own data for analysis please check/validate the output. In the design the preference was made for 'more flexable' > 'less safe'.

---

# Perl dependencies

perl scripts require the following libraries:

- **Data::Dumper**
- **Getopt::Std**
- **List::MoreUtils** qw(uniq)
- **List::Util** qw(sum min max)
- **Scalar::Util** qw /looks_like_number/ 
- **Vcf** (the vcftools library for reading vcf files)

----

# Script Descriptions

What more can I say. Read the descriptions. Also the usemessages from the scripts should be helpful.

#### AddAdBasedAnnotations.pl

- Adds Z scores based on the AD annotation fields. Needs cleanup.

#### dbNSFPAnnotator.pl

- DIY dbnsfp annotator: experimental: works similarly like SNPSifts dbnsfp annotator needs tabix and bgzip.

#### ForGerardTeMeermanVcfAnnotation.pl

- Pipeline annotator script should work out of the box on the dataset

#### GTfGet1000bpExonsBeforeTES.pl

- GTF extraction tool: extracts the 1000 bp in exon space from the genes/transcripts and dumps them to a file: this is 
made especially for 3 prime mrna-seq and derived protocols

## PlotFloatsOnInterVals.pl

- A beast of legency code used for production: Does some reformatting on Varscan2 output files to perform a basic aCGH  like analysis within R using the package 'DNAcopy', making nice plots and making a SEG file for downstream analysis

## tickerRefine.pl

- Utility script used by 'trimByBed.pl'. See 'trimByBed.pl' for description

## trimByBed.pl

- Trims sequences based on mapping location and specified bed files. Neatly trims off 'landing probes' used in the 
Nugene protocol. Wrapper script using 'tickerRefine.pl' and 'tickertape.pl' they need to be in the same folder for it to work, also 'bwa mem' , samtools and bedtools should be available in path
- Now has an experimental option to specify samfile instead of letting 'bwa mem' doing the dirty work, not everything ending with '.sam' should work though, also see version/commit messages on tested aligners.

## AddAlleleFrequenciesToSeg.pl

- Adds F scores based on the AD annotation fields. Needs cleanup.

## DumpFieldsForVariantsToTable.pl

- For easy conversion of your vcf to table format using the GATK VariantsToTable, also read the use message.

## GenerateTableDescriptionByVcfHeader.pl

- It walks trought your VCF file and prints the annotations and the filtering results in a tsv like format.

## NugeneDigitalSplitter.pl

- This crudly demultipexes data based on insertion of Sequencing info in the header of your file

## tickertape.pl

- Utility script used by 'trimByBed.pl'. See 'trimByBed.pl' for description

## VcfaddAdBasedZScores.pl

- Adds Z scores based on the AD annotation fields, needs cleanup. 

-----
