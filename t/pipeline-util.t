# Before 'make install' is performed this script should be runnable with
# 'make test'. After 'make install' it should work as 'perl pipeline-util.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;


my $exes = [
	'src/AddAdBasedAnnotations.pl',
        'src/AddAlleleFrequenciesToSeg.pl',
        'src/AddFragCountsToVariants.pl',
	'src/AddFragCountsToVariantsSamtools.pl',
        'src/AdFilter.pl',
        'src/CalleriseVcf.pl',
        'src/CollectNugeneLandingProbeMetrics.pl',
        'src/dbNSFPAnnotator.pl',
        'src/DumpFieldsForVariantsToTable.pl',
        'src/filterCombinedVariantsForGatk.pl',
        'src/FilterFreeBayes.pl',
        'src/ForGerardTeMeermanVcfAnnotation.pl',
        'src/GenerateTableDescriptionByVcfHeader.pl',
        'src/GTfGet1000bpExonsBeforeTES.pl',
        'src/multiIntersectSeg.pl',
        'src/MutectAnnotationsToSampleFormat.pl',
        'src/NugeneDigitalSplitter.pl',
        'src/PlotFloatsOnInterVals.pl',
        'src/RecoverSampleAnnotationsAfterCombineVariants.pl',
        'src/refineSam.pl',
        'src/SampleSheetAddVal.pl',
        'src/tickerRefine.pl',
        'src/tickertape.pl',
        'src/trimByBed.pl',
        'src/updateGTPLByAD.pl',
        'src/updateSimilarPos.pl',
        'src/VcfaddAdBasedZScores.pl',
        'src/VcfSnpEffAsGatk.pl',
	'src/VcfQssFix.pl',
	'src/VcfTableExportOneVariantPerSample.pl',
	'src/RenameChromosomes.pl'
];


use Test::More tests => 36;
#(6+scalar(@{$exes}));#1 pm test #5 funtional tests and other syntax tests for scripts
BEGIN { use_ok('pipeline::util') };

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

for my $exe (@{$exes}){
	ok(system('bash -c "echo > /dev/stderr && set -ex && perl -wc '.$exe.'"') == 0, "$exe syntax test");
}
# test bedtrim
ok(system( 'bash -c "echo > /dev/stderr && set -e && ' .
	'perl src/trimByBed.pl -s t/data/bedtrim_before.sam  -b t/data/bedtrim_probe.bed -o t/data/bedtrim_test -n t/data/bedtrim_test && ' . 
	'diff <(samtools view t/data/bedtrim_after.bam ) <(samtools view t/data/bedtrim_test.bam )"' ) == 0, "bedtrim minimal functional test");
# test RecoverSampleAnnotationsAfterCombineVariantsByPosWalk
ok(system('bash -c "echo > /dev/stderr && set -e -o pipefail && '.
	'diff <(export PERL5LIB=\"blib/lib/\":$PERL5LIB && perl src/RecoverSampleAnnotationsAfterCombineVariantsByPosWalk.pl complex.vcf t/data/annot.vcf t/data/annot.call1.vcf t/data/annot.call2.vcf ) '.
	' t/data/recov.vcf &>/dev/stderr"') == 0 , 'RecoverSampleAnnotationsAfterCombineVariantsByPosWalk functional test');
#multiintersectseg
ok(system('bash -c "echo > /dev/stderr && set -e -o pipefail && mkdir -p tmp && '.
        'diff <(export PERL5LIB=\"blib/lib/\":$PERL5LIB && perl src//multiIntersectSeg.pl test/data/ref.dict  ./tmp/  test/data/multiinstersectseg/*.seg &>/dev/stderr &&  cat tmp/merged.tsv && rm -rv tmp/ &> /dev/stderr ) test/data/multiinstersectseg/merged.tsv  &>/dev/stderr"') == 0 , 'MultiInstersectSeg functional test');
#RenameChromosomes
ok(system('bash -c "echo > /dev/stderr && set -e -o pipefail && '.
	'diff <(export PERL5LIB=\"blib/lib/\":$PERL5LIB && perl src/RenameChromosomes.pl 0 test/data/renamechromosomes/renamechroms.tsv test/data/renamechromosomes/renamechroms.inputbed.bed) test/data/renamechromosomes/renamechroms.outputbed.bed &>/dev/stderr"') == 0 , 'RenameChromosomes functional test');
#made with perl src/RenameChromosomes.pl 0 test/data/renamechromosomes/renamechroms.tsv test/data/renamechromosomes/renamechroms.inputbed.bed > test/data/renamechromosomes/renamechroms.outputbed.bed

#ADfilter test
#made with `cat <(export PERL5LIB="blib/lib/":$PERL5LIB && perl src/AdFilter.pl -f 0.1 -c 4 t/data/target.vcf t/data/filter.vcf 2>/dev/null )> t/data/filtered.vcf` 
ok(system('bash -c "echo > /dev/stderr && set -e -o pipefail && '.
        'diff <(export PERL5LIB=\"blib/lib/\":$PERL5LIB && perl src/AdFilter.pl -f 0.1 -c 4 t/data/target.vcf t/data/filter.vcf ) '.
	't/data/filtered.vcf '.
        ' &>/dev/stderr"') == 0 , 'AdFilter functional test');

