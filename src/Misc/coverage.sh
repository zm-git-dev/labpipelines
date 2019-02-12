####################################
## base/cpg coverage and uniformity
####################################

function __wzseq_uniformity_1M() {
  cmd='
  bedtools makewindows -w 1000000 -g '${WZSEQ_REFERENCE}'.fai | grep -v random | grep -v chrUn | grep -v hap | sortbed | bedtools coverage -a - -b <(samtools view -O BAM -q 40 '$input_bam') -sorted >uniformity/'${sname}'_1Mb.bed
'
  jobname='uniformity_1m_'$sname
}

# create coverage track, unique and nonunique mapping
function wzseq_bam_coverage {

  base=$(pwd);
  [[ -d tracks ]] || mkdir tracks
  [[ -d pbs ]] || mkdir pbs
  [[ -d qc ]] || mkdir qc
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    __wzseq_bam_coverage
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wgbs_cpgcoverage_OBSOLETE {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d cpg ]] || mkdir cpg
  for f in pileup/*.vcf.gz; do
    bfn=$(basename $f .vcf.gz)
    cmd="
cd $base

## coverage
~/tools/biscuit/master/biscuit/biscuit vcf2bed -t cg -k 0 -c $f > cpg/${bfn}_cg.bedg
bedtools intersect -a $WZSEQ_CPGBED -b cpg/${bfn}_cg.bedg -sorted -loj | awk '{if(\$11==\".\")\$11=0;print \$1,\$2,\$3,\$11}' | bedtools groupby -g 1-3 -c 4 -o sum >cpg/cpgCoverage_${bfn}.bedg
wzplot cumhist -t cpg/cpgCoverage_${bfn}.bedg -c 4 --xlabel \"Coverage\" --ylabel \"Cumulative Distribution\" -o cpg/cpgCoverage_wholeGenome_${bfn}.png
bedtools intersect -a cpg/cpgCoverage_${bfn}.bedg -b $WZSEQ_CGIBED | wzplot cumhist -t - -c 4 --xlabel \"Coverage\" --ylabel \"Cumulative Distribution\" -o cpg/cpgCoverage_CGI_${bfn}.png

## beta value distribution
wzplot hist --maxline 1000000000 -t cpg/${bfn}_cg.bedg -c 4 --xlabel \"Beta Values\" -o cpg/betaDist_wholeGenome_${bfn}.png
bedtools intersect -a cpg/${bfn}_cg.bedg -b $WZSEQ_CGIBED | wzplot hist --maxline 1000000000 -t - -c 4 --xlabel \"Beta Values\" -o cpg/betaDist_CGI_${bfn}.png

## Methylation average in CGI
awk '\$5>5' cpg/${bfn}_cg.bedg | bedtools intersect -a $WZSEQ_CGIBED -b - -wo -sorted | bedtools groupby -i - -g 1-5 -c 14,14 -o mean,count > cpg/CGImethAverage_cov5_${bfn}.bed
wzplot hist --maxline 1000000000 -t cpg/CGImethAverage_cov5_${bfn}.bed -c 6 --xlabel \"CGImethAverage\" -o cpg/CGImethAverage_cov5_${bfn}.png

# bedtools closest -a $WZSEQ_TSSBED -b cpg/$bfn.cgi.bed -d >cpg/CGImethAverage_cov5_${bfn}_tss.bed
# rm -f cpg/$bfn.cg.bedg
"
    jobname="cpgCoverage_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}


