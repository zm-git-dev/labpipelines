################################################################################
# DNA-seq, mutation calling etc.
################################################################################

function examplepipeline_dna {
  cat <<- EOF
=== pipeline 2016-02-05 ===

=> (+) wzseq_pindel => (+) wzseq_pindel_somatic

EOF
}

# samtools mpileup call mutation #
# check http://samtools.sourceforge.net/mpileup.shtml

# mpileup one sample
# usage: wzseq_mpileup1 bam/sname.bam
# one could use "vcfutils.pl varFilter -d10" to control vcf minimum depth
function wzseq_mpileup1 {

  base=$(pwd)
  bam=$1;
  sname=$(basename $bam .bam)
  [[ -d mpileup ]] || mkdir -p mpileup
  cmd="
cd $base
samtools mpileup -q 10 -uf $WZSEQ_REFERENCE $bam | bcftools call -cv -O b -o mpileup/$sname.raw.bcf
bcftools view mpileup/$sname.raw.bcf | vcf-sort -t mpileup/ | bgzip -c > mpileup/$sname.sorted.vcf.gz
"
  [[ ${!#} == "do" ]]  && qsub $pbsfn
}

# GATK best practice
function wzseq_seqtk_trimfq1 {
  [[ -d trimfq ]] || mkdir -p trimfq;
  [[ -d pbs ]] || mkdir -p pbs;
  $fq=$1;
  $fqname=${fq#fastq};
  $fqname=${fqname%.*};
  cmd="
seqtk trimfq $1 fastq/$fqname.fastq >trimfq/${fqname}_trimmed.fastq
"
  jobname="trimfq_$fqname"
  pbsfn="pbs/$jobname.pbs"
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 2 -ppn 1
  [[ ${!#} == "do" ]] && qsub $pbsfn
}

function __wzseq_bwa_mem_PE {

  cmd='
cd '$base'
mkdir -p bam;
bwa mem -M -R \"@RG\tLB:$WZSEQ_REFVERSION\tID:${sname}\tPL:Illumina\tPU:nextseq500\tSM:${sname}\" -t 28 $WZSEQ_BWA_INDEX fastq/$sread1 fastq/$sread2 | samtools sort -O bam -T bam/$sname.tmp -o bam/$sname.bam -
samtools index bam/$sname.bam
samtools flagstat bam/${sname}.bam > bam/$sname.bam.flagstat
'
  jobname="bwamem_$sname"
}

function __wzseq_bwa_mem_SE {
  cmd='
cd '$base'
mkdir -p bam;
bwa mem -M -R "@RG\tLB:'$WZSEQ_REFVERSION'\tID:'$sname'\tPL:Illumina\tPU:nextseq500\tSM:'$sname'" -t '$ppn' '$WZSEQ_BWA_INDEX' '$fastq' | samtools sort -O bam -T bam/'$sname'.tmp -o bam/'$sname'.bam -
samtools index bam/'$sname'.bam
samtools flagstat bam/'$sname'.bam > bam/'$sname'.bam.flagstat
'
  jobname="bwamem_"$sname
}

# wzseq_bwa_mem
function wzseq_bwa_mem {

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      __wzseq_bwa_mem
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# indel realignment
function wzseq_GATK_indelrealign {

  [[ -d pbs ]] || mkdir -p pbs
  [[ -d indelrealn ]] || mkdir -p indelrealn
  while read sname sread1 sread2; do
    nthreads=28;
    interval=logs/${sname}_gatk_indels.intervals
    cmd="
# create interval
logs/$(sample_id)_gatk_indels.intervals:bams/$(sample_id)_rg_dedup.bam.bai bams/$(sample_id)_rg_dedup.bam
java -Xmx60g -Djava.io.tmpdir=./indelrealn/ -jar ~/software/GATK/GATK-3.3.0/GenomeAnalysisTK.jar --fix_misencoded_quality_scores --num_threads $nthreads -I bam/$bam -R $WZSEQ_REFERENCE -T RealignerTargetCreator --filter_mismatching_base_and_quals -o $interval --known $WZSEQ_GATK_KNOWN_INDEL

# actual realignment
java -Xmx2g -Djava.io.tmpdir=./indelrealn/ -jar ~/software/GATK/GATK-3.3.0/GenomeAnalysisTK.jar --fix_misencoded_quality_scores -o indelrealn/${sname}.indelrealn.bam -I bam/$bam -R $WZSEQ_REFERENCE -T IndelRealigner --filter_mismatching_base_and_quals -rf BadCigar -targetIntervals $interval --maxReadsForRealignment 200000 -known ${WZSEQ_GATK_KNOWN_INDEL}
"
    jobname="indelrealn_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
  # -known $(broad_files)/Mills_and_1000G_gold_standard.indels.hg19.vcf -known $(broad_files)/1000G_phase1.indels.hg19.vcf
}

function wzseq_pindel {

  # call pindel_somatic directly if you want somatic indels
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d pindel ]] || mkdir pindel

  awk '/^\[/{p=0}/\[pindel\]/{p=1;next} p&&!/^$/' samples |
    while read sname bam; do
      cmd="
cd $base
[[ -d pindel/$sname ]] || mkdir pindel/$sname
echo -e $(readlink -f $bam)\"\\t250\\t\"$sname >pindel/$sname/pindel.config
~/software/pindel/default/pindel -i pindel/$sname/pindel.config -f $WZSEQ_REFERENCE -o pindel/$sname/$sname -c ALL -T 28 2>&1 > pindel/$sname/$sname.pindel_log
"
      jobname="pindel_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done

}

function wzseq_pindel_somatic {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d pindel ]] || mkdir pindel
  # sname normal tumor
  awk '/^\[/{p=0}/\[pindel_somatic\]/{p=1;next} p&&!/^$/' samples |
    while read sname sname1 sname2; do
      cmd="
cd $base
[[ -d pindel/$sname ]] || mkdir pindel/$sname
echo -e bam/$sname1.bam\"\\t250\\t\"$sname1 >pindel/$sname/pindel.config
echo -e bam/$sname2.bam\"\\t250\\t\"$sname2 >>pindel/$sname/pindel.config
~/software/pindel/default/pindel -i pindel/$sname/pindel.config -f $WZSEQ_REFERENCE -o pindel/$sname/$sname -c ALL -T 28 2>&1 > pindel/$sname/$sname.pindel_log

# filter somatic events
grep \"ChrID\" pindel/$sname/${sname}_D >pindel/$sname/head_DSI
grep \"ChrID\" pindel/$sname/${sname}_SI >>pindel/$sname/head_DSI
cat <<EOT >pindel/$sname/somatic.config
indel.filter.input = pindel/$sname/head_DSI
indel.filter.vaf = 0.1
indel.filter.cov = 5
indel.filter.hom = 2
indel.filter.pindel2vcf = /primary/vari/software/pindel/default/pindel2vcf
indel.filter.reference = $WZSEQ_REFERENCE
indel.filter.referencename = NA
indel.filter.referencedate = NA
indel.filter.output = pindel/$sname/$sname.indel.filter
EOT
perl ~/software/pindel/default/somatic_filter/somatic_indelfilter.pl pindel/$sname/somatic.config
"
      jobname="pindel_somatic_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# compute base coverage for both tumor and normal
function wzseq_basecov_somatic {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d basecov ]] || mkdir basecov
  awk '/^\[/{p=0}/\[pindel_somatic\]/{p=1;next} p&&!/^$/' samples |
    while read sname sname1 sname2; do
      cmd="
cd $base
samtools mpileup -s bam/$sname1.bam bam/$sname2.bam | perl -alne 'BEGIN{\$b=0}{\$c=0; foreach(split //, \$F[6]){if (ord(\$_)>53) {\$c++;}} \$d=0; foreach(split //, \$F[10]){if(ord(\$_)>53){\$d++;}} if(\$c>4 && \$d>4){\$b++;}}END{print \$b}' > basecov/${sname1}_vs_${sname2}.cov;
"
      jobname="basecov_somatic_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function wzseq_lofreq {

  base=$(pwd)
  [[ -d lofreq ]] || mkdir lofreq
  [[ -d pbs ]] || mkdir pbs

  for bam in bam/*.bam; do
    sname=$(basename $bam .bam)
    cmd="
cd $base
~/software/lofreq/default/bin/lofreq call -f $WZSEQ_REFERENCE -s -S $WZSEQ_DBSNP -o lofreq/$sname.vcf $bam -b 1
transvar ganno --vcf lofreq/$sname.vcf --ccds >lofreq/$sname.vcf.transvar
"
    jobname="lofreq_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 60 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wzseq_vardict {
  base=$(pwd)
  [[ -d vardict ]] || mkdir vardict
  [[ -d pbs ]] || mkdir pbs
  bedregion=""
  for bam in bam/*.bam; do
    sname=$(basename $bam .bam)
    cmd="
cd $base
export PATH=/primary/vari/software/vardict/VarDict:$PATH
## AF_THR=0.001 # minimum allele frequency
# vardict -G $WZSEQ_REFERENCE -f 0.001 -N ${sname} -b $bam -c 1 -S 2 -E 3 -g 4 $WZSEQ_EXOME_CAPTURE > vardict/${sname}_out.raw
# cat vardict/${sname}_out.raw | teststrandbias.R >vardict/${sname}_out.strandbias
# cat vardict/${sname}_out.strandbias | var2vcf_valid.pl -N sample_name -E -f 0.001 >vardict/${sname}_out.vcf
transvar ganno --vcf vardict/${sname}_out.vcf --ccds >vardict/${sname}_out.vcf.transvar
"
    jobname="vardict_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

