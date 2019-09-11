################################################################################
# Whole Genome Bisulfite Sequencing
################################################################################

function examplepipeline_wgbs {
  cat <<- EOF
=== pipeline 2015-10-01 ===
[o] wgbs_adaptor => [o] wzseq_fastqc

 => (+) wgbs_biscuit_align / (+) wgbs_bwameth / (+) wgbs_bismark_bowtie1 / (+) wgbs_bismark_bowtie2 / wgbs_bsmap

 => (+) wgbs_biscuit_align_lambdaphage => (+) wgbs_biscuit_pileup_lambdaphage (TODO: exclude human reads)

 => (+) wzseq_GATK_realign (TODO: wgbs_indel_realign) => (+) wzseq_picard_markdup (TODO: wgbs_biscuit_markdup) => (+) wzseq_clean_intermediate => (+) TODO: wgbs_basequal_recal

 => [o] wzseq_merge_bam => (+) wzseq_qualimap => (defunct) wzseq_picard_WGSmetrics => (+) wzseq_bam_coverage

 => (+) wgbs_methpipe => (+) wgbs_methpipe_methylome => (+) wgbs_methpipe_allele

 => (+) wgbs_biscuit_pileup [-nome] => (+) wgbs_vcf2tracks [-nome] => (+) wgbs_cpgcoverage => (+) wgbs_repeat => (+) wgbs_repeat_diff

 => (+) wgbs_methylKit_summary

 => (+) wgbs_diffmeth_simple => (+) wgbs_methylKit_diffmeth => (+) wgbs_methpipe_diff => (+) wgbs_metilene
EOF
}

# check adaptor conversion rate
function wgbs_adaptor() {
  local base=$(pwd);
  [[ -d adaptor ]] || mkdir adaptor;
  [[ -d pbs ]] || mkdir pbs
  cmd="
cd $base
parallel -j 4 ~/wzlib/pyutils/wzadapter.py {} '>' adaptor/{/.}.txt ::: $base/fastq/*R1*.fastq.gz
for f in adaptor/*.txt; do
  tail -2 \$f | cut -d\":\" -f2 | cut -d\" \" -f2 | awk -v fn=\$f '{a[NR]=\$1}END{print fn,a[1],a[2]}'
done | awk 'BEGIN{print \"filename\tadaptorC\tadaptorT\tadaptorC2T\"}{print \$1,\$2,\$3,\$3/\$2*100\"%\"}' > adaptor/adaptor.stats
"
  jobname="adaptor_analysis_"$(basename $base)
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 8 -ppn 4
  [[ ${!#} == "do" ]] && qsub $pbsfn
}

function __wzseq_clean_bam {

  # input: bamfn
  cmd="
set -xe
cd $base
rm -f bam/$bamfn
"
  jobname="clean_bam_$bamfn"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 1 -memG 1 -ppn 1
}

function __wzseq_clean_fastq {

  # input: sname, sread1, sread2
  cmd="
set -xe
cd $base
[[ -f fastq/$sread1 ]] && rm -f fastq/$sread1
[[ -f fastq/$sread2 ]] && rm -f fastq/$sread2
"
  jobname="clean_fastq_$sname"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 1 -memG 1 -ppn 1
}


########################
## section1: alignment
########################

# biscuit alignment
# wgbs_biscuit_index_reference
# run mm10generic to setup $WZ_BISCUIT_INDEX
# use ~/pbs
function wgbs_biscuit_index_reference {
  base=$(pwd);
  [[ -d ~/pbs ]] || mkdir ~/pbs
  cmd="
biscuit index $WZSEQ_BISCUIT_INDEX
"
  jobname="biscuit_index"
  pbsfn=~/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 20 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function __wgbs_biscuit_align_SE {
  cmd='
cd '$base'
mkdir -p bam
~/bin/biscuit align -b 3 '$WZSEQ_BISCUIT_INDEX' -t '$ppn' '$fastq' | samtools sort -T '$output_bam'_tmp -O bam -o '$output_bam'
'
  jobname="biscuit_align_"$sname"_SE"
}

function __wgbs_biscuit_align_SE_both {
  cmd='
cd '$base'
mkdir -p bam
~/bin/biscuit align '$WZSEQ_BISCUIT_INDEX' -t '$ppn' '$fastq' | samtools sort -T '$output_bam'_tmp -O bam -o '$output_bam'
'
  jobname="biscuit_align_"$sname"_SE"
}

function __wgbs_biscuit_align_PE {
  cmd='
cd '$base'
mkdir -p bam
~/bin/biscuit align '$WZSEQ_BISCUIT_INDEX' -b 1 -t '$ppn' '$fastq1' '$fastq2' | samtools sort -T '$output_bam'_tmp -O bam -o '$output_bam'
'
  jobname="biscuit_align_"$sname"_PE"
}

function __wgbs_biscuit_align_PE_both {
  cmd='
cd '$base'
mkdir -p bam
# ~/tools/biscuit/master/biscuit/biscuit
~/bin/biscuit align '$WZSEQ_BISCUIT_INDEX' -t '$ppn' '$fastq1' '$fastq2' >'$output_bam'.sam
samtools sort -T '$output_bam'_tmp -O bam -o '$output_bam' '$output_bam'.sam
#samtools index '$output_bam';
#samtools flagstat '$output_bam' > '$output_bam'.flagstat
'
  jobname="biscuit_align_"${sname}"_PE_both"
}

function __wgbs_biscuit_align_PE_Walid_lib {
  cmd='
cd '$base'
mkdir -p bam
~/bin/biscuit align '$WZSEQ_BISCUIT_INDEX' -t '$ppn' -J AGATCGGAAGAGC -K AGATCGGAAGAGC '$fastq1' '$fastq2' | samtools sort -T '$output_bam'_tmp -O bam | ~/bin/biscuit bsconv -m 10 '$WZSEQ_REFERENCE' - '$output_bam'
samtools index '$output_bam';
samtools flagstat '$output_bam' > '$output_bam'.flagstat
'
  jobname="biscuit_align_"${sname}"_PE_both"
}

function __wgbs_biscuit_align_PE_POETIC {
  cmd='
cd '$base'
mkdir -p bam
~/bin/biscuit align '$WZSEQ_BISCUIT_INDEX' -t '$ppn' -J AGATCGGAAGAGC -K AAATCAAAAAAAC '$fastq1' '$fastq2' >'$output_bam'.sam
samtools sort -T '$output_bam'_tmp -O bam -o '$output_bam' '$output_bam'.sam
# samtools index '$output_bam';
# samtools flagstat '$output_bam' > '$output_bam'.flagstat
'
  jobname="biscuit_align_"${sname}"_PE_POETIC"
}

function __wgbs_biscuit_markdup {
  cmd='
cd '$base'
~/tools/biscuit/development/biscuit/biscuit markdup '$input_bam' '$output_bam'_unsrt.bam 2>'$output_bam'_markdup_report.txt
samtools sort -T '$output_bam'_unsrt.tmp -o '$output_bam' '$output_bam'_unsrt.bam
rm -f '$output_bam'_unsrt.bam
mkdir -p multiqc/raw/biscuit/
ln -sf `readlink -f '$output_bam'_markdup_report.txt` multiqc/raw/biscuit/
'
  jobname='biscuit_markdup_'$sname
}

function __wgbs_biscuit_QC {
  cmd='
cd '$base'
~/tools/biscuit/development/biscuit/scripts/QC.sh -v '$input_vcf' '$WZSEQ_BISCUIT_QC_SETUP' '$sname' '$input_bam'
mkdir -p multiqc/raw/BISCUITqc/'$sname'
ln -sf `readlink -f BISCUITqc` multiqc/raw/BISCUITqc
'
  jobname='biscuit_QC_'$sname
}

function wgbs_biscuit_align_lambdaphage {
  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
biscuit align $WZSEQ_BISCUIT_INDEX_LAMBDAPHAGE -t 28 $base/fastq/$sread1 $base/fastq/$sread2 | samtools view -h -F 0x4 - | samtools sort -T $base/bam/${sname} -O bam -o $base/bam/${sname}_lambdaphage.bam
samtools index $base/bam/${sname}_lambdaphage.bam
samtools flagstat $base/bam/${sname}_lambdaphage.bam > $base/bam/${sname}_lambdaphage.bam.flagstat
"
      jobname="biscuit_align_lambdaphage_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# BWA-meth
function wgbs_bwameth_index_reference {
  [[ -d ~/pbs ]] || mkdir ~/pbs
  cmd="
bwameth.py index $WZSEQ_BWAMETH_INDEX
"
  jobname="bwameth_index"
  pbsfn=~/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 20 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function wgbs_bwameth() {

  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmds="
cd $base
if [[ \"$sread2\" == \".\" ]]; then
  bwameth --reference $WZSEQ_BWAMETH_INDEX fastq/$sread1 -t 28 --prefix bam/${sname}_bwameth
else
  bwameth --reference $WZSEQ_BWAMETH_INDEX fastq/$sread1 fastq/$sread2 -t 28 --prefix bam/${sname}_bwameth
fi
samtools flagstat bam/${sname}_bwameth.bam > bam/${sname}_bwameth.bam.flagstat
"
      jobname="bwameth_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmds" -name $jobname -dest $pbsfn -hour 8 -ppn 28 -memG 250
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# bismark with bowtie1
function wgbs_bismark_bowtie1_prepare_reference {
  base=$(pwd)
  cmd="
export PATH=~/tools/bismark/default:~/tools/bowtie1/default:$PATH
bismark_genome_preparation --bowtie1 $WZSEQ_BISMARK_BT1_INDEX
"
  jobname="bismark_bt1_prepare_reference"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function wgbs_bismark_bowtie1 {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
export PATH=~/tools/bismark/default:~/tools/bowtie1/default:$PATH
cd $base
bismark $WZSEQ_BISMARK_BT1_INDEX --chunkmbs 2000 -1 $base/fastq/$sread1 -2 $base/fastq/$sread2
"
      jobname="bismark_bt1_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 100 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

## TRIM GALORE!
##################
function __wzseq_trim_galore_SE {
  cmd='
set -xe
cd '$base'
[[ -d '$trim_galore_dir' ]] && rm -rf '$trim_galore_dir'
mkdir -p '$trim_galore_dir'
~/software/trim_galore/default/trim_galore --fastqc '$fastq' --gzip -o '$trim_galore_dir'
mkdir -p multiqc/raw/trim_galore/
ln -fs `readlink -f '$trim_galore_dir'` multiqc/raw/trim_galore/
'
  jobname="trim_galore_SE_"$sname
}

function __wzseq_trim_galore_PE {
  cmd='
set -xe
cd '$base'
[[ -d '$trim_galore_dir' ]] && rm -rf '$trim_galore_dir'
mkdir -p '$trim_galore_dir'
~/software/trim_galore/default/trim_galore --fastqc --paired '$fastq1' '$fastq2' -o '$trim_galore_dir'
mkdir -p multiqc/raw/trim_galore/
ln -s `readlink -f '$trim_galore_dir'` multiqc/raw/trim_galore/
'
  jobname="trim_galore_PE_"$sname
}

function __wzseq_trim_galore_PE2 {
  cmd='
set -xe
cd '$base'
[[ -d '$trim_galore_dir' ]] && rm -rf '$trim_galore_dir'
mkdir -p '$trim_galore_dir'
~/software/trim_galore/default/trim_galore --quality 28 --phred33 --fastqc --clip_R1 9 --clip_R2 9 --three_prime_clip_R1 9 --three_prime_clip_R2 9 --paired '$fastq1' '$fastq2' -o '$trim_galore_dir'
mkdir -p multiqc/raw/trim_galore/
ln -s `readlink -f '$trim_galore_dir'` multiqc/raw/trim_galore/
'
  jobname="trim_galore_PE_"$sname
}

function __wzseq_customize {
  cmd="
set -xe
cd $base
$customized_command
"
  jobname="wzseq_customize_$sname"
}




## bismark with bowtie2
###########################
function wgbs_bismark_bowtie2_prepare_reference {
  base=$(pwd)
  cmd="
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
bismark_genome_preparation --bowtie2 --verbose $WZSEQ_BISMARK_BT2_INDEX
"
  jobname="bismark_bt2_prepare_reference"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 10 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function __wgbs_bismark_bowtie2_SE {
  : '
fastq=fastq/${sname}_trim_galore/${sname}_merged.fq.gz
direction="--non_directional"
bismark_bt2_dir=bam/${sname}_bismark_bt2
bismark_bt2_bam_unsorted=bam/${sname}_bismark_bt2/${sname}_merged.fq.gz_bismark_bt2.bam
bismark_bt2_bam_final=bam/${sname}_bismark_bt2.bam
hour=200; memG=180; ppn=28
'
  cmd='
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
set -xe
cd '$base'
mkdir -p bam
rm -rf '$bismark_bt2_dir'
bismark '$direction' '$WZSEQ_BISMARK_BT2_INDEX' --bowtie2 --chunkmbs 2000 -p '$ppn' -o '$bismark_bt2_dir' '$fastq' --temp_dir '$bismark_bt2_dir'/tmp
samtools sort -O bam -o '$bismark_bt2_bam_final' -T '${bismark_bt2_bam_unsorted}'_tmp '${bismark_bt2_bam_unsorted}'
samtools index '${bismark_bt2_bam_final}'
samtools flagstat '${bismark_bt2_bam_final}' >'${bismark_bt2_bam_final}'.flagstat
mkdir -p multiqc/raw/bismark/
ln -sf `readlink -f '$bismark_bt2_dir'/[PS]E_report.txt` multiqc/raw/bismark/; done
'
  jobname="bismark_bt2_SE_"$sname
}

function __wgbs_bismark_bowtie2_PE {
  # fastq1=fastq/${sname}_trim_galore/${sname}_R1_val_1.fq.gz
  # fastq2=fastq/${sname}_trim_galore/${sname}_R2_val_2.fq.gz
  # direction="--non_directional"
  # pbat=""
  # bismark_bt2_dir=bam/${sname}_bismark_bt2
  # bismark_bt2_bam_unsorted=bam/${sname}_bismark_bt2/${sname}_pe.bam
  # bismark_bt2_bam_final=bam/${sname}_bismark_bt2.bam
  # hour=200; memG=180; ppn=28
  cmd='
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
set -xe
cd '$base'
mkdir -p bam
rm -rf '$bismark_bt2_dir'
bismark '$direction' '$pbat' '$WZSEQ_BISMARK_BT2_INDEX' --bowtie2 --chunkmbs 2000 -p '$ppn' -o '$bismark_bt2_dir' -B '$sname' -1 '$fastq1' -2 '$fastq2' -temp_dir '$bismark_bt2_dir'/tmp
samtools sort -O bam -o '$bismark_bt2_bam_final' -T '$bismark_bt2_bam_unsorted'_tmp '$bismark_bt2_bam_unsorted'
samtools index '$bismark_bt2_bam_final'
samtools flagstat '$bismark_bt2_bam_final' >'$bismark_bt2_bam_final'.flagstat
mkdir -p multiqc/raw/bismark
ln -s `readlink -f '$bismark_bt2_dir'` multiqc/raw/bismark/
'
  jobname="bismark_bt2_PE_$sname"
}

function __wgbs_bismark_deduplicate {
  # Please note that this function is very memory-consuming
  # input_bam=bam/${sname}_bismark_bt2.bam
  # hour=48; memG=50; ppn=2
  # library="-p"
  ## paired-end
  # library="-p"
  ## single-end
  # library="-s"
  cmd='
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
cd '$base'
mkdir -p multiqc/raw/bismark/
deduplicate_bismark '$library' --bam '$input_bam'
ln -sf `readlink -f '${input_bam%.bam}'.deduplication_report.txt` multiqc/raw/bismark/
'
  jobname="bismark_deduplicate_"$sname
}

function __wgbs_bismark_methylextraction {
  : '
# for paired-end
library="--no_overlap"
# for single-end
library=""

input_bam=bam/${sname}_bismark_bt2.deduplicated.bam
hour=48; memG=10; ppn=1
'
  cmd='
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
cd '$base'
mkdir -p bismark_methylextraction
bismark_methylation_extractor '$library' --gzip --bedGraph '$input_bam' -o bismark_methylextraction
zcat bismark_methylextraction/'$(basename $input_bam .bam)'.bismark.cov.gz | awk '\''{print $1,$2-1,$3,$4/100,$5+$6}'\'' | sortbed | biscuit mergecg '$WZSEQ_REFERENCE' - >bismark_methylextraction/'$sname'.cpg_methylation.bed
ln -sf `readlink -f bismark_methylextraction/'$(basename $input_bam)'_splitting_report.txt` multiqc/raw/bismark/
'
  jobname="bismark_methylation_extraction_"$sname
}

## bsmap
###########
function wgbs_bsmap {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam

  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
cd $base
if [[ \"$sread2\" == \".\" ]]; then
  ~/tools/bsmap/default/bsmap -a fastq/$sread1 -d $WZSEQ_BSMAP_INDEX -o bam/${sname}_bsmap_unsorted.bam -p 28 -s 16 -v 10 -q 2 -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
else
  ~/tools/bsmap/default/bsmap -a fastq/$sread1 -b fastq/$sread2 -d $WZSEQ_BSMAP_INDEX -o bam/${sname}_bsmap_unsorted.bam -p 28 -s 16 -v 10 -q 2 -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
fi

samtools sort -O bam -o bam/${sname}_bsmap.bam -T bam/${sname}_bsmap.tmp bam/${sname}_bsmap_unsorted.bam
rm -f bam/${sname}_bsmap_unsorted.bam
samtools index bam/${sname}_bsmap.bam
samtools flagstat bam/${sname}_bsmap.bam > bam/${sname}_bsmap.bam.flagstat
"
      jobname="bsmap_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 8 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

#### rmapbs (doesn't allow gaps) ###
# rmapbs -c hg18 -o Human_NHFF.mr Human_NHFF.fastq
# rmapbs-pe -c hg18 -o Human_ESC.mr Human_ESC_1.fastq Human_ESC_2.fastq

# methylKit summary
function wgbs_methylKit_summary {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d methylKit ]] || mkdir methylKit
  for f in pileup/*.vcf.gz; do
    bfn=$(basename $f .vcf.gz)
    cmd="
cd $base
mkdir methylKit/$bfn
biscuit vcf2bed -t cg -u -c $f | pybiscuit.py to_methylKit >methylKit/$bfn/$bfn.methylKit
~/wzlib/Rutils/bin/bioinfo/methylKit.r summary methylKit/$bfn/$bfn.methylKit -o methylKit/$bfn/
"
    jobname="methylKit_summary_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 50 -ppn 3
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

## TODO bissnp ####
function wgbs_merge_methlevelaverages {
  # usage: biscuit_merge_methlevelaverages data/methlevelaverage

  local basedir=$1;
  for f in $basedir/*.txt; do
    sed -e 's/://g' -e 's/%//' $f | awk -f wanding.awk -e 'BEGIN{split("",k);split("",n);split("",v);}(length($0)>0){k[length(k)+1]=$1;n[length(n)+1]=$2;v[length(v)+1]=$3/100;}END{for(i=1;i<=length(k);++i){kn[i]=k[i]"n"};print(join(k,1,length(k),"\t")"\t"join(kn,1,length(kn),"\t"));print(join(v,1,length(v),"\t")"\t"join(n,1,length(n),"\t"))}' > ${f%.txt}.processed
  done

  wzmanip concat -f $basedir/*.processed | awk -f wanding.awk -e '{n=NF;print;}END{repeat("NA", n, rep); print joina(rep,"\t");}' > $basedir/merge;
}

#############
## MethPIPE
#############

# note that I use my own version of to-mr which is agnostic of input bam type
# /home/wanding.zhou/software/methpipe/default/to-mr
# ~/software/methpipe/default/to-mr -o methpipe/$bfn -m bsmap 
# ~/software/methpipe/default/to-mr -o methpipe bam/Undetermined_bismark_bt2/Undetermined_L000_R1_001_val_1.fq.gz_bismark_bt2_pe.bam -m bismark
# ~/software/methpipe/default/to-mr -o methpi2 bam/Undetermined_bsmap.bam -m bsmap
# duplicate-remover -S 
function wgbs_methpipe {

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d methpipe ]] || mkdir methpipe
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)
    obfn=methpipe/$bfn
    cmd="
cd $base
pybiscuit.py to_mr -i $fn -o $obfn.mr
LC_ALL=C sort -T methpipe/ -k1,1 -k2,2n -k3,3n -k6,6 -o $obfn.mr.sorted_start $obfn.mr

echo \"[\$(date)] Remove duplicates...\"
~/tools/methpipe/default/bin/duplicate-remover -S ${obfn}_dremove_stat.txt -o $obfn.mr.dremove.untrim $obfn.mr.sorted_start

echo \"[\$(date)] Remove random contigs and haplotype...\"
cut -f1 $obfn.mr.dremove.untrim | uniq | while read chrm; do [[ -s $WZSEQ_REFERENCE_SPLIT/\$chrm.fa ]] || echo \"^\"\$chrm; done >$obfn.mr.uchr
grep -v -f $obfn.mr.uchr $obfn.mr.dremove.untrim > $obfn.mr.dremove

echo \"[\$(date)] Estimate bisulfite conversion rate...\"
~/tools/methpipe/default/bin/bsrate -c $WZSEQ_REFERENCE_SPLIT -o $obfn.bsrate $obfn.mr.dremove
LC_ALL=C sort -T methpipe/ -k1,1 -k3,3n -k2,2n -k6,6 -o $obfn.mr.sorted_end_first $obfn.mr.dremove

echo \"[\$(date)] Compute single-site methylation levels...\"
~/tools/methpipe/default/bin/methcounts -c $WZSEQ_REFERENCE_SPLIT -o $obfn.meth $obfn.mr.sorted_end_first

echo \"[\$(date)] Merge C,G in CpG context...\"
~/tools/methpipe/default/bin/symmetric-cpgs -m -o $obfn.CpG.meth $obfn.meth

echo \"[\$(date)] Get statistics of methylation levels...\"
~/tools/methpipe/default/bin/levels -o $obfn.levels $obfn.meth

echo \"[\$(date)] Clean up intermediate files...\"
rm -f $obfn.mr $obfn.mr.sorted_start $obfn.mr.dremove.untrim $obfn.mr.uchr 
# not sure if these 2 should be removed as intermediate: $obfn.mr.dremove $obfn.mr.sorted_end_first
echo \"[\$(date)] Done.\"
"
    jobname="methpipe_basic_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

## methpipe
function wgbs_methpipe_methylome {
  base=$(pwd)
  [[ -d methpipe ]] || mkdir methpipe
  for f in methpipe/*.mr; do
    bfn=$(basename $f .mr)
    obfn=methpipe/$bfn
    cmd="
cd $base

echo \"[\$(date)] Hypomethylated (primarily) and hypermethylated region...\"
~/tools/methpipe/default/bin/hmr -p $obfn.hmr.params -o $obfn.hmr $obfn.CpG.meth

# use trained parameter
# echo \"[\$(date)] Train parameters...\"
# ~/tools/methpipe/default/bin/hmr -p $obfn.hmr.params -o $obfn.hmr $obfn.CpG.meth
# echo \"[\$(date)] Use trained parameter to call HMR...\"
# ~/tools/methpipe/default/bin/hmr -P $obfn.hmr.params -o $obfn.hmr $obfn.CpG.meth

echo \"[\$(date)] Partially methylated regions (PMR)...\"
# not sure what this is
~/tools/methpipe/default/bin/hmr -partial -o $obfn.hmr $obfn.CpG.meth

# echo \"[\$(date)] Hypermethylated region in plant (HyperMR)...\"
# ~/tools/methpipe/default/bin/hypermr -o $obfn.hypermr $obfn.CpG.meth

echo \"[\$(date)] Partially methylated domain (PMD, 1kb bin level segmentation).\"
~/tools/methpipe/default/bin/pmd -o $obfn.pmd $obfn.CpG.meth

echo \"[\$(date)] Done.\"
"
    jobname="methpipe_methylome_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function wgbs_methpipe_allele {
  base=$(pwd)
  [[ -d methpipe ]] || mkdir methpipe
  for f in methpipe/*.mr.dremove; do
    obfn=${f%.mr.dremove}
    cmd="
cd $base

echo \"[\$(date)] Make read distribution, (.epiread)\"
~/tools/methpipe/default/bin/methstates -c $WZSEQ_REFERENCE_SPLIT -o ${obfn}.epiread ${obfn}.mr.dremove

echo \"[\$(date)] Single site ASM scoring (.allelic)\"
~/tools/methpipe/default/bin/allelicmeth -c $WZSEQ_REFERENCE_SPLIT -o ${obfn}.allelic ${obfn}.epiread

echo \"[\$(date)] Allelically methylated region (AMRs)\"
~/tools/methpipe/default/bin/amrfinder -o $obfn.amr -c $WZSEQ_REFERENCE_SPLIT $obfn.epiread.tmp
sort -k1,1 -k2,2n -T methpipe/ $obfn.epiread.tmp > $obfn.epiread

# echo \"[\$(date)] Test allele-specific methylation at promoter of imprinted region\"
# amrtester -o ${obfn}.amr -c hg19 target_interval.bed ${obfn}.epiread

echo \"[\$(date)] Calculate meth-entropy\"
# -v verbose
~/tools/methpipe/default/bin/methentropy -w 5 -v -o $obfn.entropy $WZSEQ_REFERENCE_SPLIT $obfn.epiread

echo \"[\$(date)] Done.\"
"
    jobname="methpipe_allele_"$(basename $obfn)
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function wgbs_methpipe_diff {
  # assuming the corresponding methpipe run (mr file) has finished
  # for each listed vcf file in the corresponding methpipe folder
  base=$(pwd)
  awk '/^\[/{p=0}/\[diffmeth\]/{p=1;next} p&&!/^$/' samples |
    while read sname vcfpath1 vcfpath2; do
      echo 1
    done
}

####################
## BISCUIT
####################
function __wgbs_biscuit_pileup {
  : '
input_bam=bam/${sname}_markdup.bam
output_vcf=pileup/${sname}.vcf.gz
hour=24; memG=100; ppn=10
pipeline_depend bamjob
'
  cmd='
cd '$base'
mkdir -p pileup
~/tools/biscuit/development/biscuit/biscuit pileup -q '$ppn' '$WZSEQ_REFERENCE' '$input_bam' >pileup/'$sname.vcf'
bgzip -f pileup/'$sname.vcf'
tabix -p vcf pileup/'$sname.vcf.gz'

# ~/tools/biscuit/development/biscuit/biscuit vcf2bed -k 1 -t cg pileup/'$sname.vcf.gz' | cut -f1-5 | LC_ALL=C sort -k1,1 -k2,2n -T pileup | ~/tools/biscuit/development/biscuit/biscuit mergecg '$WZSEQ_REFERENCE' - | cut -f1-5 | gzip -c >pileup/'$sname'_cpg.b.gz

~/bin/biscuit vcf2bed -k 1 -t cg pileup/'$sname'.vcf.gz | cut -f1-5 | LC_ALL=C sort -k1,1 -k2,2n -T pileup | ~/bin/biscuit mergecg '$WZSEQ_REFERENCE' - | bedtools intersect -a '$WZSEQ_CPGBED' -b - -sorted -loj | cut -f1-3,7,8 | bgzip -fc >pileup/'$sname'_cpg.b.gz
tabix -p bed pileup/'$sname'_cpg.b.gz

# make a bw file, arbitrarily chose cov5
# ~/tools/biscuit/development/biscuit/biscuit vcf2bed -k 5 -t cg pileup/'$sname.vcf.gz' | LC_ALL=C sort -k1,1 -k2,2n -T pileup >pileup/'$sname'_cg_cov5.bed
zcat pileup/'$sname'_cpg.b.gz | awk '\''$5!="." && $5>=5'\'' | gzip -c >pileup/'$sname'_cpg_cov5.bed.gz
# rm -f pileup/'$sname'_cpg_cov5.bed.gz

zcat pileup/'$sname'_cpg_cov5.bed.gz | cut -f1-4 >pileup/'$sname'_cpg_cov5_tmp.bedg
bedGraphToBigWig pileup/'$sname'_cpg_cov5_tmp.bedg '$WZSEQ_REFERENCE'.fai pileup/'$sname'_cpg_cov5.bw
rm -f pileup/'$sname'_cpg_cov5_tmp.bedg
'
  jobname="biscuit_pileup_"$sname
}

function __wgbs_biscuit_pileup_nome {
  # input bamfn
  cmd="
set -xe
cd $base
mkdir -p pileup
biscuit pileup -r $WZSEQ_REFERENCE -N -i bam/$bamfn.bam -o pileup/$bamfn.vcf -q 28
bgzip pileup/$bamfn.vcf
tabix -p vcf pileup/$bamfn.vcf.gz
biscuit vcf2bed -k 0 -c -t hcg pileup/${bamfn}.vcf.gz | cut -f1-5 | LC_ALL=C sort -k1,1 -k2,2n -T pileup | gzip -c >pileup/${bamfn}.hcg.bed.gz
biscuit vcf2bed -k 10 -t hcg pileup/${bamfn}.vcf.gz | LC_ALL=C sort -k1,1 -k2,2n -T pileup >pileup/${bamfn}.hcg.cov10.bedg
biscuit vcf2bed -k 10 -t gch pileup/${bamfn}.vcf.gz | LC_ALL=C sort -k1,1 -k2,2n -T pileup >pileup/${bamfn}.gch.cov10.bedg
bedGraphToBigWig pileup/${bamfn}.hcg.cov10.bedg ${WZSEQ_REFERENCE}.fai pileup/${bamfn}.hcg.cov10.bw
bedGraphToBigWig pileup/${bamfn}.gch.cov10.bedg ${WZSEQ_REFERENCE}.fai pileup/${bamfn}.gch.cov10.bw
"
  jobname="biscuit_pileup_$bamfn"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 8 -memG 200 -ppn 28
}

function wgbs_biscuit_pileup {
  # wgbs_biscuit_pileup [-nome] do

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d pileup ]] || mkdir pileup
  for bamfn0 in bam/*.bam; do
    bamfn=$(basename $bamfn0 .bam);
    # awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    # while read bamfn sread1 sread2; do
    __wgbs_biscuit_pileup
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wgbs_biscuit_pileup_lambdaphage() {

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d pileup ]] || mkdir pileup
  for f in bam/*_lambdaphage.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)
    cmd="
cd $base
biscuit pileup -r $WZSEQ_REFERENCE_LAMBDAPHAGE -i $fn -o pileup/$bfn.vcf -q 28
bgzip pileup/$bfn.vcf
tabix -p vcf pileup/$bfn.vcf.gz
"
    jobname="biscuit_lambdaphage_pileup_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 28
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function __wgbs_vcf2tracks {

  # input: bamfn, items

  cmd="
set -xe
cd $base
mkdir -p tracks

for pt in ${items[@]}; do

  echo processing pileup type \$pt >&2

  biscuit vcf2bed -k 0 -c -t \${pt} pileup/${bamfn}.vcf.gz | cut -f1-5 | LC_ALL=C sort -k1,1 -k2,2n -T tracks | gzip -c >tracks/${bamfn}.\${pt}.bed.gz

  biscuit vcf2bed -k 10 -t \${pt} pileup/${bamfn}.vcf.gz | LC_ALL=C sort -k1,1 -k2,2n -T tracks >tracks/${bamfn}.\${pt}.bedg
  bedGraphToBigWig tracks/${bamfn}.\${pt}.bedg ${WZSEQ_REFERENCE}.fai tracks/${bamfn}.\${pt}.bw

  for i in 100000,100k 1000,1k 100,100; do 
    IFS=\",\"; set \$i; 
    echo processing \$i >&2
    bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w \$1 | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ | bedtools intersect -a - -b tracks/${bamfn}.\${pt}.bedg -wo | bedtools groupby -i - -g 1-3 -c 7 -o mean >tracks/${bamfn}.\${pt}.window\$2.bedg;

    bedGraphToBigWig tracks/${bamfn}.\${pt}.window\$2.bedg ${WZSEQ_REFERENCE}.fai tracks/${bamfn}.\${pt}.window\$2.bw
#rm -f tracks/${bamfn}.\${pt}.window\$2.bedg
    unset IFS
  done
done
"
  jobname="vcf2tracks_$bamfn"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 200 -ppn 1
}

# wgbs_vcf2tracks [-nome] do
function wgbs_vcf2tracks {

  # compute 1) base-pair resolution methylation track 2) mean methylation in window track
  [[ -d tracks ]] || mkdir tracks
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ $1 == "-nome" ]] && items=(hcg gch) || items=(cg)
  for f in pileup/*.vcf.gz; do
    fn=$(readlink -f $f)
    bamfn=$(basename $f .vcf.gz)
    __wgbs_vcf2tracks
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wgbs_vcf2chtracks {

  # compute 1) base-pair resolution methylation track 2) mean methylation in window track
  [[ -d tracks ]] || mkdir tracks
  [[ -d pbs ]] || mkdir pbs
  items=(cg ch)
  base=$(pwd)
  for f in pileup/*.vcf.gz; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .vcf.gz)
    cmd="
cd $base
for pt in ${items[@]}; do
  echo processing pileup type \$pt >&2
  biscuit vcf2bed -k 5 -t \${pt} pileup/${bfn}.vcf.gz | LC_ALL=C sort -k1,1 -k2,2n -T tracks > tracks/${bfn}.\${pt}.bedg
  bedGraphToBigWig tracks/${bfn}.\${pt}.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}.\${pt}.bw

  for i in 100000,100k 1000,1k 100,100; do 
    IFS=\",\"; set \$i; 
    echo processing \$i >&2
    bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w \$1 | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ | bedtools intersect -a - -b tracks/${bfn}.\${pt}.bedg -wo | bedtools groupby -i - -g 1-3 -c 7 -o mean >tracks/${bfn}.\${pt}.window\$2.bedg;

    bedGraphToBigWig tracks/${bfn}.\${pt}.window\$2.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}.\${pt}.window\$2.bw
    rm -f tracks/${bfn}.\${pt}.window\$2.bedg
    unset IFS
  done
done
"
    jobname="vcf2chtracks_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

# obsolete make tracks from USC format beds
function wgbs_methcallbed_to_tracks {

  [[ -d tracks ]] || mkdir tracks;
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  for f in bed/*.bed; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bed)
    cmd="
cd $base
awk 'NR>1&&\$8>3{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$7}' $fn | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ >tracks/${bfn}_cov3.bedg
bedGraphToBigWig tracks/${bfn}_cov3.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}_cov3.bw

bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w 100000 | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ | bedtools intersect -a - -b tracks/${bfn}_cov3.bedg -wo | bedtools groupby -i - -grp 1-3 -c 7 -o mean >tracks/${bfn}_cov3_window100k.bedg;
bedGraphToBigWig tracks/${bfn}_cov3_window100k.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}_cov3_window100k.bw
rm -f tracks/${bfn}_cov3_window100k.bedg

bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w 100 | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ | bedtools intersect -a - -b tracks/${bfn}_cov3.bedg -sorted -wo | bedtools groupby -i - -grp 1-3 -c 7 -o mean >tracks/${bfn}_cov3_window100.bedg;
bedGraphToBigWig tracks/${bfn}_cov3_window100.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}_cov3_window100.bw
rm -f tracks/${bfn}_cov3_window100.bedg

bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w 1000 | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ | bedtools intersect -a - -b tracks/${bfn}_cov3.bedg -sorted -wo | bedtools groupby -i - -grp 1-3 -c 7 -o mean >tracks/${bfn}_cov3_window1k.bedg;
bedGraphToBigWig tracks/${bfn}_cov3_window1k.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}_cov3_window1k.bw
rm -f tracks/${bfn}_cov3_window1k.bedg
"
    jobname="methcallbed_to_tracks_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function wgbs_methylKit_diffmeth {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d methylKit ]] || mkdir methylKit
  awk '/^\[/{p=0}/\[diffmeth\]/{p=1;next} p&&!/^$/' samples |
    while read sname vcfpath1 vcfpath2; do
      echo $sname $vcfpath1
      base1=$(basename $vcfpath1 .vcf.gz) # tumor
      base2=$(basename $vcfpath2 .vcf.gz) # normal
      [[ $base1 == $base2 ]] && base2=$base2"_1"
      dir1=$(dirname $vcfpath1)
      dir2=$(dirname $vcfpath2)
      dir1=${dir1%/pileup}
      dir2=${dir2%/pileup}
      cmd="
cd $base
mkdir methylKit/$sname
biscuit vcf2bed -t cg -u -c $vcfpath1 | pybiscuit.py to_methylKit >methylKit/$sname/$base1.methylKit
biscuit vcf2bed -t cg -u -c $vcfpath2 | pybiscuit.py to_methylKit >methylKit/$sname/$base2.methylKit
~/wzlib/Rutils/bin/bioinfo/methylKit.r diff -t 1,0 -b $base1,$base2 methylKit/$sname/$base1.methylKit,methylKit/$sname/$base2.methylKit -o methylKit/$sname/${base1}_vs_${base2}.tsv -g $WZSEQ_CGIBED_METHYLKIT
rm -f methylKit/$sname/$base1.methylKit
rm -f methylKit/$sname/$base2.methylKit
"
      jobname="methylKit_diff_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 50 -ppn 3
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function wgbs_diffmeth_simple {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d diffmeth ]] || mkdir diffmeth
  awk '/^\[/{p=0}/\[diffmeth\]/{p=1;next} p&&!/^$/' samples |
    while read sname vcfpath1 vcfpath2; do
      cmd="
cd $base
[[ -d diffmeth/$sname ]] || mkdir -p diffmeth/$sname
~/wzlib/bash/wzmethdiff.sh -t $vcfpath1 -n $vcfpath2 -b $base/diffmeth/$sname -v $WZSEQ_REFVERSION -g $WZSEQ_CGIBED -r $WZSEQ_REFERENCE -s $WZSEQ_TSSBED 2>$base/diffmeth/$sname/run.log
"
      jobname="simple_methdiff_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 5 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function wgbs_metilene {

  # output format:
  # chr start stop
  # q-value mean-methylation-difference #CpGs
  # p (Mann-Whitney U) p (2D KS) mean-g1 mean-g2
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d metilene ]] || mkdir metilene
  awk '/^\[/{p=0}/\[diffmeth\]/{p=1;next} p&&!/^$/' samples |
    while read sname vcf1 vcf2; do
      sname1=$(basename $vcf1 .vcf.gz)
      sname2=$(basename $vcf2 .vcf.gz)
      cmd="
cd $base

echo \"\$(date) Converting to bedgraph\"
biscuit vcf2bed -t cg $vcf1 >metilene/$sname1.bedg
biscuit vcf2bed -t cg $vcf2 >metilene/$sname2.bedg

echo \"\$(date) Merge\"
~/software/metilene/default/metilene_input.pl -in1 metilene/$sname1.bedg -in2 metilene/$sname2.bedg --h1 $sname1 --h2 $sname2 --out metilene/$sname.metilene
~/software/metilene/default/metilene_linux64 -t 4 -a $sname1 -b $sname2 metilene/$sname.metilene >metilene/$sname.metilene.result
sort -k4,4n metilene/$sname.metilene.result | head -n 1000 > metilene/$sname.metilene.result.top1k.tsv

echo \"\$(date) Clean\"
rm -f metilene/$sname1.bedg
rm -f metilene/$sname2.bedg
rm -f metilene/$sname.metilene
"
      jobname="metilene_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 4
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

## wgbs_repeat [-nome]
function wgbs_repeat {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d rmsk ]] || mkdir rmsk
  [[ $1 == "-nome" ]] && items=(hcg gch) || items=(cg)
  for bfn0 in bam/*.bam; do
    bfn=$(basename $bfn0 .bam);
    # awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    # while read bfn; do
    cmd="
cd $base

for pt in ${items[@]}; do
  biscuit vcf2bed -t \${pt} -k 5 -c pileup/$bfn.vcf.gz > rmsk/$bfn.\${pt}.bedg
  bedtools map -a $WZSEQ_RMSK -b rmsk/$bfn.\${pt}.bedg -c 4 -o count,mean,collapse >rmsk/$bfn.\${pt}.bedg.rmsk
done
"
    jobname="rmsk_meth_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 4
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wgbs_repeat_diff {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d rmsk ]] || mkdir rmsk
  awk '/^\[/{p=0}/\[diffmeth\]/{p=1;next} p&&!/^$/' samples |
    while read sname vcf1 vcf2; do
      base=$(pwd)
      [[ -d pbs ]] || mkdir pbs
      sname1=$(basename $vcf1 .vcf.gz)
      sname2=$(basename $vcf2 .vcf.gz)
      cmd="
cd $base
mkdir -p rmsk/$sname
biscuit vcf2bed -t cg -k 5 $vcf1 > rmsk/$sname/$sname1.cg.bedg
bedtools map -a $WZSEQ_RMSK -b rmsk/$sname/$sname1.cg.bedg -c 4 -o count,mean,collapse >rmsk/$sname/$sname1.cg.bedg.rmsk
biscuit vcf2bed -t cg -k 5 $vcf2 > rmsk/$sname/$sname2.cg.bedg
bedtools map -a $WZSEQ_RMSK -b rmsk/$sname/$sname2.cg.bedg -c 4 -o count,mean,collapse >rmsk/$sname/$sname2.cg.bedg.rmsk
paste rmsk/$sname/$sname1.cg.bedg.rmsk rmsk/$sname/$sname2.cg.bedg.rmsk >rmsk/$sname/merged.rmsk
awk -f wanding.awk -e '\$9!=\".\"&&\$19!=\".\"{print joinr(1,9)\"\\t\"\$10\"\\t\"\$19\"\\t\"\$20}' rmsk/$sname/merged.rmsk | wzstats Utest -c1 10 -c2 12 - -p 1-9,11 --outfc >rmsk/$sname/$sname.tsv

# awk '\$11<0.01&&\$12>1' rmsk/$sname/$sname.tsv | sort -k12,12nr > rmsk/$sname/$sname.top_hyper.tsv
# awk '\$11<0.01&&\$12<-1' rmsk/$sname/$sname.tsv | sort -k12,12n > rmsk/$sname/$sname.top_hypo.tsv

rm -f rmsk/$sname/$sname1.cg.bedg
rm -f rmsk/$sname/$sname1.cg.bedg.rmsk
rm -f rmsk/$sname/$sname2.cg.bedg
rm -f rmsk/$sname/$sname2.cg.bedg.rmsk
rm -f rmsk/$sname/merged.rmsk
"

      # output format
      # 1-7 basic repeat info
      # 8: mean beta of sample 1
      # 9: mean beta of sample 2
      # 10: p-value
      # 11: log2(fold change)
      # chr1 3000000 3002128 - L1_Mus3 LINE L1  5   0.9028  0.8968  0.917  -0.010
      # chr1 3003152 3003994 - L1Md_F  LINE L1  8   0.90775 0.94025 0.345   0.051
      # chr1 3004270 3005001 + L1_Rod  LINE L1  1   0.895   0.859   0.317   -0.059
      # chr1 3005570 3006764 + Lx9     LINE L1  3   0.95    0.957   0.663   0.011
      jobname="rmsk_diffmeth_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}
