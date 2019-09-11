
#
################################################################################
# RNA-seq pipeline
################################################################################

function examplepipeline_rnaseq {
cat<<- EOF
=== pipeline 2015-10-01 ===
(+) rnaseq_tophat2_firststrand / rnaseq_tophat2 => (+) rnaseq_cufflinks => (+) rnaseq_cuffmerge => (+) rnaseq_cuffquant (optional) => (+) rnaseq_cuffdiff
EOF

cat<<- EOF
=== pipeline 2016-01-12 ===
 (+) wzseq_fastqc

 => (+) edit samples [alignment]; (+) editrnaseq_tophat2_firststrand / (+) rnaseq_tophat2 / (+) rnaseq_STAR / (+) rnaseq_gsnap / (+) rnaseq_subjunc / (+) rnaseq_mapsplice / (+) rnaseq_hisat2

 => (+) rnaseq_kallisto => rnaseq_kallisto_diff

 => (+) wzseq_qualimap rnaseq_se/rnaseq_pe_stranded

 => [o] rnaseq_splitallele (for ASE) => (+) wzseq_bam_coverage
 => (+) rnaseq_allelomePro

 => (+) rnaseq_cufflinks => (+) rnaseq_cuffmerge => (+) edit samples [diffexp] ; rnaseq_cuffdiff

 => (+) rnaseq_featureCounts => (+) rnaseq_edgeR => (+) rnaseq_DESeq2 => (+) rnaseq_edgeR_rmsk

 => (+) rnaseq_splitstrand_pe / rnaseq_splitstrand_se => (+) rnaseq_count_rmsk_stranded => (+) rnaseq_count_rmsk_stranded_edgeR
 => (+) rnaseq_count_rmsk_unstranded => (+) rnaseq_count_rmsk_unstranded_edgeR

 => (+) rnaseq_dexseq

EOF
}

function pipeline_20170818_rnaseq_SE_stranded {

  pipeline_start=${!#}
  base=$(pwd)
  mkdir -p pbs
  all_aln_jobs=""
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
  while read sname sread1 sread2; do

    pipeline_init

    pipeline_eval 1 __rnaseq_hisat2
    all_aln_jobs="$all_aln_jobs:$jobid"

    pipeline_eval 2 __rnaseq_splitstrand_se

    pipeline_eval 3 __rnaseq_count_rmsk_stranded

  done 

  allbams="bam/*.bam"
  stranded="-s 1"
  pairEnd=""
  depend="-W depend=afterok:$all_aln_jobs"
  echo $depend
  pipeline_eval 4 __rnaseq_featureCounts
}

#########################
## section 1: alignment
#########################
# rnaseq_tophat2
function rnaseq_tophat2() {
  # single-end reads use "." as placeholder
  if [[ ! -s samples ]]; then
    echo "No [samples] file. Abort."
    return 1;
  fi
  base=$(pwd);
  [[ -d bam ]] || mkdir -p bam
  [[ -d pbs ]] || mkdir -p pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      sfile1=$(readlink -f fastq/$sread1);
      [[ $sread2 == "." ]] && sfile2="" || sfile2=$(readlink -f fastq/$sread2); # in case of single ended.
      odir=$(readlink -f bam/$sname);
      # customize the following
      # default: no novel junction, 28 threads
      cmd="
cd $base
tophat2 -p 28 -G $WZSEQ_GTF --library-type fr-unstranded -o $odir --no-novel-juncs $WZSEQ_BOWTIE2_INDEX $sfile1 $sfile2
cd $base/bam
ln -s $sname/accepted_hits.bam $sname.bam
samtools index $sname.bam
samtools flagstat $sname.bam >$sname.bam.flagstat
"
      jobname="tophat_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# rnaseq_tophat2_firststrand samplefn
# first strand
# fr-firststrand:dUTP, NSR, NNSR Same as above except we enforce the rule that
# the right-most end of the fragment (in transcript coordinates) is the first
# sequenced (or only sequenced for single-end reads).
# F2R1 for forward transcript
function rnaseq_tophat2_firststrand() {
  if [[ ! -s samples ]]; then
    echo "file: samples missing. Abort."
    return 1;
  fi
  base=$(pwd);
  [[ -d bam ]] || mkdir -p bam
  [[ -d pbs ]] || mkdir -p pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      sfile1=$(readlink -f fastq/$sread1);
      sfile2=$(readlink -f fastq/$sread2);
      odir=$(readlink -f bam/$sname);
      cmd="
tophat2 -p 28 -G $WZSEQ_GTF --library-type fr-firststrand -o $odir --no-novel-juncs $WZSEQ_BOWTIE2_INDEX $sfile1 $sfile2
cd $base/bam
ln -s $sname/accepted_hits.bam $sname.bam
samtools index $sname.bam
samtools flagstat $sname.bam >$sname.bam.flagstat
"
      jobname="tophat_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function __rnaseq_hisat2_SE {

  # hour=24; memG=30; ppn=14; queue=shortq
  cmd='
set -xe
cd '$base'
mkdir -p bam

export PATH='$PATH':/primary/vari/software/hisat/default
hisat2 -x '$WZSEQ_HISAT2_INDEX' -U '$fastq' -p '$ppn' | samtools view -bS - | samtools sort -T bam/'$sname'.tmp -O bam -o bam/'$sname'.bam
samtools index bam/'$sname'.bam
samtools flagstat bam/'$sname'.bam >bam/'$sname'.bam.flagstat
'
  jobname="hisat_${sname}_SE"
}

function __rnaseq_hisat2_PE {

  cmd='
set -xe
cd '$base'
mkdir -p bam

export PATH=$PATH:/primary/vari/software/hisat/default
hisat2 -x '$WZSEQ_HISAT2_INDEX' -1 '$fastq1' -2 '$fastq2' -p '$ppn' | samtools view -bS - | samtools sort -T bam/'$sname'.tmp -O bam -o bam/'$sname'.bam
samtools index bam/'$sname'.bam
samtools flagstat bam/'$sname'.bam >bam/'$sname'.bam.flagstat
'
  jobname="hisat_${sname}_PE"
}

function rnaseq_hisat2 {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      fastq=fastq/$sread1
      ppn=28
      queue=shortq
      __rnaseq_hisat2_SE
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 240 -ppn $ppn
      [[ $1 == "do" ]] && qsub $pbsfn
  done  
}

# STAR
# make a genome index first
# STAR --runThreadN 28 --runMode genomeGenerate --genomeDir $WZSEQ_STAR_INDEX --genomeFastaFiles $WZSEQ_REFERENCE --sjdbGTFfile $WZSEQ_GTF
# STAR --runThreadN 28 --runMode genomeGenerate --genomeDir . --genomeFastaFiles mm10.fa --sjdbGTFfile ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming [--sjdbOverhang 49 (default 100)]
# note that the GTF must be unzipped!!
# 
# mapping quality:
# 255 = uniquely mapped reads
# 3 = read maps to 2 locations
# 2 = read maps to 3 locations
# 1 = reads maps to 4-9 locations
# 0 = reads maps to 10 or more locations
function __rnaseq_STAR {
  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  cmd='
set -xe
cd '$base'
mkdir '$base'/bam/'$sname'
~/software/STAR/default/bin/Linux_x86_64/STAR --runThreadN '$ppn' --genomeDir '$WZSEQ_STAR_INDEX' --readFilesIn '$fq1' '$fq2' --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/'$sname' --outBAMsortingThreadN '$ppn_sort'
cd bam
ln -s bam/'$sname'/'${sname}'Aligned.sortedByCoord.out.bam bam/'$sname'.bam
samtools index bam/'$sname'.bam
samtools flagstat bam/'$sname'.bam >bam/'$sname'.bam.flagstat
'
  jobname="STAR_$sname"
}

function rnaseq_STAR {

  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
cd $base
mkdir $base/bam/$sname
~/software/STAR/default/source/STAR --runThreadN 28 --genomeDir $WZSEQ_STAR_INDEX --readFilesIn $base/fastq/$sread1 $base/fastq/$sread2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $base/bam/$sname/$sname --outBAMsortingThreadN 6
cd bam
ln -s $sname/${sname}Aligned.sortedByCoord.out.bam $sname.bam
samtools index $sname.bam
samtools flagstat $sname.bam >$sname.bam.flagstat
"
      jobname="STAR_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# gsnap require some prefix setting at compilation, e.g., ./configure --prefix /home/wanding.zhou/software/gsnap/gmap-2015-12-31/
# indexing: ~/software/gsnap/default/util/gmap_build -D ~/references/mm10/gsnap -k 15 -d mm10 ~/references/mm10/bowtie1/*.fa
# note: -d gives the name of the reference, -D gives the location where one looks for the reference with the name (supplied from -d)
# 
# make splice sites and intron sites
# cat ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming | ~/software/gsnap/default/bin/gtf_splicesites > ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gsnap.splicesites
# cat ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming | ~/software/gsnap/default/bin/gtf_introns > ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gsnap.introns
# 
# make splice sites and intron sites from UCSC refGene
# gunzip -c refGene.txt.gz | psl_splicesites -s 1 > foo.splicesites # -s 1 means skip 1 column from the left
# gunzip -c refGene.txt.gz | psl_introns -s 1 > foo.introns

# make .iit files
# cat ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gsnap.introns ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gsnap.splicesites | ~/software/gsnap/default/bin/iit_store -o ~/references/mm10/gsnap/mm10/mm10.maps/Mus_musculus.GRCm38.82.gtf.gsnap.splicesites.iit
# make sure those files are stored within the directory of mm10/gsnap/mm10/mm10.maps/
# WZSEQ_GSNAP_SPLICE=Mus_musculus.GRCm38.82.gtf.gsnap.splicesites

function rnaseq_gsnap {
  # GMAP, the original program is designed for mapping cDNA and EST sequences.
  # GSNAP is the derivative for mapping shot-gun sequencing reads
  # recommend nthreads not exceed 4, otherwise performance is compromised. program design problem
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      
      cmd="
cd $base
~/software/gsnap/default/src/gsnap -A sam -t 4 --gunzip --novelsplicing=1 --use-splicing=$WZSEQ_GSNAP_SPLICE -D $WZSEQ_GSNAP_INDEX -d $WZSEQ_REFVERSION $base/fastq/$sread1 $base/fastq/$sread2 -o bam/$sname.sam
samtools sort -O bam -T bam/$sname.tmp -o bam/$sname.bam bam/$sname.sam
samtools index bam/$sname.bam
samtools flagstat bam/$sname.bam >bam/$sname.bam.flagstat
rm -f bam/$sname.sam
"
      jobname="GSNAP_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 80 -memG 50 -ppn 4
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# mapsplice
# input fastq must be uncompressed and with no space
function rnaseq_mapsplice() {

  if [[ ! -s samples ]]; then
    echo "file: samples missing. Abort"
    return 1
  fi
  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
zcat $base/fastq/$sread1 | sed 's/ /_/g' >$base/fastq/_${sread1}_tmp
zcat $base/fastq/$sread2 | sed 's/ /_/g' >$base/fastq/_${sread2}_tmp
python /home/wanding.zhou/tools/mapsplice/MapSplice-v2.2.0/mapsplice.py -p 28 -o $base/bam/${sname}_mapsplice --bam -c $WZSEQ_BOWTIE1_INDEX -x $WZSEQ_BOWTIE1_INDEX/mm10 -1 $base/fastq/_${sread1}_tmp -2 $base/fastq/_${sread2}_tmp --gene-gtf $WZSEQ_GTF_ENSEMBL
rm -f $base/fastq/_${sread1}_tmp $base/fastq/_${sread2}_tmp
"
      jobname="mapsplice_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# subjunc, subread (if only differential gene expression is concerned, faster), here I
# only use subjunc, which is supposed to be more optimal
# need to build subread index first (for subjunc)
# bin/subread-buildindex -o ~/references/hg19/subread/hg19 ~/references/hg19/hg19.fa
# other options: -S fr -d <min_insert_size> -D <max_insert_size>
function rnaseq_subjunc() {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam

  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do

      if [[ -s bam/$sname.bam ]]; then
        echo "Bam bam/$sname.bam exists. skip"
        continue
      fi
      
      cmd="
cd $base
[[ -d bam/$sname ]] || mkdir bam/$sname
/primary/home/wanding.zhou/tools/subread/subread-1.4.6-p5-Linux-x86_64/bin/subjunc -T 28 -I 16 -i $WZSEQ_SUBREAD_INDEX -r fastq/$sread1 -R fastq/$sread2 --gzFASTQinput -o bam/$sname/$sname.unsrtd.bam --BAMoutput
samtools sort -O bam -o bam/${sname}.bam -T bam/$sname/${sname}.tmp bam/$sname/${sname}.unsrtd.bam
samtools index bam/${sname}.bam
samtools flagstat bam/${sname}.bam >bam/${sname}.flagstat
#rm -f bam/$sname/$sname.unsrtd.bam
"
      jobname="subjunc_${sname}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

#################################################
## section 2: differential expression
#################################################

# rnaseq_cufflinks <do>
function rnaseq_cufflinks() {
  [[ -d cufflinks ]] || mkdir -p cufflinks;
  base=$(pwd);
  for sample in bam/*.bam; do
    sample=$(basename $sample .bam)
    # [[ -s bam/$sample.bam ]] || continue
    cmd="cufflinks -p 28 -g $WZSEQ_GTF -o $base/cufflinks/$sample $base/bam/$sample.bam -q"
    jobname="cufflinks_$sample"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 100 -ppn 28
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

# rnaseq_cuffmerge <do>
# require: cuffmerge, gtf_to_sam, cuffcompare
#
# example of cuffmerge/assemblies.txt
# cufflinks/82_SL121326/transcripts.gtf
# cufflinks/83_SL121327/transcripts.gtf
# cufflinks/84_SL121328/transcripts.gtf
# cufflinks/85_SL121329/transcripts.gtf
# cufflinks/BC/transcripts.gtf
function rnaseq_cuffmerge() {
  [[ -d cuffmerge ]] || mkdir -p cuffmerge;
  base=$(pwd)
  if [[ ! -s cuffmerge/assemblies.txt ]]; then
    echo "cuffmerge/assemblies.txt nonexistent. Create new."
    :>cuffmerge/assemblies.txt
    for f in $base/cufflinks/*; do
      [[ -s $f/transcripts.gtf ]] && echo $f/transcripts.gtf >> cuffmerge/assemblies.txt;
    done
  fi
  cmd="
cd $base/cuffmerge/
cuffmerge -g $WZSEQ_GTF -s $WZSEQ_REFERENCE -p 10 $base/cuffmerge/assemblies.txt
"
  jobname='cuffmerge_'$(basename base)
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 1 -memG 2 -ppn 10
  [[ $1 == "do" ]] && qsub $pbsfn
}

# run cuffquant before this
# cuffquant <do>
function rnaseq_cuffquant() {

  base=$(pwd)
  gtf=$(readlink -f cuffmerge/merged_asm/merged.gtf)
  if [[ ! -s $gtf ]]; then
    echo "cuffmerge/merged_asm/merged.gtf missing. fall back to $WZSEQ_GTF"
    return 1
  fi

  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bf=$(basename $f .bam)
    cmd="
cuffquant $gtf $fn -o $base/cuffmerge/cuffquant_$bf -p 8 -q
"
    jobname="cuffquant_$bf"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 8
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

# require [diffexp] section in "samples"
# format: cond1\tcond2\tbams1\tbams2
# bams1 and bams2 can be comma separated if containing multiple samples
# e.g.
# [diffexp]
# aza     PBS     bam/aza.bam     bam/PBS.bam
# 
# for cuffdiff >2.2 bams can be replaced by cxb files.
function rnaseq_cuffdiff() {
  base=$(pwd);
  [[ -d cuffdiff ]] || mkdir -p cuffdiff;

  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # use merged gtf if available otherwise, use back-up gtf
      gtf=$base/cuffmerge/merged_asm/merged.gtf
      [[ -s $gtf ]] || gtf=$WZSEQ_GTF
      
      cmd="
cd $base
cuffdiff -o $base/cuffdiff/${cond1}_vs_${cond2} -q -b $WZSEQ_REFERENCE -p 8 -L $cond1,$cond2 -u $gtf $bams1 $bams2
"
      jobname="cuffdiff_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 8
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# TODO: cuffnorm
function rnaseq_cuffnorm() {
  return 1
}

###### kallisto
###### index for mm10
# make transcript fasta using the tophat's gtf_to_fasta and some postprocessing
# ~/software/tophat2/default/gtf_to_fasta ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming ~/references/mm10/mm10.fa ~/references/mm10/kallisto/mm10.transcripts.fa
# awk 'match($0,/>(\S+)\s+(\w+)/,a){print ">"a[2];next;}1' ~/references/mm10/kallisto/mm10.transcripts.fa | gzip -c >~/references/mm10/kallisto/mm10.transcripts.fa.gz
# rm ~/references/mm10/kallisto/mm10.transcripts.fa
# ~/software/kallisto/kallisto_linux-v0.42.4/kallisto index ~/references/mm10/kallisto/mm10.transcripts.fa.gz -i ~/references/mm10/kallisto/mm10.kallisto
# 
###### index for human
# cd ~/references/hg19
# ~/software/tophat2/default/gtf_to_fasta gtf/Homo_sapiens.GRCh37.75.gtf.UCSCnaming hg19.fa kallisto/hg19.transcripts.fa
# awk 'match($0,/>(\S+)\s+(\w+)/,a){print ">"a[2];next;}1' kallisto/hg19.transcripts.fa | gzip -c >kallisto/hg19.transcripts.fa.gz
# ~/software/kallisto/kallisto_linux-v0.42.4/kallisto index kallisto/hg19.transcripts.fa.gz -i kallisto/hg19.kallisto
function rnaseq_kallisto {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d kallisto ]] || mkdir kallisto

  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples | while read sname fastq1 fastq2; do
    if [[ "$fastq2" == "." ]]; then
      input2=""
      ## the following is requried, read length is 75 sd 10. this needs to be tuned for each library
      additional_option="--single -l 75 -s 10"
    else
      input2="fastq/$fastq2"
      additional_option=""
    fi
    cmd="
cd $base
[[ -d kallisto/$sname ]] || mkdir kallisto/$sname
~/software/kallisto/kallisto_linux-v0.42.4/kallisto quant -i $WZSEQ_KALLISTO_INDEX $additional_option -o kallisto/$sname fastq/$fastq1 $input2 -t 4
~/wzlib/pyutils/wzseqtk.py ensembl2name --transcript -g $WZSEQ_GTF_ENSEMBL_UCSCNAMING -i kallisto/$sname/abundance.tsv -o kallisto/$sname/abundance.anno.tsv
"
    jobname="kallisto_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 4
    [[ $1 == "do" ]] && qsub $pbsfn
  done 
}

function rnaseq_kallisto_diff {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      snames1=""
      n1=0
      for f in $(echo $bams1 | tr ',' ' '); do
        [[ -d kallisto/$(basename $f .bam) ]] || continue;
        [[ ${#snames1} == 0 ]] || snames1=$snames1",";
        snames1=$snames1$(basename $f .bam);
        ((n1++))
      done
      
      snames2=""
      n2=0
      for f in $(echo $bams2 | tr ',' ' '); do
        [[ -d kallisto/$(basename $f .bam) ]] || continue;
        [[ ${#snames2} == 0 ]] || snames2=$snames2",";
        snames2=$snames2$(basename $f .bam);
        ((n2++))
      done

      if [[ $n1 -gt 1 ]] && [[ $n2 -gt 1 ]]; then # voom requires more than 1 replica
        cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/kallisto_voom.r -a $cond1 -b $cond2 -A $snames1 -B $snames2 -o kallisto/diff_${cond1}_vs_${cond2}.tsv
~/wzlib/pyutils/wzseqtk.py ensembl2name --transcript -g $WZSEQ_GTF_ENSEMBL_UCSCNAMING -i kallisto/diff_${cond1}_vs_${cond2}.tsv -o kallisto/diff_${cond1}_vs_${cond2}.anno.tsv
"
        jobname="kallisto_diff_${cond1}_vs_${cond2}"
        pbsfn=$base/pbs/$jobname.pbs
        pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 3 -ppn 1
        [[ $1 == "do" ]] && qsub $pbsfn
      fi
    done
}

# TODO: sailfish
# TODO: salmon

function __rnaseq_featureCounts {

  # input: base, pairEnd, stranded, allbams
  prog=~/tools/subread/default/bin/featureCounts

  # -p paired-end
  # -g goup by gene_id in GTF
  # -t exon, together with -g, it means count exon (as feature) and gene_id (as meta_feature)
  # -T: number of threads
  # -Q: min mapping quality 20

  # optional:
  # -P -d 50 -D 600 set insert size range to [50,600]
  # -O allowMultiOverlap
  # -M allow multi-mapping
  # -F SAF or GTF, format of annotation file

  cmd="
set -xe
cd $base
mkdir -p featureCounts

# count genes
$prog $pairEnd $stranded -T 5 -t exon -g gene_id -a $WZSEQ_GTF_ENSEMBL_UCSCNAMING -o featureCounts/genes.tsv $allbams --primary -Q 20 --ignoreDup

~/wzlib/pyutils/wzseqtk.py cnt2rpkm -i featureCounts/genes.tsv | ~/wzlib/pyutils/wzseqtk.py ensembl2name --gene -g $WZSEQ_GTF_ENSEMBL_UCSCNAMING -H >featureCounts/genes.rpkm.tsv

# count repeat loci, -f suppresses meta-feature counts
$prog $pairEnd $stranded -T 5 -t exon -g gene_id -f -a $WZSEQ_RMSK_GTF -o featureCounts/rmsk_loci.tsv $allbams --primary -Q 20 --ignoreDup

~/wzlib/pyutils/wzseqtk.py cnt2rpkm -i featureCounts/rmsk_loci.tsv >featureCounts/rmsk_loci.rpkm.tsv

awk -f wanding.awk -e 'NR==FNR{\$2+=1; a[\$1\":\"\$2\"-\"\$3]=joinr(4,NF)}NR!=FNR{if(\$1==\"ID\") print \"strand\tfam1\tfam2\tfam3\t\"\$0; else {ak=\$2\":\"\$3\"-\"\$4; print a[ak],\$0;}}' $WZSEQ_RMSK featureCounts/rmsk_loci.rpkm.tsv >featureCounts/rmsk_loci.rpkm.fam.tsv
"
  jobname="featurecounts_"$(basename $base)
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 5
}

# R-based methods
# this doesn't work very well for repeats, since they lost the family information
# the most generic version
function rnaseq_featureCounts {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  grep '\[experiment\] single-end' samples && pairEnd="" || pairEnd="-P"
  # grep '\[experiment\] unstranded' samples && stranded="-s 0" || stranded="-s 1"
  # stranded="-s 0" # unstranded
  # stranded="-s 1" # for forward-stranded
  # stranded="-s 2" # for reverse-stranded (most of Minmin's library)
  [[ -z "$pairEnd" ]] && stranded=""

  allbams="bam/*.bam"
  [[ -d bam_allele ]] && allbams=$allbams" bam_allele/*.bam"
  __rnaseq_featureCounts
  [[ ${!#} == "do" ]] && qsub $pbsfn
}

function rnaseq_featureCounts_SE_stranded {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  pairEnd=""
  stranded="-s 1"
  allbams="bam/*.bam"
  [[ -d bam_allele ]] && allbams=$allbams" bam_allele/*.bam"
  __rnaseq_featureCounts
  [[ ${!#} == "do" ]] && qsub $pbsfn
}


function rnaseq_featureCounts_SE_revstranded {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  pairEnd=""
  stranded="-s 2"
  allbams="bam/*.bam"
  [[ -d bam_allele ]] && allbams=$allbams" bam_allele/*.bam"
  __rnaseq_featureCounts
  [[ ${!#} == "do" ]] && qsub $pbsfn
}

function rnaseq_edgeR {
  # 2016-02-16 update:
  # switch from GenomicAlignment to featureCounts
  ## original command: ~/wzlib/Rutils/bin/bioinfo/edgeR.r -g $WZSEQ_REFVERSION -G $WZSEQ_GTF_ENSEMBL -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o edgeR/${cond1}_vs_${cond2}_diffexp.tsv 2> edgeR/${cond1}_vs_${cond2}_diffexp.log
  base=$(pwd)
  [[ -d edgeR ]] || mkdir edgeR
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # use merged gtf if available otherwise, use back-up gtf
      gtf=$base/cuffmerge/merged_asm/merged.gtf
      [[ -s $gtf ]] || gtf=$WZSEQ_GTF
      cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/edgeR_featureCounts.r -G $WZSEQ_GTF_ENSEMBL_UCSCNAMING -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o edgeR/${cond1}_vs_${cond2}_diffexp.tsv 2> edgeR/${cond1}_vs_${cond2}_diffexp.log
~/wzlib/pyutils/wzseqtk.py ensembl2name --gene -g $WZSEQ_GTF_ENSEMBL_UCSCNAMING -i edgeR/${cond1}_vs_${cond2}_diffexp.tsv -o edgeR/${cond1}_vs_${cond2}_diffexp.anno.tsv
"
      jobname="edgeR_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 20 -ppn 2
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function rnaseq_edgeR_rmsk {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d featureCounts ]] || mkdir featureCounts
  [[ -d edgeR/rmsk ]] || mkdir -p edgeR/rmsk
  [[ -d edgeR/rmsk_categories ]] || mkdir -p edgeR/rmsk_categories
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/edgeR_featureCounts.r -s 7 -F featureCounts/rmsk_loci.tsv -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o edgeR/rmsk/${cond1}_vs_${cond2}.tsv 2> edgeR/rmsk/${cond1}_vs_${cond2}.log
~/wzlib/Rutils/bin/bioinfo/edgeR_featureCounts.r -s 3 -F featureCounts/rmsk_categories.tsv -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o edgeR/rmsk_categories/${cond1}_vs_${cond2}.tsv 2> edgeR/rmsk_categories/${cond1}_vs_${cond2}.log
"
      jobname="edgeR_rmsk_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 20 -ppn 2
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function rnaseq_DESeq2 {
  base=$(pwd)
  [[ -d DESeq2 ]] || mkdir DESeq2
  [[ -d pbs ]] || mkdir pbs
  grep '\[experiment\] stranded' samples && stranded="" || stranded="--ignoreStrand"
  grep '\[experiment\] single-end' samples && singleEnd="--singleEnd" || singleEnd=""
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # use merged gtf if available otherwise, use back-up gtf
      gtf=$base/cuffmerge/merged_asm/merged.gtf
      [[ -s $gtf ]] || gtf=$WZSEQ_GTF
      cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/DESeq2.r $stranded $singleEnd -g $WZSEQ_REFVERSION -G $WZSEQ_GTF_ENSEMBL -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o DESeq2/${cond1}_vs_${cond2}_diffexp.tsv 2> DESeq2/${cond1}_vs_${cond2}_diffexp.log
"
      jobname="DESeq2_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 36 -memG 100 -ppn 14
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# TODO: EBSeq, Voom
##############################
## section 3: repeats
##############################
function rnaseq_splitstrand_pe {
  # split stranded RNAseq into 2 strands
  # that only works for paired-end
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d stranded ]] || mkdir stranded
  for bam in bam/*.bam; do
    sname=$(basename $bam .bam)
    cmd="
cd $base
minmapq=10

# first read, positive strand
samtools view -H bam/${sname}.bam > stranded/${sname}_p.sam
samtools view -q \$minmapq -f 0x50 -F 0x120 bam/${sname}.bam >> stranded/${sname}_p.sam
# second read, reverse strand
samtools view -q \$minmapq -f 0xa0 -F 0x110 bam/${sname}.bam >> stranded/${sname}_p.sam
samtools view -b stranded/${sname}_p.sam | samtools sort -o stranded/${sname}_p.bam -O bam -T stranded/${sname}_tmp

# first read, reverse strand
samtools view -H bam/${sname}.bam > stranded/${sname}_r.sam
samtools view -q \$minmapq -f 0x60 -F 0x110 bam/${sname}.bam >> stranded/${sname}_r.sam
# second read, positive strand
samtools view -q \$minmapq -f 0x90 -F 0x120 bam/${sname}.bam >> stranded/${sname}_r.sam
samtools view -b stranded/${sname}_r.sam | samtools sort -o stranded/${sname}_r.bam -O bam -T stranded/${sname}_tmp

bedtools genomecov -ibam stranded/${sname}_p.bam -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_COLLATE=C sort -k1,1 -k2,2n -T stranded/ >stranded/${sname}_p.bedg
bedGraphToBigWig stranded/${sname}_p.bedg ${WZSEQ_REFERENCE}.fai stranded/${sname}_p.bw

bedtools genomecov -ibam stranded/${sname}_r.bam -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_COLLATE=C sort -k1,1 -k2,2n -T stranded/ | awk -F\"\\t\" -v OFS=\"\t\" '{print \$1,\$2,\$3,-\$4}' >stranded/${sname}_r.bedg
bedGraphToBigWig stranded/${sname}_r.bedg ${WZSEQ_REFERENCE}.fai stranded/${sname}_r.bw
rm -f stranded/${sname}_p.sam stranded/${sname}_r.sam
rm -f stranded/${sname}_p.bam stranded/${sname}_r.bam
"
    jobname="splitstrand_${sname}"
    pbsfn=pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function __rnaseq_splitstrand_se {
  # hour=12; memG=10; ppn=1; queue=shortq
  # input: base, sname
  cmd='
set -xe
cd '$base'
minmapq=10
mkdir -p stranded

# first read, positive strand
samtools view -H bam/'${sname}'.bam > stranded/'${sname}'_p.sam
samtools view -q $minmapq -F 0x110 bam/'${sname}'.bam >> stranded/'${sname}'_p.sam
samtools view -b stranded/'${sname}'_p.sam | samtools sort -o stranded/'${sname}'_p.bam -O bam -T stranded/'${sname}'_tmp

# first read, reverse strand
samtools view -H bam/'${sname}'.bam > stranded/'${sname}'_r.sam
samtools view -q $minmapq -f 0x10 -F 0x100 bam/'${sname}'.bam >> stranded/'${sname}'_r.sam
samtools view -b stranded/'${sname}'_r.sam | samtools sort -o stranded/'${sname}'_r.bam -O bam -T stranded/'${sname}'_tmp

totalreads=$(awk -F" " '\''NR==1{print $1}'\'' bam/'${sname}'.bam.flagstat)

bedtools genomecov -ibam stranded/'${sname}'_p.bam -g '${WZSEQ_REFERENCE}'.fai -bga -split | LC_COLLATE=C sort -k1,1 -k2,2n -T stranded/ >stranded/'${sname}'_p.bedg
awk -v totalreads=$totalreads '\''{print $1,$2,$3,$4/totalreads*10000000}'\'' stranded/'${sname}'_p.bedg >stranded/'${sname}'_p_norm.bedg
bedGraphToBigWig stranded/'${sname}'_p_norm.bedg '${WZSEQ_REFERENCE}'.fai stranded/'${sname}'_p_norm.bw
bedGraphToBigWig stranded/'${sname}'_p.bedg '${WZSEQ_REFERENCE}'.fai stranded/'${sname}'_p.bw

## negative strand always have negative counts
bedtools genomecov -ibam stranded/'${sname}'_r.bam -g '${WZSEQ_REFERENCE}'.fai -bga -split | LC_COLLATE=C sort -k1,1 -k2,2n -T stranded/ | awk -F"\t" -v OFS="\t" '\''{print $1,$2,$3,-$4}'\'' >stranded/'${sname}'_r.bedg
awk -v totalreads=$totalreads '\''{print $1,$2,$3,$4/totalreads*10000000}'\'' stranded/'${sname}'_r.bedg >stranded/'${sname}'_r_norm.bedg
bedGraphToBigWig stranded/'${sname}'_r_norm.bedg '${WZSEQ_REFERENCE}'.fai stranded/'${sname}'_r_norm.bw
bedGraphToBigWig stranded/'${sname}'_r.bedg '${WZSEQ_REFERENCE}'.fai stranded/'${sname}'_r.bw

rm -f stranded/'${sname}'_p.sam stranded/'${sname}'_r.sam
rm -f stranded/'${sname}'_p.bam stranded/'${sname}'_r.bam
'
  jobname="splitstrand_${sname}_se"
}

function rnaseq_splitstrand_se {
  # split stranded RNAseq into 2 strands
  # that only works for single-end
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  for bam in bam/*.bam; do
    sname=$(basename $bam .bam)
    __rnaseq_splitstrand_se
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

## the following restrict to intergenic ERVs
## for f in rmsk/*.tsv; do echo $f; bedtools intersect -a <(sed -n '2,$p' $f) -b ERVIntergenic.bed -sorted -wo | cut -f1-13 | cat <(echo -e "chrm\tbeg\tend\tstr\tclass1\tclass2\tclass3\ttlen\tBaseNeg\tBasePos\tRPKMneg\tRPKMpos\tRPKM") - >rmsk/$(basename $f).LTRintergenic.bed; done
## for f in rmsk/*.tsv; do echo $f; bedtools intersect -a <(sed -n '2,$p' $f) -b ERVIntergenic.bed -sorted -wo | cut -f1-13 | awk -f wanding.awk -e '$9<0 && $10>0 && min(-$9,$10)/max(-$9,$10)>=0.5' | cat <(echo -e "chrm\tbeg\tend\tstr\tclass1\tclass2\tclass3\ttlen\tBaseNeg\tBasePos\tRPKMneg\tRPKMpos\tRPKM") - >rmsk/$(basename $f).LTRintergenic.doubleStranded.bed; done

function __rnaseq_count_rmsk_stranded {

  # input: sname, base
  # hour=24; memG=20; ppn=2;
  cmd='
set -xe
cd '$base'
mkdir -p rmsk

bedtools intersect -a '$WZSEQ_RMSK' -b stranded/'${sname}'_r.bedg -wao -sorted | awk -f wanding.awk -e '\''{print joinr(1,7)"\t"$11*$12;}'\'' | bedtools groupby -g 1-7 -c 8 -o sum > rmsk/'${sname}'_r.tsv

## find all cummulative base counts
## all bases mapped to positive strand
all_p=$(awk '\''{a+=($3-$2)*$4}END{print a}'\'' stranded/'${sname}'_p.bedg)
## all bases mapped to negative strand
all_r=$(awk '\''{a+=($3-$2)*$4}END{print a}'\'' stranded/'${sname}'_r.bedg)

nmap=$(($all_p-$all_r))

echo "cumulative base cnt: $nmap"

bedtools intersect -a '$WZSEQ_RMSK' -b stranded/'${sname}'_p.bedg -wao -sorted | awk -f wanding.awk -e '\''{print joinr(1,7)"\t"$11*$12;}'\'' | bedtools groupby -g 1-7 -c 8 -o sum > rmsk/'${sname}'_p.tsv

paste rmsk/'${sname}'_p.tsv rmsk/'${sname}'_r.tsv | awk '\''$2==$10 && $5==$13'\'' | cut -f1-7,8,16 | awk -v alln=$nmap '\''BEGIN{print "chrm\tbeg\tend\tstrand\tcat1\tcat2\tcat3\ttlen\tposBaseCnt\tnegBaseCnt\tposRPKM\tnegRPKM\tRPKM"}{tlen=$3-$2; pp=$8/tlen*1000/alln*1000000; rr=$9/tlen*1000/alln*1000000; print $0"\t"tlen"\t"pp"\t"rr"\t"pp-rr}'\'' >rmsk/'$sname'.tsv

## count category
awk -v alln=$nmap '\''{n=$8-$9; a[$5]+=n; b[$6]+=n; c[$7]+=n;} END{print "Genome\t0\t"alln"\t1.0"; for (i in a) {print i"\t1\t"a[i]"\t"a[i]/alln} for(i in b){print i"\t2\t"b[i]"\t"b[i]/alln} for(i in c){print i"\t3\t"c[i]"\t"c[i]/alln}}'\'' rmsk/'$sname'.tsv | sort -k2,2n -k1,1 >rmsk/'$sname'.tsv.categories

rm -f rmsk/'${sname}'_p.tsv rmsk/'${sname}'_r.tsv
'
  # for locating double-strand transcription
  # awk -v OFS="\t" -f wanding.awk -e 'max($14,-$15)>0 && min($14,-$15)/max($14,-$15)>0.5 && min($14,-$15)>100{print $_"\t"($14-$15)/($8-$9+10)}' merged.rmsk.bed | sort -k16,16nr >merged.rmsk.bed.double.up
  jobname="rmsk_${sname}"
}

function rnaseq_count_rmsk_stranded {
  # require stranded/, this actually gives the cumulative
  # base counts, not the read counts
  # equivalent to RPKM because base counts are normalized by total base counts
  # both numerator and denominator are scaled by the read length
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  for bedg_r in stranded/*_r.bedg; do
    sname=$(basename $bedg_r _r.bedg)
    __rnaseq_count_rmsk_stranded
    pbsfn=$base/pbs/$jobname.pbs
    hour=24; memG=20; ppn=2
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour $hour -memG $memG -ppn $ppn
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function rnaseq_count_rmsk_stranded_edgeR {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d rmsk ]] || mkdir rmsk
  [[ -d rmsk/diff ]] || mkdir rmsk/diff
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # skip non-replicated designs
      # Remark: At non-replicated design, we would sometimes have
      # "f() values at end points not of opposite sign" error
      # for "deviance" estimate of common dispersion
      # othertimes, it may work (usually for genes but not rmsk)
      # [[ -n $(grep -o "," <<< "$bams1") ]] || continue
      # [[ -n $(grep -o "," <<< "$bams2") ]] || continue

      # make sure the tsv all exists
      allexist=1
      tsv1=""
      for i in ${bams1//,/ }; do
        _tsv=rmsk/$(basename $i .bam).tsv 
        if [[ ! -e $_tsv ]]; then
          allexist=0;
          break;
        fi
        [[ -n $tsv1 ]] && tsv1=$tsv1","
        tsv1=$tsv1""$_tsv;
      done

      tsv2=""
      for i in ${bams2//,/ }; do
        _tsv=rmsk/$(basename $i .bam).tsv
	      if [[ ! -e $_tsv ]]; then
          allexist=0;
          break;
        fi
        [[ -n $tsv2 ]] && tsv2=$tsv2","
        tsv2=$tsv2""$_tsv;
      done

      [[ $allexist == 0 ]] && continue;

      # compare rmsk
      cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/edgeR_rmsk.r -a $cond1 -b $cond2 -A $tsv1 -B $tsv2 -o rmsk/diff/${cond1}_vs_${cond2}.diff.tsv
"
      jobname="rmsk_diff_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn

    done
}

function rnaseq_count_rmsk_unstranded {

  # this count the cumulative base coverage
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d rmsk ]] || mkdir rmsk
  for bam in bam/*.bam; do
    sname=$(basename $bam .bam)
    cmd="
cd $base
minmapq=10
samtools view -q \$minmapq -b $bam | bedtools genomecov -ibam - -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_COLLATE=C sort -k1,1 -k2,2n -T rmsk/ >rmsk/${sname}.bedg

bedGraphToBigWig rmsk/${sname}.bedg ${WZSEQ_REFERENCE}.fai rmsk/${sname}.bw

all=\$(awk '{a+=(\$3-\$2)*\$4}END{print a}' rmsk/${sname}.bedg)

bedtools intersect -a $WZSEQ_RMSK -b rmsk/${sname}.bedg -wao -sorted | awk -f wanding.awk -e '{print joinr(1,7)\"\\t\"\$11*\$12;}' | bedtools groupby -g 1-7 -c 8 -o sum | awk -v alln=\$all 'BEGIN{print \"chrm\tbeg\tend\tstrand\tcat1\tcat2\tcat3\ttlen\tbaseCnt\tRPKM\"}{tlen=\$3-\$2; rpkm=\$8/tlen*1000/alln*1000000; print \$0\"\t\"tlen\"\t\"rpkm}' > rmsk/${sname}.tsv

## count category 
awk -v all=\$all '{a[\$5]+=\$8; b[\$6]+=\$8; c[\$7]+=\$8;} END{print \"Genome\t0\t\"all\"\t1.0\"; for (i in a) {print i\"\t1\t\"a[i]\"\t\"a[i]/all} for(i in b){print i\"\t2\t\"b[i]\"\t\"b[i]/all} for(i in c){print i\"\t3\t\"c[i]\"\t\"c[i]/all}}' rmsk/$sname.tsv | sort -k2,2n -k1,1 >rmsk/$sname.tsv.categories

rm -f rmsk/${sname}.bedg
"
    # The following only count reads
    # bedtools coverage -a $WZSEQ_RMSK -b bam/${bam}.bam -sorted -split -counts > repeatmasker/$bam.rmsk
    # awk '{a[\$5]+=\$8;b[\$6]+=\$8;c[\$7]+=\$8}END{for(i in a){print "1\t"i"\t"a[i]}; for(i in b){print "2\t"i"\t"b[i]}; for(i in c){print "3\t"i"\t"c[i]}}' repeatmasker/${bam}.rmsk | sort -k1,1 -k2,2nr > repeatmasker/${bam}.rmsk.cat
    jobname="rmsk_${sname}"
    pbsfn=pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 2 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function rnaseq_count_rmsk_unstranded_edgeR {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d rmsk ]] || mkdir rmsk
  [[ -d rmsk/diff ]] || mkdir rmsk/diff
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # skip non-replicated designs
      # Remark: At non-replicated design, we would sometimes have
      # "f() values at end points not of opposite sign" error
      # for "deviance" estimate of common dispersion
      # othertimes, it may work (usually for genes but not rmsk)
      # [[ -n $(grep -o "," <<< "$bams1") ]] || continue
      # [[ -n $(grep -o "," <<< "$bams2") ]] || continue

      # make sure the tsv all exists
      allexist=1
      tsv1=""
      for i in ${bams1//,/ }; do
        _tsv=rmsk/$(basename $i .bam).tsv 
        if [[ ! -e $_tsv ]]; then
          allexist=0;
          break;
        fi
        [[ -n $tsv1 ]] && tsv1=$tsv1","
        tsv1=$tsv1""$_tsv;
      done

      tsv2=""
      for i in ${bams2//,/ }; do
        _tsv=rmsk/$(basename $i .bam).tsv
	      if [[ ! -e $_tsv ]]; then
          allexist=0;
          break;
        fi
        [[ -n $tsv2 ]] && tsv2=$tsv2","
        tsv2=$tsv2""$_tsv;
      done

      [[ $allexist == 0 ]] && continue;

      # compare rmsk
      cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/edgeR_rmsk.r -U -a $cond1 -b $cond2 -A $tsv1 -B $tsv2 -o rmsk/diff/${cond1}_vs_${cond2}.diff.tsv
"
      jobname="rmsk_diff_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn

    done
}

###################################
# section 4: alternative splicing
###################################

function rnaseq_dexseq {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d DEXSeq ]] || mkdir DEXSeq
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # skip non-replicated designs
      # Remark: I haven't figured out how to make DEX-seq
      # work without replica
      [[ -n $(grep -o "," <<< "$bams1") ]] || continue
      [[ -n $(grep -o "," <<< "$bams2") ]] || continue
      outdir=DEXSeq/${cond1}_vs_${cond2}
      cmd="
cd $base
[[ -d $outdir ]] || mkdir $outdir

quant1=\"\"
for b in $(echo $bams1 | sed 's/,/ /g'); do
  bb=\$(basename \$b .bam)
  python ~/.Renv/versions/3.2.3/lib64/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam $WZSEQ_GTF_DEXSEQ \$b $outdir/\$bb.quant.txt
  [[ -z \"\$quant1\" ]] || quant1=\$quant1\",\"
  quant1=\$quant1\"$outdir/\$bb.quant.txt\"
done

quant2=\"\"
for b in $(echo $bams2 | sed 's/,/ /g'); do
  bb=\$(basename \$b .bam)
  python ~/.Renv/versions/3.2.3/lib64/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam $WZSEQ_GTF_DEXSEQ \$b $outdir/\$bb.quant.txt
  [[ -z \"\$quant2\" ]] || quant2=\$quant2\",\"
  quant2=\$quant2\"$outdir/\$bb.quant.txt\"
done

echo \"input quants:\" \$quant1
echo \"input quants:\" \$quant2

~/wzlib/Rutils/bin/bioinfo/DEXSeq.r -G $WZSEQ_GTF_DEXSEQ -a mut -b wt -A \$quant1 -B \$quant2 -o $outdir/${cond1}_vs_${cond2}.tsv

~/wzlib/pyutils/wzseqtk.py ensembl2name --gene -g $WZSEQ_GTF_ENSEMBL_UCSCNAMING -c 2 -i $outdir/${cond1}_vs_${cond2}.tsv -o $outdir/${cond1}_vs_${cond2}.anno.tsv
"
      jobname="dexseq_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# TODO: MATS 

# TODO: MISO (pretty old)

# TODO: HTseq-count

# TODO: RSEM
# ~/tools/rsem/rsem-1.2.22/rsem-prepare-reference
# rsem-calculate-expression -p 20 --calc-ci --ci-memory 12294 --bowtie-chunkmbs 2000 --paired-end --bowtie-path $WZSEQ_BOWTIE1 --rsem-index RSEM

# TODO: ALEXA-Seq http://www.alexaplatform.org/alexa_seq/

# TODO: Trinity for RNA-seq
# de novo assembly

##########################################
# section 5: allele-specific expression
##########################################

# wzbam splitallele -bam bam/Dot1LE10.bam -snp ~/references/mm10/snp/pairwise/BL6_vs_JF1.snp.mm10.bed -name BL6,JF1 -o testsplit

function rnaseq_splitallele {
  # [split_allele]
  # PL_female_wt.bam        /primary/vari/genomicdata/genomes/mm10/snp/pairwise/BL6_vs_JF1.snp.mm10.bed     BL6,JF1
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam_allele ]] || mkdir bam_allele
  awk '/^\[/{p=0}/\[split_allele\]/{p=1;next} p&&!/^$/' samples |
    while read bamfile snpfile background; do
      bfn=$(basename $bamfile .bam)
      cmd="
cd $base
wzbam splitallele -bam bam/$bamfile -snp $snpfile -name $background -o bam_allele/$bfn 2>bam_allele/$bfn.splitallele.log
IFS=',' read -r -a samples <<< \"$background\"
sub1=bam_allele/$bfn.\${samples[0]}.bam
sub2=bam_allele/$bfn.\${samples[1]}.bam
samtools index \$sub1
samtools index \$sub2
samtools flagstat \$sub1 >\$sub1.flagstat
samtools flagstat \$sub2 >\$sub2.flagstat

# make coverage tracks
minmapq=10
for sub in \$sub1 \$sub2; do
  subbase=\$(basename \$sub .bam)
  bedtools genomecov -ibam \$sub -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_ALL=C sort -k1,1 -k2,2n -T bam_allele/  >bam_allele/\${subbase}.coverage.bedg
  bedGraphToBigWig bam_allele/\${subbase}.coverage.bedg ${WZSEQ_REFERENCE}.fai bam_allele/\${subbase}.coverage.bw
  rm -f bam_allele/\${subbase}.coverage.bedg

  samtools view -q \$minmapq -b \$sub | bedtools genomecov -ibam stdin -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_ALL=C sort -k1,1 -k2,2n -T bam_allele/  >bam_allele/\${subbase}.coverage.q10.bedg
  bedGraphToBigWig bam_allele/\${subbase}.coverage.q10.bedg ${WZSEQ_REFERENCE}.fai bam_allele/\${subbase}.coverage.q\$minmapq.bw
  rm -f bam_allele/\${subbase}.coverage.q10.bedg
done
"
      jobname="split_allele_$bfn"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 3 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function rnaseq_allelomePro {

  # see wzlib/other/allelomePro.config for example of config file
  # abbreviations:
  # MAT: maternal bias
  # PAT: paternal bias
  # strain1 bias, strain2 bias
  # BAE: Biallelic expression
  # NI: non-informative (low SNP coverage)
  # NS: no SNP
  # I_score: imprinted score, S_score: strain bias score
  
  # RPSM: (10^6*A)/(B*C), A: number of mappable reads at the given single nucleotide position
  # B: number of all mappable reads in the sample, C: read length
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[allelomePro\]/{p=1;next} p&&!/^$/' samples |
    while read config; do
      cmd="
cd $base
mkdir -p allelomePro/$(basename $config .config)
~/tools/allelomepro/Allelome_PRO/allelome_pro.sh -c allelomePro/$config
"
      jobname="allelomePro_${config}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}


###########################
# RSeQC
###########################
# need to define WZSEQ_RSEQ_GENE_BED
function rnaseq_rseqc() {
  base=$(pwd);
  [[ -d rseqc ]] || mkdir rseqc
  [[ -d pbs ]] || mkdir pbs
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)

    cmd="
cd $base

# compute and plot SAM quality score
echo `date` 'running read_quality.py ...' 1>&2
read_quality.py -i $f -o rseqc/${bfn}_read_quality

# computes nucleotide composition tables and make plots
echo `date` 'running read_NVC.py ...' 1>&2
read_NVC.py -i $f -o rseqc/${bfn}_nucleotide_composition

# number of reads mapped to each genomic feature, CDS_Exons, 3'UTR_Exons, etc
echo `date` 'running read_distribution.py ...' 1>&2
read_distribution.py -i $fn -r $WZSEQ_RSEQC_GENE_BED >rseqc/${bfn}_read_distribution

# 5'->3' gene coverage metrics for each sample (coverage uniformity)
echo `date` 'running geneBody_coverage.py ...' 1>&2
geneBody_coverage.py -i $fn -r $WZSEQ_RSEQC_GENE_BED -o rseqc/${bfn}_genebody_coverage
"
    jobname="rseqc_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}
