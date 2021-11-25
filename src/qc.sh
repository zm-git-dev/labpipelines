##########
## fastqc
##########

function __wzseq_fastqc() {
  # hour=12; memG=5; ppn=1
  # note that fastq_sname include _R1/R2
  cmd='
set -xe
cd '$base'
mkdir -p fastqc/'$fastq_sname'
fastqc -f fastq '$fastq' -o fastqc/'$fastq_sname'
mkdir -p multiqc/raw/fastqc/
ln -sf `readlink -f fastqc` multiqc/raw/fastqc/
'
  jobname='fastqc_'$fastq_sname
}

function __fastqc20200718 {
  cmd='
set -xe
cd '$base'
mkdir -p fastqc
fastqc -f fastq -t 10 -o fastqc '$fastq'
'
  jobname='fastqc_'$sname
}

#####################
## Qualimap
#####################

function __wzseq_qualimap_bamqc() {
  # NOTE: be sure to set memG correctly, otherwise, it would be very very slow and weird
  # hour=24; memG=200; ppn=28; queue=longq
  # input_bam=
  # output_sname=
  cmd='
cd '$base'
qualimap --java-mem-size='$memG'G bamqc -nt '$ppn' -bam '$input_bam' -outdir qualimap/'$output_sname' -c
mkdir -p multiqc/raw/qualimap
ln -sf `readlink -f qualimap/'$sname'` multiqc/raw/qualimap/
'
  jobname='qualimap_bamqc_'$sname
}

function __wzseq_qualimap_rnaseqSE() {
  # hour=24; memG=50; ppn=10
  cmd='
cd '$base'
qualimap --java-mem-size=10G rnaseq -bam '$input_bam' -gtf '$WZSEQ_GTF' -outdir qualimap/rnaseq_'$sname'
'
  jobname='qualimap_rnaseqSE_'$sname
}

function __wzseq_qualimap_rnaseqPEstranded() {
  # hour=24; memG=50; ppn=10
  cmd='
cd '$base'
qualimap --java-mem-size=10G rnaseq -bam '$input_bam' -gtf '$WZSEQ_GTF' -outdir qualimap/rnaseq_'$sname' --paired -p strand-specific-forward
'
  jobname='qualimap_rnaseqPEst_'$sname
}

