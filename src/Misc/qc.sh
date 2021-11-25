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
