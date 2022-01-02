source $SLURM_ENTRY
pipeline_prepare

function runBwa(){
cmd='
bwa mem -M -t '${ppn}' ref/primer.fa fastq/'${fastq}' | samtools sort -T bam/'${sname}' -o bam/'${sname}'.sorted.bam
samtools index -@ '${ppn}' bam/'${sname}'.sorted.bam
samtools flagstat bam/'${sname}'.sorted.bam > bam/'${sname}'.sorted.bam.flagstat
'
jobname=${sname}'_bwa'
}

for fastq in `ls fastq/`; do
  sname=${fastq//.fastq.gz/};
  echo ${sname};
  pipeline_dependlevel
  
  ppn=24; memG=100;
  
  pipeline_eval 1 runBwa;
  
  #memG=20; ppn=10
  #fastq=fastq/${sname}.fastq.gz
  #fastq=TrimGalore/${sname}_trimmed.fq.gz
  #pipeline_eval 2 __bwa_aln_se_20200716;

  #bam="bam/"$sname".bam"
  #pipeline_eval 3 __markDupAndFilter_20200719;

  #species="mm"
  #pipeline_eval 4 bam2bigwig;
done;
