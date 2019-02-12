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

