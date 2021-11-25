source $WZSEQ_ENTRY
wzref_hg19
pipeline_prepare

while read sname fast5; do
  jump_comments

  pipeline_dependlevel

  memG=20; ppn=12;
  pipeline_eval 1 GuppyBaseCalling_20201023;
  
  pipeline_depend 1
  fastq=BaseCalling/${sname}.fastq
  pipeline_eval 2 Minimap2Aligning_20201023;
  
  pipeline_depend 2
  bam=Alignment/${sname}.sorted.bam
  pipeline_eval 3 NanopolishCallMeth_20201023;

done <<EOF
# 20201117_CH_HepG2gDNA  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201117_CH_HepG2gDNA/20201117_CH_HepG2gDNA_nolabeling_rep1/20201117_1919_MN29373_FAO74220_019bf225/fast5
# 20201119_CH_HepG2gDNA_Meth_invivo_rep1  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201119_CH_HepG2gDNA_Meth_invivo_rep1
SRR11925984_1  /mnt/isilon/zhou_lab/projects/20210127_check_methylation_NA12878/
SRR11925985_1  /mnt/isilon/zhou_lab/projects/20210127_check_methylation_NA12878/
EOF
