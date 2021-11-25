source $SLURM_ENTRY
#!/bin/bash hostname
# WZSEQ_REFERENCE="/mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/puc19.fasta"
wzref_hg38 # for human genome
pipeline_prepare

while read sname fast5; do
  jump_comments

  pipeline_dependlevel

  memG=8; ppn=12;
  # pipeline_eval 1 GuppyBaseCalling_fast_20210219; 
  pipeline_eval 1 GuppyBaseCalling_rna_fast_20210415;
  
  pipeline_depend 1
  fastq=BaseCalling/${sname}.fastq
  # pipeline_eval 2 Minimap2Aligning_20201023; #DNA
  pipeline_eval 2 Minimap2Aligning_RNA_20210412; #RNA
  
  pipeline_depend 2
  bam=Alignment/${sname}.sorted.bam
  # bam=Alignment/${sname}.sorted.bam
  pipeline_eval 3 NanopolishCallMeth_20201023;

  pipeline_depend 2
  pipeline_eval 4 NanopolishEventsToReference_scale;

  pipeline_depend 2
  pipeline_eval 5 NanopolishEventsToReference;

done <<EOF
# SRR11925984_1  /mnt/isilon/zhou_lab/projects/20210127_check_methylation_NA12878/un_methylated/fast5/NA12878_Unmethylated/
# SRR11925985_1  /mnt/isilon/zhou_lab/projects/20210127_check_methylation_NA12878/methylated/fast5/NA12878_CpGGpC/
# 20201119_CH_HepG2gDNA_Meth_invivo_rep1  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201119_CH_HepG2gDNA_Meth_invivo_rep1
# 20201119_CH_HepG2gDNA_N3_invivo_rep1  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201119_CH_HepG2gDNA_N3_invivo_rep1
# 20201117_CH_HepG2gDNA  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201117_CH_HepG2gDNA/20201117_CH_HepG2gDNA_nolabeling_rep1
20201204CH_HepG2_RNA_ctrl_rep1  /mnt/isilon/zhou_lab/projects/20201207_nanopore_RNA/20201204CH_HepG2_RNA_ctrl_rep1
20201204CH_HepG2_RNA_N3_invivo_rep1  /mnt/isilon/zhou_lab/projects/20201207_nanopore_RNA/20201204CH_HepG2_RNA_N3_invivo_rep1
EOF
