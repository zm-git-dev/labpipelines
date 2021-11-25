source $SLURM_ENTRY
# WZSEQ_REFERENCE="/mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/puc19.fasta"
#!/bin/bash hostname
wzref_hg38 # for human genome
pipeline_prepare

while read sname fast5; do
  jump_comments

  pipeline_dependlevel

  memG=12; ppn=8;
  pipeline_eval 1 GuppyBaseCalling_fast_20210219;
  
  pipeline_depend 1
  fastq=BaseCalling/${sname}.fastq
  pipeline_eval 2 Minimap2Aligning_20201023;
  
  pipeline_depend 2
  bam=Alignment/${sname}.sorted.bam
  pipeline_eval 3 NanopolishCallMeth_20201023;

  pipeline_depend 2
  pipeline_eval 4 NanopolishEventsToReference_scale;

  pipeline_depend 2
  pipeline_eval 5 NanopolishEventsToReference;

done <<EOF
# NB07_accurate  /mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/NB07
# 20201119_CH_HepG2gDNA_Meth_invivo_rep1  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201119_CH_HepG2gDNA_Meth_invivo_rep1
# N2a_LREM_seq /mnt/isilon/zhou_lab/projects/20210517_Nanopore_LR_EM_seq/20210517_N2a_LR_EM_seq_Sun_et_al/fast5
# GM12878_rad004 /mnt/isilon/zhou_lab/projects/20210521_GM12878_gDNA_RAD004/GM12878_rad004/20210521_1406_MN36407_FAP29302_9df97345/fast5
N2a_LREM_seq_small /mnt/isilon/zhou_lab/projects/20210517_Nanopore_LR_EM_seq/20210517_N2a_LR_EM_seq_Sun_et_al/small_set_check/fast5
EOF
