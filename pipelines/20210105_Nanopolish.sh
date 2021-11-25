source $WZSEQ_ENTRY
# WZSEQ_REFERENCE="/mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/puc19.fasta"
wzref_hg19 # for human genome
pipeline_prepare

while read sname fast5; do
  jump_comments

  pipeline_dependlevel

  memG=20; ppn=20;
  pipeline_eval 1 GuppyBaseCalling_20201023;
  
  pipeline_depend 1
  fastq=BaseCalling/${sname}.fastq
  pipeline_eval 2 Minimap2Aligning_20201023;
  
  pipeline_depend 2
  bam=Alignment_timp/${sname}.sorted.bam
  pipeline_eval 3 NanopolishCallMeth_20201023;

done <<EOF
NA12_small  /mnt/isilon/zhou_lab/projects/20210127_check_methylation_NA12878/FAST5_small
# NB08  /mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/NB08
# 20201117_CH_HepG2gDNA  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201117_CH_HepG2gDNA/20201117_CH_HepG2gDNA_nolabeling_rep1
# 20201119_CH_HepG2gDNA_Meth_invivo_rep1  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201119_CH_HepG2gDNA_Meth_invivo_rep1
# 20201119_CH_HepG2gDNA_N3_invivo_rep1  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201119_CH_HepG2gDNA_N3_invivo_rep1
EOF
