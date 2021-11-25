source $WZSEQ_ENTRY
WZSEQ_REFERENCE="/mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/puc19.fasta"
# wzref_hg19 # for human genome
pipeline_prepare

while read sname fast5; do
  jump_comments

  pipeline_dependlevel

  memG=20; ppn=10;
  pipeline_eval 1 GuppyBaseCalling_fast_20210219;
  
  pipeline_depend 1
  fastq=BaseCalling/${sname}.fastq
  pipeline_eval 2 Minimap2Aligning_20201023;
  
  pipeline_depend 2
  bam=Alignment/${sname}.sorted.bam
  pipeline_eval 3 NanopolishCallMeth_20201023;

  pipeline_depend 3
  pipeline_eval 4 NnopolishEventsToReference_scale;

  pipeline_depend 2
  pipeline_eval 5 NanopolishEventsToReference;

done <<EOF
NB07_fast  /mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/NB07
NB08_fast  /mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/NB08
# NB07_check  /mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/NB07
# NB08_check  /mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/NB08
#20201119_CH_HepG2gDNA_Meth_invivo_rep1  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201119_CH_HepG2gDNA_Meth_invivo_rep1
#20201119_CH_HepG2gDNA_N3_invivo_rep1  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201119_CH_HepG2gDNA_N3_invivo_rep1
#20201117_CH_HepG2gDNA  /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201117_CH_HepG2gDNA/20201117_CH_HepG2gDNA_nolabeling_rep1
EOF
