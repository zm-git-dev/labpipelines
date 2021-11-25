source $WZSEQ_ENTRY
# WZSEQ_REFERENCE="/mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/puc19.fa.txt"
wzref_hg38
pipeline_prepare

while read sname fast5; do
  jump_comments

  pipeline_dependlevel

  memG=20; ppn=24;
  pipeline_eval 1 GuppyBaseCalling_hac_fast5_out_20210223;
  
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
20201117_CH_HepG2gDNA /mnt/isilon/zhou_lab/projects/20201120_nanopore/20201117_CH_HepG2gDNA/20201117_CH_HepG2gDNA_nolabeling_rep1/20201117_1919_MN29373_FAO74220_019bf225/fast5
# cont  /mnt/isilon/zhou_lab/projects/20201120_nanopore/test/cont
# N3  /mnt/isilon/zhou_lab/projects/20201120_nanopore/test/N3
# NB07  /mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/NB07
# NB08  /mnt/isilon/zhou_lab/projects/20201220_pUC19_plasmid/NB08
# 20210107_1Dnorepair_HepG2gDNA_MGO_vivo /mnt/isilon/zhou_lab/projects/20210107_nanopore_1D_no_repair
# 20210107_1Dnorepair_HepG2gDNA_noLabeling /mnt/isilon/zhou_lab/projects/20210107_nanopore_1D_no_repair
# 20210107_1Dnorepair_HepG2gDNA_N3-Ke_vivo /mnt/isilon/zhou_lab/projects/20210107_nanopore_1D_no_repair
# M  /mnt/isilon/zhou_lab/projects/20210127_check_methylation_NA12878/methylated/fast5
# U  /mnt/isilon/zhou_lab/projects/20210127_check_methylation_NA12878/un_methylated/fast5
EOF
