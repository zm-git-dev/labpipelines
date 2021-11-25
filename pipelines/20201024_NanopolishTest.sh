source $WZSEQ_ENTRY
wzref_hg19
pipeline_prepare

while read fast5 sname; do
  jump_comments

  pipeline_dependlevel

  memG=20; ppn=20;
  pipeline_eval 1 GuppyBaseCalling_20201023;
  
  pipeline_depend 1
  fastq=BaseCalling/${sname}.fastq
  pipeline_eval 2 Minimap2Aligning_20201023;
  
  pipeline_depend 2
  bam=Alignment/${sname}.sorted.bam
  pipeline_eval 3 NanopolishCallMeth_20201023;

done <<EOF
/mnt/isilon/zhou_lab/projects/20201010_nanopore/Test/CallMethylation/methylation_example/fast5_files	methylation_example
EOF
