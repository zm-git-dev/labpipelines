source $WZSEQ_ENTRY
wzref_hg38
pipeline_prepare

# cd /mnt/isilon/zhou_lab/projects/20191212_GEO_datasets/GSE128731_TruSeqEPIC
while read sname srr_ids; do
jump_comments

srr_ids=${srr_ids//,/ };
pipeline_dependlevel
hour=48; memG=10; ppn=1;
pipeline_eval 1 __wzseq_fastq_dump_PE;

done <<EOF
whole_blood_GSM3683954  SRR9888304
whole_blood_GSM3683960  SRR9888310
whole_blood_GSM3683961  SRR9888311
whole_blood_GSM3683967  SRR9888317
whole_blood_GSM3683968  SRR9888318
whole_blood_GSM3683974  SRR9888324
whole_blood_GSM3683975  SRR9888325
CD4_GSM3683979  SRR9888329
Neutrophils_GSM3683983  SRR9888333
Neutrophils_GSM3683987  SRR9888337
CD4_GSM3683991  SRR9888341
EOF
