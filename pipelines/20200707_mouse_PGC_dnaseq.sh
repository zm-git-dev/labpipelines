source $WZSEQ_ENTRY
wzref_mm10
pipeline_prepare

# cd /mnt/isilon/zhou_lab/projects/20200706_PGC_DNaseseq
while read srr_ids sname; do
jump_comments

srr_ids=${srr_ids//,/ };
pipeline_dependlevel
memG=20; ppn=1;
pipeline_eval 1 __sra_fasterq_dump_SE_20200706;

done <<EOF
SRR6519360      E9.5_PGC_DNase-seq
SRR6519361      E9.5_PGC_DNase-seq
SRR6519362      E10.5_PGC_DNase-seq
SRR6519363      E10.5_PGC_DNase-seq
SRR6519364      E12.5_female_PGC_DNase-seq
SRR6519365      E12.5_female_PGC_DNase-seq
SRR6519366      E12.5_female_soma_DNase-seq
SRR6519367      E12.5_male_PGC_DNase-seq
SRR6519368      E12.5_male_PGC_DNase-seq
SRR6519369      E12.5_male_soma_DNase-seq
SRR6519370      E13.5_female_PGC_DNase-seq
SRR6519371      E13.5_female_PGC_DNase-seq
SRR6519372      E13.5_male_PGC_DNase-seq
SRR6519373      E13.5_male_PGC_DNase-seq
SRR6519374      E14.5_female_PGC_DNase-seq
SRR6519375      E14.5_female_PGC_DNase-seq
SRR6519376      E14.5_male_PGC_DNase-seq
SRR6519377      E14.5_male_PGC_DNase-seq
SRR6519378      E16.5_female_PGC_DNase-seq
SRR6519379      E16.5_female_PGC_DNase-seq
SRR6519380      E16.5_male_PGC_DNase-seq
SRR6519381      E16.5_male_PGC_DNase-seq
SRR6519382      E16.5_male_soma_DNase-seq
SRR6519383      E16.5_male_soma_DNase-seq
SRR6519384      BVSC_ESC_DNase-seq
SRR6519385      BVSC_ESC_DNase-seq
SRR6519386      ESC_DNase-seq
SRR6519387      ESC_DNase-seq
SRR6519388      EpiLC_DNase-seq
SRR6519389      EpiLC_DNase-seq
EOF
