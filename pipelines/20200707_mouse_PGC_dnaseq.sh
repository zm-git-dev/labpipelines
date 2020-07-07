source $WZSEQ_ENTRY
wzref_mm10
pipeline_prepare

# cd /mnt/isilon/zhou_lab/projects/20200706_PGC_DNaseseq
while read srr_ids gsm sname; do
jump_comments

srr_ids=${srr_ids//,/ };
sname=$sname"_"$gsm
pipeline_dependlevel
memG=20; ppn=1;
pipeline_eval 1 __sra_fasterq_dump_SE_20200706;

done <<EOF
SRR6519360      GSM2966930      E9.5_PGC_DNase-seq
SRR6519361      GSM2966931      E9.5_PGC_DNase-seq
SRR6519362      GSM2966932      E10.5_PGC_DNase-seq
SRR6519363      GSM2966933      E10.5_PGC_DNase-seq
SRR6519364      GSM2966934      E12.5_female_PGC_DNase-seq
SRR6519365      GSM2966935      E12.5_female_PGC_DNase-seq
SRR6519366      GSM2966936      E12.5_female_soma_DNase-seq
SRR6519367      GSM2966937      E12.5_male_PGC_DNase-seq
SRR6519368      GSM2966938      E12.5_male_PGC_DNase-seq
SRR6519369      GSM2966939      E12.5_male_soma_DNase-seq
SRR6519370      GSM2966940      E13.5_female_PGC_DNase-seq
SRR6519371      GSM2966941      E13.5_female_PGC_DNase-seq
SRR6519372      GSM2966942      E13.5_male_PGC_DNase-seq
SRR6519373      GSM2966943      E13.5_male_PGC_DNase-seq
SRR6519374      GSM2966944      E14.5_female_PGC_DNase-seq
SRR6519375      GSM2966945      E14.5_female_PGC_DNase-seq
SRR6519376      GSM2966946      E14.5_male_PGC_DNase-seq
SRR6519377      GSM2966947      E14.5_male_PGC_DNase-seq
SRR6519378      GSM2966948      E16.5_female_PGC_DNase-seq
SRR6519379      GSM2966949      E16.5_female_PGC_DNase-seq
SRR6519380      GSM2966950      E16.5_male_PGC_DNase-seq
SRR6519381      GSM2966951      E16.5_male_PGC_DNase-seq
SRR6519382      GSM2966952      E16.5_male_soma_DNase-seq
SRR6519383      GSM2966953      E16.5_male_soma_DNase-seq
SRR6519384      GSM2966954      BVSC_ESC_DNase-seq
SRR6519385      GSM2966955      BVSC_ESC_DNase-seq
SRR6519386      GSM2966956      ESC_DNase-seq
SRR6519387      GSM2966957      ESC_DNase-seq
SRR6519388      GSM2966958      EpiLC_DNase-seq
SRR6519389      GSM2966959      EpiLC_DNase-seq
EOF
