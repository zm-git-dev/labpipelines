#!/bin/bash

source ~/wzlib/bash/wzseq.sh
pipeline_prepare

while read sname pe_or_se srr_ids; do
  jump_comments
  srr_ids=${srr_ids//,/ };
  pipeline_dependlevel
  hour=48; memG=10; ppn=1
  if [[ $pe_or_se == "SINGLE" ]]; then
    pipeline_eval 1 __wzseq_fastq_dump_SE
  else
    pipeline_eval 1 __wzseq_fastq_dump_PE
  fi
done << EOM
GSM497246_Apis_mellifera_PAIRED	PAIRED	SRR042617,SRR042618,SRR042619,SRR042620
GSM497248_Bombyx_mori_PAIRED	PAIRED	SRR063324,SRR063326
GSM497248_Bombyx_mori_SINGLE	SINGLE	SRR042623,SRR042624,SRR042625
GSM497250_Chlorella_sp.NC64A_SINGLE	SINGLE	SRR042626
GSM497251_Ciona_intestinalis_PAIRED	PAIRED	SRR042627,SRR042628,SRR060811
GSM497253_Coprinopsis_cinerea_SINGLE	SINGLE	SRR042630
GSM497255_Drosophila_melanogaster_SINGLE	SINGLE	SRR042631,SRR060807
GSM497256_Laccaria_bicolor_SINGLE	SINGLE	SRR042632,SRR042633
GSM497258_Nematostella_vectensis_PAIRED	PAIRED	SRR042634,SRR042635,SRR042636,SRR060812,SRR063322
GSM497260_Oryza_sativa_PAIRED	PAIRED	SRR042638
GSM497260_Oryza_sativa_SINGLE	SINGLE	SRR042639,SRR042640,SRR042641,SRR042642
GSM497262_Phycomyces_blakesleeanus_SINGLE	SINGLE	SRR042643
GSM497264_Physcomitrella_pattens_PAIRED	PAIRED	SRR042644
GSM497264_Physcomitrella_pattens_SINGLE	SINGLE	SRR042645,SRR042646,SRR042647
GSM497266_Postia_placenta_SINGLE	SINGLE	SRR042648,SRR042649
GSM497268_Selaginella_moellendorffii_SINGLE	SINGLE	SRR042650,SRR042651
GSM497270_Tetraodon_nigroviridis_PAIRED	PAIRED	SRR042652,SRR042653
GSM497270_Tetraodon_nigroviridis_SINGLE	SINGLE	SRR042654,SRR042655
GSM497273_Tribolium_castaneum_PAIRED	PAIRED	SRR042656
GSM497274_Uncinocarpus_reesi_SINGLE	SINGLE	SRR042657
GSM497276_Volvox_carteri_PAIRED	PAIRED	SRR042658
GSM497276_Volvox_carteri_SINGLE	SINGLE	SRR042659,SRR042660
EOM

while read sname design species; do
  jump_comments
  pipeline_dependlevel

  WZSEQ_BISCUIT_INDEX=~/genomes/clean/$species/biscuit/$species.fa
  output_bam=bam/${sname}.bam
  hour=48; memG=150; ppn=20;
  if [[ $design == "PAIRED" ]]; then
    fastq1=fastq/${sname}_R1.fastq.gz
    fastq2=fastq/${sname}_R2.fastq.gz
    pipeline_eval 2 __wgbs_biscuit_align_PE_both
  else
    fastq=fastq/${sname}.fastq.gz
    pipeline_eval 2 __wgbs_biscuit_align_SE_both
  fi

  WZSEQ_REFERENCE=~/genomes/clean/$species/$species.fa
  input_bam=bam/${sname}.bam
  hour=24; memG=10; ppn=1; queue=longq
  pipeline_eval 5 __wgbs_biscuit_pileup

  input_bam=bam/${sname}.bam
  input_vcf=pileup/${sname}.vcf.gz
  hour=48; memG=20; ppn=2
  pipeline_eval 6 __wgbs_biscuit_QC
done << EOM
# GSM497246_Apis_mellifera_PAIRED PAIRED  apis_mellifera_Amel_4.5
GSM497248_Bombyx_mori_PAIRED    PAIRED  bombyx_mori_ASM15162v1
GSM497248_Bombyx_mori_SINGLE    SINGLE  bombyx_mori_ASM15162v1
# GSM497250_Chlorella_sp.NC64A_SINGLE     SINGLE  ChlNC64A_1
# GSM497251_Ciona_intestinalis_PAIRED     PAIRED  ciona_intestinalis_KH
# GSM497253_Coprinopsis_cinerea_SINGLE    SINGLE  Copci1
# GSM497255_Drosophila_melanogaster_SINGLE        SINGLE  drosophila_melanogaster_BDGP6
# GSM497256_Laccaria_bicolor_SINGLE       SINGLE  Lacbi2
# GSM497258_Nematostella_vectensis_PAIRED PAIRED  nematostella_vectensis_ASM20922v1
# GSM497260_Oryza_sativa_PAIRED   PAIRED  oryza_sativa_IRGSP-1.0
# GSM497260_Oryza_sativa_SINGLE   SINGLE  oryza_sativa_IRGSP-1.0
# GSM497262_Phycomyces_blakesleeanus_SINGLE       SINGLE  Phybl2
# GSM497264_Physcomitrella_pattens_PAIRED PAIRED  physcomitrella_patens_ASM242v1
# GSM497264_Physcomitrella_pattens_SINGLE SINGLE  physcomitrella_patens_ASM242v1
# GSM497266_Postia_placenta_SINGLE        SINGLE  Pospl1
# GSM497268_Selaginella_moellendorffii_SINGLE     SINGLE  selaginella_moellendorffii_v1.0
# GSM497270_Tetraodon_nigroviridis_PAIRED PAIRED  tetraodon_nigroviridis_TETRAODON8
# GSM497270_Tetraodon_nigroviridis_SINGLE SINGLE  tetraodon_nigroviridis_TETRAODON8
# GSM497273_Tribolium_castaneum_PAIRED    PAIRED  tribolium_castaneum_Tcas5.2
# GSM497274_Uncinocarpus_reesi_SINGLE     SINGLE  Uncre1
# GSM497276_Volvox_carteri_PAIRED PAIRED  Vcarteri_317_v2
# GSM497276_Volvox_carteri_SINGLE SINGLE  Vcarteri_317_v2
EOM
