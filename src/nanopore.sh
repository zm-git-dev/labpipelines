#!/bin/bash

function tomboResquiggle() {
  guppy_dir="/mnt/isilon/zhoulab/labsoftware/guppy/ont-guppy-cpu"
  nanopolish_dir="/mnt/isilon/zhoulab/labsoftware/nanopolish/nanopolish-0.13.2"
  minimap_dir="/mnt/isilon/zhoulab/labsoftware/minimap2"
  cfg=${guppy_dir}/data/dna_r9.4.1_450bps_hac.cfg #${guppy_dir}/data/dna_r9.4.1_450bps_fast.cfg
  multi_to_single_fast5="/mnt/isilon/zhoulab/labsoftware/anaconda/anaconda3_2020/bin/multi_to_single_fast5"
  #cfg=${guppy_dir}/data/dna_r9.4.1_450bps_hac.cfg
  
  cmd='
module load gcc8/8.4.0
base=$("pwd")
cd '${base}'
mkdir -p Tombo
#Note: tombo only accept single_fast5 as input.
'${multi_to_single_fast5}' --recursive -t 1 -i BaseCalling/fast5_test/workspace -s Tombo/'${sname}'
#tombo preprocess annotate_raw_with_fastqs --overwrite --processes '${ppn}' --fast5-basedir Tombo/'${sname}' --fastq-filenames '${fastq}'
tombo resquiggle Tombo/'${sname}' '${WZSEQ_REFERENCE}' --processes '${ppn}' --num-most-common-errors 5 --overwrite --ignore-read-locks
'
  jobname="resquiggle_"$sname
}


guppy_dir="/mnt/isilon/zhoulab/labsoftware/guppy/ont-guppy-cpu"
guppy_gpu_dir="/mnt/isilon/zhoulab/labsoftware/guppy/ont-guppy-gpu"
nanopolish_dir="/mnt/isilon/zhoulab/labsoftware/nanopolish/nanopolish"
minimap_dir="/mnt/isilon/zhoulab/labsoftware/minimap2"
NanoMod_dir="/mnt/isilon/zhoulab/labsoftware/NanoMod/NanoMod/bin"
export LD_LIBRARY_PATH=/mnt/isilon/zhoulab/labsoftware/shared_Renv/versions/4.0.4/lib64/R/lib/:${LD_LIBRARY_PATH}
multi_to_single_fast5="/mnt/isilon/zhoulab/labsoftware/anaconda/anaconda3_2020/bin/multi_to_single_fast5"
#cfg=${guppy_dir}/data/dna_r9.4.1_450bps_hac.cfg #${guppy_dir}/data/dna_r9.4.1_450bps_fast.cfg

function GuppyBaseCalling_20201023() {
  #This function will take ${fast5} as input and generate BaseCalling/'${sname}' in the current work directory.
  #Needed variates: base, ppn, fast5, sname.
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p BaseCalling
'${guppy_dir}'/bin/guppy_basecaller -c '${guppy_dir}'/data/dna_r9.4.1_450bps_hac.cfg --num_callers '${ppn}' -r --compress_fastq -i '${fast5}' -s BaseCalling/'${sname}'
#--chunks_per_runner 768 --chunk_size 500 --fast5_out
zcat BaseCalling/'${sname}'/*.fastq.gz > BaseCalling/'${sname}'.fastq
'${nanopolish_dir}'/nanopolish index -d '${fast5}' BaseCalling/'${sname}'.fastq
'
  jobname="guppy_basecaller_"$sname
}

function GuppyBaseCalling_fast_20210219() {
  #This function will take ${fast5} as input and generate BaseCalling/'${sname}' in the current work directory.
  #Needed variates: base, ppn, fast5, sname.
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p BaseCalling
'${guppy_dir}'/bin/guppy_basecaller -c '${guppy_dir}'/data/dna_r9.4.1_450bps_fast.cfg --num_callers '${ppn}' -r --compress_fastq -i '${fast5}' -s BaseCalling/'${sname}'
#--chunks_per_runner 768 --chunk_size 500 
zcat BaseCalling/'${sname}'/*.fastq.gz > BaseCalling/'${sname}'.fastq
'${nanopolish_dir}'/nanopolish index -d '${fast5}' BaseCalling/'${sname}'.fastq
'
  jobname="guppy_basecaller_"$sname
}

function GuppyBaseCalling_rna_fast_20210415() {
#This function will take ${fast5} as input and generate BaseCalling/'${sname}' in the current work directory.
#Needed variates: base, ppn, fast5, sname.
cmd='
base=$("pwd")
cd '${base}'
mkdir -p BaseCalling
'${guppy_dir}'/bin/guppy_basecaller -c '${guppy_dir}'/data/rna_r9.4.1_70bps_fast.cfg --num_callers '${ppn}' -r --compress_fastq -i '${fast5}' -s BaseCalling/'${sname}'
#--chunks_per_runner 768 --chunk_size 500
zcat BaseCalling/'${sname}'/*.fastq.gz > BaseCalling/'${sname}'.fastq
'${nanopolish_dir}'/nanopolish index -d '${fast5}' BaseCalling/'${sname}'.fastq
'

  jobname="guppy_basecaller_"$sname
}

function NanopolishIndex20210406() {
cmd='
base=$("pwd")
cd '${base}'
if [ ! -f '${fastq}' ]; then
    gunzip '${fastq}'.gz
fi
'${nanopolish_dir}'/nanopolish index -d '${fast5}' '${fastq}'
'
  jobname="nanopolish_index_"$sname
}

function GuppyBaseCalling_hac_fast5_out_20210223() {
  #This function will take ${fast5} as input and generate BaseCalling/'${sname}' in the current work directory.
  #Needed variates: base, ppn, fast5, sname.
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p BaseCalling
'${guppy_dir}'/bin/guppy_basecaller -c '${guppy_dir}'/data/dna_r9.4.1_450bps_hac.cfg --num_callers '${ppn}' -r --compress_fastq --fast5_out -i '${fast5}' -s BaseCalling/'${sname}'
#output fast5 should be located at BaseCalling/sname/workspace
zcat BaseCalling/'${sname}'/*.fastq.gz > BaseCalling/'${sname}'.fastq
'${nanopolish_dir}'/nanopolish index -d '${fast5}' BaseCalling/'${sname}'.fastq
'
  jobname="guppy_basecaller_"$sname
}

function GuppyBaseCalling_fast_fast5_out_20210711() {
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p BaseCalling
'${guppy_dir}'/bin/guppy_basecaller -c '${guppy_dir}'/data/rna_r9.4.1_70bps_fast.cfg --num_callers '${ppn}' -r --compress_fastq --fast5_out -i '${fast5}' -s BaseCalling/'${sname}'
#output fast5 should be located at BaseCalling/sname/workspace
zcat BaseCalling/'${sname}'/*.fastq.gz > BaseCalling/'${sname}'.fastq
'${nanopolish_dir}'/nanopolish index -d '${fast5}' BaseCalling/'${sname}'.fastq
'
  jobname="guppy_basecaller_"$sname
}

function GuppyBaseCalling_hac_gpu_20210227() {
  #This function will take ${fast5} as input and generate BaseCalling/'${sname}' in the current work directory.
  #Needed variates: base, ppn, fast5, sname.
  #https://esr-nz.github.io/gpu_basecalling_testing/gpu_benchmarking.html
  #https://gist.github.com/disulfidebond/00ff5a6f84a0a81057c6e5817c540569
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p BaseCalling
'${guppy_gpu_dir}'/bin/guppy_basecaller -c '${guppy_gpu_dir}'/data/dna_r9.4.1_450bps_hac.cfg \
			-x "auto" --gpu_runners_per_device 8 --num_callers '${ppn}' \
			 --chunks_per_runner 512 --chunk_size 1000 \
			 -r --compress_fastq -i '${fast5}' -s BaseCalling/'${sname}'
#-x "cuda:0"
zcat BaseCalling/'${sname}'/*.fastq.gz > BaseCalling/'${sname}'.fastq
rm -rf BaseCalling/'${sname}'
'${nanopolish_dir}'/nanopolish index -d '${fast5}' BaseCalling/'${sname}'.fastq
'
  jobname="guppy_basecaller_"$sname
}

function Minimap2Aligning_20201023() {
  #This function will take BaseCalling/'${sname}'.fastq, '${WZSEQ_REFERENCE}' as input and generate Alignment/'${sname}'.sorted.bam.
  #Needed variates: base, fastq, ppn,sname and WZSEQ_REFERENCE.
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p Alignment
'${minimap_dir}'/minimap2 -a -x map-ont -t '${ppn}' '${WZSEQ_REFERENCE}' '${fastq}' | samtools sort -T Alignment/'${sname}' -o Alignment/'${sname}'.sorted.bam
samtools index -@ '${ppn}' Alignment/'${sname}'.sorted.bam
samtools flagstat Alignment/'${sname}'.sorted.bam > Alignment/'${sname}'.sorted.bam.flagstat
'
  jobname="minimap2_"$sname
}

function Minimap2Aligning_RNA_20210412() {
#This function will take BaseCalling/'${sname}'.fastq, '${WZSEQ_REFERENCE}' as input and generate Alignment/'${sname}'.sorted.bam.
#Needed variates: base, fastq, ppn,sname and WZSEQ_REFERENCE.
cmd='
base=$("pwd")
cd '${base}'
mkdir -p Alignment
'${minimap_dir}'/minimap2 -a -x splice -uf -k14 -t '${ppn}' '${WZSEQ_REFERENCE}' '${fastq}' | samtools sort -T Alignment/'${sname}' -o Alignment/'${sname}'.sorted.bam
samtools index -@ '${ppn}' Alignment/'${sname}'.sorted.bam
samtools flagstat Alignment/'${sname}'.sorted.bam > Alignment/'${sname}'.sorted.bam.flagstat
'
jobname="minimap2_"$sname
}


function Minimap2Aligning_for_timp() {
  #This function will take BaseCalling/'${sname}'.fastq, '${WZSEQ_REFERENCE}' as input and generate Alignment/'${sname}'.sorted.bam.
  #Needed variates: base, fastq, ppn,sname and WZSEQ_REFERENCE.
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p Alignment_timp
'${minimap_dir}'/minimap2 -d '${WZSEQ_REFERENCE}'.mmi '${WZSEQ_REFERENCE}'
'${minimap_dir}'/minimap2 --MD -a -x map-ont -t '${ppn}' '${WZSEQ_REFERENCE}'.mmi -f '${WZSEQ_REFERENCE}' '${fastq}' | samtools sort -T tmp -o Alignment_timp/'${sname}'.sorted.bam
samtools index -@ '${ppn}' Alignment_timp/'${sname}'.sorted.bam
samtools flagstat Alignment_timp/'${sname}'.sorted.bam > Alignment_timp/'${sname}'.sorted.bam.flagstat
'
  jobname="minimap2_"$sname
}

function NanopolishCallMeth_20201023() {
  #Runing nanopolish to call methylation.
  #Needed variates: base, fastq, bam, ppn, sname and WZSEQ_REFERENCE.
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p CallMethylation
if [ ! -f '${fastq}'.index ]; then
    echo "Run nanopolish index.."
    '${nanopolish_dir}'/nanopolish index -d '${fast5}' '${fastq}'
fi
'${nanopolish_dir}'/nanopolish call-methylation -t '${ppn}' -r '${fastq}' -b '${bam}' -g '${WZSEQ_REFERENCE}' > CallMethylation/'${sname}'.tsv
'${nanopolish_dir}'/scripts/calculate_methylation_frequency.py -s CallMethylation/'${sname}'.tsv > CallMethylation/'${sname}'.frequency.tsv
'
  jobname="nanopolish_"$sname
}

function NanopolishCallMeth_for_timp() {
  #Runing nanopolish to call methylation.
  #Needed variates: base, fastq, bam, ppn, sname and WZSEQ_REFERENCE.
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p CallMethylation_timp
'${nanopolish_dir}'/nanopolish call-methylation -t '${ppn}' -r '${fastq}' -b '${bam}' -g '${WZSEQ_REFERENCE}' > CallMethylation_timp/'${sname}'.tsv
'${nanopolish_dir}'/scripts/calculate_methylation_frequency.py -s CallMethylation_timp/'${sname}'.tsv > CallMethylation_timp/'${sname}'.frequency.tsv
'
  jobname="nanopolish_"$sname
}



function NanopolishEventsToReference() {
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p EventsToReference
if [ ! -f EventsToReference/'${sname}'.fasta ]; then
  seqtk seq -a '${fastq}' >  EventsToReference/'${sname}'.fasta
  #paste - - - - < '${fastq}' | cut -f 1,2 | sed '"'"'s/^@/>/'"'"' | tr ''"''\t''"'' ''"''\n''"'' > EventsToReference/'${sname}'.fasta
fi
'${nanopolish_dir}'/nanopolish index -d '${fast5}' EventsToReference/'${sname}'.fasta
#'${nanopolish_dir}'/nanopolish eventalign -n -t '${ppn}' --reads EventsToReference/'${sname}'.fasta --bam '${bam}' --genome '${WZSEQ_REFERENCE}' > EventsToReference/'${sname}'.eventalign.txt
#To make eventalign.txt smaller, get_minimum_sl.py is a custom python script to select one record with the smallest abs(standardized_level) for one position, it was put in the same directory with nanopolish.
'${nanopolish_dir}'/nanopolish eventalign -t '${ppn}' --reads EventsToReference/'${sname}'.fasta --bam '${bam}' --genome '${WZSEQ_REFERENCE}' | awk '"'"'BEGIN{OFS="\t"}; {if (NR>1) {print($1,$2,$2+6,$10,$7,$8,$9,$13)}}'"'"' | '${nanopolish_dir}'/get_minimum_sl.py | bgzip  > EventsToReference/'${sname}'.eventalign.txt.gz
#output columns: "chr","start","end","model_kmer","event_level_mean","event_stdv","event_length","standardized_level"
'
  jobname="eventalign_"$sname
}

function NanopolishEventsToReference_scale() {
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p EventsToReference_scale
if [ ! -f EventsToReference_scale/'${sname}'.fasta ]; then
  seqtk seq -a '${fastq}' >  EventsToReference_scale/'${sname}'.fasta
  #paste - - - - < '${fastq}' | cut -f 1,2 | sed '"'"'s/^@/>/'"'"' | tr ''"''\t''"'' ''"''\n''"'' > EventsToReference_scale/'${sname}'.fasta
fi
'${nanopolish_dir}'/nanopolish index -d '${fast5}' EventsToReference_scale/'${sname}'.fasta
#'${nanopolish_dir}'/nanopolish eventalign -n -t '${ppn}' --scale-events --reads EventsToReference_scale/'${sname}'.fasta --bam '${bam}' --genome '${WZSEQ_REFERENCE}' > EventsToReference_scale/'${sname}'.eventalign.txt
#To make eventalign.txt smaller, get_minimum_sl.py is a custom python script to select one record with the smallest abs(standardized_level) for one position, it was put in the same directory with nanopolish.
'${nanopolish_dir}'/nanopolish eventalign -n -t '${ppn}' --scale-events --reads EventsToReference_scale/'${sname}'.fasta --bam '${bam}' --genome '${WZSEQ_REFERENCE}' | awk '"'"'BEGIN{OFS="\t"}; {if (NR>1) {print($1,$2,$2+6,$10,$7,$8,$9,$13)}}'"'"' | '${nanopolish_dir}'/get_minimum_sl.py | bgzip  > EventsToReference_scale/'${sname}'.eventalign.txt.gz
#output columns: "chr","start","end","model_kmer","event_level_mean","event_stdv","event_length","standardized_level"
'
  jobname="eventalign_"$sname
}

function NanoModAnnotate_20210625() {
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p SingleFast5/'${sname}'
'${multi_to_single_fast5}' --recursive -t '${ppn}'  -i BaseCalling/'${sname}'/workspace -s SingleFast5/'${sname}'
#Input for NanoMod Annotate should be basecalled fast5.
'${NanoMod_dir}'/NanoMod.py Annotate --wrkBase1 SingleFast5/'${sname}' --Ref '${WZSEQ_REFERENCE}' --threads '${ppn}' 
# How to use parameter kmer_model_file
'
  jobname="NanoModAnnotate_"$sname
}

function NanoModDetect_20210630() {
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p Detect
'${NanoMod_dir}'/NanoMod.py detect --wrkBase1 SingleFast5/'${sname1}' --wrkBase2 SingleFast5/'${sname2}' --outFolder Detect
'
  jobname="NanoModDetect_"$sname
}
