#!/bin/bash
guppy_dir="/mnt/isilon/zhoulab/labsoftware/guppy/ont-guppy-cpu"
nanopolish_dir="/mnt/isilon/zhoulab/labsoftware/nanopolish/nanopolish-0.13.2"
minimap_dir="/mnt/isilon/zhoulab/labsoftware/minimap2"
cfg=${guppy_dir}/data/dna_r9.4.1_450bps_hac.cfg #${guppy_dir}/data/dna_r9.4.1_450bps_fast.cfg
multi_to_single_fast5="/mnt/isilon/zhoulab/labsoftware/anaconda/anaconda3_2020/bin/multi_to_single_fast5"
#cfg=${guppy_dir}/data/dna_r9.4.1_450bps_hac.cfg

function tomboResquiggle() {
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

