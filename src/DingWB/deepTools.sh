#!/bin/bash

function bam2bigwig()
{
  cmd='
base=$("pwd")
cd ${base}
mkdir -p bigwig

if [[ '$species' == "hs" ]]; then
	gsize=2790000000
elif [[ '$species' == "mm" ]]; then
	gsize=1870000000
else
	gsize=2000000000
fi

if [ ! -f filtered_bam/'$sname'.filtered.bam.bai ]; then
  samtools index filtered_bam/'$sname'.filtered.bam
fi

bamCoverage -b filtered_bam/'$sname'.filtered.bam -o bigwig/'$sname'.bw -of bigwig -p '$ppn' -bs 20 --effectiveGenomeSize ${gsize} --normalizeUsing RPKM --minMappingQuality 20 -bl /home/dingw1/Mytools/data/mm10-blacklist.v2.bed
#--ignoreDuplicates
'
	jobname="bam2bigwig_"$sname
}

#example: 