#!/bin/bash

function macs2mergeNoControl() {
  cmd='
base=$("pwd")
cd ${base}
mkdir -p macs2Merged
macs2 callpeak -t filtered_bam/'$sname'*.filtered.bam -g mm -q 0.05 -n '$sname' --outdir macs2Merged -B
#ls filtered_bam/${arg1}*.filtered.bam
'
  jobname="macs2mergeNoControl_"$sname
}


function runMacs2() {
  cmd='
base=$("pwd")
cd ${base}
mkdir -p macs2
macs2 callpeak -t filtered_bam/'$sname'*.filtered.bam -g mm -q 0.05 -n '$sname' --outdir macs2 -B --broad
#ls filtered_bam/${arg1}*.filtered.bam
'
  jobname="runMacs2_"$sname
}

function runMacs2WithControl() {
  cmd='
base=$("pwd")
cd ${base}
mkdir -p macs2
macs2 callpeak -t filtered_bam/'$sname'.filtered.bam -c filtered_bam/'$control' -g mm -q 0.05 -n '$sname' --outdir macs2 -B --broad
'
  jobname="macs2WithControl_"$sname
}
