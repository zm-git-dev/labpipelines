#########
## others
#########

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

## assume chrm, beg, end, beta, coverage
function wzseq_cov5 {
  f=$1
  zcat $f | awk '$5>=5' | gzip -c >${f%.bed.gz}.cov5.bed.gz
  echo `zcat ${f%.bed.gz}.cov5.bed.gz | wc -l` "CpGs covered 5X"
}

## assume chrm, beg, end, beta, coverage
function wzseq_cov10 {
  f=$1
  zcat $f | awk '$5>=10' | gzip -c >${f%.bed.gz}.cov10.bed.gz
  echo `zcat ${f%.bed.gz}.cov10.bed.gz | wc -l` "CpGs covered 10X"
}

function wzseq_liftbw {
  # liftOver bigwig file
  # Usage: wzseq_liftbw input.bigWig ~/tools/liftover/mm9ToMm10.over.chain.gz output.bigWig ~/references/mm10/mm10.fa.fai
  input=$1
  chain=$2
  output=$3
  chromsize=$4
  echo "[$(date)] Converting bigwig to bedgraph.."
  bigWigToBedGraph $input $input.bedg.tmp
  echo "[$(date)] Lifting over.."
  liftOver $input.bedg.tmp $chain $output.bedg.tmp $output.tmp.unmapped
  echo "  Mapped:   $(wc -l $output.bedg.tmp)"
  echo "  Unmapped: $(wc -l $output.tmp.unmapped)"
  echo "[$(date)] Sorting bedGraph and skip overlapping ..."
  sortbed $output.bedg.tmp | wzbedtools deoverlap -i - -o $output.bedg.tmp.sorted
  echo "  Before skipping: "$(wc -l $output.bedg.tmp)
  echo "  After skipping:  "$(wc -l $output.bedg.tmp.sorted)
  echo "[$(date)] Converting bedGraph to bigWig .."
  bedGraphToBigWig $output.bedg.tmp.sorted $chromsize $output
  echo "[$(date)] Cleaning"
  rm -f $input.bedg.tmp $output.tmp.unmapped $output.bedg.tmp $output.bedg.tmp.sorted
  echo "[$(date)] Done."
}

# fetch all SRR using SRX 
# wzseq_srx2srr SRX306253
function wzseq_srx2srr {
  curl "https://www.ncbi.nlm.nih.gov/sra?term=$1" | grep -o "SRR[[:digit:]]*" | sort | uniq | paste -d, -s
}

# wzseq_merge_fastq target1.fastq.gz source1.fq.gz,source2.fq.gz target2.fq.gz source3.fq.gz,source4.fq.gz ... <do>
function wzseq_merge_fastq {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[merge_fastq\]/{p=1;next} p&&!/^$/' samples |
    while read target sources; do
      sources=${sources//,/ }
      cmd="
cd $base/fastq
zcat $sources | pigz -p 4 -c > $target
"
      jobname="merge_fastq_${target//\//_}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 4
      [[ ${!#} == "do" ]] && qsub $pbsfn
    done
}

function wzseq_merge_bam {

  # "samtools merge -r" only add read group to each record, not to
  # the header, reheader helped insert read groups to header
  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[merge\]/{p=1;next} p&&!/^$/' samples | 
    while read merged sourcebams; do
      cmd="
cd $base
samtools merge -r bam/$merged.tmp.bam \$(echo $sourcebams | tr ',' ' ')
~/wzlib/pyutils/addRGtoSAMHeader.py -H -i bam/$merged.tmp.bam -o - | samtools reheader - bam/$merged.tmp.bam >bam/$merged.bam
rm -f bam/$merged.tmp.bam
samtools index bam/$merged.bam
samtools flagstat bam/$merged.bam >bam/$merged.bam.flagstat
"
      jobname="bam_merge_$merged"
      pbsfn=pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function __wzseq_index_bam() {
  : '
bam=bam/${sname}.bam
hour=5; memG=5; ppn=1
pipeline_eval x __wzseq_index_bam
'
  cmd='
samtools index '$bam';
samtools flagstat '$bam' > '$bam'.flagstat
mkdir -p multiqc/raw/flagstats/
ln -sf `readlink -f '$bam'.flagstat` multiqc/raw/flagstats/
'
  jobname="bamindex_"${bam//\//_}
}

# basic summary of the folder
function wzseq_basic_summary() {
  # for RNAseq
  for f in bam/*.bam; do b=${f%.bam}; b=${b#bam/}; sec=$(cat ${f%.bam}/accepted_hits.bam.flagstat | grep 'secondary' | awk '{match($1, /([0-9]*) \+/, a); print a[1];}'); tot=$(cat ${f%.bam}/accepted_hits.bam.flagstat | grep 'total (' | awk '{match($1, /([0-9]*) \+/, a); print a[1];}'); echo ${b%.bam} $(($tot-$sec)); done
  echo -e "\nraw read counts (paired-end, single-end)"
  for f in fastq/*.fastq.gz; do c=$(zcat $f | lc); echo $f $(($c / 2)) $(($c / 4)); done
}

function __wzseq_bam_mapq() {
  cmd='
cd '$base'
mkdir -p multiqc/raw/mapq
samtools view -F 0x100 -f 0x4 '$input_bam' | wc -l | cat <(echo -ne "unmapped\t") - >'$input_bam'_mapq_table
samtools view -F 0x104 '$input_bam' | awk '\''{cnt[$5]+=1}END{for(mapq in cnt) {print mapq"\t"cnt[mapq];}}'\'' | sort -k1,1n >>'$input_bam'_mapq_table
mkdir -p multiqc/raw/mapq/
ln -fs `readlink -f '$input_bam'_mapq_table` multiqc/raw/mapq/
'
  jobname='mapqdist_'$sname
}

