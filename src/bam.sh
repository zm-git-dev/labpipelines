function __zlab_PicardMarkdup_20211125 {
  cmd='
  cd '$base'
  [[ -d bam/picard_mdup ]] || mkdir bam/picard_mdup
  java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar ~/zhoulab/labsoftware/picard/picard-2.23.2.jar MarkDuplicates\
       CREATE_INDEX=true ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT\
       METRICS_FILE=bam/picard_mdup/${input}.mdup.stats READ_NAME_REGEX=null\
       INPUT='${input}' OUTPUT=bam/picard_mdup/'${input}' TMP_DIR=tmp
  samtools flagstat bam/picard_mdup/'${input}' > bam/picard_mdup/'${input}'.flagstat'
  jobname='_PicardMarkdup_'$sname
}

function __zlab_PicardMarkdup_20200719 {
  cmd='
cd '$base'
[[ -d pbs ]] || mkdir pbs
outdir="filtered_bam"
mkdir -p ${outdir}
outBasename=${outdir}/'$sname'
map_thresh=20
if [[ ! -e ${outBasename}.markDup.bam ]]; then
  java -jar /mnt/isilon/zhoulab/labsoftware/picard/picard-2.23.2.jar MarkDuplicates\
  I='$bam' O=${outBasename}.markDup.bam METRICS_FILE=${outBasename}_dup_qc.txt\
  ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true REMOVE_DUPLICATES=false
fi
'
  jobname='_PicardMarkdup_'$sname
}

function __zlab_bam2fastq_20211125 {
  # reverse bam to fastq
  # wzseq_bam2fastq
  # unmark secondary (some bam have only secondary mapping but primary mapping missing)
  # samtools view -h in.bam | awk '!/^@/{$2=and($2,compl(0x100)); print $0}/^@/' | samtools view -bo tmp0.bam
  cmd='
set -xe
# group reads by read names, collate is faster than sort
cd '$base';
mkdir -p bam/collate;
mkdir -p fastq
i=1;

## clean existing fastq files
rm -f fastq/'$sname'.pe1.fq.gz
rm -f fastq/'$sname'.pe2.fq.gz
rm -f fastq/'$sname'.se.fq.gz
rm -f fastq/'$sname'.paired_nolabel.fq.gz
for sourcebam in '${sourcebams//,/ }'; do

  ## NOTE: I extracted only primary mapping
  samtools view -h $sourcebam | awk -F"\t" -v OFS="\t" '\''!/^@/{$2=and($2,compl(0x100)); print $0}/^@/'\'' > bam/collate/'$sname'_${i}_tmp1.sam
  samtools collate -u bam/collate/'$sname'_${i}_tmp1.sam bam/collate/'$sname'_${i}_tmp2
  ## NOTE: use the following to skip primary mapping extraction
  ## samtools collate -u $sourcebam collate/${sample}_\$i.bam;
  
  ## NOTE: sometimes I used -O too. But -O can cause seg-fault on some malformed bams
  samtools fastq -n -0 fastq/'$sname'_$i.paired_nolabel.fq -1 fastq/'$sname'_$i.pe1.fq -2 fastq/'$sname'_$i.pe2.fq -s fastq/'$sname'_$i.se.fq bam/collate/'$sname'_${i}_tmp2.bam;
  gzip -c fastq/'$sname'_$i.paired_nolabel.fq >>fastq/'$sname'.paired_nolabel.fq.gz
  gzip -c fastq/'$sname'_$i.pe1.fq >>fastq/'$sname'.pe1.fq.gz
  gzip -c fastq/'$sname'_$i.pe2.fq >>fastq/'$sname'.pe2.fq.gz
  gzip -c fastq/'$sname'_$i.se.fq >>fastq/'$sname'.se.fq.gz 
  rm -f fastq/'$sname'_$i.paired_nolabel.fq fastq/'$sname'_$i.pe1.fq fastq/'$sname'_$i.pe2.fq fastq/'$sname'_$i.se.fq
  rm -f bam/collate/'$sname'_$i_tmp{1,2}*;
  i=$((i+1))
done
'
  jobname='bam2fastq_'$sname
}

function __zlab_indexBAM_20211125 {
  cmd='
samtools index '$bam';
samtools flagstat '$bam' > '$bam'.flagstat
mkdir -p multiqc/raw/flagstats/
ln -sf `readlink -f '$bam'.flagstat` multiqc/raw/flagstats/
'
  jobname="bamindex_"${bam//\//_}
}

function _zlab_downsampleBAM_20211125 {
  [[ -z ${n_replicates+x} ]] && n_replicates=50;
  [[ -z ${n_reads+x} ]] && n_reads=100
  cmd='
cd '$base'
mkdir subBAM
frac=$(samtools idxstats ${input_bam}| awk '\''{n+=$3} END {f='${n_reads}'/n; if(f>1) {f=1} print f;}'\'' | sed '\''s/^0*//'\'')
  
echo "Downsampling based on fraction: '${n_reads}', ${frac}"
for j in {1..'$n_replicates'}; do
  samtools view -b -@ $ppn -s ${j}${frac} ${input_bam} >subBAM/${sname}_down_${n_reads}_${j}.bam
  echo "$n_reads, $frac"
done
'
}

# function bam2bigwig() {
#   # need deepTools?
#   cmd='
# base=$("pwd")
# cd ${base}
# mkdir -p bigwig
# if [[ '$species' == "hs" ]]; then
# 	gsize=2790000000
# elif [[ '$species' == "mm" ]]; then
# 	gsize=1870000000
# else
# 	gsize=2000000000
# fi
# if [ ! -f filtered_bam/'$sname'.filtered.bam.bai ]; then
#   samtools index filtered_bam/'$sname'.filtered.bam
# fi
# bamCoverage -b filtered_bam/'$sname'.filtered.bam -o bigwig/'$sname'.bw -of bigwig -p '$ppn' -bs 20 --effectiveGenomeSize ${gsize} --normalizeUsing RPKM --minMappingQuality 20 -bl /home/dingw1/Mytools/data/mm10-blacklist.v2.bed
# #--ignoreDuplicates
# '
# 	jobname="bam2bigwig_"$sname
# }

# ## assume chrm, beg, end, beta, coverage
# function wzseq_cov5 {
#   f=$1
#   zcat $f | awk '$5>=5' | gzip -c >${f%.bed.gz}.cov5.bed.gz
#   echo `zcat ${f%.bed.gz}.cov5.bed.gz | wc -l` "CpGs covered 5X"
# }

# ## assume chrm, beg, end, beta, coverage
# function wzseq_cov10 {
#   f=$1
#   zcat $f | awk '$5>=10' | gzip -c >${f%.bed.gz}.cov10.bed.gz
#   echo `zcat ${f%.bed.gz}.cov10.bed.gz | wc -l` "CpGs covered 10X"
# }

# function wzseq_liftbw {
#   # liftOver bigwig file
#   # Usage: wzseq_liftbw input.bigWig ~/tools/liftover/mm9ToMm10.over.chain.gz output.bigWig ~/references/mm10/mm10.fa.fai
#   input=$1
#   chain=$2
#   output=$3
#   chromsize=$4
#   echo "[$(date)] Converting bigwig to bedgraph.."
#   bigWigToBedGraph $input $input.bedg.tmp
#   echo "[$(date)] Lifting over.."
#   liftOver $input.bedg.tmp $chain $output.bedg.tmp $output.tmp.unmapped
#   echo "  Mapped:   $(wc -l $output.bedg.tmp)"
#   echo "  Unmapped: $(wc -l $output.tmp.unmapped)"
#   echo "[$(date)] Sorting bedGraph and skip overlapping ..."
#   sortbed $output.bedg.tmp | wzbedtools deoverlap -i - -o $output.bedg.tmp.sorted
#   echo "  Before skipping: "$(wc -l $output.bedg.tmp)
#   echo "  After skipping:  "$(wc -l $output.bedg.tmp.sorted)
#   echo "[$(date)] Converting bedGraph to bigWig .."
#   bedGraphToBigWig $output.bedg.tmp.sorted $chromsize $output
#   echo "[$(date)] Cleaning"
#   rm -f $input.bedg.tmp $output.tmp.unmapped $output.bedg.tmp $output.bedg.tmp.sorted
#   echo "[$(date)] Done."
# }

# # fetch all SRR using SRX 
# # wzseq_srx2srr SRX306253
# function wzseq_srx2srr {
#   curl "https://www.ncbi.nlm.nih.gov/sra?term=$1" | grep -o "SRR[[:digit:]]*" | sort | uniq | paste -d, -s
# }

# # wzseq_merge_fastq target1.fastq.gz source1.fq.gz,source2.fq.gz target2.fq.gz source3.fq.gz,source4.fq.gz ... <do>
# function wzseq_merge_fastq {
#   base=$(pwd)
#   [[ -d pbs ]] || mkdir pbs
#   awk '/^\[/{p=0}/\[merge_fastq\]/{p=1;next} p&&!/^$/' samples |
#     while read target sources; do
#       sources=${sources//,/ }
#       cmd="
# cd $base/fastq
# zcat $sources | pigz -p 4 -c > $target
# "
#       jobname="merge_fastq_${target//\//_}"
#       pbsfn=$base/pbs/$jobname.pbs
#       pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 4
#       [[ ${!#} == "do" ]] && qsub $pbsfn
#     done
# }

# function wzseq_merge_bam {

#   # "samtools merge -r" only add read group to each record, not to
#   # the header, reheader helped insert read groups to header
#   base=$(pwd)
#   [[ -d bam ]] || mkdir bam
#   [[ -d pbs ]] || mkdir pbs
#   awk '/^\[/{p=0}/\[merge\]/{p=1;next} p&&!/^$/' samples | 
#     while read merged sourcebams; do
#       cmd="
# cd $base
# samtools merge -r bam/$merged.tmp.bam \$(echo $sourcebams | tr ',' ' ')
# ~/wzlib/pyutils/addRGtoSAMHeader.py -H -i bam/$merged.tmp.bam -o - | samtools reheader - bam/$merged.tmp.bam >bam/$merged.bam
# rm -f bam/$merged.tmp.bam
# samtools index bam/$merged.bam
# samtools flagstat bam/$merged.bam >bam/$merged.bam.flagstat
# "
#       jobname="bam_merge_$merged"
#       pbsfn=pbs/$jobname.pbs
#       pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
#       [[ $1 == "do" ]] && qsub $pbsfn
#     done
# }

# function __wzseq_index_bam() {
#   # Example:
#   # bam=bam/${sname}.bam
#   # hour=5; memG=5; ppn=1
#   # pipeline_eval x __wzseq_index_bam

#   cmd='
# samtools index '$bam';
# samtools flagstat '$bam' > '$bam'.flagstat
# mkdir -p multiqc/raw/flagstats/
# ln -sf `readlink -f '$bam'.flagstat` multiqc/raw/flagstats/
# '
#   jobname="bamindex_"${bam//\//_}
# }

# # basic summary of the folder
# function wzseq_basic_summary() {
#   # for RNAseq
#   for f in bam/*.bam; do b=${f%.bam}; b=${b#bam/}; sec=$(cat ${f%.bam}/accepted_hits.bam.flagstat | grep 'secondary' | awk '{match($1, /([0-9]*) \+/, a); print a[1];}'); tot=$(cat ${f%.bam}/accepted_hits.bam.flagstat | grep 'total (' | awk '{match($1, /([0-9]*) \+/, a); print a[1];}'); echo ${b%.bam} $(($tot-$sec)); done
#   echo -e "\nraw read counts (paired-end, single-end)"
#   for f in fastq/*.fastq.gz; do c=$(zcat $f | lc); echo $f $(($c / 2)) $(($c / 4)); done
# }

# function __wzseq_bam_mapq() {
#   cmd='
# cd '$base'
# mkdir -p multiqc/raw/mapq
# samtools view -F 0x100 -f 0x4 '$input_bam' | wc -l | cat <(echo -ne "unmapped\t") - >'$input_bam'_mapq_table
# samtools view -F 0x104 '$input_bam' | awk '\''{cnt[$5]+=1}END{for(mapq in cnt) {print mapq"\t"cnt[mapq];}}'\'' | sort -k1,1n >>'$input_bam'_mapq_table
# mkdir -p multiqc/raw/mapq/
# ln -fs `readlink -f '$input_bam'_mapq_table` multiqc/raw/mapq/
# '
#   jobname='mapqdist_'$sname
# }

# function __filterMAPQ_20200719 {
#   # filter reads given MAPQ threshold
#   cmd='
# cd '$base'
# if [[ ! -e '${outBasename}'.filtered.bam ]]; then
#   #remove reads with QCFAIL (flag 512)
#   samtools view -F 512 -q '${map_thresh}' -b '${outBasename}'.markDup.bam > '${outBasename}'.filtered.bam
#   samtools flagstat '${outBasename}'.filtered.bam > '${outBasename}'.filtered.bam.flagstat
#   samtools stats '${outBasename}'.filtered.bam > '${outBasename}'.filtered.bam.stats
# fi 
# '
# }

# function wzseq_picard_index_fasta {
#   dictfn=${$WZSEQ_REFERENCE%.fa}
#   cmd="
# java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar CreateSequenceDictionary REFERENCE=$WZSEQ_REFERENCE OUTPUT=${WZSEQ_REFERENCE}.dict
# "
#   jobname="picard_index_reference"
#   pbsfn=~/pbs/pbs/$jobname.pbs
#   pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
#   [[ $1 == "do" ]] && qsub $pbsfn
# }

# function wzseq_picard_WGSmetrics {
#   # right now, I keep getting this error
#   # Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: 14945
#   base=$(pwd);
#   [[ -d qc ]] || mkdir qc
#   for f in bam/*.bam; do
#     bfn=$(basename $f .bam)
#     cmd="
# cd $base
# java -Xmx5g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar CollectWgsMetrics INPUT=$fn OUTPUT=qc/$bfn.wgsmetrics REFERENCE_SEQUENCE=$WZSEQ_REFERENCE
# "
#     jobname="picard_WGSmetrics_$bfn"
#     pbsfn=$base/pbs/$jobname.pbs
#     pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
#     [[ ${!#} == "do" ]] && qsub $pbsfn
#   done
# }

# function wzseq_picard_markdup {

#   base=$(pwd)
#   [[ -d pbs ]] || mkdir pbs
#   [[ -d bam/before_mdup ]] || mkdir -p bam/before_mdup
#   [[ -d tmp ]] || mkdir tmp
#   awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
#     while read bfn sread1 sread2; do
#       f=bam/$bfn.bam
#       cmd="
# cd $base
# if [[ ! -e bam/before_mdup/$bfn.bam ]]; then
#   mv $f bam/before_mdup/$bfn.bam
#   [[ -e $f.bai ]] && mv $f.bai bam/before_mdup/$bfn.bam.bai
#   [[ -e $f.flagstat ]] && mv $f.flagstat bam/before_mdup/$bfn.bam.flagstat
# fi
# [[ -d bam/picard_mdup ]] || mkdir bam/picard_mdup
# java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=bam/picard_mdup/$bfn.mdup.stats READ_NAME_REGEX=null INPUT=bam/before_mdup/$bfn.bam OUTPUT=bam/picard_mdup/$bfn.bam TMP_DIR=tmp
# samtools flagstat bam/picard_mdup/$bfn.bam >bam/picard_mdup/$bfn.bam.flagstat

# cd bam; 
# [[ -h $bfn.bam ]] || [[ ! -e $bfn.bam ]] && ln -sf picard_mdup/$bfn.bam .
# [[ -h $bfn.bam.bai ]] || [[ ! -e $bfn.bam.bai ]] && ln -sf picard_mdup/$bfn.bai $bfn.bam.bai
# [[ -h $bfn.bam.flagstat ]] || [[ ! -e $bfn.bam.flagstat ]] && ln -sf picard_mdup/$bfn.bam.flagstat .
# <F11>
# " # other options: REMOVE_DUPLICATES=true
#     jobname="picard_markdup_$bfn"
#     pbsfn=$base/pbs/$jobname.pbs
#     pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 60 -memG 50 -ppn 5
#     [[ $1 == "do" ]] && qsub $pbsfn
#   done
# }

# function wzseq_merge_bam_picard {

#   # this properly handles read groups, though not sure how useful read groups are...
#   # this add RG to original source bam and output temporary files, merge works on
#   # those temporary files
#   base=$(pwd)
#   [[ -d bam ]] || mkdir bam
#   [[ -d pbs ]] || mkdir pbs
#   awk '/^\[/{p=0}/\[merge\]/{p=1;next} p&&!/^$/' samples |
#     while read merged sourcebams; do
#       cmd="
# cd $base
# rm -rf bam/tmp
# mkdir -p bam/tmp
# sourceid=0
# mergeinput=\"\"
# for sourcebam in \$(echo $sourcebams | tr ',' ' '); do
#   sourceid=\$((\$sourceid + 1))
#   rgname=\$(basename \$sourcebam .bam)
#   java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar AddOrReplaceReadGroups I=\$sourcebam O=bam/tmp/\${sourceid}.bam RGID=\$rgname RGLB=NA RGPL=illumina RGPU=NA RGSM=\$rgname
#   mergeinput=\$mergeinput\" I=bam/tmp/\${sourceid}.bam\"
# done

# java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar MergeSamFiles \$mergeinput O=bam/$merged.bam ASSUME_SORTED=true CREATE_INDEX=true
# samtools flagstat bam/$merged.bam >bam/$merged.bam.flagstat
# rm -rf bam/tmp
# "
#       jobname="picard_bam_merge_$merged"
#       pbsfn=pbs/$jobname.pbs
#       pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
#       [[ $1 == "do" ]] && qsub $pbsfn
#     done
# }

# function __markDupAndFilter_20200719 {
#   cmd='
# cd '$base'
# [[ -d pbs ]] || mkdir pbs
# outdir="filtered_bam"
# mkdir -p ${outdir}
# outBasename=${outdir}/'$sname'
# map_thresh=20
# if [[ ! -e ${outBasename}.markDup.bam ]]; then
#   java -jar /mnt/isilon/zhoulab/labsoftware/picard/picard-2.23.2.jar MarkDuplicates\
#   I='$bam' O=${outBasename}.markDup.bam METRICS_FILE=${outBasename}_dup_qc.txt\
#   ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true REMOVE_DUPLICATES=false
# fi

# #filter bam with given mapQ
# if [[ ! -e ${outBasename}.filtered.bam ]]; then
#   #remove reads with QCFAIL flag (512)
#   samtools view -F 512 -q ${map_thresh} -b ${outBasename}.markDup.bam > ${outBasename}.filtered.bam
#   samtools flagstat ${outBasename}.filtered.bam > ${outBasename}.filtered.bam.flagstat
#   samtools stats ${outBasename}.filtered.bam > ${outBasename}.filtered.bam.stats
# fi
# '
#   jobname="markdupAndFilter_${sname}"
# }

# function __wzseq_bam2fastq {

#   # reverse bam to fastq
#   #
#   # wzseq_bam2fastq
#   # unmark secondary (some bam have only secondary mapping but primary mapping missing)
#   # samtools view -h in.bam | awk '!/^@/{$2=and($2,compl(0x100)); print $0}/^@/' | samtools view -bo tmp0.bam
  
#   cmd='
# set -xe
# # group reads by read names, collate is faster than sort
# cd '$base';
# mkdir -p bam/collate;
# mkdir -p fastq
# i=1;

# ## clean existing fastq files
# rm -f fastq/'$sname'.pe1.fq.gz
# rm -f fastq/'$sname'.pe2.fq.gz
# rm -f fastq/'$sname'.se.fq.gz
# rm -f fastq/'$sname'.paired_nolabel.fq.gz
# for sourcebam in '${sourcebams//,/ }'; do

#   ## NOTE: I extracted only primary mapping
#   samtools view -h $sourcebam | awk -F"\t" -v OFS="\t" '\''!/^@/{$2=and($2,compl(0x100)); print $0}/^@/'\'' > bam/collate/'$sname'_${i}_tmp1.sam
#   samtools collate -u bam/collate/'$sname'_${i}_tmp1.sam bam/collate/'$sname'_${i}_tmp2
#   ## NOTE: use the following to skip primary mapping extraction
#   ## samtools collate -u $sourcebam collate/${sample}_\$i.bam;
  
#   ## NOTE: sometimes I used -O too. But -O can cause seg-fault on some malformed bams
#   samtools fastq -n -0 fastq/'$sname'_$i.paired_nolabel.fq -1 fastq/'$sname'_$i.pe1.fq -2 fastq/'$sname'_$i.pe2.fq -s fastq/'$sname'_$i.se.fq bam/collate/'$sname'_${i}_tmp2.bam;
#   gzip -c fastq/'$sname'_$i.paired_nolabel.fq >>fastq/'$sname'.paired_nolabel.fq.gz
#   gzip -c fastq/'$sname'_$i.pe1.fq >>fastq/'$sname'.pe1.fq.gz
#   gzip -c fastq/'$sname'_$i.pe2.fq >>fastq/'$sname'.pe2.fq.gz
#   gzip -c fastq/'$sname'_$i.se.fq >>fastq/'$sname'.se.fq.gz 
#   rm -f fastq/'$sname'_$i.paired_nolabel.fq fastq/'$sname'_$i.pe1.fq fastq/'$sname'_$i.pe2.fq fastq/'$sname'_$i.se.fq
#   rm -f bam/collate/'$sname'_$i_tmp{1,2}*;
#   i=$((i+1))
# done
# '
#   jobname='bam2fastq_'$sname
# }

# remove suffix _1, _2, mark read info to flag and
#  samtools view -h bam/SRR1029055.bam chr19 | awk '!/^@/{inpair=substr($1,length($1),1);$1=substr($1,1,length($1)-2);if(inpair==1) {$2=or($2,0x40);} else {$2=or($2,0x80);} $2=or($2,0x1); if (!(and($2, 0x100))) print $0}/^@/' | samtools collate -uO - SRR1029055tmp >SRR1029055.collate.bam
#  samtools view -h SRR1029055.collate.bam | awk 'BEGIN{key=""; line=""}!/^@/{if (key==$1) {print line; print $0;} key=$1; line=$0;}/^@/' | samtools fastq - -1 fastq_chr19/read1.fastq -2 fastq_chr19/read2.fastq -0 fastq_chr19/unpaired.fastq
# function wzseq_bam2fastq {
#   base=$(pwd);
#   [[ -d bam ]] || mkdir bam
#   [[ -d pbs ]] || mkdir pbs
#   awk '/^\[/{p=0}/\[bam2fastq\]/{p=1;next} p&&!/^$/' samples |
#     while read sample sourcebams; do
#       cmd="
# set -xe
# # group reads by read names, collate is faster than sort
# cd $base;
# mkdir -p bam/collate;
# i=1;
# for sourcebam in ${sourcebams//,/ }; do
#   samtools view -h \$sourcebam | awk -F\"\\t\" -v OFS=\"\t\" '!/^@/{\$2=and(\$2,compl(0x100)); print \$0}/^@/' > bam/collate/${sample}_\${i}_tmp1.sam
#   samtools collate -u bam/collate/${sample}_\${i}_tmp1.sam bam/collate/${sample}_\${i}_tmp2
#   # samtools collate -u \$sourcebam collate/${sample}_\$i.bam;
#   samtools fastq -On -0 fastq/${sample}_\$i.paired_nolabel.fq -1 fastq/${sample}_\$i.pe1.fq -2 fastq/${sample}_\$i.pe2.fq -s fastq/${sample}_\$i.se.fq bam/collate/${sample}_\${i}_tmp2.bam;
#   gzip -c fastq/${sample}_\$i.paired_nolabel.fq >>fastq/${sample}.paired_nolabel.fq.gz
#   gzip -c fastq/${sample}_\$i.pe1.fq >>fastq/${sample}.pe1.fq.gz
#   gzip -c fastq/${sample}_\$i.pe2.fq >>fastq/${sample}.pe2.fq.gz
#   gzip -c fastq/${sample}_\$i.se.fq >>fastq/${sample}.se.fq.gz 
#   rm -f fastq/${sample}_\$i.paired_nolabel.fq fastq/${sample}_\$i.pe1.fq fastq/${sample}_\$i.pe2.fq fastq/${sample}_\$i.se.fq
#   rm -f bam/collate/${sample}_\${i}_tmp{1,2}*;
#   i=\$((i+1))
# done
# "
# ## the following cleaning creates piping error, better check the size of output fastq manually
# # [[ -z \$(gunzip -c fastq/${sample}.paired_nolabel.fq.gz | head -c1) ]] && rm -f fastq/${sample}.paired_nolabel.fq.gz
# # [[ -z \$(gunzip -c fastq/${sample}.pe1.fq.gz | head -c1) ]] && rm -f fastq/${sample}.pe1.fq.gz
# # [[ -z \$(gunzip -c fastq/${sample}.pe2.fq.gz | head -c1) ]] && rm -f fastq/${sample}.pe2.fq.gz
# # [[ -z \$(gunzip -c fastq/${sample}.se.fq.gz | head -c1) ]] && rm -f fastq/${sample}.se.fq.gz
#     jobname="bam2fastq_$sample"
#     pbsfn=$base/pbs/$jobname.pbs
#     pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 5 -ppn 1
#     [[ ${!#} == "do" ]] && qsub $pbsfn
#   done
# }


# function __wzseq_uniformity_1M() {
#   cmd='
#   bedtools makewindows -w 1000000 -g '${WZSEQ_REFERENCE}'.fai | grep -v random | grep -v chrUn | grep -v hap | sortbed | bedtools coverage -a - -b <(samtools view -O BAM -q 40 '$input_bam') -sorted >uniformity/'${sname}'_1Mb.bed
# '
#   jobname='uniformity_1m_'$sname
# }

# # create coverage track, unique and nonunique mapping
# function wzseq_bam_coverage {

#   base=$(pwd);
#   [[ -d tracks ]] || mkdir tracks
#   [[ -d pbs ]] || mkdir pbs
#   [[ -d qc ]] || mkdir qc
#   for f in bam/*.bam; do
#     fn=$(readlink -f $f)
#     __wzseq_bam_coverage
#     [[ ${!#} == "do" ]] && qsub $pbsfn
#   done
# }

# function wgbs_cpgcoverage_OBSOLETE {

#   base=$(pwd)
#   [[ -d pbs ]] || mkdir pbs
#   [[ -d cpg ]] || mkdir cpg
#   for f in pileup/*.vcf.gz; do
#     bfn=$(basename $f .vcf.gz)
#     cmd="
# cd $base

# ## coverage
# ~/tools/biscuit/master/biscuit/biscuit vcf2bed -t cg -k 0 -c $f > cpg/${bfn}_cg.bedg
# bedtools intersect -a $WZSEQ_CPGBED -b cpg/${bfn}_cg.bedg -sorted -loj | awk '{if(\$11==\".\")\$11=0;print \$1,\$2,\$3,\$11}' | bedtools groupby -g 1-3 -c 4 -o sum >cpg/cpgCoverage_${bfn}.bedg
# wzplot cumhist -t cpg/cpgCoverage_${bfn}.bedg -c 4 --xlabel \"Coverage\" --ylabel \"Cumulative Distribution\" -o cpg/cpgCoverage_wholeGenome_${bfn}.png
# bedtools intersect -a cpg/cpgCoverage_${bfn}.bedg -b $WZSEQ_CGIBED | wzplot cumhist -t - -c 4 --xlabel \"Coverage\" --ylabel \"Cumulative Distribution\" -o cpg/cpgCoverage_CGI_${bfn}.png

# ## beta value distribution
# wzplot hist --maxline 1000000000 -t cpg/${bfn}_cg.bedg -c 4 --xlabel \"Beta Values\" -o cpg/betaDist_wholeGenome_${bfn}.png
# bedtools intersect -a cpg/${bfn}_cg.bedg -b $WZSEQ_CGIBED | wzplot hist --maxline 1000000000 -t - -c 4 --xlabel \"Beta Values\" -o cpg/betaDist_CGI_${bfn}.png

# ## Methylation average in CGI
# awk '\$5>5' cpg/${bfn}_cg.bedg | bedtools intersect -a $WZSEQ_CGIBED -b - -wo -sorted | bedtools groupby -i - -g 1-5 -c 14,14 -o mean,count > cpg/CGImethAverage_cov5_${bfn}.bed
# wzplot hist --maxline 1000000000 -t cpg/CGImethAverage_cov5_${bfn}.bed -c 6 --xlabel \"CGImethAverage\" -o cpg/CGImethAverage_cov5_${bfn}.png

# # bedtools closest -a $WZSEQ_TSSBED -b cpg/$bfn.cgi.bed -d >cpg/CGImethAverage_cov5_${bfn}_tss.bed
# # rm -f cpg/$bfn.cg.bedg
# "
#     jobname="cpgCoverage_$bfn"
#     pbsfn=$base/pbs/$jobname.pbs
#     pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
#     [[ ${!#} == "do" ]] && qsub $pbsfn
#   done
# }


