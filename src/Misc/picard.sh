############
## Picard
############

function wzseq_picard_index_fasta {
  dictfn=${$WZSEQ_REFERENCE%.fa}
  cmd="
java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar CreateSequenceDictionary REFERENCE=$WZSEQ_REFERENCE OUTPUT=${WZSEQ_REFERENCE}.dict
"
  jobname="picard_index_reference"
  pbsfn=~/pbs/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function wzseq_picard_WGSmetrics {
  # right now, I keep getting this error
  # Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: 14945
  base=$(pwd);
  [[ -d qc ]] || mkdir qc
  for f in bam/*.bam; do
    bfn=$(basename $f .bam)
    cmd="
cd $base
java -Xmx5g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar CollectWgsMetrics INPUT=$fn OUTPUT=qc/$bfn.wgsmetrics REFERENCE_SEQUENCE=$WZSEQ_REFERENCE
"
    jobname="picard_WGSmetrics_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wzseq_picard_markdup {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam/before_mdup ]] || mkdir -p bam/before_mdup
  [[ -d tmp ]] || mkdir tmp
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read bfn sread1 sread2; do
      f=bam/$bfn.bam
      cmd="
cd $base
if [[ ! -e bam/before_mdup/$bfn.bam ]]; then
  mv $f bam/before_mdup/$bfn.bam
  [[ -e $f.bai ]] && mv $f.bai bam/before_mdup/$bfn.bam.bai
  [[ -e $f.flagstat ]] && mv $f.flagstat bam/before_mdup/$bfn.bam.flagstat
fi
[[ -d bam/picard_mdup ]] || mkdir bam/picard_mdup
java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=bam/picard_mdup/$bfn.mdup.stats READ_NAME_REGEX=null INPUT=bam/before_mdup/$bfn.bam OUTPUT=bam/picard_mdup/$bfn.bam TMP_DIR=tmp
samtools flagstat bam/picard_mdup/$bfn.bam >bam/picard_mdup/$bfn.bam.flagstat

cd bam; 
[[ -h $bfn.bam ]] || [[ ! -e $bfn.bam ]] && ln -sf picard_mdup/$bfn.bam .
[[ -h $bfn.bam.bai ]] || [[ ! -e $bfn.bam.bai ]] && ln -sf picard_mdup/$bfn.bai $bfn.bam.bai
[[ -h $bfn.bam.flagstat ]] || [[ ! -e $bfn.bam.flagstat ]] && ln -sf picard_mdup/$bfn.bam.flagstat .
<F11>
" # other options: REMOVE_DUPLICATES=true
    jobname="picard_markdup_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 60 -memG 50 -ppn 5
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function wzseq_picard_markdup2 {

base=$(pwd)

cd $base
[[ -d bam/picard_mdup ]] || mkdir bam/picard_mdup
java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar ~/zhoulab/labsoftware/picard/picard-2.23.2.jar MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=bam/picard_mdup/${input}.mdup.stats READ_NAME_REGEX=null INPUT=${input} OUTPUT=bam/picard_mdup/${input} TMP_DIR=tmp
samtools flagstat bam/picard_mdup/${input} > bam/picard_mdup/${input}.flagstat
done
}


function wzseq_merge_bam_picard {

  # this properly handles read groups, though not sure how useful read groups are...
  # this add RG to original source bam and output temporary files, merge works on
  # those temporary files
  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[merge\]/{p=1;next} p&&!/^$/' samples |
    while read merged sourcebams; do
      cmd="
cd $base
rm -rf bam/tmp
mkdir -p bam/tmp
sourceid=0
mergeinput=\"\"
for sourcebam in \$(echo $sourcebams | tr ',' ' '); do
  sourceid=\$((\$sourceid + 1))
  rgname=\$(basename \$sourcebam .bam)
  java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar AddOrReplaceReadGroups I=\$sourcebam O=bam/tmp/\${sourceid}.bam RGID=\$rgname RGLB=NA RGPL=illumina RGPU=NA RGSM=\$rgname
  mergeinput=\$mergeinput\" I=bam/tmp/\${sourceid}.bam\"
done

java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar MergeSamFiles \$mergeinput O=bam/$merged.bam ASSUME_SORTED=true CREATE_INDEX=true
samtools flagstat bam/$merged.bam >bam/$merged.bam.flagstat
rm -rf bam/tmp
"
      jobname="picard_bam_merge_$merged"
      pbsfn=pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function __markDupAndFilter_20200719 {
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

#filter bam with given mapQ
if [[ ! -e ${outBasename}.filtered.bam ]]; then
  #remove reads with QCFAIL flag (512)
  samtools view -F 512 -q ${map_thresh} -b ${outBasename}.markDup.bam > ${outBasename}.filtered.bam
  samtools flagstat ${outBasename}.filtered.bam > ${outBasename}.filtered.bam.flagstat
  samtools stats ${outBasename}.filtered.bam > ${outBasename}.filtered.bam.stats
fi
'
  jobname="markdupAndFilter_${sname}"
}
