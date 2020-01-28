
##########################
# reverse bam to fastq
##########################
#
# bam2fastq
#
# wzseq_bam2fastq
# unmark secondary (some bam have only secondary mapping but primary mapping missing)
# samtools view -h in.bam | awk '!/^@/{$2=and($2,compl(0x100)); print $0}/^@/' | samtools view -bo tmp0.bam

function __wzseq_bam2fastq {

  cmd='
set -xe
# group reads by read names, collate is faster than sort
cd '$base';
mkdir -p bam/collate;
mkdir -p fastq
i=1;

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

# remove suffix _1, _2, mark read info to flag and
#  samtools view -h bam/SRR1029055.bam chr19 | awk '!/^@/{inpair=substr($1,length($1),1);$1=substr($1,1,length($1)-2);if(inpair==1) {$2=or($2,0x40);} else {$2=or($2,0x80);} $2=or($2,0x1); if (!(and($2, 0x100))) print $0}/^@/' | samtools collate -uO - SRR1029055tmp >SRR1029055.collate.bam
#  samtools view -h SRR1029055.collate.bam | awk 'BEGIN{key=""; line=""}!/^@/{if (key==$1) {print line; print $0;} key=$1; line=$0;}/^@/' | samtools fastq - -1 fastq_chr19/read1.fastq -2 fastq_chr19/read2.fastq -0 fastq_chr19/unpaired.fastq
function wzseq_bam2fastq {
  base=$(pwd);
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[bam2fastq\]/{p=1;next} p&&!/^$/' samples |
    while read sample sourcebams; do
      cmd="
set -xe
# group reads by read names, collate is faster than sort
cd $base;
mkdir -p bam/collate;
i=1;
for sourcebam in ${sourcebams//,/ }; do
  samtools view -h \$sourcebam | awk -F\"\\t\" -v OFS=\"\t\" '!/^@/{\$2=and(\$2,compl(0x100)); print \$0}/^@/' > bam/collate/${sample}_\${i}_tmp1.sam
  samtools collate -u bam/collate/${sample}_\${i}_tmp1.sam bam/collate/${sample}_\${i}_tmp2
  # samtools collate -u \$sourcebam collate/${sample}_\$i.bam;
  samtools fastq -On -0 fastq/${sample}_\$i.paired_nolabel.fq -1 fastq/${sample}_\$i.pe1.fq -2 fastq/${sample}_\$i.pe2.fq -s fastq/${sample}_\$i.se.fq bam/collate/${sample}_\${i}_tmp2.bam;
  gzip -c fastq/${sample}_\$i.paired_nolabel.fq >>fastq/${sample}.paired_nolabel.fq.gz
  gzip -c fastq/${sample}_\$i.pe1.fq >>fastq/${sample}.pe1.fq.gz
  gzip -c fastq/${sample}_\$i.pe2.fq >>fastq/${sample}.pe2.fq.gz
  gzip -c fastq/${sample}_\$i.se.fq >>fastq/${sample}.se.fq.gz 
  rm -f fastq/${sample}_\$i.paired_nolabel.fq fastq/${sample}_\$i.pe1.fq fastq/${sample}_\$i.pe2.fq fastq/${sample}_\$i.se.fq
  rm -f bam/collate/${sample}_\${i}_tmp{1,2}*;
  i=\$((i+1))
done
"
## the following cleaning creates piping error, better check the size of output fastq manually
# [[ -z \$(gunzip -c fastq/${sample}.paired_nolabel.fq.gz | head -c1) ]] && rm -f fastq/${sample}.paired_nolabel.fq.gz
# [[ -z \$(gunzip -c fastq/${sample}.pe1.fq.gz | head -c1) ]] && rm -f fastq/${sample}.pe1.fq.gz
# [[ -z \$(gunzip -c fastq/${sample}.pe2.fq.gz | head -c1) ]] && rm -f fastq/${sample}.pe2.fq.gz
# [[ -z \$(gunzip -c fastq/${sample}.se.fq.gz | head -c1) ]] && rm -f fastq/${sample}.se.fq.gz
    jobname="bam2fastq_$sample"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 5 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}



