################################################################################
# ChIP-seq pipeline
################################################################################

function examplepipeline_chipseq {
  cat <<- EOF
=== pipeline 2015-01-28 ===
 (+) wzseq_bwa_aln_se / (+) wzseq_bwa_mem =>
 => (+) wzseq_bam_coverage

 (+) chipseq_bcp => (+) chipseq_macs2 => (+) chipseq_gem => (+) chipseq_sissrs
 (+) chipseq_HOMER => (+) chipseq_HOMER_peak

  => (+) chipseq_diffbind_window_with_internal_control
  => (+) chipseq_bam2track_with_internal_control
EOF
}

function __bwa_aln_se_20200716 {
  cmd='
cd '$base'
mkdir -p bam
if [ ! -f bam/'${sname}'.bam ]; then
  bwa aln -t '$ppn' '$WZSEQ_BWA_INDEX' '$fastq' >bam/'${sname}'.sai
  bwa samse '$WZSEQ_BWA_INDEX' bam/'${sname}'.sai '$fastq' | samtools view -bS - | samtools sort -T '$base'/bam/'${sname}' -O bam -o bam/'${sname}'.bam
  samtools index bam/'${sname}'.bam
  samtools flagstat bam/'${sname}'.bam > bam/'${sname}'.bam.flagstat
  rm -f bam/'${sname}'.sai
else
  echo bam/'${sname}'.bam existed
fi
'
  jobname="bwa_align_"${sname}"_SE"
}

# BWA-aln single-ended
function wzseq_bwa_aln_se {

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam;
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread _junk_; do
      sfile=$(readlink -f fastq/$sread);
      cmd="
cd $base
bwa aln -t 10 $WZSEQ_BWA_INDEX $sfile >bam/${sname}.sai
bwa samse $WZSEQ_BWA_INDEX bam/${sname}.sai $sfile | samtools view -bS - | samtools sort -T $base/bam/$sname -O bam -o bam/${sname}.bam
samtools index bam/${sname}.bam
samtools flagstat bam/${sname}.bam > bam/${sname}.bam.flagstat
rm -f bam/${sname}.sai
"
      jobname="bwa_aln_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 10
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# http://homer.salk.edu/homer/
# got to install seqlogo from webLogo
# http://weblogo.berkeley.edu/
# the sequences can be obtained from
# perl configureHomer.pl -install mm10
# sequence and annotation gets downloaded to where you install HOMER
function chipseq_HOMER {
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam;
  [[ -d HOMER ]] || mkdir HOMER;
  for f in bam/*.bam; do
    bfn=$(basename $f .bam)
    cmd="
cd $base
echo Creating Tags for HOMER/$bfn
~/software/HOMER/default/bin/makeTagDirectory HOMER/$bfn $f
"
    jobname="HOMER_$(basename $base)"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 2 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function chipseq_HOMER_peak {
  ## make sure you have seqlogo, gs, blat and samtools available from command line.
  ## targettype: factor, histone, super, groseq, tss, dnase, mC
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      if [[ $targettype == "factor" ]]; then
        targetsize=100;
      else
        targetsize="given"
      fi

      cmd="
cd $base
~/software/HOMER/default/bin/findPeaks HOMER/$targetname -style $targettype -o HOMER/$targetname/peaks.txt -i HOMER/$controlname
export PATH=$PATH:~/software/HOMER/default/bin/:~/software/webLogo/weblogo/
~/software/HOMER/default/bin/findMotifsGenome.pl HOMER/$targetname/peaks.txt $WZSEQ_REFVERSION HOMER/${targetname}/Motif -size $targetsize # -mask? 
~/software/HOMER/default/bin/annotatePeaks.pl HOMER/$targetname/peaks.txt $WZSEQ_REFVERSION > HOMER/$targetname/peaks.txt.annotation
"
      # TODO: there are some more plotting methods to add
      # see http://homer.salk.edu/homer/ngs/quantification.html
      jobname="HOMER_peak_$targetname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 4 -ppn 2
      [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function chipseq_bcp() {

  base=$(pwd);
  [[ -d bcp ]] || mkdir bcp
  [[ -d pbs ]] || mkdir pbs

  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      tfile=$(readlink -f bam/$targetname);
      cfile=$(readlink -f bam/$controlname);

      mkdir bcp/$targetname;
      tbed=bcp/$targetname/$targetname.bcp.bed
      cbed=bcp/$targetname/$controlname.bcp.bed

      if [[ $targettype == "factor" ]]; then
        runcmd1="~/tools/bcp/BCP_v1.1/BCP_TF -1 $tbed -2 $cbed -3 bcp/$targetname/peaks.tfbs"
      fi
      if [[ $targettype == "histone" ]]; then
        runcmd1="~/tools/bcp/BCP_v1.1/BCP_HM -1 $tbed -2 $cbed -3 bcp/$targetname/peaks.hm";
      fi
      
      cmd="
cd $base
wzbam bed6 -bam bam/$targetname.bam -o $tbed
wzbam bed6 -bam bam/$controlname.bam -o $cbed

$runcmd1

rm -f $tbed $cbed
"
      jobname="bcp_$targetname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function __chipseq_macs2_noCtrl {
  : '
target_bam=bam/
hour=12; memG=10; ppn=1
pipeline_eval x __chipseq_macs2
'
  cmd='
cd '$base'
mkdir -p macs2
macs2 callpeak -t '$target_bam' -f BAM -g '$WZSEQ_MACS_SHORT' -n macs2/'$sname' -B --broad
'
  jobname="macs2_"$sname
}

function __dnase_macs2_20200722 {
  cmd='
cd '$base'
[[ -d pbs ]] || mkdir pbs
mkdir -p macs2
/mnt/isilon/zhoulab/labsoftware/anaconda/anaconda3_2020/bin/macs2 callpeak -t filtered_bam/'$sname'.filtered.bam -g '$WZSEQ_MACS_SHORT' -q 0.05 -n '$sname' --outdir macs2 -B
'
  jobname="macs2_"$sname

}

function __chipseq_macs2 {
  : '
target_bam=bam/
control_bam=bam/
hour=12; memG=10; ppn=1
pipeline_eval x __chipseq_macs2
'
  cmd='
cd '$base'
mkdir -p macs2
macs2 callpeak -t '$target_bam' -c '$control_bam' -f BAM -g '$WZSEQ_MACS_SHORT' -n macs2/'$sname' -B --broad
'
  jobname="macs2_"$sname
}

function chipseq_macs2 {

  base=$(pwd);
  [[ -s macs2 ]] || mkdir macs2
  [[ -d pbs ]] || mkdir pbs

  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      cmd="
cd $base
macs2 callpeak -t bam/$targetname.bam -c bam/$controlname.bam -f BAM -g $WZSEQ_MACS_SHORT -n macs2/$targetname -B --broad
"
      jobname="macs2_$targetname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function chipseq_gem {

  base=$(pwd);
  [[ -d gem ]] || mkdir -p gem
  [[ -d pbs ]] || mkdir -p pbs

  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      [[ -d gem/$targetname ]] || mkdir gem/$targetname;
      tbed=gem/$targetname/$targetname.bed
      cbed=gem/$targetname/$controlname.bed
      cmd="
cd $base
wzbam bed6 -bam bam/$targetname.bam -o $tbed
wzbam bed6 -bam bam/$controlname.bam -o $cbed

java -Xmx100G -jar ~/tools/gem/gem/gem.jar --d ~/tools/gem/gem/Read_Distribution_default.txt --g $WZSEQ_REFERENCE.fai --genome $WZSEQ_REFERENCE_SPLIT --s 2000000000 --expt $tbed --ctrl $cbed --f BED --out gem/$targetname --k_min 6 --k_max 13

rm -f $tbed $cbed
"
      jobname="gem_$targetname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 50 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# http://dir.nhlbi.nih.gov/papers/lmi/epigenomes/sissrs/SISSRs-Manual.pdf
# from Keji Zhao lab, allows both with and without background
# WZSEQ_REFERENCE_SIZE can come from the last row of .fai
function chipseq_sissrs {

  [[ -d sissrs ]] || mkdir sissrs
  [[ -d pbs ]] || mkdir pbs
  base=$(pwd)
  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      [[ -d sissrs/$targetname ]] || mkdir sissrs/$targetname;
      tbed=sissrs/$targetname/$targetname.bed
      cbed=sissrs/$targetname/$controlname.bed
      cmd="
cd $base
wzbam bed6 -bam bam/$targetname.bam -o $tbed
wzbam bed6 -bam bam/$controlname.bam -o $cbed
~/software/SISSERs/v1.4/sissrs.pl -i $tbed -b $cbed -o sissrs/$targetname/$targetname.bsites -s $WZSEQ_REFERENCE_SIZE
rm -f $tbed $cbed
"
      jobname="sissrs_$targetname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 20 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function chipseq_diffbind_window_with_internal_control {
  [[ -d pbs ]] || mkdir pbs
  [[ -d diffbind ]] || mkdir diffbind

  # read in internal control
  declare -A sname2control
  while read sname controlval; do
    sname2control[$sname]=$controlval
  done <<EOF
$(awk '/^\[/{p=0}/\[internalcontrol\]/{p=1;next} p&&!/^$/' samples)
EOF

  base=$(pwd)
  awk '/^\[/{p=0}/\[diffbind\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 sname1 sname2; do
      odir=diffbind/${cond1}_vs_${cond2}
      # bam is usually unordered, would be very memory-consuming
      # here bam is first converted to bam and sorted
      [[ ${sname2control[$sname1]+x} ]] && control1=${sname2control[$sname1]} || control1=1
      [[ ${sname2control[$sname2]+x} ]] && control2=${sname2control[$sname2]} || control2=1
      cmd="
cd $base
[[ -d $odir ]] || mkdir $odir

bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w 100 -s 50 | sort -k1,1 -k2,2n -T $odir | bedtools intersect -a - -b <(bedtools bamtobed -i bam/$sname1.bam | awk '\$5>=20' | sort -k1,1 -k2,2n -T $odir) -wao -sorted | awk '{if (\$10==\".\") \$10=0; print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$10}' | bedtools groupby -i - -g 1-3 -c 4 -o sum >$odir/$sname1.bed
wc -l $odir/$sname1.bed
bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w 100 -s 50 | sort -k1,1 -k2,2n -T $odir | bedtools intersect -a - -b <(bedtools bamtobed -i bam/$sname2.bam | awk '\$5>=20' | sort -k1,1 -k2,2n -T $odir) -wao -sorted | awk '{if (\$10==\".\") \$10=0; print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$10}' | bedtools groupby -i - -g 1-3 -c 4 -o sum >$odir/$sname2.bed
wc -l $odir/$sname2.bed

# control fold change and maximum of the two
paste $odir/$sname1.bed $odir/$sname2.bed | awk -f wanding.awk -e '{rawma=max(\$4,\$8); n1=(\$4+1); n2=(\$8+1)*$control1/$control2; ma=max(n1,n2); mi=min(n1,n2); strand=n1<n2?\"+\":\"-\"; if (rawma>1000 && (ma/mi)>5) print \$1\"\t\"\$2\"\t\"\$3\"\t.\t0\t\"strand\"\t\"log(n2/n1)/log(2)\"\t\"\$4\"\t\"\$8}' | sort -k1,1 -k2,2n -T $odir | bedtools merge -s -i - -delim \";\" -c 7,7,8,9 -o mean,count,mean,mean >$odir/${cond1}_vs_${cond2}_diffbind.tsv
"
      jobname="diffbind_${cond1}_vs_${cond2}_with_internal_control"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 80 -ppn 7
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

##################################################
## normalize chipseq signal with internal control
##################################################
function chipseq_bam2track_with_internal_control {

  [[ -d pbs ]] || mkdir pbs
  [[ -d tracks ]] || mkdir tracks

  # read in internal control
  declare -A sname2control
  while read sname controlval; do
    sname2control[$sname]=$controlval
  done <<EOF
$(awk '/^\[/{p=0}/\[internalcontrol\]/{p=1;next} p&&!/^$/' samples)
EOF

  minmapq=10

  base=$(pwd)
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname fastq1 fastq2; do
      # bam is usually unordered, would be very memory-consuming
      # here bam is first converted to bam and sorted
      [[ ${sname2control[$sname]+x} ]] && control=${sname2control[$sname]} || control=1
      cmd="
cd $base
[[ -d tracks ]] || mkdir tracks

bedtools genomecov -ibam bam/$sname.bam -g ${WZSEQ_REFERENCE}.fai -bga -split | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4/$control}' | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ >tracks/${sname}.coverage.bedg
bedGraphToBigWig tracks/${sname}.coverage.bedg ${WZSEQ_REFERENCE}.fai tracks/${sname}.coverage.bw

# rm -f tracks/${sname}.coverage.bedg

samtools view -q $minmapq -b bam/$sname.bam | bedtools genomecov -ibam stdin -g ${WZSEQ_REFERENCE}.fai -bga -split | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4/$control}' | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ >tracks/${sname}.coverage.q10.bedg
bedGraphToBigWig tracks/${sname}.coverage.q10.bedg ${WZSEQ_REFERENCE}.fai tracks/${sname}.coverage.q${minmapq}.bw

# rm -f tracks/${sname}.coverage.q10.bedg
"
      jobname="chipseq_bam2track_with_internal_control_${sname}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 5 -ppn 2
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}
