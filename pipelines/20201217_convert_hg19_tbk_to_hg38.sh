source $WZSEQ_ENTRY
pipeline_prepare
indir="/mnt/isilon/zhou_lab/projects/20200106_human_WGBS/bed"
out_bed_dir="/mnt/isilon/zhou_lab/projects/20200106_human_WGBS/bed_hg38"
out_tbk_dir="/mnt/isilon/zhou_lab/projects/20200106_human_WGBS/tbk_hg38"
mkdir -p ${out_bed_dir}
mkdir -p ${out_tbk_dir}

function hg19_to_hg38_20201217() {
  cmd='
base=$("pwd")
cd '${base}'
#Convert hg19 .bed.gz into hg38 .bed.gz with customed python script (use pyliftover).
if [ ! -f '${outfile1}' ]; then
	python /mnt/isilon/zhou_lab/projects/20200928_EPIREJ/Wubin/soloWCGW/bed_hg19_to_hg38.py '${infile}' '${outfile1}'
else
	echo '${outfile1}'" exists"
fi

if [ ! -f '${outfile2}' ]; then
	#zcat '${idx}' | bedtools intersect -a stdin -b '${outfile1}' -loj -F 1 |cut -f 1-3,8-9 | sed "s/\t\./\t-1/g" | bedtools groupby -g 1-3 -c 4,5 -o mean | tbmate pack -s float.int -m '${idx}' - '${outfile2}'
	zcat '${idx}' | bedtools intersect -a stdin -b '${outfile1}' -loj -F 1 |cut -f 1-3,8-9 | tbmate pack -s float.int -m '${idx}' - '${outfile2}'
else
	echo '${outfile2}'" exists"
fi
'
  jobname="tbk_hg19_to_hg38_"$sname
}


for file in `ls ${indir}`; do
  sname=${file//.bed.gz/};
  pipeline_dependlevel
  #echo $sname

  ppn=20;memG=20;
  idx="/mnt/isilon/zhou_lab/projects/20191221_references/hg38/annotation/cpg/idx.gz"
  infile=${indir}/${sname}.bed.gz
  outfile1=${out_bed_dir}/${sname}.bed.gz
  outfile2=${out_tbk_dir}/${sname}.tbk
  pipeline_eval 1 hg19_to_hg38_20201217;
done;
#sh 20201217_convert_hg19_tbk_to_hg38.sh 1 do
