source $WZSEQ_ENTRY
pipeline_prepare
tbk_dir="/mnt/isilon/zhou_lab/projects/20200106_human_WGBS/tbk_hg38"
hg38_100kb_bin_path="/mnt/isilon/zhou_lab/projects/20200928_EPIREJ/Wubin/MethPredictor/data/hg38_100kb_bin.bed"
outdir1="100kb_bin_avg_beta"
outdir2="100kb_bin_avg_depth"

function calculate_100kb_avg_beta() {
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p '${outdir1}'

if [ ! -f '${outfile}' ]; then
	# tbmate view -N ''"''"'' -d '${infile}' | awk '"'"'$4!=""'"'"'| bedtools intersect -a '${hg38_100kb_bin_path}' -b stdin -loj | awk '"'"'$5!="." && $8!=""'"'"' | bedtools groupby -g 1-4 -c 8 -o mean >  '${outfile}'
	tbmate view -N ''"''"'' -d '${infile}' | awk '"'"'$4!=""'"'"'|bedtools intersect -a stdin -b '${solowcgw}' -wa | bedtools intersect -a '${hg38_100kb_bin_path}' -b stdin -loj | awk '"'"'$5!="." && $8!=""'"'"' | bedtools groupby -g 1-4 -c 8 -o mean >  '${outfile}'
else
        echo '${outfile}'" exists"
fi

'
  jobname="calculate_100kb_avg_beta_"$sname
}

function calculate_100kb_avg_depth() {
  cmd='
base=$("pwd")
cd '${base}'
mkdir -p '${outdir2}'

if [ ! -f '${outfile}' ]; then
	#tbmate view -N ''"''"'' -bd '${infile}' | cut -f 1-3,5| bedtools intersect -a '${hg38_100kb_bin_path}' -b stdin -loj | awk '"'"'$5!="." && $8!=-1'"'"' | bedtools groupby -g 1-4 -c 8 -o mean >  '${outfile}'
	tbmate view -N ''"''"'' -bd '${infile}' | cut -f 1-3,5 |bedtools intersect -a stdin -b '${solowcgw}' -wa | bedtools intersect -a '${hg38_100kb_bin_path}' -b stdin -loj | awk '"'"'$5!="." && $8!=-1'"'"' | bedtools groupby -g 1-4 -c 8 -o mean >  '${outfile}'
else
	echo '${outfile}'" exists"
fi

'
  jobname="calculate_100kb_avg_depth_"$sname
}

solowcgw="/mnt/isilon/zhou_lab/projects/20200928_EPIREJ/Wubin/MethPredictor/data/solowcgw.bed"

if [ ! -f ${solowcgw} ]; then
	zcat /mnt/isilon/zhou_lab/projects/20200928_EPIREJ/Wubin/MethPredictor/data/cpgs.info.gz | sed '1d' | awk '$8==1 && $9==""' |cut -f 1-3 > ${solowcgw}
fi

for tbk_file in `ls ${tbk_dir}`; do
  sname=${tbk_file//.tbk/};
  pipeline_dependlevel
  #echo $sname

  ppn=1;memG=10;
  infile=${tbk_dir}/${sname}.tbk
  outfile=${outdir1}/${sname}.bed
  pipeline_eval 1 calculate_100kb_avg_beta;
  
  outfile=${outdir2}/${sname}.bed
  pipeline_eval 2 calculate_100kb_avg_depth;
done;
#sh ~/labpipeline/pipelines/20210104_calculate_100kb_avg_beta.sh 1 do
