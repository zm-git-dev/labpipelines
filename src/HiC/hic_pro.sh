###############################################################################
# HiC pipeline
################################################################################
# TODO: http://homer.salk.edu/homer/interactions/HiCtagDirectory.html

function examplepipeline_hic {
cat<<- EOF
=== pipeline ===
(+) hic_hicpro
EOF
}

## generate fragments Hind-III A^AGCTT, Mbol: ^GATC
## utils/digest_genome.py -r A^AGCTT -o mm9_hindiii.bed mm9.fasta
## [~/software/HiC_Pro/default/annotation]$ ../bin/utils/digest_genome.py -r ^GATC -o Mbol_resfrag_hg19.bed ~/references/hg19_noContig/hg19_noContig.fa
## 
## java -Xmx8g -jar ~/software/juicer/juicer_tools_linux_0.8.jar pre hicpro/O3_57_CTCF_HiCHP-41521568_output/hic_results/data/O3_57_CTCF_HiCHP-41521568/O3_57_CTCF_HiCHP-41521568_allValidPairs.forjuicer hicpro/O3_57_CTCF_HiCHP-41521568_output/hic_results/data/O3_57_CTCF_HiCHP-41521568/O3_57_CTCF_HiCHP-41521568_allValidPairs.forjuicer.hic hg19
function hic_hicpro {

  base=$(pwd)
  [[ -d pbs ]] || mkdir -p pbs;
  [[ -d hicpro ]] || mkdir -p hicpro;
  awk '/^\[/{p=0}/\[hicpro\]/{p=1;next} p&&!/^$/' samples |
    while read sname fastqs; do
      ## need config file with name $sname.hicpro
      mkdir -p hicpro/${sname}_input/$sname;

      ## check whether to overwrite existing output
      if [[ -d hicpro/${sname}_output ]]; then
        read -p "Do you wish to replace existing output: hicpro/${sname}_output [yn]? " yn </dev/tty
        case $yn in
          [Yy]* ) rm -rf hicpro/${sname}_output; ;;
          [Nn]* ) echo "Skip $sname."; continue;;
              * ) echo "Please answer yes or no."; exit;;
        esac
      fi

      ## hicpro_config file, sample-specific config file has higher precedence
      for fastq in ${fastqs//,/ }; do
        ln -sf `rf fastq/$fastq` hicpro/${sname}_input/$sname/$fastq;
      done
      configfile=hicpro/${sname}_hicpro_config
      if [[ ! -e $configfile ]]; then configfile=hicpro/hicpro_config; fi;

      ## write pbs file
    cmd="
cd $base
/home/wanding.zhou/software/HiC_Pro/default/bin/HiC-Pro -c $configfile -i hicpro/${sname}_input -o hicpro/${sname}_output
mkdir -p hicpro/hic
awk '{gsub(/chr/,\"\",\$2); gsub(/chr/,\"\",\$5); print 0,\$2,\$3,0,0,\$5,\$6,1,60}' hicpro/${sname}_output/hic_results/data/${sname}/${sname}_allValidPairs | sort -k1,1 -k6,6 >hicpro/hic/${sname}.juicer
java -Xmx8g -jar ~/software/juicer/juicer_tools_linux_0.8.jar pre hicpro/hic/${sname}.juicer hicpro/hic/${sname}.hic $WZSEQ_REFVERSION
rm -f hicpro/hic/${sname}.juicer
"
      jobname="hicpro_${sname}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 100 -ppn 28
      [[ ${!#} == "do" ]] && qsub $pbsfn
    done
}

function __hic_hicpro {

  ## sname, fastq
  ## need config file with name $sname.hicpro
  ## note that in the config file we specify N_CPU=20, this is hard-coded
  ## hour=48; memG=100; ppn=20; queue=shortq
  mkdir -p hicpro/${sname}_input/$sname;

  sed "s/WZREPLACE_N_CPU/"$ppn"/" $WZSEQ_HICPRO_CONFIG >hicpro/config
  
  ## check whether to overwrite existing output
  rm -rf hicpro/${sname}_output
  ## hicpro_config file, sample-specific config file has higher precedence
  for fastq in ${fastqs//,/ }; do
    ln -sf `readlink -f fastq/$fastq` hicpro/${sname}_input/$sname/
  done
  
  ## write pbs file
  cmd='
cd '$base'
/home/wanding.zhou/software/HiC_Pro/default/bin/HiC-Pro -c hicpro/config -i hicpro/'${sname}'_input -o hicpro/'${sname}'_output
mkdir -p hicpro/hic
awk '\''{gsub(/chr/,"",$2); gsub(/chr/,"",$5); print 0,$2,$3,0,0,$5,$6,1,60}'\'' hicpro/'${sname}'_output/hic_results/data/'${sname}'/'${sname}'_allValidPairs | sort -k1,1 -k6,6 >hicpro/hic/'${sname}'.juicer
java -Xmx8g -jar ~/software/juicer/juicer_tools_linux_0.8.jar pre hicpro/hic/'${sname}'.juicer hicpro/hic/'${sname}'.hic '$WZSEQ_REFVERSION'
rm -f hicpro/hic/'${sname}'.juicer
'
  jobname="hicpro_"${sname}
  pbsfn=$base/pbs/$jobname.pbs
}
