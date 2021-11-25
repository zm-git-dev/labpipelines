#################
## SRA-toolkits
#################

## I have spent a lot of time finagling with fastq-dump but it gave me the timeout error randomly.
## I think this issue is specific to respublica.
## Please use the following (prefetch + fasterq-dump)
function __sra_fasterq_dump_PE_20200706 {
  cmd='
set -xe
mkdir -p '$base'/fastq
cd '$base'/fastq
rm -f '${sname}'_R1.fastq.gz '${sname}'_R2.fastq.gz
for f in '$srr_ids'; do

  /mnt/isilon/zhoulab/labsoftware/sra-toolkit/sratoolkit.2.10.8-centos_linux64/bin/prefetch --resume yes --max-size 1T -O sra/ ${f}
  /mnt/isilon/zhoulab/labsoftware/sra-toolkit/sratoolkit.2.10.8-centos_linux64/bin/fasterq-dump --force --temp sra/${f} --split-3 --split-files sra/${f}.sra
  
  pigz -p '$ppn' -c ${f}.sra_1.fastq >>'${sname}'_R1.fastq.gz
  [[ -e ${f}.sra_2.fastq ]] && pigz -p '$ppn' -c ${f}.sra_2.fastq >>'${sname}'_R2.fastq.gz

  # cleaning, delete .fastq only in case .sra needs to be used later
  rm -f ${f}.sra_1.fastq ${f}.sra_2.fastq
  # rm -f sra/${f}.sra
done
'
  jobname="fasterqdump_"$sname
}

function __sra_fasterq_dump_SE_20200706 {
  cmd='
set -xe
mkdir -p '$base'/fastq
cd '$base'/fastq
rm -f '${sname}'_R1.fastq.gz
for f in '$srr_ids'; do

  /mnt/isilon/zhoulab/labsoftware/sra-toolkit/sratoolkit.2.10.8-centos_linux64/bin/prefetch --resume yes --max-size 1T -O sra/ ${f}
  /mnt/isilon/zhoulab/labsoftware/sra-toolkit/sratoolkit.2.10.8-centos_linux64/bin/fasterq-dump --force --temp sra/${f} --split-3 --split-files sra/${f}.sra
  
  pigz -p '$ppn' -c ${f}.sra.fastq >>'${sname}'_R1.fastq.gz

  # cleaning, delete .fastq only in case .sra needs to be used later
  rm -f ${f}.sra.fastq 
  # rm -f sra/${f}.sra
done
'
  jobname="fasterqdump_"$sname
}

function __wzseq_fastq_dump_SE {
  # input: base, sname, srr_ids
  cmd='
set -xe
mkdir -p '$base'/fastq
cd '$base'/fastq
rm -f '${sname}'.fastq.gz
for f in '$srr_ids'; do
  /mnt/isilon/zhoulab/labsoftware/sra-toolkit/default/bin/fastq-dump -S -3 --gzip $f;
  cat ${f}.fastq.gz >>'${sname}'.fastq.gz
  rm -f ${f}.fastq.gz
done
'
  jobname='fastqdump_'$sname
}

function __wzseq_fastq_dump_PE {
  cmd='
set -xe
mkdir -p '$base'/fastq
cd '$base'/fastq
rm -f '${sname}'_R1.fastq.gz '${sname}'_R2.fastq.gz
for f in '$srr_ids'; do
  /mnt/isilon/zhoulab/labsoftware/sra-toolkit/default/bin/fastq-dump --split-files --gzip $f;
  cat ${f}_1.fastq.gz >>'${sname}'_R1.fastq.gz
  [[ -e ${f}_2.fastq.gz ]] && cat ${f}_2.fastq.gz >>'${sname}'_R2.fastq.gz
  rm -f ${f}_1.fastq.gz ${f}_2.fastq.gz
done
'
  jobname="fastqdump_"$sname
}

# # download using fastq-dump with accession numbers
# # section
# # [sra]
# # HUES64_derived_CD56_Mesoderm  PE	SRR1067566,SRR1067568,SRR1067569,SRR1067570
# function wzseq_fastq_dump() {
#   base=$(pwd);
#   [[ -d pbs ]] || mkdir pbs
#   awk '/^\[/{p=0}/\[sra\]/{p=1;next} p&&!/^$/' samples |
#     while read sname pe_or_se srr_ids; do
#       srr_ids=${srr_ids//,/ };
#       pipeline_init
#       hour=48; memG=10; ppn=1
#       if [[ $pe_or_se == "PE" ]]; then
#         pipeline_eval 1 __wzseq_fastq_dump_PE
#       else
#         pipeline_eval 1 __wzseq_fastq_dump_SE
#       fi
#     done
# }

# # convert downloaded sra to fastq
# function wzseq_sra_to_fastq() {
#   base=$(pwd);
#   [[ -d pbs ]] || mkdir pbs
#   [[ -d fastq ]] || mkdir fastq
#   for f in sra/*.sra; do
#     fn=$(readlink -f $f)
#     bfn=$(basename $f .sra)
#     # if you want to append 1/2 to read name use -I. But this
#     # usually causes trouble
#     cmd="
# fastq-dump --split-files $fn -O $base/fastq --gzip
# "
#     jobname="sra2fastq_$bfn"
#     pbsfn=$base/pbs/$jobname.pbs
#     pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 10 -memG 10 -ppn 1
#     [[ $1 == "do" ]] && qsub $pbsfn
#   done
# }
