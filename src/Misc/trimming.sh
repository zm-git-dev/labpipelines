################
## trimmomatic
################
# trimming of adaptor and quality
# see http://www.usadellab.org/cms/?page=trimmomatic for explanation
function __wzseq_trimmomatic_SE {
  # hour=12; memG=10; ppn=10
  cmd='
mkdir trimmomatic
java -jar ~/tools/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads '$ppn' -phred33 '$input_fastq' '$output_fastq' ILLUMINACLIP:~/tools/trimmomatic/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>'$input_fastq'_trimmomatic_report
mkdir -p multiqc/raw/trimmomatic/
ln -sf `readlink -f '$input_fastq'_trimmomatic_report` multiqc/raw/trimmomatic/
'
  jobname='trimmomatic_SE_'$sname
}

function __wzseq_trimmomatic_PE {
  # hour=12; memG=10; ppn=10
  cmd='
mkdir -p multiqc/raw/trimmomatic/
java -jar ~/tools/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads '$ppn' -phred33 '$input_fastq1' '$input_fastq2' '$output_fastq1' '$output_fastq1'_unpaired.fastq.gz '$output_fastq2' '$output_fastq2'_unpaired.fastq.gz ILLUMINACLIP:~/tools/trimmomatic/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>multiqc/raw/trimmomatic/'$sname'_trimming_report
'
  jobname='trimmomatic_PE_'$sname
}


