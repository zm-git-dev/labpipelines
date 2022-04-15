################################################################################
# reference configuration
# on a new machine, you just need to update this section
################################################################################


##########
## mouse
##########
function wzref_mm10LambdaT4 {
  refbase=/mnt/isilon/zhou_lab/projects/20191221_references
  export WZSEQ_BISMARK_BT2_INDEX=$refbase/mm10LambdaT4/bismark_bt2
}


###### mouse mm10 #####
function wzref_mm10 {
  refbase=/mnt/isilon/zhou_lab/projects/20191221_references
  export WZSEQ_REFVERSION=mm10
  export WZSEQ_BSSEQ_FEAT=$refbase/mm10/features/
  export WZSEQ_REFERENCE=$refbase/mm10/mm10.fa
  export WZSEQ_REFERENCE_SIZE=2785373478

  export WZSEQ_REFERENCE_LAMBDAPHAGE=$refbase/lambdaphage/biscuit/NC_001416.fa
  export WZSEQ_GTF=$refbase/mm10/tophat/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf
  export WZSEQ_GTF_ENSEMBL=$refbase/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gz

  # $ zcat ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gz | awk '$1~/^#/{print}$1~/^[0-9XYM]/{if ($1=="MT"){$1="chrM"}else{$1="chr"$1};print $0}' | gzip -c >~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming.gz
  export WZSEQ_GTF_ENSEMBL_UCSCNAMING=$refbase/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming

  # [~/references/mm10/gtf]$ python ~/.Renv/versions/3.2.3/lib64/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py Mus_musculus.GRCm38.82.gtf.UCSCnaming Mus_musculus.GRCm38.82.gtf.UCSCnaming.DEXSeq.gff
  export WZSEQ_GTF_DEXSEQ=$refbase/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming.DEXSeq.gff

  export WZSEQ_BOWTIE2_INDEX=$refbase/mm10/bowtie2/mm10
  export WZSEQ_BOWTIE1_INDEX=$refbase/mm10/bowtie1/mm10
  export WZSEQ_BWA_INDEX=$refbase/mm10/bwa/mm10.fa
  export WZSEQ_STAR_INDEX=$refbase/mm10/STAR
  export WZSEQ_GSNAP_INDEX=$refbase/mm10/gsnap
  # needs $WZSEQ_GSNAP_INDEX/$WZSEQ_GSNAP_SPLICE.iit
  export WZSEQ_GSNAP_SPLICE=Mus_musculus.GRCm38.82.gtf.gsnap.splicesites
  export WZSEQ_SUBREAD_INDEX=$refbase/mm10/subread/mm10

  export WZSEQ_KALLISTO_INDEX=$refbase/mm10/kallisto/mm10.kallisto

  # WGBS indices
  export WZSEQ_BISCUIT_INDEX_LAMBDAPHAGE=$refbase/lambdaphage/biscuit/NC_001416.fa
  export WZSEQ_BISCUIT_INDEX=$refbase/mm10/biscuit/mm10.fa
  export WZSEQ_BISCUIT_QC_SETUP=/home/zhouw3/tools/biscuit/development/biscuit/test_shen/QC_assets/mm10_QC_assets/setup.sh
  export WZSEQ_BWAMETH_INDEX=$refbase/mm10/bwameth/mm10.fa

  export WZSEQ_BISMARK_BT1_INDEX=$refbase/mm10/bismark_bt1
  export WZSEQ_BISMARK_BT2_INDEX=$refbase/mm10/bismark_bt2
  export WZSEQ_BSMAP_INDEX=$refbase/mm10/bsmap/mm10.fa

  export WZSEQ_REFERENCE_SPLIT=$refbase/mm10/bowtie1/

  export WZSEQ_MACS_SHORT=mm
  export WZSEQ_CPGBED=$refbase/mm10/annotation/cpg/cpg.bed.gz
  export WZSEQ_CGIBED=$refbase/mm10/annotation/cpgisland/cpgIslandExt.bed
  export WZSEQ_CGIBED_METHYLKIT=$refbase/mm10/annotation/cpgisland/cpgIslandExt.methylKit.bed
  export WZSEQ_TSSBED=$refbase/mm10/annotation/TSS/mm10.refseq.tss.bed

  # UCSC table
  export WZSEQ_UCSC_REFSEQ=$refbase/mm10/annotation/mm10_RefSeq_FromUCSC.bed
  export WZSEQ_UCSC_CGIBED=$refbase/mm10/annotation/mm10_cpgIsland_FromUCSC.bed

  # rmsk
  export WZSEQ_RMSK=$refbase/mm10/annotation/rmsk/rmsk.bed
  # build the following using UCSC table builder
  export WZSEQ_RMSK_GTF=$refbase/mm10/annotation/rmsk/rmsk.mm10.gtf
}

#############
## Human
#############

###### human hg19 ####
function wzref_hg19 {
  refbase=/mnt/isilon/zhou_lab/projects/20191221_references
  export WZSEQ_REFVERSION=hg19
  export WZSEQ_REFERENCE=$refbase/hg19/hg19.fa
  export WZSEQ_REFERENCE_SIZE=3199901561

  export WZSEQ_GTF=$refbase/hg19/tophat/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
  export WZSEQ_EXON=$refbase/hg19/annotation/gtf_exon_merged.bed
  export WZSEQ_BWA_INDEX=$refbase/hg19/bwa/hg19.fa
  export WZSEQ_BOWTIE2_INDEX=$refbase/hg19/bowtie2/hg19
  export WZSEQ_BOWTIE1_INDEX=$refbase/hg19/bowtie1/hg19
  export WZSEQ_STAR_INDEX=$refbase/hg19/STAR
  export WZSEQ_RSEQC_GENE_BED=$refbase/hg19/rseqc/hg19_GENCODE_GENE_V19_comprehensive.bed

  # awk '!/^#/{if($1~/^[0-9XY]*$/) $1="chr"$1; if($1=="MT") $1="chrM"; print $0}/^#/' gtf/Homo_sapiens.GRCh37.75.gtf >gtf/Homo_sapiens.GRCh37.75.gtf.UCSCnaming
  export WZSEQ_GTF_ENSEMBL_UCSCNAMING=$refbase/hg19/gtf/Homo_sapiens.GRCh37.75.gtf.UCSCnaming
  export WZSEQ_BISMARK_BT2_INDEX=$refbase/hg19/bismark_bt2
  export WZSEQ_BWAMETH_INDEX=$refbase/hg19/bwameth/hg19.fa
  export WZSEQ_SUBREAD_INDEX=$refbase/hg19/subread/hg19
  export WZSEQ_REFERENCE_SPLIT=$refbase/hg19/tophat/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes
  # export WZSEQ_HISAT2_INDEX=$refbase/hg19/hisat/genome
  export WZSEQ_HISAT2_INDEX=$refbase/hg19_noContig/hisat2/genome
  export WZSEQ_EXOME_CAPTURE=$refbase/hg19/annotation/hg19.exomes.bed
  export WZSEQ_KALLISTO_INDEX=$refbase/hg19/kallisto/hg19.kallisto

  # bedtools makewindows -w 100 -g hg19.fa.fai >windows/hg19.windows100bp.bed
  # seqtk comp -r windows/hg19.windows100bp.bed ~/references/hg19/hg19.fa | awk '$4+$5+$6+$7>0{print $1,$2,$3,($5+$6)/($4+$5+$6+$7)}' >windows/hg19.windows100bp.gc_content.bed
  # sort -k4,4n windows/hg19.windows100bp.gc_content.bed >windows/hg19.windows100bp.gc_content.srted.bed
  # head -3000000 windows/hg19.windows100bp.gc_content.srted.bed | sortbed >windows/hg19.windows100bp.gc_content.top10p.bed
  # tail -3000000 windows/hg19.windows100bp.gc_content.srted.bed | sortbed >windows/hg19.windows100bp.gc_content.bot10p.bed
  export WZSEQ_TOPGC_BED=$refbase/hg19/windows/hg19.windows100bp.gc_content.top10p.bed
  export WZSEQ_BOTGC_BED=$refbase/hg19/windows/hg19.windows100bp.gc_content.bot10p.bed

  # WGBS
  export WZSEQ_BISCUIT_INDEX=$refbase/hg19/biscuit/hg19.fa
  export WZSEQ_BSMAP_INDEX=$refbase/hg19/bsmap/hg19.fa
  export WZSEQ_CPGBED=$refbase/hg19/annotation/cpg/cpg.bed.gz
  export WZSEQ_CGIBED=$refbase/hg19/annotation/cpgisland/cpgIslandExt.bed
  export WZSEQ_TSSBED=$refbase/hg19/annotation/hg19.refseq.tss.bed
  export WZSEQ_MACS_SHORT=hs
  export WZSEQ_BISCUIT_QC_ASSETS=$refbase/hg19/biscuit/QC_assets

  # rmsk
  export WZSEQ_RMSK=$refbase/hg19/annotation/rmsk/rmsk.txt.bed
  # build the following using the UCSC table builder
  export WZSEQ_RMSK_GTF=$refbase/hg19/annotation/rmsk/rmsk.hg19.gtf

  export WZSEQ_HICPRO_CONFIG=/home/zhouw3/wzprojects/2017/2017_05_13_CTCF_EZH2_hicpro_config_template
}

###### human hg19_noContig ####
function wzref_hg19_noContig {
  refbase=/mnt/isilon/zhou_lab/projects/20191221_references
  export WZSEQ_REFVERSION=hg19
  export WZSEQ_REFERENCE=$refbase/hg19_noContig/hg19_noContig.fa

  export WZSEQ_BISMARK_BT2_INDEX=$refbase/hg19_noContig/bismark_bt2
  export WZSEQ_BWAMETH_INDEX=$refbase/hg19_noContig/bwameth/hg19_noContig.fa

  # WGBS
  export WZSEQ_BISCUIT_INDEX=$refbase/hg19_noContig/biscuit/hg19_noContig.fa
  export WZSEQ_BSMAP_INDEX=$refbase/hg19_noContig/bsmap/hg19_noContig.fa
}


function wzref_hg38 {
  refbase=/mnt/isilon/zhou_lab/projects/20191221_references
  export WZSEQ_REFVERSION=hg38
  export WZSEQ_REFERENCE=$refbase/hg38/hg38.fa
  # export WZSEQ_REFERENCE_SIZE=3199901561

  # export WZSEQ_GTF=$refbase/hg19/tophat/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
  # export WZSEQ_BWA_INDEX=$refbase/hg19/bwa/hg19.fa
  # export WZSEQ_BOWTIE2_INDEX=$refbase/hg19/bowtie2/hg19
  # export WZSEQ_BOWTIE1_INDEX=$refbase/hg19/bowtie1/hg19
  # export WZSEQ_STAR_INDEX=$refbase/hg19/STAR
  export WZSEQ_STAR_INDEX=$refbase/hg38/STAR
  # export WZSEQ_RSEQC_GENE_BED=$refbase/hg19/rseqc/hg19_GENCODE_GENE_V19_comprehensive.bed
  export WZSEQ_BISMARK_BT2_INDEX=$refbase/hg38/bismark_bt2
  export WZSEQ_BISCUIT_INDEX=$refbase/hg38/biscuit/hg38.fa
  export WZSEQ_BISCUIT_QC_SETUP=$refbase/hg38/biscuit/QC/
  export WZSEQ_GTF_ENSEMBL_UCSCNAMING=$refbase/hg38/annotation/gencode.v28.annotation.gtf
  # export WZSEQ_SUBREAD_INDEX=$refbase/hg19/subread/hg19
  # export WZSEQ_REFERENCE_SPLIT=$refbase/hg19/tophat/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes
  
  # WGBS
  # export WZSEQ_BISCUIT_INDEX=$refbase/hg19/biscuit/hg19.fa
  # export WZSEQ_CGIBED=$refbase/hg19/annotation/cpgisland/cpgIslandExt.bed
  # export WZSEQ_MACS_SHORT=hs
  export WZSEQ_CPGBED=$refbase/hg38/annotation/cpg/cpg_noDecoy.bed.gz

  # rmsk
  # export WZSEQ_RMSK=$refbase/hg19/annotation/rmsk/rmsk.txt.bed
  export WZSEQ_RMSK_GTF=$refbase/hg38/annotation/rmsk/rmsk_hg38.gtf
}

###### human hg19 rCRS ####
function wzref_hg19rCRS {
  refbase=/mnt/isilon/zhou_lab/projects/20191221_references
  export WZSEQ_REFVERSION=hg19
  export WZSEQ_REFERENCE=$refbase/hg19-rCRS/hg19_rCRS.fa
  export WZSEQ_DBSNP=$refbase/hg19-rCRS/dbsnp_137.hg19.vcf
  export WZSEQ_EXOME_CAPTURE=$refbase/hg19/annotation/hg19.exomes.bed
}


