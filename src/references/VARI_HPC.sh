################################################################################
# reference configuration
# on a new machine, you just need to update this section
################################################################################

##########
## mouse
##########
###### mouse mm10 #####
function wzref_mm10 {
  export WZSEQ_REFVERSION=mm10
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/mm10/mm10.fa
  export WZSEQ_REFERENCE_SIZE=2785373478

  export WZSEQ_REFERENCE_LAMBDAPHAGE=/home/wanding.zhou/references/lambdaphage/biscuit/NC_001416.fa
  export WZSEQ_GTF=/primary/vari/genomicdata/genomes/mm10/tophat/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf
  export WZSEQ_GTF_ENSEMBL=/primary/vari/genomicdata/genomes/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gz

  # $ zcat ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gz | awk '$1~/^#/{print}$1~/^[0-9XYM]/{if ($1=="MT"){$1="chrM"}else{$1="chr"$1};print $0}' | gzip -c >~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming.gz
  export WZSEQ_GTF_ENSEMBL_UCSCNAMING=/primary/vari/genomicdata/genomes/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming

  # [~/references/mm10/gtf]$ python ~/.Renv/versions/3.2.3/lib64/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py Mus_musculus.GRCm38.82.gtf.UCSCnaming Mus_musculus.GRCm38.82.gtf.UCSCnaming.DEXSeq.gff
  export WZSEQ_GTF_DEXSEQ=/home/wanding.zhou/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming.DEXSeq.gff

  export WZSEQ_BOWTIE2_INDEX=/home/wanding.zhou/references/mm10/bowtie2/mm10
  export WZSEQ_BOWTIE1_INDEX=/primary/vari/genomicdata/genomes/mm10/bowtie1/mm10
  export WZSEQ_BWA_INDEX=/home/wanding.zhou/references/mm10/bwa/mm10.fa
  export WZSEQ_STAR_INDEX=/primary/vari/genomicdata/genomes/mm10/STAR
  export WZSEQ_GSNAP_INDEX=/primary/vari/genomicdata/genomes/mm10/gsnap
  # needs $WZSEQ_GSNAP_INDEX/$WZSEQ_GSNAP_SPLICE.iit
  export WZSEQ_GSNAP_SPLICE=Mus_musculus.GRCm38.82.gtf.gsnap.splicesites
  export WZSEQ_SUBREAD_INDEX=/primary/vari/genomicdata/genomes/mm10/subread/mm10

  export WZSEQ_KALLISTO_INDEX=/primary/vari/genomicdata/genomes/mm10/kallisto/mm10.kallisto

  # WGBS indices
  export WZSEQ_BISCUIT_INDEX_LAMBDAPHAGE=/home/wanding.zhou/references/lambdaphage/biscuit/NC_001416.fa
  export WZSEQ_BISCUIT_INDEX=/home/wanding.zhou/references/mm10/biscuit/mm10.fa
  export WZSEQ_BISCUIT_QC_SETUP=/home/wanding.zhou/tools/biscuit/development/biscuit/test/QC_assets/mm10_QC_assets/setup.sh
  export WZSEQ_BWAMETH_INDEX=/home/wanding.zhou/references/mm10/bwameth/mm10.fa

  export WZSEQ_BISMARK_BT1_INDEX=/home/wanding.zhou/references/mm10/bismark_bt1
  export WZSEQ_BISMARK_BT2_INDEX=/home/wanding.zhou/references/mm10/bismark_bt2
  export WZSEQ_BSMAP_INDEX=/home/wanding.zhou/references/mm10/bsmap/mm10.fa

  export WZSEQ_REFERENCE_SPLIT=/primary/vari/genomicdata/genomes/mm10/bowtie1/

  export WZSEQ_MACS_SHORT=mm
  export WZSEQ_CPGBED=/primary/vari/genomicdata/genomes/mm10/annotation/cpg.bed
  export WZSEQ_CGIBED=/home/wanding.zhou/references/mm10/annotation/cpgisland/cpgIslandExt.bed
  export WZSEQ_CGIBED_METHYLKIT=/home/wanding.zhou/references/mm10/annotation/cpgisland/cpgIslandExt.methylKit.bed
  export WZSEQ_TSSBED=/home/wanding.zhou/references/mm10/annotation/TSS/mm10.refseq.tss.bed

  # UCSC table
  export WZSEQ_UCSC_REFSEQ=/home/wanding.zhou/references/mm10/annotation/mm10_RefSeq_FromUCSC.bed
  export WZSEQ_UCSC_CGIBED=/home/wanding.zhou/references/mm10/annotation/mm10_cpgIsland_FromUCSC.bed

  # rmsk
  export WZSEQ_RMSK=/primary/vari/genomicdata/genomes/mm10/annotation/rmsk/rmsk.bed
  # build the following using UCSC table builder
  export WZSEQ_RMSK_GTF=/home/wanding.zhou/references/mm10/annotation/rmsk/rmsk.mm10.gtf
}

#############
## Human
#############

###### human hg19 ####
function wzref_hg19 {
  export WZSEQ_REFVERSION=hg19
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/hg19/hg19.fa
  export WZSEQ_REFERENCE_SIZE=3199901561

  export WZSEQ_GTF=/primary/vari/genomicdata/genomes/hg19/tophat/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
  export WZSEQ_EXON=/primary/vari/genomicdata/genomes/hg19/annotation/gtf_exon_merged.bed
  export WZSEQ_BWA_INDEX=/home/wanding.zhou/references/hg19/bwa/hg19.fa
  export WZSEQ_BOWTIE2_INDEX=/home/wanding.zhou/references/hg19/bowtie2/hg19
  export WZSEQ_BOWTIE1_INDEX=/home/wanding.zhou/references/hg19/bowtie1/hg19
  export WZSEQ_STAR_INDEX=/primary/vari/genomicdata/genomes/hg19/STAR
  export WZSEQ_RSEQC_GENE_BED=/primary/vari/genomicdata/genomes/hg19/rseqc/hg19_GENCODE_GENE_V19_comprehensive.bed

  # awk '!/^#/{if($1~/^[0-9XY]*$/) $1="chr"$1; if($1=="MT") $1="chrM"; print $0}/^#/' gtf/Homo_sapiens.GRCh37.75.gtf >gtf/Homo_sapiens.GRCh37.75.gtf.UCSCnaming
  export WZSEQ_GTF_ENSEMBL_UCSCNAMING=/primary/vari/genomicdata/genomes/hg19/gtf/Homo_sapiens.GRCh37.75.gtf.UCSCnaming
  export WZSEQ_BISMARK_BT2_INDEX=/home/wanding.zhou/references/hg19/bismark_bt2
  export WZSEQ_BWAMETH_INDEX=/home/wanding.zhou/references/hg19/bwameth/hg19.fa
  export WZSEQ_SUBREAD_INDEX=/primary/vari/genomicdata/genomes/hg19/subread/hg19
  export WZSEQ_REFERENCE_SPLIT=/primary/vari/genomicdata/genomes/hg19/tophat/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes
  # export WZSEQ_HISAT2_INDEX=/primary/vari/genomicdata/genomes/hg19/hisat/genome
  export WZSEQ_HISAT2_INDEX=/primary/vari/genomicdata/genomes/hg19_noContig/hisat2/genome
  export WZSEQ_EXOME_CAPTURE=/primary/vari/genomicdata/genomes/hg19/annotation/hg19.exomes.bed
  export WZSEQ_KALLISTO_INDEX=/primary/vari/genomicdata/genomes/hg19/kallisto/hg19.kallisto

  # bedtools makewindows -w 100 -g hg19.fa.fai >windows/hg19.windows100bp.bed
  # seqtk comp -r windows/hg19.windows100bp.bed ~/references/hg19/hg19.fa | awk '$4+$5+$6+$7>0{print $1,$2,$3,($5+$6)/($4+$5+$6+$7)}' >windows/hg19.windows100bp.gc_content.bed
  # sort -k4,4n windows/hg19.windows100bp.gc_content.bed >windows/hg19.windows100bp.gc_content.srted.bed
  # head -3000000 windows/hg19.windows100bp.gc_content.srted.bed | sortbed >windows/hg19.windows100bp.gc_content.top10p.bed
  # tail -3000000 windows/hg19.windows100bp.gc_content.srted.bed | sortbed >windows/hg19.windows100bp.gc_content.bot10p.bed
  export WZSEQ_TOPGC_BED=/primary/vari/genomicdata/genomes/hg19/windows/hg19.windows100bp.gc_content.top10p.bed
  export WZSEQ_BOTGC_BED=/primary/vari/genomicdata/genomes/hg19/windows/hg19.windows100bp.gc_content.bot10p.bed

  # WGBS
  export WZSEQ_BISCUIT_INDEX=/home/wanding.zhou/references/hg19/biscuit/hg19.fa
  export WZSEQ_BSMAP_INDEX=/home/wanding.zhou/references/hg19/bsmap/hg19.fa
  export WZSEQ_CPGBED=/primary/vari/genomicdata/genomes/hg19/annotation/cpg.bed
  export WZSEQ_CGIBED=/primary/vari/genomicdata/genomes/hg19/annotation/cpgisland/cpgIslandExt.bed
  export WZSEQ_TSSBED=/primary/vari/genomicdata/genomes/hg19/annotation/hg19.refseq.tss.bed
  export WZSEQ_MACS_SHORT=hs
  export WZSEQ_BISCUIT_QC_SETUP=/home/wanding.zhou/tools/biscuit/development/biscuit/test/QC_assets/hg19_QC_assets/setup.sh

  # rmsk
  export WZSEQ_RMSK=/primary/vari/genomicdata/genomes/hg19/annotation/rmsk/rmsk.txt.bed
  # build the following using the UCSC table builder
  export WZSEQ_RMSK_GTF=/home/wanding.zhou/references/hg19/annotation/rmsk/rmsk.hg19.gtf

  export WZSEQ_HICPRO_CONFIG=/home/wanding.zhou/wzprojects/2017/2017_05_13_CTCF_EZH2_hicpro_config_template
}

###### human hg19_noContig ####
function wzref_hg19_noContig {
  export WZSEQ_REFVERSION=hg19
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/hg19_noContig/hg19_noContig.fa

  export WZSEQ_BISMARK_BT2_INDEX=/home/wanding.zhou/references/hg19_noContig/bismark_bt2
  export WZSEQ_BWAMETH_INDEX=/home/wanding.zhou/references/hg19_noContig/bwameth/hg19_noContig.fa

  # WGBS
  export WZSEQ_BISCUIT_INDEX=/home/wanding.zhou/references/hg19_noContig/biscuit/hg19_noContig.fa
  export WZSEQ_BSMAP_INDEX=/home/wanding.zhou/references/hg19_noContig/bsmap/hg19_noContig.fa
}


function wzref_hg38 {
  export WZSEQ_REFVERSION=hg38
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/hg38/hg38.fa
  # export WZSEQ_REFERENCE_SIZE=3199901561

  # export WZSEQ_GTF=/primary/vari/genomicdata/genomes/hg19/tophat/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
  # export WZSEQ_BWA_INDEX=/home/wanding.zhou/references/hg19/bwa/hg19.fa
  # export WZSEQ_BOWTIE2_INDEX=/home/wanding.zhou/references/hg19/bowtie2/hg19
  # export WZSEQ_BOWTIE1_INDEX=/home/wanding.zhou/references/hg19/bowtie1/hg19
  # export WZSEQ_STAR_INDEX=/primary/vari/genomicdata/genomes/hg19/STAR
  # export WZSEQ_RSEQC_GENE_BED=/primary/vari/genomicdata/genomes/hg19/rseqc/hg19_GENCODE_GENE_V19_comprehensive.bed
  export WZSEQ_BISMARK_BT2_INDEX=/home/wanding.zhou/references/hg38/bismark_bt2
  export WZSEQ_BISCUIT_INDEX=/primary/vari/genomicdata/genomes/hg38/biscuit/hg38.fa
  export WZSEQ_BISCUIT_QC_SETUP=/home/wanding.zhou/tools/biscuit/development/biscuit/test/QC_assets/hg38_QC_assets/setup.sh
  # export WZSEQ_SUBREAD_INDEX=/primary/vari/genomicdata/genomes/hg19/subread/hg19
  # export WZSEQ_REFERENCE_SPLIT=/primary/vari/genomicdata/genomes/hg19/tophat/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes

  # WGBS
  # export WZSEQ_BISCUIT_INDEX=/home/wanding.zhou/references/hg19/biscuit/hg19.fa
  # export WZSEQ_CGIBED=/primary/vari/genomicdata/genomes/hg19/annotation/cpgisland/cpgIslandExt.bed
  # export WZSEQ_MACS_SHORT=hs

  # rmsk
  # export WZSEQ_RMSK=/primary/vari/genomicdata/genomes/hg19/annotation/rmsk/rmsk.txt.bed
}

###### human hg19 rCRS ####
function wzref_hg19rCRS {
  export WZSEQ_REFVERSION=hg19
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/hg19-rCRS/hg19_rCRS.fa
  export WZSEQ_DBSNP=/primary/vari/genomicdata/genomes/hg19-rCRS/dbsnp_137.hg19.vcf
  export WZSEQ_EXOME_CAPTURE=/primary/vari/genomicdata/genomes/hg19/annotation/hg19.exomes.bed
}


##################
## lancelet
##################

function wzref_braFlo1 {
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/braFlo1/braFlo1.fa
  export WZSEQ_BISCUIT_INDEX=/primary/vari/genomicdata/genomes/braFlo1/biscuit/braFlo1.fa
}

function wzref_braFlo2 {
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/braFlo2/braFlo2.fa
  export WZSEQ_BISCUIT_INDEX=/primary/vari/genomicdata/genomes/braFlo2/biscuit/braFlo2.fa
}

##################
## Elephant Shark
##################

function wzref_calMil1 {
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/calMil1/calMil1.fa
  export WZSEQ_BISCUIT_INDEX=/primary/vari/genomicdata/genomes/calMil1/biscuit/calMil1.fa
}

############
## Lamprey
############

function wzref_petMar2 {
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/petMar2/petMar2.fa
  export WZSEQ_BISCUIT_INDEX=/primary/vari/genomicdata/genomes/petMar2/biscuit/petMar2.fa
}

##############
## Sea Squirt
##############

function wzref_ci2 {
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/ci2/ci2.fa
}