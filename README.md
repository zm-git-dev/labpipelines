# Zhou Lab's Pipeline Repository

## Usage

Clone this repository point an environment variable __SLURM_ENTRY__ to your installed `repo/entry/slurm.sh`. Then you should be able to use scripts in the __pipelines__ folder. The modules should stay in the __modules__ folder.

## Example pipelines

[pipelines/20211125_Waleed_WGBS.sh](20211125_Waleed_WGBS.sh)

## Available modules

Note that new module code should have the naming convension `__zlab_function_date`. Other modules noncompliant with this convention are obsolete.
All modules must provide project base name.

### BAM utilities - [src/bam.sh](src/bam.sh)

```
__zlab_PicardMarkdup_20200719
__zlab_PicardMarkdup_20211125
__zlab_sortIndexBam_20211126
__zlab_bam2fastq_20211125
__zlab_indexBAM_20211125

```

### Bisulfite-seq

```
__zlab_trimGaloreSE_20211125 
__zlab_trimGalorePE_20211125
__zlab_BismarkBt2SE_20211125
__zlab_BismarkMethExtract_20211125

```

## SRA

```
__zlab_fasterqDumpSE_20200706
__zlab_fastqDumpSE_20211126
```

## Obsolete module functions

```
############
## bam.sh ##
############
bam2bigwig
__filterMAPQ_20200719
wzseq_picard_index_fasta
wzseq_picard_WGSmetrics
wzseq_picard_markdup
wzseq_merge_bam_picard
__markDupAndFilter_20200719
wzseq_cov5
wzseq_cov10
wzseq_liftbw
wzseq_merge_fastq
wzseq_srx2srr
wzseq_merge_bam
__wzseq_index_bam
wzseq_basic_summary
__wzseq_bam_mapq
wzseq_bam2fastq
__wzseq_uniformity_1M
wzseq_bam_coverage
wgbs_cpgcoverage_OBSOLETE

#######################
## bisulfite_seq.sh ##
#######################
wgbs_adaptor
wzseq_fastqc
wgbs_biscuit_align
wgbs_bwameth
wgbs_bismark_bowtie1
wgbs_bismark_bowtie2
wgbs_bsmap
wgbs_biscuit_align_lambdaphage
wgbs_biscuit_pileup_lambdaphage (TODO: exclude human reads)
zseq_GATK_realign
wzseq_picard_markdup 
wzseq_clean_intermediate
wgbs_basequal_recal

wzseq_merge_bam
wzseq_qualimap 
wzseq_picard_WGSmetrics
wzseq_bam_coverage

wgbs_methpipe
wgbs_methpipe_methylome
wgbs_methpipe_allele

wgbs_biscuit_pileup [-nome]
wgbs_vcf2tracks [-nome] 
wgbs_cpgcoverage 
wgbs_repeat => (+) wgbs_repeat_diff
wgbs_methylKit_summary
wgbs_diffmeth_simple
wgbs_methylKit_diffmeth
wgbs_methpipe_diff 
wgbs_metilene

```

## To be added
### PacBio SMRT (Single Molecule Real Time) Sequencing Preprocessing

1.  preliminary processing using SMRTlink 5.0 software (Pacific Biosciences)
2. generate Circular consensus sequence (CCS)
3. The non-chimeric reads, which include non-full length and full-length transcripts, were then clustered by isoform level clustering (ICE) algorithm.
4. The produced clusters were finally polished using ARROW software (Pacific Biosciences)
5. Additional nucleotide errors in consensus reads were corrected using the Illumina RNA-seq data by the LoRDEC software
6. Consensus reads were aligned to reference annotations (Ensemble 38 release 91) using GMAP with the following parameters ```--no-chimeras --cross-species --expand-offsets 1 -B 5 -K 50000 -f samse -n 1```
7. gene model were downloaded from [here](ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapi ens.GRCh38.dna.primary_assembly.fa.gz)
8. [HTSeq v0.6.1](https://htseq.readthedocs.io/en/release_0.9.1/) was used to calculate FPKM.
9. Unmapped transcripts and novel gene transcripts were scanned and annotated by Diamond BLASTX with parameter `e value '1e-5'` in the following protein/peptide database: NR (NCBI non-redundant protein sequences), KOG/COG (Clusters of Orthologous Groups of proteins), Swiss-Prot (a manually annotated and reviewed protein sequence database), KEGG ortholog database.
10. novel transcripts were also searched against Pfam using [Hmmscan](http://hmmer.org/download.html)

### Alternative Splicing

1. Alternative splicing were analyzed using [SUPPA](https://github.com/comprna/SUPPA)

### Coding potential
1. [CNCI](https://github.com/www-bioinfo-org/CNCI)
2. [PLEK](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-311)
3. [COME](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5224497/)

### Hi-C

1. [AFC](https://www.nature.com/articles/nature12644)
2. [Fit-Hi-C](https://noble.gs.washington.edu/proj/fit-hi-c/)
3. [HiCNorm](https://academic.oup.com/bioinformatics/article/28/23/3131/192582)
4. [HUGIn](https://academic.oup.com/bioinformatics/article/33/23/3793/3861336) a web browser for HiC Count data, developed jointly by Yun Li, Ming Hu and Bing Ren.
5. [FastHiC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013904/) cannot find code
6. [HMRF](https://academic.oup.com/bioinformatics/article/32/5/650/1744391) cannot find code

