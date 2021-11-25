# Pipeline Repository

## Usage

Just clone this repository point an environment variable __WZSEQ_ENTRY__ to your installed `repo/entry/wzseq.sh`. Then you should be able to use scripts in the __pipelines__ folder.



## PacBio SMRT (Single Molecule Real Time) Sequencing Preprocessing

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

## Alternative Splicing

1. Alternative splicing were analyzed using [SUPPA](https://github.com/comprna/SUPPA)

## Coding potential
1. [CNCI](https://github.com/www-bioinfo-org/CNCI)
2. [PLEK](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-311)
3. [COME](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5224497/)


# Hi-C Methods

1. [AFC](https://www.nature.com/articles/nature12644)
2. [Fit-Hi-C](https://noble.gs.washington.edu/proj/fit-hi-c/)
3. [HiCNorm](https://academic.oup.com/bioinformatics/article/28/23/3131/192582)
4. [HUGIn](https://academic.oup.com/bioinformatics/article/33/23/3793/3861336) a web browser for HiC Count data, developed jointly by Yun Li, Ming Hu and Bing Ren.
5. [FastHiC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013904/) cannot find code
6. [HMRF](https://academic.oup.com/bioinformatics/article/32/5/650/1744391) cannot find code

