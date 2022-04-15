source $SLURM_ENTRY
wzref_mm10
pipeline_prepare

base=/mnt/isilon/zhoulab/labprojects/20220312_JoeZhou/
cd $base

while read sname bismark_meth_file; do
  jump_comments
  pipeline_depend none

  bismark_meth_file=download/$bismark_meth_file
  output=bed/$sname.bed.gz
  days=2; memG=10; ppn=2      # 1-2h
  pipeline_eval 1 __zlab_parseBismarkCG_20220322

  input=bed/$sname.bed.gz
  output=tmp/CGFULL/$sname
  days=2; memG=10; ppn=2      # 5min
  pipeline_eval 2 __zlab_cgBedToFull_20220322
  ## paste tmp/CGFULL/*_B >tmp/CG_B.txt
  ## paste tmp/CGFULL/*_C >tmp/CG_C.txt

  input=bed/$sname.bed.gz
  outdir=bigwig
  days=2; memG=10; ppn=2      # 5min
  pipeline_eval 3 __zlab_convertBigwigCG_20220219

  input=bed/$sname.bed.gz
  outdir=tmp/FeatureCG/
  group=1
  days=2; memG=10; ppn=2
  pipeline_eval 11 __zlab_featureMean_20220322

  input=bed/$sname.bed.gz
  outdir=tmp/FeatureCG/
  group=2
  days=2; memG=10; ppn=2
  pipeline_eval 12 __zlab_featureMean_20220322

  input=bed/$sname.bed.gz
  outdir=tmp/FeatureCG/
  group=3
  days=2; memG=10; ppn=2
  pipeline_eval 13 __zlab_featureMean_20220322

  input=bed/$sname.bed.gz
  outdir=tmp/FeatureCG/
  group=5
  days=2; memG=10; ppn=2
  pipeline_eval 15 __zlab_featureMean_20220322
  
done << EOM
R106W_bACE_R1	R106W_bACEseq-CAGATC_nonCG.filtered.mm10.134bp.se.bismark.bt2_CpG.meth.gz
R106W_bACE_R2	R106W_bACEseq-CTTGTA_nonCG.filtered.mm10.134bp.se.bismark.bt2_CpG.meth.gz
R106W_BS_R1	R106W_BSseq-ACTTGA_nonCG.filtered.mm10.134bp.se.bismark.bt2_CpG.meth.gz
R106W_BS_R2	R106W_BSseq-TAGCTT_nonCG.filtered.mm10.134bp.se.bismark.bt2_CpG.meth.gz
TKO_bACE	TKO_bACEseq_nonCG.filtered.mm10.134bp.se.bismark.bt2_CpG.meth.gz
WT_bACE_R1	WT_bACEseq-ATCACG_nonCG.filtered.mm10.134bp.se.bismark.bt2_CpG.meth.gz
WT_bACE_R2	WT_bACEseq-TGACCA_nonCG.filtered.mm10.134bp.se.bismark.bt2_CpG.meth.gz
WT_BS_R1	WT_BSseq-CGATGT_nonCG.filtered.mm10.134bp.se.bismark.bt2_CpG.meth.gz
WT_BS_R2	WT_BSseq-GCCAAT_nonCG.filtered.mm10.134bp.se.bismark.bt2_CpG.meth.gz
EOM




while read sname bismark_meth_file; do
  jump_comments
  pipeline_depend none

  bismark_meth_file=download/$bismark_meth_file
  output=tmp/CH/$sname.bed.gz
  days=2; memG=10; ppn=2      # 1-2h
  pipeline_eval 32 __zlab_parseBismarkCH_20220322

  input=tmp/CH/$sname.bed.gz
  outdir=bigwig
  days=2; memG=10; ppn=2
  pipeline_eval 33 __zlab_convertBigwigCH_20220414

  input=tmp/CH/$sname.bed.gz
  outdir=tmp/FeatureCH/
  group=100
  days=2; memG=10; ppn=2
  pipeline_eval 34 __zlab_featureMean_20220322
  
done << EOM
R106W_bACE_R1	R106W_bACEseq-CAGATC_nonCG.filtered.mm10.134bp.se.bismark.bt2_nonCG.meth.gz
R106W_bACE_R2	R106W_bACEseq-CTTGTA_nonCG.filtered.mm10.134bp.se.bismark.bt2_nonCG.meth.gz
R106W_BS_R1	R106W_BSseq-ACTTGA_nonCG.filtered.mm10.134bp.se.bismark.bt2_nonCG.meth.gz
R106W_BS_R2	R106W_BSseq-TAGCTT_nonCG.filtered.mm10.134bp.se.bismark.bt2_nonCG.meth.gz
TKO_bACE	TKO_bACEseq_nonCG.filtered.mm10.134bp.se.bismark.bt2_nonCG.meth.gz
WT_bACE_R1	WT_bACEseq-ATCACG_nonCG.filtered.mm10.134bp.se.bismark.bt2_nonCG.meth.gz
WT_bACE_R2	WT_bACEseq-TGACCA_nonCG.filtered.mm10.134bp.se.bismark.bt2_nonCG.meth.gz
WT_BS_R1	WT_BSseq-CGATGT_nonCG.filtered.mm10.134bp.se.bismark.bt2_nonCG.meth.gz
WT_BS_R2	WT_BSseq-GCCAAT_nonCG.filtered.mm10.134bp.se.bismark.bt2_nonCG.meth.gz
EOM

while read sname; do
  jump_comments
  pipeline_depend none

  # combine is done separately
  input=bed/$sname.bed.gz
  output=tmp/CGFULL/$sname
  days=2; memG=10; ppn=2      # 5min
  pipeline_eval 50 __zlab_cgBedToFull_20220322
  ## paste tmp/CGFULL/*_B >tmp/CG_B.txt
  ## paste tmp/CGFULL/*_C >tmp/CG_C.txt

  input=bed/$sname.bed.gz
  outdir=bigwig
  days=2; memG=10; ppn=2      # 5min
  pipeline_eval 51 __zlab_convertBigwigCG_20220219
  
done << EOM
#R106W_SUB_combined
#WT_SUB_combined
WT_BS_combined
WT_bACE_combined
R106W_bACE_combined
R106W_BS_combined
EOM
