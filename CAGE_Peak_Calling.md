# Plan for peak calling on CAGE data 


## Run strand-aware peak-calling methods on CAGE data

Methods
- DBSCAN
- DBSCAN + sliding window 
- scipy find peaks 
- MACS3
- CAGEr

*Inputs* from bam format :

### BAM 1 
- /private/groups/brookslab/gabai/projects/yeastMeth/data/rna/cage/bam/alpha_YPD/DRR524315.sorted.bam
### BAM 2 
- /private/groups/brookslab/gabai/projects/yeastMeth/data/rna/cage/bam/alpha_YPD/DRR524316.sorted.bam

*Outputs* 
- bed files predicting promoter regions 

## QC

- Length dist. of called regions 

### Prepare Long-read data 

Input : *Long-read teloprime aligned bam file*
/private/groups/brookslab/gabai/projects/yeastMeth/data/rna/teloprime/tx_alignment/250609_teloprime_ys18_rep1.bam

Output : *TSS's for long reads* 

### Run bedtools closest

*activate bedtools_env*

Input : *Long-read TSS's* , *CAGE peak bed file*

Output :

- Precision/Recall of long-reads in CAGE peaks 
- Dist from long-read to annotated CAGE peak (% within 50 bp of CAGE peak)

### Prepare Annotation 

Input : *Annotation File *
/private/groups/brookslab/gabai/projects/yeastMeth/data/ref/sacCer3_ares_v13_sorted.bed

```bash
$ head /private/groups/brookslab/gabai/projects/yeastMeth/data/ref/sacCer3_ares_v13_sorted.bed
chrI    334     649     YAL069W 1000    +       334     649     0       1       315,    0,
chrI    537     792     YAL068W-A       1000    +       537     792     0       1       255,    0,
chrI    1806    2169    PAU8    1000    -       1806    2169    0       1       363,    0,
chrI    2479    2707    YAL067W-A       1000    +       2479    2707    0       1       228,    0,
chrI    7234    9016    SEO1    1000    -       7234    9016    0       1       1782,   0,
chrI    10090   10399   YAL066W 1000    +       10090   10399   0       1       309,    0,
chrI    11564   11951   YAL065C 1000    -       11564   11951   0       1       387,    0,
chrI    12045   12426   YAL064W-B       1000    +       12045   12426   0       1       381,    0,
chrI    13362   13743   YAL064C-A       1000    -       13362   13743   0       1       381,    0,
chrI    21565   21850   YAL064W 1000    +       21565   21850   0       1       285,    0,
```

Output : TSV of start sites (chr  |  Pos   | Strand)

### Run bedtools closest 

*activate bedtools_env*

Input : *Annotation TSS's* , *CAGE peak bed file*

Output :
- Number of peaks x Annotated Start Sites baseline 
- Dist from annotated start-site to annotated CAGE peak (% within 50 bp of CAGE peak)

### Summarization box plots and bars 

One plot for each CAGE bam file 

Input : 
- Length dist. of called regions 

Output : 
- Box plots with deviations, seperated by method 

Input : TSV with - 
- Precision/Recall of long-reads in CAGE peaks 
- Dist from long-read to annotated CAGE peak (% within 50 bp of CAGE peak)
- Number of peaks x Annotated Start Sites baseline 
- Dist from annotated start-site to annotated CAGE peak (% within 50 bp of CAGE peak)

Output :
Bar plots for each metric 