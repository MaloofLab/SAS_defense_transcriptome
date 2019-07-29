# SAS_defense_transcriptome
## Please download data from following places:
## For mapping reads
* raw fastq files: 
* Arabidopsis thaliana cDNA reference sequences

## For differential gene expression analysis
* Saved edgeR dge object ("dge.Rdata")
* Arabidopsis thaliana gene annotationta

## explanation of DEG file names
Please find DEG files under DEGs_with_descritipn folder (KN; where the description has been added in the script?)

RNAseq data was splitted by time points (1h or 49h).

* Interaction model (gt*trt) RNAseq data analysis.Rmd
    + DEgenes.Col.1h.rH.csv: Col, 1h, shade responsive genes (reference is H).
    + DEgenes.Col.49h.4H.csv: Col, 49h, shade responsive genes (reference is H)
    + DEgenes.myc234.1h.rCol.rH.csv: myc234, 1h, (double coefficients, c("gtmyc234", "gtmyc234:trtL") )

* subsetted by genotype further?
    + DEgenes.myc234.1h.rH.csv: myc234, 1h, (KN, what is the difference between "DEgenes.myc234.1h.rCol.rH.csv" and "DEgenes.myc234.1h.rH.csv"? )


## GO term over-representation analysis (GOseq_analysis_CL.Rmd)
* Arabidopsis thaliana GOslim 

## gene clustering (DE_genes_heatmap_scale_center_f.Rmd)

## Schweizer (2013) myc234 transcriptome data reanalized (2013_data_reanalized.Rmd)

## Detailed analysis of PIFtargets and myc234 misregulated genes (MYC234_PIF4target.Rmd, .md, .html, MYC234_PIF4target_files) (by Kazu)

## All output files are stoerd in output folder not included in this repository.
  |-input
  
  |-output
  
  |-SAS_defense_transcriptome (this repository)
 
 



