---
title: "MYC234_PIFtargets"
author: "Kazu"
date: "5/22/2019"
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---


###package needed

```r
library(tidyverse)
```

```
## ── Attaching packages ──────────────────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ ggplot2 3.1.1       ✔ purrr   0.3.2  
## ✔ tibble  2.1.1       ✔ dplyr   0.8.0.1
## ✔ tidyr   0.8.3       ✔ stringr 1.4.0  
## ✔ readr   1.3.1       ✔ forcats 0.4.0
```

```
## Warning: package 'ggplot2' was built under R version 3.5.2
```

```
## Warning: package 'tibble' was built under R version 3.5.2
```

```
## Warning: package 'tidyr' was built under R version 3.5.2
```

```
## Warning: package 'purrr' was built under R version 3.5.2
```

```
## Warning: package 'dplyr' was built under R version 3.5.2
```

```
## Warning: package 'stringr' was built under R version 3.5.2
```

```
## Warning: package 'forcats' was built under R version 3.5.2
```

```
## ── Conflicts ─────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library(ggdendro)
library(plotly)
```

```
## 
## Attaching package: 'plotly'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     last_plot
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```
## The following object is masked from 'package:graphics':
## 
##     layout
```

```r
library(edgeR)
```

```
## Loading required package: limma
```

```r
library(readxl)
```

```
## Warning: package 'readxl' was built under R version 3.5.2
```

```r
library(goseq) # goseq
```

```
## Loading required package: BiasedUrn
```

```
## Loading required package: geneLenDataBase
```

```
## 
```
# read custom category list

```r
#load("../input/At.hormone.responsive.list.Rdata")
# load custom categories
load(file.path("..","input","hormone.responsiveness6.DF.s.Rdata")) # using SAup_combined & SAdown_combined
#hormone.responsiveness6.DF.s %>% dplyr::select(PIFtarget)
#At.hormone.responsive.list<-read_xlsx(path=file.path("..","input","Supplemental_Dataset3_source_of_custom_categories.xlsx"),sheet=1,skip=15)
#At.hormone.responsive.list %>% dplyr::select(PIFtarget)
# from tpc120857_SupplementalTable1.pdf (Leivar and Monte, 2014)
PIFup<-c("AT5G39860","AT5G15160","AT4G16780","AT2G44910","AT5G47370","AT3G60390","AT1G18400","AT1G69010","AT1G03790","AT1G02340","AT1G14920","AT2G42870","AT2G46970","AT3G62090","AT2G43060","AT5G02200","AT2G37678","AT1G52830","AT1G04250","AT1G09570","AT4G28720","AT1G04180","AT1G70560","AT4G39950","AT4G27260","AT1G70940","AT3G14370","AT5G18010","AT4G37770","AT5G65800","AT1G02400","AT5G07010","AT5G07000","AT1G75450","AT4G10240","AT4G32280","AT3G15540","AT5G04190","AT4G25470")
PIFdown<-c("AT1G06040","AT5G44190","AT4G26150","AT5G56860","AT5G13630","AT1G70700","AT1G05010","AT5G67030","AT1G01060","AT2G46830")
PIFcomplex<-c("AT5G11260","AT5G02760")
At.PIFtarget.list<-list(PIFup=PIFup,PIFdown=PIFdown,PIFcomplex=PIFcomplex)
```
# convert into custom category list
# prep for GOseq with custom category

```r
 # bias.data vector must have the same length as DEgenes vector! 
getwd()
```

```
## [1] "/Volumes/data_work/Data8/NGS_related/Arabidopsis_analysis/SAS_defense_transcriptome_project/SAS_defense_transcriptome"
```

```r
At_cdna<-Biostrings::readDNAStringSet("https://www.arabidopsis.org/download_files/Sequences/TAIR10_blastsets/TAIR10_cdna_20110103_representative_gene_model_updated")
At_cdna
```

```
##   A DNAStringSet instance of length 33602
##         width seq                                      names               
##     [1]  2394 AGAAAACAGTCGACCGTCA...TTGGTAATTTTTTGAGTC AT1G50920.1 | Sym...
##     [2]   546 ATGACTCGTTTGTTGCCTT...GTTGATTCTGGTACATAG AT1G36960.1 | Sym...
##     [3]  2314 ATGGATTCAGAGTCAGAGT...GGTGCATTGTGTTTCTCC AT1G44020.1 | Sym...
##     [4]  1658 TCGTTTCGTCGTCGATCAG...GATTACATGCTACATTTT AT1G15970.1 | Sym...
##     [5]  1453 ATTGAAAAGAAAACACATC...CACCAAAATCTTCTCATA AT1G73440.1 | Sym...
##     ...   ... ...
## [33598]    87 GGATGGATGTCTGAGCGGT...CGAATCCCTCTCCATCCG ATMG00420.1 | Sym...
## [33599]   384 ATGCTCCCCGCCGGTTGTT...CGATACTTAACTATATAA ATMG01330.1 | Sym...
## [33600]   573 ATGGATAACCAATTCATTT...CAGCGTAGCGACGGATAA ATMG00070.1 | Sym...
## [33601]   366 ATGGCATCAAAAATCCGCA...CCTTCTGCATACGCATAA ATMG00130.1 | Sym...
## [33602]    74 GCGCTCTTAGTTCAGTTCG...CAAATCCTACAGAGCGTG ATMG00930.1 | Sym...
```
# format At.PIFtarget.list2 into goseq() compatible list (run once)

```r
# extract AGI from At_cdna
At.gene.name<-tibble(name=names(At_cdna)) %>% separate(name,into=c("name2","Symbol","description"),sep=" \\|",extra="drop",fill="left") %>%  mutate(AGI=str_remove(name2,pattern="\\.[[:digit:]]+")) 
# make presence/absense hormone responsive gene table (dataframe)
  for(i in 1:3) {
  genes<-At.PIFtarget.list[[i]] %>% as_vector()
  temp<-data.frame(AGI=genes,category=names(At.PIFtarget.list)[i])
  #At.gene.name %>% left_join(temp,by=c("name"="AGI")) -> At.gene.name
At.gene.name %>% left_join(temp,by="AGI") -> At.gene.name 
  }
names(At.gene.name)[5:7] <- names(At.PIFtarget.list)
At.gene.name %>% filter(str_detect(AGI,pattern="AT1G|AT2G|AT3G|AT4G|AT5G|ATC|ATM")) %>% unite(category,5:7,sep=",")->test2
test2[1:10,] %>%  mutate(category=gsub("NA","",category)) %>% mutate(category=gsub(",","",category))
# only select genes with AGI name and concatenate categories
At.gene.name %>% filter(str_detect(AGI,pattern="AT1G|AT2G|AT3G|AT4G|AT5G|ATC|ATM")) %>% unite(category,5:7,sep=",") %>%  mutate(category=gsub("NA,","",category)) %>% mutate(category=gsub(",NA","",category)) %>% mutate(category=gsub("NA","",category)) ->test3
# convert into list object
temp.list<-list()
for(i in 1:dim(test3)[1]) {
  temp.list[[i]]<-test3 %>% slice(i) %>% pull(category)
}
names(temp.list)<-test3 %>% pull(AGI)
# split concatanated categories in each gene
At.PIFtarget.list2<-lapply(temp.list,function(x) unlist(strsplit(x, split=",")))
# save
save(At.PIFtarget.list2,file=file.path("..","output","At.PIFtarget.list2.Rdata"))
```


# prep for GOseq with custom category

```r
#  # bias.data vector must have the same length as DEgenes vector! 
# getwd()
# At_cdna<-Biostrings::readDNAStringSet("https://www.arabidopsis.org/download_files/Sequences/TAIR10_blastsets/TAIR10_cdna_20110103_representative_gene_model_updated")
# At_cdna

## if you want to test expressed genes as background, add background in this function
load(file.path("..","output","At.PIFtarget.list2.Rdata"))
GOseq.At.customcat.ORA<-function(genelist,padjust=0.05,custom.category.list=At.hormone.responsive.list,cdna=At_cdna,background="") { # return GO enrichment table, padjus, padjust=0.05. New BLAST2GO results with only BP Brgo.v2.5.BP.list is used (042518). Brgo.v2.5.BP.list is based on Hajar's Blast2go, which should be replaced with Brgo.Atgoslim.BP.list and giving a new name to this function (080818). cdna is DNAstring object. background is a vector of genes.
  #bias<-nchar(cdna)
  bias<-Biostrings::width(cdna)
  # add name to "bias". cdna
  #names(bias)<-names(cdna)
 names(bias) <- tibble(name=names(cdna)) %>% separate(name,into=c("name2","Symbol","description"),sep=" \\|",extra="drop",fill="left") %>%  mutate(AGI=gsub(pattern="\\.[[:digit:]]+", "", name2)) %>% pull(AGI)
 if(background=="") {print("Use all genes in genome as background")} else {
 # select only expressed geens 
  bias <- bias[background$gene_id]
  print("Expessed genes are used as background.")
  print(paste("length of bias is ",length(bias)))
 }
  TF<-(names(bias) %in% genelist)*1
  names(TF)<-names(bias)
  #print(TF)
  pwf<-nullp(TF,bias.data=bias)
  #print(pwf$DEgenes)
  GO.pval <- goseq(pwf,gene2cat=custom.category.list,use_genes_without_cat=TRUE) # format became different in new goseq version (021111). Does not work (042716)
  # calculate p ajust by BH
  GO.pval$over_represented_padjust<-p.adjust(GO.pval$over_represented_pvalue,method="BH")
  if(GO.pval$over_represented_padjust[1]>padjust) return("no enriched GO")
  else {
    enriched.GO<-GO.pval[GO.pval$over_represented_padjust<padjust,] 
    print("enriched.GO is")
    print(enriched.GO)
    return(enriched.GO)
  }
}
```

# Importing enrichment result table and have summary table: Table of DEGs (under construction)

```r
# Ding (2018) Cell DEGs
list.files(path=file.path("..","output","DEGs_with_description"),pattern=".csv") # finding files in different blogs did not work.. Did work this time... why?
```

```
##  [1] "DEgenes.Col.1h.rH.csv"          "DEgenes.Col.49h.rH.csv"        
##  [3] "DEgenes.myc234.1h.rCol.rH.csv"  "DEgenes.myc234.1h.rH.csv"      
##  [5] "DEgenes.myc234.49h.rCol.rH.csv" "DEgenes.myc234.49h.rH.csv"     
##  [7] "DEgenes.npr.1h.rCol.rH.csv"     "DEgenes.npr.1h.rH.csv"         
##  [9] "DEgenes.npr.49h.rCol.rH.csv"    "DEgenes.npr.49h.rH.csv"        
## [11] "DEgenes.sid.1h.rCol.rH.csv"     "DEgenes.sid.1h.rH.csv"         
## [13] "DEgenes.sid.49h.rCol.rH.csv"    "DEgenes.sid.49h.rH.csv"
```

```r
getwd()
```

```
## [1] "/Volumes/data_work/Data8/NGS_related/Arabidopsis_analysis/SAS_defense_transcriptome_project/SAS_defense_transcriptome"
```

```r
DEG.objs<-list.files(path=file.path("..","output","DEGs_with_description"),pattern="\\.csv$") # under construction
DEG.objs
```

```
##  [1] "DEgenes.Col.1h.rH.csv"          "DEgenes.Col.49h.rH.csv"        
##  [3] "DEgenes.myc234.1h.rCol.rH.csv"  "DEgenes.myc234.1h.rH.csv"      
##  [5] "DEgenes.myc234.49h.rCol.rH.csv" "DEgenes.myc234.49h.rH.csv"     
##  [7] "DEgenes.npr.1h.rCol.rH.csv"     "DEgenes.npr.1h.rH.csv"         
##  [9] "DEgenes.npr.49h.rCol.rH.csv"    "DEgenes.npr.49h.rH.csv"        
## [11] "DEgenes.sid.1h.rCol.rH.csv"     "DEgenes.sid.1h.rH.csv"         
## [13] "DEgenes.sid.49h.rCol.rH.csv"    "DEgenes.sid.49h.rH.csv"
```

```r
# read csv file
DEG.list<-lapply(DEG.objs, function(x) read_csv(paste(file.path("..","output","DEGs_with_description"),"/",x,sep="")))
```

```
## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]

## Warning: Missing column names filled in: 'X1' [1]
```

```r
names(DEG.list)<-gsub(".csv","",DEG.objs)
DEG.list
```

```
## $DEgenes.Col.1h.rH
## # A tibble: 151 x 9
##       X1 gene_id   logFC logCPM    LR   PValue      FDR description   name 
##    <dbl> <chr>     <dbl>  <dbl> <dbl>    <dbl>    <dbl> <chr>         <chr>
##  1     1 AT3G154… -1.40    5.34  87.1 1.01e-20 1.75e-16 Aluminium in… <NA> 
##  2     2 AT1G809… -1.13    4.75  61.7 3.95e-15 3.42e-11 Chaperone Dn… AtJ8…
##  3     3 AT5G202… -1.48    4.87  60.8 6.17e-15 3.56e-11 Raffinose sy… DIN1…
##  4     4 AT4G376… -2.32    2.56  60.0 9.45e-15 3.87e-11 BTB and TAZ … BT5  
##  5     5 AT2G158… -0.947   6.05  59.7 1.12e-14 3.87e-11 maternal eff… CBP1…
##  6     6 AT2G076…  4.20    3.16  58.0 2.64e-14 7.64e-11 Cytochrome c… <NA> 
##  7     7 AT4G398…  0.892   6.23  51.6 6.75e-13 1.67e- 9 myo-inositol… ATIP…
##  8     8 AT3G255…  6.69    1.09  47.7 4.86e-12 1.05e- 8 Adenosylmeth… SAMD…
##  9     9 AT1G050… -6.92    1.25  47.2 6.31e-12 1.14e- 8 dentin sialo… <NA> 
## 10    10 AT3G100… -1.00    6.36  47.1 6.59e-12 1.14e- 8 unknown prot… <NA> 
## # … with 141 more rows
## 
## $DEgenes.Col.49h.rH
## # A tibble: 196 x 9
##       X1 gene_id   logFC logCPM    LR   PValue      FDR description   name 
##    <dbl> <chr>     <dbl>  <dbl> <dbl>    <dbl>    <dbl> <chr>         <chr>
##  1     1 AT2G106… -4.80    2.30  89.8 2.66e-21 4.61e-17 transposable… <NA> 
##  2     2 AT5G139…  0.961   6.60  51.8 6.03e-13 5.23e- 9 Chalcone and… ATCH…
##  3     3 AT3G154… -0.950   5.31  43.9 3.51e-11 2.03e- 7 Aluminium in… <NA> 
##  4     4 AT1G757…  1.60    3.42  40.4 2.08e-10 9.01e- 7 GAST1 protei… GASA1
##  5     5 AT1G268… -1.16    4.08  38.5 5.59e-10 1.94e- 6 RING/U-box s… MPSR1
##  6     6 ATCG002… -0.885   7.91  36.9 1.26e- 9 3.65e- 6 photosystem … PSBD 
##  7     7 AT5G541… -1.27    3.63  36.4 1.58e- 9 3.90e- 6 protochlorop… PORA 
##  8     8 AT5G428…  1.68    3.65  35.6 2.41e- 9 5.23e- 6 dihydroflavo… DFR;…
##  9     9 AT5G523… -0.775   6.58  34.2 4.91e- 9 9.45e- 6 low-temperat… COR7…
## 10    10 AT2G425… -0.932   6.89  32.6 1.16e- 8 1.86e- 5 cold regulat… COR1…
## # … with 186 more rows
## 
## $DEgenes.myc234.1h.rCol.rH
## # A tibble: 518 x 10
##       X1 gene_id logFC.gtmyc234 logFC.gtmyc234.… logCPM    LR    PValue
##    <dbl> <chr>            <dbl>            <dbl>  <dbl> <dbl>     <dbl>
##  1     1 AT3G01…         -10.2            1.58     6.74 3576. 0.       
##  2     2 AT4G31…          -5.32          -1.47     4.68  732. 1.36e-159
##  3     3 AT1G16…          -6.48           0.0320   6.06  604. 6.61e-132
##  4     4 AT2G30…          -3.27          -0.153    8.29  585. 9.95e-128
##  5     5 AT3G19…          -6.71           0.480    5.49  559. 3.40e-122
##  6     6 AT3G02…          -6.24          -0.813    4.88  528. 1.73e-115
##  7     7 AT3G58…          -5.35          -0.883    5.03  488. 8.73e-107
##  8     8 AT5G14…          -4.21           0.161    6.08  454. 2.63e- 99
##  9     9 AT4G39…          -5.15          -0.199    4.76  424. 9.44e- 93
## 10    10 AT2G43…          -4.98          -0.443    6.16  418. 1.98e- 91
## # … with 508 more rows, and 3 more variables: FDR <dbl>,
## #   description <chr>, name <chr>
## 
## $DEgenes.myc234.1h.rH
## # A tibble: 75 x 9
##       X1 gene_id   logFC logCPM    LR   PValue      FDR description   name 
##    <dbl> <chr>     <dbl>  <dbl> <dbl>    <dbl>    <dbl> <chr>         <chr>
##  1     1 AT3G154… -1.19   5.06   76.4 2.30e-18 3.99e-14 Aluminium in… <NA> 
##  2     2 AT1G233… -0.972  6.68   72.4 1.72e-17 1.49e-13 Kelch repeat… KFB  
##  3     3 AT5G493… -0.923  4.88   54.1 1.94e-13 1.12e- 9 beta-xylosid… ATBX…
##  4     4 AT5G202… -1.13   4.26   48.3 3.69e-12 1.60e- 8 Raffinose sy… DIN1…
##  5     5 AT2G158… -0.831  5.10   41.2 1.39e-10 4.80e- 7 Leucine-rich… LRX10
##  6     6 AT5G219… -1.15   3.76   40.3 2.21e-10 6.38e- 7 unknown prot… <NA> 
##  7     7 AT5G615… -0.953  4.27   39.7 2.92e-10 7.23e- 7 Integrase-ty… DEWA…
##  8     8 AT2G259… -0.530  7.38   36.7 1.41e- 9 3.06e- 6 Zinc finger … ATCT…
##  9     9 AT4G109…  5.84   0.544  33.4 7.49e- 9 1.27e- 5 UDP-D-glucos… UGE5 
## 10    10 AT2G158… -0.739  5.75   33.3 7.73e- 9 1.27e- 5 maternal eff… CBP1…
## # … with 65 more rows
## 
## $DEgenes.myc234.49h.rCol.rH
## # A tibble: 893 x 10
##       X1 gene_id logFC.gtmyc234 logFC.gtmyc234.… logCPM    LR    PValue
##    <dbl> <chr>            <dbl>            <dbl>  <dbl> <dbl>     <dbl>
##  1     1 AT3G01…          -8.38           0.449    6.28 4591. 0.       
##  2     2 AT1G16…          -7.63          -0.0606   6.85 1179. 1.09e-256
##  3     3 AT2G43…          -5.17          -0.612    6.88 1160. 1.02e-252
##  4     4 AT5G23…          -4.02          -1.00     7.05 1015. 4.02e-221
##  5     5 AT3G19…          -7.23           0.661    5.75  982. 7.28e-214
##  6     6 AT5G14…          -4.40          -0.659    6.63  921. 1.01e-200
##  7     7 AT4G13…          -4.49          -0.486    6.63  877. 4.37e-191
##  8     8 AT2G30…          -3.26          -0.454    8.47  855. 2.28e-186
##  9     9 AT3G02…          -6.70           0.620    5.18  824. 1.22e-179
## 10    10 AT3G58…          -6.73          -0.496    5.52  723. 9.17e-158
## # … with 883 more rows, and 3 more variables: FDR <dbl>,
## #   description <chr>, name <chr>
## 
## $DEgenes.myc234.49h.rH
## # A tibble: 195 x 9
##       X1 gene_id   logFC logCPM    LR   PValue      FDR description   name 
##    <dbl> <chr>     <dbl>  <dbl> <dbl>    <dbl>    <dbl> <chr>         <chr>
##  1     1 AT5G366… -8.17   2.07  124.  7.49e-29 1.30e-24 RING/FYVE/PH… <NA> 
##  2     2 AT1G299… -4.86   1.98   86.8 1.22e-20 1.06e-16 conserved pe… <NA> 
##  3     3 AT1G233… -0.705  6.92   75.1 4.47e-18 2.58e-14 Kelch repeat… KFB  
##  4     4 AT5G404… -0.613  8.26   57.3 3.71e-14 1.61e-10 unknown prot… RBB1 
##  5     5 AT3G100… -0.632  6.19   47.3 6.07e-12 2.10e- 8 unknown prot… <NA> 
##  6     6 AT5G139…  0.831  6.52   45.6 1.46e-11 4.21e- 8 Chalcone and… ATCH…
##  7     7 AT5G505… -6.76   0.838  44.1 3.10e-11 7.67e- 8 SUMO-activat… AT-S…
##  8     8 AT4G109… -6.70   0.832  43.6 3.95e-11 8.56e- 8 UDP-D-glucos… UGE5 
##  9     9 AT5G497… -0.557  7.26   42.8 6.21e-11 1.20e- 7 ferric reduc… ATFR…
## 10    10 AT1G180… -6.71   0.867  42.1 8.81e-11 1.53e- 7 Major facili… <NA> 
## # … with 185 more rows
## 
## $DEgenes.npr.1h.rCol.rH
## # A tibble: 1,581 x 10
##       X1 gene_id logFC.gtnpr logFC.gtnpr.trtL logCPM    LR    PValue
##    <dbl> <chr>         <dbl>            <dbl>  <dbl> <dbl>     <dbl>
##  1     1 AT2G48…       -7.06          -0.900    6.02 1531. 0.       
##  2     2 AT3G01…       -6.58          -1.56     6.74 2983. 0.       
##  3     3 AT3G42…        8.22           0.156    4.29 1100. 1.20e-239
##  4     4 AT1G63…       -1.70          -0.0279   7.26  514. 3.02e-112
##  5     5 AT5G35…        6.95          -0.990    2.67  454. 2.28e- 99
##  6     6 AT3G55…       -7.02           0.132    3.38  428. 1.28e- 93
##  7     7 AT1G53…       -3.04           0.253    5.26  332. 6.64e- 73
##  8     8 AT5G46…        2.71           0.373    3.75  306. 3.11e- 67
##  9     9 AT5G28…        6.43           0.954    2.33  253. 9.21e- 56
## 10    10 AT1G18…        2.76          -0.245    4.90  238. 1.66e- 52
## # … with 1,571 more rows, and 3 more variables: FDR <dbl>,
## #   description <chr>, name <chr>
## 
## $DEgenes.npr.1h.rH
## # A tibble: 60 x 9
##       X1 gene_id   logFC logCPM    LR   PValue      FDR description   name 
##    <dbl> <chr>     <dbl>  <dbl> <dbl>    <dbl>    <dbl> <chr>         <chr>
##  1     1 ATCG010…  2.35    3.52  74.5 6.18e-18 1.07e-13 NADH:ubiquin… NDHG 
##  2     2 ATCG010…  1.74    4.49  56.9 4.69e-14 4.07e-10 NADH-Ubiquin… NDHD 
##  3     3 AT1G233… -1.07    5.82  54.1 1.90e-13 1.10e- 9 Kelch repeat… KFB  
##  4     4 AT2G076…  2.02    3.30  45.8 1.32e-11 5.70e- 8 Cytochrome c… <NA> 
##  5     5 AT2G158… -1.51    4.90  41.5 1.19e-10 4.12e- 7 Leucine-rich… LRX10
##  6     6 AT5G229… -1.30    4.29  39.6 3.19e-10 9.20e- 7 CHY-type/CTC… AtRZ…
##  7     7 AT5G247…  1.03    4.88  35.9 2.06e- 9 4.92e- 6 vegetative s… ATVS…
##  8     8 AT5G202… -1.36    4.25  35.7 2.27e- 9 4.92e- 6 Raffinose sy… DIN1…
##  9     9 AT3G154… -1.56    5.19  33.9 5.71e- 9 1.10e- 5 Aluminium in… <NA> 
## 10    10 AT3G197…  0.767   6.03  31.3 2.19e- 8 3.49e- 5 branched-cha… BCAT4
## # … with 50 more rows
## 
## $DEgenes.npr.49h.rCol.rH
## # A tibble: 2,904 x 10
##       X1 gene_id logFC.gtnpr logFC.gtnpr.trtL logCPM    LR    PValue
##    <dbl> <chr>         <dbl>            <dbl>  <dbl> <dbl>     <dbl>
##  1     1 AT3G01…       -7.52           0.863    6.28 4092. 0.       
##  2     2 AT2G48…       -7.06           1.12     5.99 2573. 0.       
##  3     3 AT3G42…        6.84           4.01     3.82 1027. 1.03e-223
##  4     4 AT1G63…       -1.78          -0.410    7.12  728. 9.82e-159
##  5     5 AT1G53…       -3.86           1.80     4.76  449. 2.84e- 98
##  6     6 AT3G55…       -5.74          -0.542    3.18  392. 9.31e- 86
##  7     7 AT3G14…        1.57          -0.0471   7.71  316. 2.75e- 69
##  8     8 AT5G35…        5.77           0.893    2.62  245. 7.38e- 54
##  9     9 AT5G15…        2.49           0.501    3.62  212. 1.15e- 46
## 10    10 AT5G46…        2.76          -0.424    3.69  196. 2.15e- 43
## # … with 2,894 more rows, and 3 more variables: FDR <dbl>,
## #   description <chr>, name <chr>
## 
## $DEgenes.npr.49h.rH
## # A tibble: 69 x 9
##       X1 gene_id   logFC logCPM    LR   PValue      FDR description   name 
##    <dbl> <chr>     <dbl>  <dbl> <dbl>    <dbl>    <dbl> <chr>         <chr>
##  1     1 AT5G506…  7.56   1.68   71.6 2.67e-17 4.63e-13 Major facili… <NA> 
##  2     2 AT2G400…  5.66   2.19   63.8 1.41e-15 1.22e-11 unknown prot… <NA> 
##  3     3 AT2G052… -3.11   2.82   49.9 1.64e-12 9.49e- 9 DNAJ heat sh… <NA> 
##  4     4 AT2G431…  0.807  7.40   47.4 5.83e-12 2.52e- 8 isopropylmal… ATLE…
##  5     5 AT1G233… -0.943  6.29   31.2 2.37e- 8 8.21e- 5 Kelch repeat… KFB  
##  6     6 AT1G056… -1.55   3.50   29.9 4.47e- 8 1.12e- 4 Uridine diph… UGT7…
##  7     7 AT5G078… -3.62   1.90   29.9 4.51e- 8 1.12e- 4 conserved pe… <NA> 
##  8     8 AT1G136… -0.747  6.53   29.4 5.83e- 8 1.21e- 4 BEST Arabido… <NA> 
##  9     9 AT5G148…  3.67   0.952  29.3 6.30e- 8 1.21e- 4 transposable… <NA> 
## 10    10 AT5G159…  2.25   2.05   28.5 9.24e- 8 1.60e- 4 Adenosylmeth… SAMD…
## # … with 59 more rows
## 
## $DEgenes.sid.1h.rCol.rH
## # A tibble: 25 x 10
##       X1 gene_id logFC.gtsid logFC.gtsid.trtL logCPM     LR   PValue
##    <dbl> <chr>         <dbl>            <dbl>  <dbl>  <dbl>    <dbl>
##  1     1 AT3G01…      -8.09           -0.747   6.74  3308.  0.      
##  2     2 AT5G28…       5.30           -0.0535  2.33   118.  2.57e-26
##  3     3 AT3G14…       0.999          -0.0102  7.26   106.  7.55e-24
##  4     4 AT4G21…      -6.26           -0.588   0.442   83.2 8.66e-19
##  5     5 AT5G50…      -0.900          -0.178   4.71    71.8 2.55e-16
##  6     6 AT3G29…      -1.67           -0.182   2.00    42.5 5.77e-10
##  7     7 AT5G50…      -1.30            0.474   4.41    35.8 1.69e- 8
##  8     8 AT5G50…      -1.64            0.679   3.79    35.2 2.33e- 8
##  9     9 AT4G13…      -0.703           0.336   7.56    34.0 4.13e- 8
## 10    10 AT1G48…       0.800           0.0330  6.13    32.0 1.11e- 7
## # … with 15 more rows, and 3 more variables: FDR <dbl>, description <chr>,
## #   name <chr>
## 
## $DEgenes.sid.1h.rH
## # A tibble: 42 x 9
##       X1 gene_id   logFC logCPM    LR   PValue      FDR description   name 
##    <dbl> <chr>     <dbl>  <dbl> <dbl>    <dbl>    <dbl> <chr>         <chr>
##  1     1 AT2G249… -2.42   5.86  121.  3.25e-28 5.64e-24 unknown prot… <NA> 
##  2     2 AT5G366…  8.05   2.27  114.  1.10e-26 9.53e-23 RING/FYVE/PH… <NA> 
##  3     3 AT1G299…  6.14   2.11   51.0 9.47e-13 5.47e- 9 conserved pe… <NA> 
##  4     4 AT3G154… -1.13   5.84   41.2 1.40e-10 6.05e- 7 Aluminium in… <NA> 
##  5     5 AT5G505… -6.42   1.19   37.1 1.12e- 9 3.87e- 6 Major facili… <NA> 
##  6     6 AT2G338… -1.37   4.33   34.8 3.65e- 9 1.05e- 5 Dormancy/aux… AtDR…
##  7     7 AT2G158… -1.05   5.54   31.4 2.14e- 8 5.30e- 5 Leucine-rich… LRX10
##  8     8 AT5G159…  4.19   0.837  28.0 1.23e- 7 2.67e- 4 Adenosylmeth… SAMD…
##  9     9 ATCG010…  1.46   3.52   25.3 4.79e- 7 8.45e- 4 NADH:ubiquin… NDHG 
## 10    10 AT1G809… -0.870  4.96   25.3 4.88e- 7 8.45e- 4 Chaperone Dn… AtJ8…
## # … with 32 more rows
## 
## $DEgenes.sid.49h.rCol.rH
## # A tibble: 25 x 10
##       X1 gene_id logFC.gtsid logFC.gtsid.trtL  logCPM     LR   PValue
##    <dbl> <chr>         <dbl>            <dbl>   <dbl>  <dbl>    <dbl>
##  1     1 AT3G01…      -7.91           -0.647   6.28   4135.  0.      
##  2     2 AT3G14…       1.02            0.0773  7.71    148.  6.23e-33
##  3     3 AT5G28…       7.53           -2.91    2.19     69.3 8.81e-16
##  4     4 AT3G56…      -0.844           0.368   4.71     59.1 1.46e-13
##  5     5 AT5G50…      -0.972           0.0828  4.54     56.0 6.94e-13
##  6     6 AT4G21…      -4.70           -1.13   -0.0786   51.4 7.02e-12
##  7     7 AT4G04…      -2.45           -4.00    2.29     44.4 2.30e-10
##  8     8 AT4G07…      -0.687          -0.403   3.09     39.1 3.19e- 9
##  9     9 AT3G55…      -0.882           0.169   4.66     34.7 2.92e- 8
## 10    10 AT1G48…       0.671          -0.0609  6.75     29.1 4.82e- 7
## # … with 15 more rows, and 3 more variables: FDR <dbl>, description <chr>,
## #   name <chr>
## 
## $DEgenes.sid.49h.rH
## # A tibble: 116 x 9
##       X1 gene_id   logFC logCPM    LR   PValue      FDR description   name 
##    <dbl> <chr>     <dbl>  <dbl> <dbl>    <dbl>    <dbl> <chr>         <chr>
##  1     1 AT2G055… -1.90   3.53   73.9 8.39e-18 1.45e-13 Glycine-rich… <NA> 
##  2     2 AT5G505…  6.91   1.61   56.8 4.88e-14 4.23e-10 Major facili… <NA> 
##  3     3 AT1G757…  1.27   3.84   45.2 1.77e-11 1.02e- 7 GAST1 protei… GASA1
##  4     4 AT4G228…  0.983  5.14   37.6 8.62e-10 3.74e- 6 leucoanthocy… ANS;…
##  5     5 AT5G139…  0.765  6.74   35.0 3.23e- 9 1.02e- 5 Chalcone and… ATCH…
##  6     6 AT1G247… -5.76   0.534  34.9 3.52e- 9 1.02e- 5 UDP-3-O-acyl… AtLp…
##  7     7 AT3G575…  0.746  5.78   33.8 6.14e- 9 1.52e- 5 seed imbibit… AtSI…
##  8     8 AT3G228…  1.49   3.01   33.4 7.53e- 9 1.63e- 5 Chlorophyll … ELIP…
##  9     9 AT2G423… -0.723  5.22   29.3 6.05e- 8 1.17e- 4 Basic-leucin… ATBZ…
## 10    10 AT3G100… -0.667  6.41   28.4 9.78e- 8 1.69e- 4 unknown prot… <NA> 
## # … with 106 more rows
```

```r
# GOseq.Atgoslim.BP.list.ORA(genelist=genes.shade1h.up$AGI)
```

# overlap table

```r
# expressed genes... what is the best way?
#expressed.genes<- s
#  Col.49h %>% dplyr::select(gene_id) # only 151 genes?
#  gsub("(.)(\\.[[:digit:]]+)","\\1",as_vector(dge.nolow$genes)) # needs to be modified

# testingGOseq.At.customcat 
GOseq.At.customcat.ORA(genelist=At.PIFtarget.list[[1]],custom.category.list=At.PIFtarget.list2)
```

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

```
## Using manually entered categories.
```

```
## Calculating the p-values...
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```
## [1] "enriched.GO is"
##   category over_represented_pvalue under_represented_pvalue numDEInCat
## 3    PIFup                       0                        1         39
##   numInCat over_represented_padjust
## 3       39                        0
```

```
##   category over_represented_pvalue under_represented_pvalue numDEInCat
## 3    PIFup                       0                        1         39
##   numInCat over_represented_padjust
## 3       39                        0
```
# GO.ORA for each DEG list (loop)

```r
names(DEG.list)
```

```
##  [1] "DEgenes.Col.1h.rH"          "DEgenes.Col.49h.rH"        
##  [3] "DEgenes.myc234.1h.rCol.rH"  "DEgenes.myc234.1h.rH"      
##  [5] "DEgenes.myc234.49h.rCol.rH" "DEgenes.myc234.49h.rH"     
##  [7] "DEgenes.npr.1h.rCol.rH"     "DEgenes.npr.1h.rH"         
##  [9] "DEgenes.npr.49h.rCol.rH"    "DEgenes.npr.49h.rH"        
## [11] "DEgenes.sid.1h.rCol.rH"     "DEgenes.sid.1h.rH"         
## [13] "DEgenes.sid.49h.rCol.rH"    "DEgenes.sid.49h.rH"
```

```r
# make special output directory
dir.create(path=file.path("..","output","doubleCoef"))
```

```
## Warning in dir.create(path = file.path("..", "output", "doubleCoef")): '../
## output/doubleCoef' already exists
```

```r
#
for(n in 1:14) {
  genelist.all<-DEG.list[[n]]
# DEG.list[[n]]  %>%  dplyr::filter((as.name(names(DEG.list[[n]])[2]))>0)  # does not work
genelist.upup<-DEG.list[[n]] %>% dplyr::filter((UQ(as.name(names(DEG.list[[n]])[2])))>0) %>% dplyr::filter((UQ(as.name(names(DEG.list[[n]])[3])))>0)
    genelist.updown<-DEG.list[[n]] %>% dplyr::filter((UQ(as.name(names(DEG.list[[n]])[2])))>0) %>% dplyr::filter((UQ(as.name(names(DEG.list[[n]])[3])))<0)
        genelist.downup<-DEG.list[[n]] %>% dplyr::filter((UQ(as.name(names(DEG.list[[n]])[2])))<0) %>% dplyr::filter((UQ(as.name(names(DEG.list[[n]])[3])))>0)
    genelist.downdown<-  DEG.list[[n]] %>% dplyr::filter((UQ(as.name(names(DEG.list[[n]])[2])))<0) %>% dplyr::filter((UQ(as.name(names(DEG.list[[n]])[3])))<0)
#    if(as_vector(genelist.all[,"gene_id"]) ***) { # under construciton. Do not use this form this time.
# if genelist.all[,"gene_id"] has "AT***.1" pattern (with splicing variants)
#GO.ORA.temp.all<-GOseq.At.customcat.ORA(str_remove(as_vector(genelist.all[,"gene_id"]),"\\.[[:digit:]]+"),custom.category.list=At.PIFtarget.list2)
#    GO.ORA.temp.up_up<-GOseq.At.customcat.ORA(str_remove(as_vector(genelist.upup[,"gene_id"]),"\\.[[:digit:]]+"),custom.category.list=At.PIFtarget.list2)
#    GO.ORA.temp.up_down<-GOseq.At.customcat.ORA(str_remove(as_vector(genelist.updown[,"gene_id"]),"\\.[[:digit:]]+"),custom.category.list=At.PIFtarget.list2)
#    GO.ORA.temp.down_up<-GOseq.At.customcat.ORA(str_remove(as_vector(genelist.downup[,"gene_id"]),"\\.[[:digit:]]+"),custom.category.list=At.PIFtarget.list2)
#    GO.ORA.temp.down_down<-GOseq.At.customcat.ORA(str_remove(as_vector(genelist.downdown[,"gene_id"]),"\\.[[:digit:]]+"),custom.category.list=At.PIFtarget.list2)
#    } else {
# if genelist.all[,"gene_id"] has no spricing variant version (only AGI name)
    GO.ORA.temp.all<-GOseq.At.customcat.ORA(as_vector(genelist.all[,"gene_id"]),custom.category.list=At.PIFtarget.list2)
    GO.ORA.temp.up_up<-GOseq.At.customcat.ORA(as_vector(genelist.upup[,"gene_id"]),custom.category.list=At.PIFtarget.list2)
    GO.ORA.temp.up_down<-GOseq.At.customcat.ORA(as_vector(genelist.updown[,"gene_id"]),custom.category.list=At.PIFtarget.list2)
    GO.ORA.temp.down_up<-GOseq.At.customcat.ORA(as_vector(genelist.downup[,"gene_id"]),custom.category.list=At.PIFtarget.list2)
    GO.ORA.temp.down_down<-GOseq.At.customcat.ORA(as_vector(genelist.downdown[,"gene_id"]),custom.category.list=At.PIFtarget.list2)
 #   }
    
    # handling "no enriched GO" 
    # genelist.names<-c("GO.ORA.temp.up_down","GO.ORA.temp.down_up") # test
    x<-list(GO.ORA.temp.all=GO.ORA.temp.all,
      GO.ORA.temp.up_up=GO.ORA.temp.up_up,
            GO.ORA.temp.up_down=GO.ORA.temp.up_down,
            GO.ORA.temp.down_up=GO.ORA.temp.down_up,
            GO.ORA.temp.down_down=GO.ORA.temp.down_down) # list
    print(x[x=="no enriched GO"])
    x<-x[!x=="no enriched GO"] # reove "no enriched GO" result
    ## add sample info and FC info and save GO.ORA result (this is not working well; 072819)
    ## avoid error with x is empty
    if(length(x)==0) {next} else {
  for (i in 1:length(x)) {
      GO.ORA.result<-x[[i]] %>% mutate(FC = gsub("(GO.ORA.temp.)(.)","\\2",names(x)[i]),sample=DEG.objs[n])
    save(GO.ORA.result,file=file.path("..","output","doubleCoef",paste(gsub(".csv","",DEG.objs[n]),gsub("(GO.ORA.temp.)(.)","\\2",names(x)[i]),"enrich.Rdata",sep=".")))
    rm(GO.ORA.result)
  }
    } 
}
```

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

```
## [1] "enriched.GO is"
##   category over_represented_pvalue under_represented_pvalue numDEInCat
## 3    PIFup             1.22375e-06                        1          5
##   numInCat over_represented_padjust
## 3       39             3.671249e-06
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```
## [1] "enriched.GO is"
##   category over_represented_pvalue under_represented_pvalue numDEInCat
## 3    PIFup            5.209622e-09                        1          5
##   numInCat over_represented_padjust
## 3       39             1.562887e-08
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-4.png)<!-- -->

```
## [1] "enriched.GO is"
##     category over_represented_pvalue under_represented_pvalue numDEInCat
## 1 PIFcomplex              0.01227929                0.9999623          1
##   numInCat over_represented_padjust
## 1        2               0.03683786
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-5.png)<!-- -->

```
## [1] "enriched.GO is"
##     category over_represented_pvalue under_represented_pvalue numDEInCat
## 1 PIFcomplex             0.003782456                0.9999965          1
##   numInCat over_represented_padjust
## 1        2               0.01134737
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-6.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-7.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-8.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-9.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-10.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-11.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.all
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-12.png)<!-- -->![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-13.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-14.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-15.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.all
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-16.png)<!-- -->

```
## [1] "enriched.GO is"
##   category over_represented_pvalue under_represented_pvalue numDEInCat
## 3    PIFup            0.0007763874                0.9998955          6
## 2  PIFdown            0.0031048381                0.9998282          3
##   numInCat over_represented_padjust
## 3       39              0.002329162
## 2       10              0.004657257
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-17.png)<!-- -->

```
## [1] "enriched.GO is"
##     category over_represented_pvalue under_represented_pvalue numDEInCat
## 2    PIFdown             0.006068713                0.9998067          2
## 1 PIFcomplex             0.015796101                0.9999373          1
##   numInCat over_represented_padjust
## 2       10               0.01820614
## 1        2               0.02369415
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-18.png)<!-- -->

```
## [1] "enriched.GO is"
##   category over_represented_pvalue under_represented_pvalue numDEInCat
## 3    PIFup            0.0006628999                0.9999325          5
##   numInCat over_represented_padjust
## 3       39                0.0019887
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-19.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-20.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-21.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-22.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-23.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.all
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-24.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-25.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-26.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-27.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.all
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-28.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-29.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-30.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-31.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.all
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-32.png)<!-- -->

```
## [1] "enriched.GO is"
##   category over_represented_pvalue under_represented_pvalue numDEInCat
## 2  PIFdown              0.01395899                0.9981697          4
##   numInCat over_represented_padjust
## 2       10               0.04187696
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-33.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-34.png)<!-- -->

```
## [1] "enriched.GO is"
##   category over_represented_pvalue under_represented_pvalue numDEInCat
## 2  PIFdown               0.0140555                0.9986469          3
##   numInCat over_represented_padjust
## 2       10               0.04216649
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-35.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.up_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-36.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-37.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-38.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-39.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.all
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-40.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-41.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-42.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-43.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.all
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-44.png)<!-- -->

```
## [1] "enriched.GO is"
##   category over_represented_pvalue under_represented_pvalue numDEInCat
## 3    PIFup             0.001157632                0.9999826          2
##   numInCat over_represented_padjust
## 3       39              0.003472897
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-45.png)<!-- -->![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-46.png)<!-- -->

```
## [1] "enriched.GO is"
##   category over_represented_pvalue under_represented_pvalue numDEInCat
## 3    PIFup            0.0002804182                 0.999998          2
##   numInCat over_represented_padjust
## 3       39             0.0008412547
## [1] "Use all genes in genome as background"
```

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-47.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-48.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-49.png)<!-- -->![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-50.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-51.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-52.png)<!-- -->

```
## $GO.ORA.temp.all
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
## 
## [1] "Use all genes in genome as background"
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-53.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## [1] "Use all genes in genome as background"
```

```
## Warning in pcls(G): initial point very close to some inequality constraints
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-54.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-55.png)<!-- -->

```
## [1] "Use all genes in genome as background"
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

![](MYC234_PIF4target_files/figure-html/unnamed-chunk-8-56.png)<!-- -->

```
## $GO.ORA.temp.all
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.up_down
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_up
## [1] "no enriched GO"
## 
## $GO.ORA.temp.down_down
## [1] "no enriched GO"
```

# read GOseq result table

```r
eGOseqs<-list.files(pattern="enrich.Rdata",path=file.path("..","output","doubleCoef"))
eGOseqs.list2<-sapply(file.path("..","output","doubleCoef",eGOseqs),function(x) mget(load(x))) # mget will return the value of the object(or objects) in a list. see https://stackoverflow.com/questions/29398630/load-data-frames-into-list
#names(eGOseqs.list2)
eGOseqs.list2.summary<-do.call("rbind",eGOseqs.list2) 
#head(eGOseqs.list2.summary) # make sure those are file names
rownames(eGOseqs.list2.summary)<-1:nrow(eGOseqs.list2.summary)
View(eGOseqs.list2.summary)
# write in csv file
write_csv(eGOseqs.list2.summary,path=file.path("..","output","doubleCoef","PIFtarget.eGOseqs.list2.summary.csv"))
```

# plotting expression pattern

```r
###import dge.data and calculate cpm
# fo Kazu
 if(Sys.info()["user"]=="nozue") { 
load(file.path("~","Box","Chunmei_Myc_Paper","manuscripts","output/dge.Rdata"))
 }
```

```
## Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file
## '/Users/nozue/Box/Chunmei_Myc_Paper/manuscripts/output/dge.Rdata', probable
## reason 'No such file or directory'
```

```
## Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
```

```r
# for Chunmei
 if(Sys.info()["user"]=="LCM") { 
load(file.path("C:","Users","LCM","Box","Chunmei_Myc_Paper","manuscripts","output/dge.Rdata"))
 }

log2_cpm <- log2(cpm(dge.data)+1) # copied and pasted from "combined_scripts_for_figures.Rmd". Why log2(cpm(dge.data)+1)? not log2(cpm(dge.data))? (052319)
```

```
## Error in cpm(dge.data): object 'dge.data' not found
```

```r
head(log2_cpm)
```

```
## Error in head(log2_cpm): object 'log2_cpm' not found
```

```r
rownames(log2_cpm) <- gsub("\\.\\d", "", rownames(log2_cpm)) #remove ".digital" after each AGI number
```

```
## Error in rownames(log2_cpm): object 'log2_cpm' not found
```

```r
#subset to keep Col and myc234 data

log2_cpm <- log2_cpm[, 1:32]
```

```
## Error in eval(expr, envir, enclos): object 'log2_cpm' not found
```

```r
head(log2_cpm)
```

```
## Error in head(log2_cpm): object 'log2_cpm' not found
```

```r
#
cpm.wide <- bind_cols(tibble(gene_id=rownames(log2_cpm)),as_tibble(log2_cpm)) 
```

```
## Error in rownames(log2_cpm): object 'log2_cpm' not found
```

```r
cpm.wide
```

```
## Error in eval(expr, envir, enclos): object 'cpm.wide' not found
```
# plot ("up" in myc324 (logFC.gtmyc234>0))

```r
names(DEG.list)
```

```
##  [1] "DEgenes.Col.1h.rH"          "DEgenes.Col.49h.rH"        
##  [3] "DEgenes.myc234.1h.rCol.rH"  "DEgenes.myc234.1h.rH"      
##  [5] "DEgenes.myc234.49h.rCol.rH" "DEgenes.myc234.49h.rH"     
##  [7] "DEgenes.npr.1h.rCol.rH"     "DEgenes.npr.1h.rH"         
##  [9] "DEgenes.npr.49h.rCol.rH"    "DEgenes.npr.49h.rH"        
## [11] "DEgenes.sid.1h.rCol.rH"     "DEgenes.sid.1h.rH"         
## [13] "DEgenes.sid.49h.rCol.rH"    "DEgenes.sid.49h.rH"
```

```r
# npr1-1 misregulated genes 
plot.data<-DEG.list[["DEgenes.myc234.49h.rCol.rH"]] %>% left_join(cpm.wide,by="gene_id") %>% unite(AGI_desc,c("gene_id","name")) %>% dplyr::select(-LR,-PValue,-X1,-logCPM) 
```

```
## Error in tbl_vars(y): object 'cpm.wide' not found
```

```r
# plot: use labeller = label_wrap_gen(width=10) for multiple line in each gene name
p.up<-plot.data %>% filter(logFC.gtmyc234>0 & FDR<1e-10) %>% dplyr::select(-logFC.gtmyc234,-logFC.gtmyc234.trtL,-FDR,-description) %>%  gather(sample,value,-1) %>%
  separate(sample, into = c("genotype", "time_point", "treatment", "rep"), sep = "_") %>% ggplot(aes(x=fct_relevel(treatment,"H"),y=value,color=genotype,shape=rep)) + geom_jitter() + facet_wrap(AGI_desc~.,scale="free",ncol=5,labeller = label_wrap_gen(width=30))+theme(strip.text = element_text(size=6)) + labs(title="up (logFC.gtmyc234>0)")
```

```
## Error in eval(lhs, parent, parent): object 'plot.data' not found
```

```r
p.up
```

```
## Error in eval(expr, envir, enclos): object 'p.up' not found
```
# plot ("down" in myc324 (logFC.gtmyc234>0))

```r
names(DEG.list)
```

```
##  [1] "DEgenes.Col.1h.rH"          "DEgenes.Col.49h.rH"        
##  [3] "DEgenes.myc234.1h.rCol.rH"  "DEgenes.myc234.1h.rH"      
##  [5] "DEgenes.myc234.49h.rCol.rH" "DEgenes.myc234.49h.rH"     
##  [7] "DEgenes.npr.1h.rCol.rH"     "DEgenes.npr.1h.rH"         
##  [9] "DEgenes.npr.49h.rCol.rH"    "DEgenes.npr.49h.rH"        
## [11] "DEgenes.sid.1h.rCol.rH"     "DEgenes.sid.1h.rH"         
## [13] "DEgenes.sid.49h.rCol.rH"    "DEgenes.sid.49h.rH"
```

```r
# npr1-1 misregulated genes 
plot.data<-DEG.list[["DEgenes.myc234.49h.rCol.rH"]] %>% left_join(cpm.wide,by="gene_id") %>% unite(AGI_desc,c("gene_id","name")) %>% dplyr::select(-LR,-PValue,-X1,-logCPM) 
```

```
## Error in tbl_vars(y): object 'cpm.wide' not found
```

```r
# plot: use labeller = label_wrap_gen(width=10) for multiple line in each gene name
p.down<-plot.data %>% filter(logFC.gtmyc234<0 & FDR<1e-100) %>% dplyr::select(-logFC.gtmyc234,-logFC.gtmyc234.trtL,-FDR,-description) %>%  gather(sample,value,-1) %>%
  separate(sample, into = c("genotype", "time_point", "treatment", "rep"), sep = "_") %>% ggplot(aes(x=fct_relevel(treatment,"H"),y=value,color=genotype,shape=rep)) + geom_jitter() + facet_wrap(AGI_desc~.,scale="free",ncol=5,labeller = label_wrap_gen(width=30))+theme(strip.text = element_text(size=6)) + labs(title="down (logFC.gtmyc234<0)")
```

```
## Error in eval(lhs, parent, parent): object 'plot.data' not found
```

```r
p.down
```

```
## Error in eval(expr, envir, enclos): object 'p.down' not found
```


# Session info

```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS  10.14.5
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] goseq_1.32.0           geneLenDataBase_1.16.0 BiasedUrn_1.07        
##  [4] readxl_1.3.1           edgeR_3.22.5           limma_3.36.5          
##  [7] plotly_4.8.0           ggdendro_0.1-20        forcats_0.4.0         
## [10] stringr_1.4.0          dplyr_0.8.0.1          purrr_0.3.2           
## [13] readr_1.3.1            tidyr_0.8.3            tibble_2.1.1          
## [16] ggplot2_3.1.1          tidyverse_1.2.1       
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-137                bitops_1.0-6               
##  [3] matrixStats_0.54.0          lubridate_1.7.4            
##  [5] bit64_0.9-7                 progress_1.2.0             
##  [7] httr_1.4.0                  GenomeInfoDb_1.16.0        
##  [9] tools_3.5.1                 backports_1.1.4            
## [11] utf8_1.1.4                  R6_2.4.0                   
## [13] mgcv_1.8-25                 DBI_1.0.0                  
## [15] lazyeval_0.2.2              BiocGenerics_0.26.0        
## [17] colorspace_1.4-1            withr_2.1.2                
## [19] tidyselect_0.2.5            prettyunits_1.0.2          
## [21] bit_1.1-14                  compiler_3.5.1             
## [23] cli_1.1.0                   rvest_0.3.2                
## [25] Biobase_2.40.0              xml2_1.2.0                 
## [27] DelayedArray_0.6.6          rtracklayer_1.40.6         
## [29] scales_1.0.0                digest_0.6.18              
## [31] Rsamtools_1.32.3            rmarkdown_1.12             
## [33] XVector_0.20.0              pkgconfig_2.0.2            
## [35] htmltools_0.3.6             htmlwidgets_1.3            
## [37] rlang_0.3.4                 rstudioapi_0.10            
## [39] RSQLite_2.1.1               generics_0.0.2             
## [41] jsonlite_1.6                BiocParallel_1.14.2        
## [43] RCurl_1.95-4.12             magrittr_1.5               
## [45] GO.db_3.6.0                 GenomeInfoDbData_1.1.0     
## [47] Matrix_1.2-17               fansi_0.4.0                
## [49] Rcpp_1.0.1                  munsell_0.5.0              
## [51] S4Vectors_0.18.3            stringi_1.4.3              
## [53] yaml_2.2.0                  MASS_7.3-51.3              
## [55] SummarizedExperiment_1.10.1 zlibbioc_1.26.0            
## [57] plyr_1.8.4                  grid_3.5.1                 
## [59] blob_1.1.1                  parallel_3.5.1             
## [61] crayon_1.3.4                lattice_0.20-38            
## [63] Biostrings_2.48.0           haven_2.1.0                
## [65] GenomicFeatures_1.32.3      hms_0.4.2                  
## [67] locfit_1.5-9.1              knitr_1.22                 
## [69] pillar_1.3.1                GenomicRanges_1.32.7       
## [71] biomaRt_2.36.1              stats4_3.5.1               
## [73] XML_3.98-1.19               glue_1.3.1                 
## [75] evaluate_0.13               data.table_1.12.0          
## [77] modelr_0.1.4                cellranger_1.1.0           
## [79] gtable_0.3.0                assertthat_0.2.1           
## [81] xfun_0.5                    broom_0.5.1                
## [83] viridisLite_0.3.0           GenomicAlignments_1.16.0   
## [85] AnnotationDbi_1.42.1        memoise_1.1.0              
## [87] IRanges_2.14.12
```
