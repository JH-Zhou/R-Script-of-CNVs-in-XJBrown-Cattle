# Introduction
The documents inculde the R Script, CNV lists and map files used in the post-analysis of CNVs in XJBrown-Cattle.

To replicate our results, we first need to install and load the following R packages

```{r}
library(HandyCNV)
library(data.table)
library(dplyr)
```

# 1.Coverting and comparing map fies between UMD and ARS versions
```{r results='hide'}
#UMD Vs ARS
convert_map(default_map = "1_map/XJB_CNV_UMD.map",
            target_map = "1_map/XJB_CNV_ARS.map",
            defMap_title = "UMD 3.1", tarMap_title = "Offical ARS 1.2", 
            col_1 = "steelblue1", col_2 = "tomato", col_3 = "steelblue1", 
            col_4 = "steelblue1", col_5 = "tomato", col_6 = "tomato")
```

# 2. Summary and process CNVs results
## 2.1 Summrary of Part-UMD results
```{r}
#clean cnv
clean_cnv_part_umd <- cnv_clean(cnvpartition = "2_cnv_results/XJB_Part_CNV_Report_UMD.txt", 
                                drop_length = 5, 
                                folder = "2_cnv_results/1_part_umd")
#cnv summary plot
cnv_summary_plot(clean_cnv = clean_cnv_part_umd, 
                 plot_sum_1 = TRUE, plot_sum_2 = TRUE, 
                 folder = "2_cnv_results/1_part_umd/summary_plot",
                 height_sum1 = 26, width_sum1 = 20,
                 col_0 = "tomato", col_1 = "deeppink4", 
                 col_3 = "turquoise", col_4 = "royalblue")
#generate CNVR
cnvr_part_umd <- call_cnvr(clean_cnv = clean_cnv_part_umd, 
                           chr_set = 29, 
                           folder = "2_cnv_results/1_part_umd/call_cnvr")
#cnvr distribution map
cnvr_plot(cnvr = cnvr_part_umd, assembly = "UMD", 
          gain_col = "deeppink1", loss_col = "deepskyblue3", mixed_col = "springgreen3", 
          folder = "2_cnv_results/1_part_umd/cnvr_map")
#get reference gene
get_refgene()
refgene_UMD <- get_refgene(gene_version = "Cow_UMD_UCSC")
#annotate gene from UCSC reference gene list
call_gene(refgene = refgene_UMD, 
          interval = cnvr_part_umd, 
          clean_cnv = clean_cnv_part_umd, 
          folder = "2_cnv_results/1_part_umd/call_gene")
#plot high frequent CNVR with genes
cnvr_plot(cnvr = cnvr_part_umd, 
          clean_cnv = clean_cnv_part_umd, 
          sample_size = 365, gene_font_size = 2.2, common_cnv_threshold = 0.01, 
          refgene = refgene_UMD, col_gene = "yellow2", 
          folder = "2_cnv_results/1_part_umd/high_freq_cnvr")
#plot high frequent CNVR with related information
#Note: This step not available because of signal and genotype too large to share
# plot_cnvr_panorama(cnvr = cnvr_part_umd, 
#                    cnv_annotation = "2_cnv_results/1_part_umd/call_gene/cnv_annotation.txt", 
#                    intensity = "final_403_XJB_CNV_FinalReport.txt", #This file too large to share
#                    map = "1_map/XJB_CNV_UMD.map", 
#                    prefix_bed = "PLINK_150820_0913/final_403_XJB_CNV", #This file too large to share
#                    ld_heat = TRUE, sample_size = 365, common_cnv_threshold = 0.01, 
#                    col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
#                    col_3 = "turquoise", col_4 = "royalblue", col_gene = "red",
#                    width_1 = 12, height_1 = 25,
#                    folder = "2_cnv_results/1_part_umd/cnvr_panorama_heat")
```

## 2.2 Summrary of Part-ARS results
```{r}
#clean cnv
clean_cnv_part_ars <- cnv_clean(cnvpartition = "2_cnv_results/XJB_Part_CNV_Report_ARS.txt",
                                drop_length = 5,
                                folder = "2_cnv_results/2_part_ars")
#cnv summary plot
cnv_summary_plot(clean_cnv = clean_cnv_part_ars, 
                 plot_sum_1 = TRUE, 
                 plot_sum_2 = TRUE,
                 folder = "2_cnv_results/2_part_ars/summary_plot",
                 height_sum1 = 26, width_sum1 = 20,
                 col_0 = "tomato", col_1 = "deeppink4", col_3 = "turquoise", col_4 = "royalblue")
#call cnvr
cnvr_part_ars <- call_cnvr(clean_cnv = clean_cnv_part_ars, 
                           chr_set = 29, 
                           folder = "2_cnv_results/2_part_ars/call_cnvr")
#cnvr map
cnvr_plot(cnvr = cnvr_part_ars, 
          assembly = "ARS", 
          gain_col = "deeppink1", 
          loss_col = "deepskyblue3", 
          mixed_col = "springgreen3", 
          folder = "2_cnv_results/2_part_ars/cnvr_map")
#get reference gene
get_refgene()
refgene_ARS <- get_refgene(gene_version = "Cow_ARS_UCSC")
#annotate gene from UCSC reference gene list
call_gene(refgene = refgene_ARS, 
          interval = cnvr_part_ars, 
          clean_cnv = clean_cnv_part_ars, 
          folder = "2_cnv_results/2_part_ars/call_gene")
#plot high frequency cnvr
# #Note: This step not available because of signal and genotype too large to share
# plot_cnvr_panorama(cnvr = cnvr_part_ars, 
#                    cnv_annotation = "2_cnv_results/2_part_ars/call_gene/cnv_annotation.txt", #This file too large to share
#                    intensity = "cnv_result/part_o_ars/XJB_ARS_403_GStudio_FinalReport.txt", 
#                    map = "1_map/XJB_CNV_ARS.map", 
#                    prefix_bed = "cnv_result/part_o_ars/PLINK_050221_0136/official_ars_xjb", #This file too large to share
#                    ld_heat = TRUE, 
#                    sample_size = 378, common_cnv_threshold = 0.01, 
#                    col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
#                    col_3 = "turquoise", col_4 = "royalblue", col_gene = "red", 
#                    width_1 = 12, height_1 = 25, 
#                    folder = "2_cnv_results/2_part_ars/cnvr_panorama_heat")
```


## 2.3 Summrary of Penn-UMD results
```{r}
#clean cnv
clean_cnv_penn_umd <- cnv_clean(penncnv = "2_cnv_results/XJB_UMD_403_29_GC_Merge.goodcnv", 
          penn_id_sep = "cnv/", 
          drop_length = 5, 
          folder = "2_cnv_results/3_penn_umd")
#cnv summary plot
cnv_summary_plot(clean_cnv = clean_cnv_penn_umd, 
                 plot_sum_1 = TRUE, 
                 plot_sum_2 = TRUE,
                 folder = "2_cnv_results/3_penn_umd/summary_plot",
                 height_sum1 = 26, width_sum1 = 20,
                 col_0 = "tomato", col_1 = "deeppink4", col_3 = "turquoise", col_4 = "royalblue")
#call cnvr
cnvr_penn_umd <- call_cnvr(clean_cnv = clean_cnv_penn_umd, 
                           chr_set = 29, 
                           folder = "2_cnv_results/3_penn_umd/call_cnvr")
#cnvr map
cnvr_plot(cnvr = cnvr_penn_umd, 
          assembly = "UMD", 
          gain_col = "deeppink1", 
          loss_col = "deepskyblue3", 
          mixed_col = "springgreen3", 
          folder = "2_cnv_results/3_penn_umd/cnvr_map")
#annotate gene from UCSC reference gene list
call_gene(refgene = refgene_UMD, 
          interval = cnvr_penn_umd, 
          clean_cnv = clean_cnv_penn_umd, 
          folder = "2_cnv_results/3_penn_umd/call_gene")
#plot high frequent CNVR
# #Note: This step not available because of signal and genotype too large to share
# plot_cnvr_panorama(cnvr = cnvr_penn_umd, 
#                    cnv_annotation = "cnv_result/penn_umd/call_gene_ucsc/cnv_annotation.txt", 
#                    intensity = "cnv_result/part_umd/final_403_XJB_CNV_FinalReport.txt", 
#                    map = "1_map/XJB_CNV_UMD.map", 
#                    prefix_bed = "cnv_result/part_umd/PLINK_150820_0913/final_403_XJB_CNV", 
#                    ld_heat = TRUE, sample_size = 387, common_cnv_threshold = 0.01, 
#                    col_0 = "tomato", col_1 = "deeppink4", col_2 = "NA", 
#                    col_3 = "turquoise", col_4 = "royalblue", col_gene = "red", 
#                    width_1 = 12, height_1 = 25, 
#                    folder = "2_cnv_results/3_penn_umd/cnvr_panorama")
```



## 2.4 Summrary of Penn-ARS results
```{r}
#clean cnv
clean_cnv_penn_ars <- cnv_clean(penncnv = "2_cnv_results/OARS_XJB_GC_Merge.goodcnv",
                                penn_id_sep = "sity/",
                                drop_length = 5, 
                                folder = "2_cnv_results/4_penn_ars")
#cnv summary plot
cnv_summary_plot(clean_cnv = clean_cnv_penn_ars, 
                 plot_sum_1 = TRUE, 
                 plot_sum_2 = TRUE,
                 folder = "2_cnv_results/4_penn_ars/summary_plot",
                 height_sum1 = 26, width_sum1 = 20,
                 col_0 = "tomato", col_1 = "deeppink4", col_3 = "turquoise", col_4 = "royalblue")
#call cnvr
cnvr_penn_ars <- call_cnvr(clean_cnv = clean_cnv_penn_ars, 
                           chr_set = 29, 
                           folder = "2_cnv_results/4_penn_ars/call_cnvr")
#cnvr map
cnvr_plot(cnvr = cnvr_penn_ars, 
          assembly = "ARS", 
          gain_col = "deeppink1", 
          loss_col = "deepskyblue3", 
          mixed_col = "springgreen3", 
          folder = "2_cnv_results/4_penn_ars/cnvr_map")
#annotate gene from USCS reference gene list
call_gene(refgene = refgene_ARS, 
          interval = cnvr_penn_ars, 
          clean_cnv = clean_cnv_penn_ars, 
          folder = "2_cnv_results/4_penn_ars/call_gene")
#plot high frequency CNVR
#Note: This step not available because of signal and genotype too large to share
# plot_cnvr_panorama(cnvr = cnvr_penn_ars, 
#                    cnv_annotation = "2_cnv_results/4_penn_ars/call_gene/cnv_annotation.txt", 
#                    intensity = "cnv_result/part_o_ars/XJB_ARS_403_GStudio_FinalReport.txt", 
#                    map = "1_map/XJB_CNV_ARS.map", 
#                    prefix_bed = "cnv_result/part_o_ars/PLINK_050221_0136/official_ars_xjb", 
#                    ld_heat = TRUE, 
#                    sample_size = 388, common_cnv_threshold = 0.04, 
#                    col_0 = "tomato", col_1 = "deeppink4", col_2 = "NA", 
#                    col_3 = "turquoise", col_4 = "royalblue", col_gene = "red", 
#                    width_1 = 12, height_1 = 25, 
#                    folder = "2_cnv_results/4_penn_ars/cnvr_panorama")
```


# 3. Comaprison of CNVs
```{r}
#Penn-UMD Vs Part-UMD
compare_cnv(cnv_def = clean_cnv_penn_umd, 
            cnv_tar = clean_cnv_part_umd, 
            width_1 = 16, height_1 = 12, 
            plot_caption = FALSE, col_1 = "turquoise", col_2 = "plum", 
            folder = "2_cnv_results/5_compare_cnv_Penn_UMD_Vs_Part_UMD")
#Part-o-ARS Vs Penn-P-ARS
compare_cnv(cnv_def = clean_cnv_part_ars, 
            cnv_tar = clean_cnv_penn_ars, 
            width_1 = 16, height_1 = 12, 
            plot_caption = FALSE, 
            col_1 = "turquoise", col_2 = "plum", 
            folder = "2_cnv_results/6_compare_cnv_Part_OARS_Vs_Penn_OARS")
```


# 4. Comparison of CNVRs
```{r}
#Penn-UMD Vs Part-UMD
compare_cnvr(cnvr_def = cnvr_penn_umd, 
             cnvr_tar = cnvr_part_umd, 
             width_1 = 8, height_1 = 10, hjust_prop = 0.2, 
             hjust_num = 1.6, 
             col_1 = "darkorchid", col_2 = "turquoise", 
             folder = "2_cnv_results/7_compare_cnvr_Penn_UMD_Vs_Part_UMD")
#Part-o-ARS Vs Penn-P-ARS
compare_cnvr(cnvr_def = cnvr_part_ars, 
             cnvr_tar = cnvr_penn_ars, 
             width_1 = 8, height_1 = 10, hjust_prop = 0.2, 
             hjust_num = 1.6, 
             col_1 = "darkorchid", col_2 = "turquoise", 
             folder = "2_cnv_results/8_compare_cnvr_Part_OARS_Vs_Penn_OARS")
```


# 5. Customise visualiszing CNVs and CNVRs
## 5.1 Visualisze CNVs in all chromosomes
```{r}
#Part-UMD
cnv_visual(clean_cnv = clean_cnv_part_umd, max_chr = 29, 
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue",
           folder = "2_cnv_results/1_part_umd/cnv_visual")
#Part-O-ARS
cnv_visual(clean_cnv = clean_cnv_part_ars, max_chr = 29, 
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue",
           folder = "2_cnv_results/2_part_ars/cnv_visual")
#Penn-UMD
cnv_visual(clean_cnv = clean_cnv_penn_umd, max_chr = 29, 
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue",
           folder = "2_cnv_results/3_penn_umd/cnv_visual")
#Penn-P_ARS
cnv_visual(clean_cnv = clean_cnv_penn_ars, max_chr = 29,
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue", 
           folder = "2_cnv_results/4_penn_ars/cnv_visual")
```

## 5.2 Customise CNVs visualizition in some interested regions
```{r}
#Part-UMD
cnv_visual(clean_cnv = clean_cnv_part_umd, chr_id  = 28, 
           height_1 = 8, width_1 = 13,
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue",
           folder = "2_cnv_results/1_part_umd/visual_interested_chr")
cnv_visual(clean_cnv = clean_cnv_part_umd, chr_id  = 12, 
           height_1 = 8, width_1 = 13,
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue",
           folder = "2_cnv_results/1_part_umd/visual_interested_chr")
cnv_visual(clean_cnv = clean_cnv_part_umd, chr_id  = 15, 
           height_1 = 8, width_1 = 13,
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue",
           folder = "2_cnv_results/1_part_umd/visual_interested_chr")
#Part-O-ARS
cnv_visual(clean_cnv = clean_cnv_part_ars, chr_id  = 28, 
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue", 
           height_1 = 8, width_1 = 13,
           folder = "2_cnv_results/2_part_ars/visual_interested_chr")
cnv_visual(clean_cnv = clean_cnv_part_ars, chr_id  = 12, 
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue", 
           height_1 = 8, width_1 = 13,
           folder = "2_cnv_results/2_part_ars/visual_interested_chr")
cnv_visual(clean_cnv = clean_cnv_part_ars, chr_id  = 15, 
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue", 
           height_1 = 8, width_1 = 13,
           folder = "2_cnv_results/2_part_ars/visual_interested_chr")
#Penn-P_ARS
cnv_visual(clean_cnv = clean_cnv_penn_ars, chr_id  = 28,
           height_1 = 8, width_1 = 13,
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue", 
           folder = "2_cnv_results/4_penn_ars/visual_interested_chr")
cnv_visual(clean_cnv = clean_cnv_penn_ars, chr_id  = 12,
           height_1 = 8, width_1 = 13,
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue", 
           folder = "2_cnv_results/4_penn_ars/visual_interested_chr")
cnv_visual(clean_cnv = clean_cnv_penn_ars, chr_id  = 15,
           height_1 = 8, width_1 = 13,
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue", 
           folder = "2_cnv_results/4_penn_ars/visual_interested_chr")
#Penn-UMD
cnv_visual(clean_cnv = clean_cnv_penn_umd, chr_id  = 28, 
           height_1 = 8, width_1 = 13,
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue",
           folder = "2_cnv_results/3_penn_umd/visual_interested_chr")
cnv_visual(clean_cnv = clean_cnv_penn_umd, chr_id  = 12, 
           height_1 = 8, width_1 = 13,
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue",
           folder = "2_cnv_results/3_penn_umd/visual_interested_chr")
cnv_visual(clean_cnv = clean_cnv_penn_umd, chr_id  = 15, 
           height_1 = 8, width_1 = 13,
           col_0 = "tomato", col_1 = "deeppink4", col_2 = "snow1", 
           col_3 = "turquoise", col_4 = "royalblue",
           folder = "2_cnv_results/3_penn_umd/visual_interested_chr")
```


# 6. Compare the difference of specific CNVRs
```{r}
#compare cnvr and cnv
file_1 <- data.table::fread("2_cnv_results/1_part_umd/call_cnvr/individual_cnv_cnvr.txt", header = TRUE)
file_2 <- data.table::fread("2_cnv_results/2_part_ars/call_cnvr/individual_cnv_cnvr.txt", header = TRUE)
file_3 <- data.table::fread("2_cnv_results/3_penn_umd/call_cnvr/individual_cnv_cnvr.txt", header = TRUE)
file_4 <- data.table::fread("2_cnv_results/4_penn_ars/call_cnvr/individual_cnv_cnvr.txt", header = TRUE)
#get the snp for each interval
ars_map <- data.table::fread("1_map/XJB_CNV_ARS.map")
umd_map <- data.table::fread("1_map/XJB_CNV_UMD.map")
#' Title
#' Function to get the unique sample size in all possible combination results
#' @param id_1 CNVR ID from Penn-PARS 
#' @param id_2 CNVR ID from Penn_UMD
#' @param id_3 CNVR_ID from Part_OARS
#' @param id_4 CNVR_ID from Part_UMD
#' @return Unique sample size for each lists
#' @export get_unique_sample
get_unique_sample <- function(id_1, id_2, id_3, id_4){
  # 1. comparison of chr15:78-81
id_1 <- file_1 %>%
        filter(CNVR_ID == id_1) %>%
        mutate(Version = "Penn-ARS") %>%
        select(Chr, Start, End, CNVR_ID, Sample_ID, CNV_Start, CNV_End, Length, CNV_Value, Version)
id_2 <- file_2 %>%
        filter(CNVR_ID == id_2) %>%
        mutate(Version = "Penn_UMD") %>%
        select(Chr, Start, End, CNVR_ID, Sample_ID, CNV_Start, CNV_End, Length, CNV_Value, Version)
id_3 <- file_3 %>%
        filter(CNVR_ID == id_3) %>%
        mutate(Version = "Part-ARS") %>%
        select(Chr, Start, End, CNVR_ID, Sample_ID, CNV_Start, CNV_End, Length, CNV_Value, Version)
id_4 <- file_4 %>%
        filter(CNVR_ID == id_4) %>%
        mutate(Version = "Part-UMD") %>%
        select(Chr, Start, End, CNVR_ID, Sample_ID, CNV_Start, CNV_End, Length, CNV_Value, Version)
#check how many unique sample in each result
cat(paste0("Unique sample size in first result are ", length(unique(id_1$Sample_ID)), "\n"))
cat(paste0("Unique sample size in second result are ", length(unique(id_2$Sample_ID)), "\n"))
cat(paste0("Unique sample size in third results are ", length(unique(id_3$Sample_ID)), "\n"))
cat(paste0("Unique sample size in fourth results are ", length(unique(id_4$Sample_ID)), "\n"))
#merge all results
id_all <- rbind(id_1, id_2, id_3, id_4)
id_all_unique <- select(id_all,Sample_ID) %>%
                 unique(.)
cat(paste0("Unique sample size in all four results are ", length(unique(id_all$Sample_ID)), "\n"))
#View(id_all)
#merge ARS only
ars_1 <- rbind(id_1, id_3)
cat(paste0("Unique sample size in two ARS results are ", length(unique(ars_1$Sample_ID)), "\n"))
#merge UMD only
umd_1 <- rbind(id_2, id_4)
cat(paste0("Unique sample size in two UMD results are ", length(unique(umd_1$Sample_ID)), "\n"))
return(id_all_unique)
}
#small function to get SNPs for CNVR
get_snp <- function(map_plink, chr, start, end){
  target <- map_plink %>%
            arrange(V1, V4) %>%
            filter(V1 == chr & V4 >= start & V4 <= end)
  return(target)
}
# 1.unique samples on Chr15 region
get_unique_sample(id_1 = "CNVR_622", id_2 = "CNVR_561", id_3 = "CNVR_180", id_4 = "CNVR_156")
#SNPs within CNVR_618 in Penn-PARS
cnvr_622 <- get_snp(map_plink = ars_map, chr = 15, start = 78985623, end = 79575418)
#SNPs within CNVR_561 in Penn_UMD
cnvr_561 <- get_snp(map_plink = umd_map, chr = 15, start = 80280657, end = 80842722)
#SNPs within CNVR_180 in Part-OARS
cnvr_180 <- get_snp(map_plink = ars_map, chr = 15, start = 79048912, end = 79611716)
#SNPs within CNVR_156 in Part-UMD
cnvr_156 <- get_snp(map_plink = umd_map, chr = 15, start = 80329011, end = 80935724)
# 2.comparison of Chr3:54-55 Mb
get_unique_sample(id_1 = "CNVR_134", id_2 = "CNVR_111", id_3 = "CNVR_34", id_4 = "CNVR_27")
#results shown same
cnvr_134 <- get_snp(map_plink = ars_map, chr = 3, start = 54211974, end = 54781990)
# 3.comparison of Chr2:27-28 Mb
get_unique_sample(id_1 = "CNVR_80", id_2 = "CNVR_63", id_3 = "CNVR_19", id_4 = "CNVR_17")
#results shown same
# 4.comparison of Chr4:65-67
get_unique_sample(id_1 = "CNVR_179", id_2 = "CNVR_151", id_3 = "CNVR_48", id_4 = "CNVR_39")
```


# 7. Get consensus CNVRs

## 7.1 Get consensus CNVRs for UMD version
```{r}
#combine CNVs for UMD version
common_colnames <- c("Sample_ID", "Chr", "Start", "End", "Length", "CNV_Value")
cnv_penn_umd <- dplyr::select(clean_cnv_penn_umd, all_of(common_colnames))
cnv_part_umd <- dplyr::select(clean_cnv_part_umd, all_of(common_colnames))
cnv_umd_all <- base::rbind(cnv_penn_umd, cnv_part_umd)
cnvr_umd_all <- call_cnvr(clean_cnv = cnv_umd_all, folder = "2_cnv_results/9_cnvr_combine_part_penn_umd")
write.table(cnv_umd_all, file = "2_cnv_results/9_cnvr_combine_part_penn_umd/cnv_umd_all.txt", quote = F, col.names = T, sep = "\t", row.names = F)
cnvr_plot(cnvr = cnvr_umd_all, assembly = "UMD", folder = "2_cnv_results/9_cnvr_combine_part_penn_umd/cnvr_map")
unique_sample_size <- length(unique(cnv_umd_all$Sample_ID))
#annotate genes for united CNV regions
call_gene(refgene = refgene_UMD, interval = cnvr_umd_all, clean_cnv = cnv_umd_all, folder = "2_cnv_results/9_cnvr_combine_part_penn_umd/call_gene")
#a. get consensus CNVR
consensus_cnvr <- cnvr_umd_all %>%
                  filter(n_Sample > 394 * 0.05)
write.table(consensus_cnvr, file = "2_cnv_results/9_cnvr_combine_part_penn_umd/consensus_cnvr_umd.txt", quote = F, sep = "\t", col.names = T, row.names = F)
#plot all CNVRs and label the common regions 
consensus_cnvr_UMD_plot <- consensus_cnvr %>%
                           select(Chr, Start, End)
setDT(consensus_cnvr_UMD_plot)
cnvr_plot(cnvr = cnvr_umd_all, assembly = "UMD", 
          sample_size =  394, common_cnv_threshold = 0.00, 
          overlap_cnvr = consensus_cnvr_UMD_plot, 
          gain_col = "black", loss_col = "blue", mixed_col = "tomato", chr_col = "gray", 
          folder = "2_cnv_results/9_cnvr_combine_part_penn_umd/cnvr_map_common")
#b. get annotated CNVRs
cnvr_gene <- read.table(file = "2_cnv_results/9_cnvr_combine_part_penn_umd/call_gene/interval_gene_summarise_table.txt", header = T)
#c. extract consensus CNVRs from annotated list
consensus_cnvr_gene <- cnvr_gene %>%
                       filter(ID %in% consensus_cnvr$CNVR_ID)
#d. count the number of CNVRs without gene annotated
n_cnvr_without_gene <- consensus_cnvr_gene %>%
                       filter(is.na(gene_name))
n_cnvr_with_gene <- consensus_cnvr_gene %>%
                       filter(!is.na(gene_name))
write.table(n_cnvr_with_gene, file = "2_cnv_results/9_cnvr_combine_part_penn_umd/n_consensus_cnvr_with_gene_umd.txt", 
            quote = F, sep = "\t", col.names = T, row.names = F)
#calculate the statistics of Common CNVRs generated in UMD version
common_cnvr_umd <- read.table(file ="2_cnv_results/9_cnvr_combine_part_penn_umd/cnvr_map_common/cnvr_plot.txt", header = T, sep ="\t")
common_cnvr_umd_summary <- common_cnvr_umd %>%
                           group_by(Type) %>%
                           summarise(Number = n(),
                                     Average = mean(Length, na.rm = T), 
                                     Max = max(Length, na.rm = T),
                                     Min = min(Length, na.rm = T),
                                     Total_length = sum(Length, na.rm = T))
summary_chr_umd <- common_cnvr_umd %>%
                   group_by(Chr) %>%
                   summarise(N = n(),
                         Total_Length = sum(Length, na.rm = T))
write.table(common_cnvr_umd_summary, file = "2_cnv_results/9_cnvr_combine_part_penn_umd/cnvr_map_common/common_cnvr_summary.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)
common_cnvr_umd_plot <- common_cnvr_umd %>% 
                        select(Chr, Start, End)
write.table(common_cnvr_umd_plot, file = "2_cnv_results/9_cnvr_combine_part_penn_umd/cnvr_map_common/common_cnvr_plot.txt", 
            quote = F, sep = "\t", col.names = T, row.names = F)
```


## 7.2 Check weather the CNVRs generated by combined CNV list same as the manually generated? 
```{r}
#manually generated the unique CNVR list in two UMD results(PennCNV and CNVPartition)
cnvr_umd <- base::rbind(cnvr_penn_umd, cnvr_part_umd)
names(cnvr_umd)
#function to merge overlapped CNVRs
merge_cnvr <- function(cnv) {
    if (nrow(cnv) == 1) {
      return(cnv)
    }
    cnv <- cnv[order(cnv$Chr, cnv$Start),]
    cnvr_union = cnv[1, ]
    for (i in 2:nrow(cnv)) {
      rest_cnv <- cnv[i, ]
      if (cnvr_union$End[nrow(cnvr_union)] < rest_cnv$Start) {
        cnvr_union <- bind_rows(cnvr_union, rest_cnv)
      } else if (cnvr_union$End[nrow(cnvr_union)] == rest_cnv$Start) {
        cnvr_union$End[nrow(cnvr_union)] <- rest_cnv$End
      }
      
      if (rest_cnv$End > cnvr_union$End[nrow(cnvr_union)]) {
        cnvr_union$End[nrow(cnvr_union)] <- rest_cnv$End
      }
    }
    return(cnvr_union)
}
get_cnvr <- function(interval){
  
  interval <- interval
  max_chr <- max(interval$Chr)
  cnvr <- data.frame()
  
  for (i in 1:max_chr){
    cnv_chr <- interval[which(interval$Chr == i), ]
    cnvr_chr <- merge_cnvr(cnv = cnv_chr)
    cnvr <- rbind(cnvr, cnvr_chr)
    print(paste0("Chromosome ", i, " has been processed."))
  }
  
  return(cnvr)
}
#Overlapped CNVR
ocnvr_umd_penn <- read.table(file = "2_cnv_results/7_compare_cnvr_Penn_UMD_Vs_Part_UMD/overlap_cnvr_def.popu", sep = "\t", header = TRUE)
ocnvr_umd_part <- read.table(file = "2_cnv_results/7_compare_cnvr_Penn_UMD_Vs_Part_UMD/overlap_cnvr_tar.popu", sep = "\t", header = TRUE)
names(ocnvr_umd_part)
ocnvr_umd_penn <- ocnvr_umd_penn %>%
                  select(Chr_DEF, Start_DEF, End_DEF) %>%
                  rename(Chr = Chr_DEF, Start = Start_DEF, End = End_DEF)
ocnvr_umd_part <- ocnvr_umd_part %>%
                  select(Chr_TAR, Start_TAR, End_TAR) %>%
                  rename(Chr = Chr_TAR, Start = Start_TAR, End = End_TAR)
ocnvr_umd_all <- base::rbind(ocnvr_umd_penn, ocnvr_umd_part)
ocnvr_all <- data.frame()
 for (i in 1:29){
    cnv_chr <- ocnvr_umd_all[which(ocnvr_umd_all$Chr == i), ]
    cnvr_chr <- merge_cnvr(cnv = cnv_chr)
    ocnvr_all <- rbind(ocnvr_all, cnvr_chr)
    print(paste0("Chromosome ", i, " has been processed."))
 }
ocnvr_all_1 <- get_cnvr(interval = ocnvr_umd_all)
#Non-overlapped CNVR
nocnvr_umd_penn <- read.table(file = "2_cnv_results/7_compare_cnvr_Penn_UMD_Vs_Part_UMD/non_overlap_cnvr_def.popu", sep = "\t", header = TRUE)
nocnvr_umd_part <- read.table(file = "2_cnv_results/7_compare_cnvr_Penn_UMD_Vs_Part_UMD/non_overlap_cnvr_tar.popu", sep = "\t", header = TRUE)
nocnvr_umd_penn <- nocnvr_umd_penn %>%
                   select(Chr_DEF, Start_DEF, End_DEF) %>%
                   rename(Chr = Chr_DEF, Start = Start_DEF, End = End_DEF)
nocnvr_umd_part <- nocnvr_umd_part %>%
                   select(Chr_TAR, Start_TAR, End_TAR) %>%
                   rename(Chr = Chr_TAR, Start = Start_TAR, End = End_TAR)
#unite the CNVR list
nocnvr_umd_all <- base::rbind(nocnvr_umd_penn, nocnvr_umd_part, ocnvr_all)
#check if the manual generated CNVR list are competently same as the CNVR list generated by unite CNV results
#check chr, start and end position 
all(nocnvr_umd_all$Start %in% cnvr_umd_all$Start)
all(nocnvr_umd_all$End %in% cnvr_umd_all$End)
nocnvr_umd_all <- tidyr::unite(nocnvr_umd_all, col = "region", c("Chr", "Start", "End"), sep = "", remove = FALSE) %>%
                  arrange(region)
cnvr_umd_all <- tidyr::unite(cnvr_umd_all, col = "region", c("Chr", "Start", "End"), sep = "", remove = FALSE) %>%
                arrange(region)
all(nocnvr_umd_all$region == cnvr_umd_all$region)
#all results are returned TRUE value which indicate the the two method get the same results
```


## 7.3 Get consenuses CNVRs for ARS version 
```{r}
#combine CNVs for ARS version
cnv_penn_ars <- dplyr::select(clean_cnv_penn_ars, all_of(common_colnames))
cnv_part_oars <- dplyr::select(clean_cnv_part_ars, all_of(common_colnames))
cnv_ars_all <- base::rbind(cnv_penn_ars, cnv_part_oars)
cnvr_ars_all <- call_cnvr(cnv_ars_all, folder = "2_cnv_results/10_cnvr_combine_part_penn_ars")
write.table(cnv_ars_all, file = "2_cnv_results/10_cnvr_combine_part_penn_ars/cnv_ars_all.txt", 
            quote = F, sep = "\t", col.names = T, row.names = F)
unique_sample_size <- length(unique(cnv_ars_all$Sample_ID))
#annotate genes
call_gene(refgene = refgene_ARS, interval = cnvr_ars_all, clean_cnv = cnv_ars_all, 
          folder = "2_cnv_results/10_cnvr_combine_part_penn_ars/call_gene")
#a. get consensus cnvr
consensus_cnvr_ars <- cnvr_ars_all %>%
                      filter(n_Sample > unique_sample_size * 0.05)
write.table(consensus_cnvr_ars, file = "2_cnv_results/10_cnvr_combine_part_penn_ars/consensus_cnvr_ars.txt", 
            quote = F, sep = "\t", col.names = T, row.names = F)
#plot all CNVRs and label the common regions 
consensus_cnvr_ars_plot <- consensus_cnvr_ars %>%
                           select(Chr, Start, End)
setDT(consensus_cnvr_ars_plot)
cnvr_plot(cnvr = cnvr_ars_all, assembly = "ARS", 
          overlap_cnvr = consensus_cnvr_ars_plot, 
          sample_size = 393, common_cnv_threshold = 0.00, 
          gain_col = "deeppink1", loss_col = "deepskyblue3", mixed_col = "springgreen3",chr_col = "gray",
          folder = "2_cnv_results/10_cnvr_combine_part_penn_ars/cnvr_map_common")
#b. get annotated CNVRs
cnvr_gene_ars <- read.table(file = "2_cnv_results/10_cnvr_combine_part_penn_ars/call_gene/interval_gene_summarise_table.txt", header = T)
#c. extract consensus CNVRs from annotated list
consensus_cnvr_gene_ars <- cnvr_gene_ars %>%
                           filter(ID %in% consensus_cnvr_ars$CNVR_ID)
#d. count the number of CNVRs without gene annotated
n_cnvr_without_gene_ars <- consensus_cnvr_gene_ars %>%
                       filter(is.na(gene_name))
n_cnvr_with_gene_ars <- consensus_cnvr_gene_ars %>%
                       filter(!is.na(gene_name))
write.table(n_cnvr_with_gene_ars, file = "2_cnv_results/10_cnvr_combine_part_penn_ars/n_consensus_cnvr_with_gene_ars.txt", 
            quote = F, sep = "\t", col.names = T, row.names = F)
#calculate the statistics of Common CNVRs for ARS version
common_cnvr_ars <- read.table(file ="2_cnv_results/10_cnvr_combine_part_penn_ars/cnvr_map_common/cnvr_plot.txt", header = T, sep ="\t")
common_cnvr_ars_summary <- common_cnvr_ars %>%
                           group_by(Type) %>%
                           summarise(Number = n(),
                                     Average = mean(Length, na.rm = T), 
                                     Max = max(Length, na.rm = T),
                                     Min = min(Length, na.rm = T),
                                     Total_length = sum(Length, na.rm = T))
summary_chr_ars <- common_cnvr_ars %>%
                   group_by(Chr) %>%
                   summarise(N = n(),
                         Total_Length = sum(Length, na.rm = T))
write.table(common_cnvr_ars_summary, file = "2_cnv_results/10_cnvr_combine_part_penn_ars/cnvr_map_common/common_cnvr_summary.txt", 
            quote = F, sep = "\t", col.names = T, row.names = F)
common_cnvr_ars_plot <- common_cnvr_ars %>% 
                        select(Chr, Start, End)
write.table(common_cnvr_ars_plot, file = "2_cnv_results/10_cnvr_combine_part_penn_ars/cnvr_map_common/common_cnvr_plot.txt", 
            quote = F, sep = "\t", col.names = T, row.names = F)
```


# 8. Comparison of CNV-Genes and get the consensus genes
```{r}
compare_gene(gene_freq_1 = "2_cnv_results/4_penn_ars/call_gene/gene_freq_cnv.txt", title_1 = "Penn ARS",
             gene_freq_2 = "2_cnv_results/2_part_ars/call_gene/gene_freq_cnv.txt", title_2 = "Part ARS",
             gene_freq_3 = "2_cnv_results/3_penn_umd/call_gene/gene_freq_cnv.txt", title_3 = "Penn UMD",
             gene_freq_4 = "2_cnv_results/1_part_umd/call_gene/gene_freq_cnv.txt", title_4 = "Part UMD", 
             common_gene_threshold = 18, height_1 = 18, width_1 = 18, color_label = T,
             folder = "2_cnv_results/11_compare_gene")
#because we used the combine CNV list, the frequency of CNVRs should replace by n_Sample
gene_freq_cnv_ars <- fread("2_cnv_results/10_cnvr_combine_part_penn_ars/call_gene/gene_freq_cnv.txt")
cnv_annotation_ars <- fread("2_cnv_results/10_cnvr_combine_part_penn_ars/call_gene/cnv_annotation.txt")
gene_freq_cnv_ars_sample <- cnv_annotation_ars %>% 
                            group_by(name2) %>% 
                            summarise(Frequency = length(unique(Sample_ID))) %>%
                            left_join(gene_freq_cnv_ars[, c("name2", "ID","Chr", "CNVR_Start", "CNVR_End")], by = "name2") %>%
                            filter(Frequency < 350)
fwrite(gene_freq_cnv_ars_sample, "2_cnv_results/10_cnvr_combine_part_penn_ars/call_gene/gene_freq_cnv_ars_sample.txt", sep = "\t", quote = F)
gene_freq_cnv_umd <- fread("2_cnv_results/9_cnvr_combine_part_penn_umd/call_gene/gene_freq_cnv.txt")
cnv_annotation_umd <- fread("2_cnv_results/9_cnvr_combine_part_penn_umd/call_gene/cnv_annotation.txt")
gene_freq_cnv_umd_sample <- cnv_annotation_umd %>% 
                            group_by(name2) %>% 
                            summarise(Frequency = length(unique(Sample_ID))) %>%
                            left_join(gene_freq_cnv_umd[, c("name2", "ID","Chr", "CNVR_Start", "CNVR_End")], by = "name2") %>%
                            filter(Frequency < 350)
fwrite(gene_freq_cnv_umd_sample, "2_cnv_results/9_cnvr_combine_part_penn_umd/call_gene/gene_freq_cnv_umd_sample.txt", sep = "\t", quote = F)
compare_gene(gene_freq_1 = gene_freq_cnv_ars_sample, title_1 = "Gene frequency in union set of ARS version",
             gene_freq_2 = gene_freq_cnv_umd_sample, title_2 = "Gene frequency in union of UMD version",
             common_gene_threshold = 20, height_1 = 13, width_1 = 16, color_label = T,
             folder = "2_cnv_results/11_compare_gene/")
consensus_gene <- fread("2_cnv_results/11_compare_gene/two_lists_comparison.txt") 
ars_consensus_gene <- consensus_gene %>%
                      filter(!Common_Gene == "Common_Low") %>%
                      group_by(ID_1) %>%
                      summarise(Gene = paste(unique(name2_1), collapse = ",")) %>%
                      left_join(unique(consensus_gene[, c("ID_1", "Chr_1", "CNVR_Start_1", "CNVR_End_1")]), by = "ID_1") %>%
                      mutate(num_gene = lengths(strsplit(Gene, split = ",")))
fwrite(ars_consensus_gene, file = "2_cnv_results/11_compare_gene/ars_consensus_gene.txt", sep ="\t", quote = F)
umd_consensus_gene <- consensus_gene %>%
                      filter(!Common_Gene == "Common_Low") %>%
                      group_by(ID_2) %>%
                      summarise(Gene = paste(unique(name2_1), collapse = ",")) %>%
                      left_join(unique(consensus_gene[, c("ID_2", "Chr_2", "CNVR_Start_2", "CNVR_End_2")]), by = "ID_2") %>%
                      mutate(num_gene = lengths(strsplit(Gene, split = ",")))
fwrite(umd_consensus_gene, file = "2_cnv_results/11_compare_gene/umd_consensus_gene.txt", sep ="\t", quote = F)
```



# 9. Plot consensus genes
```{r}
consensus_gene <- fread("2_cnv_results/11_compare_gene/four_gene_comparison.txt") %>% 
                  filter(Common_check == "Common_high" | Common_check == "High_and_Low")
for(i in 1:nrow(consensus_gene)){
  try(cnv_visual(clean_cnv = "2_cnv_results/10_cnvr_combine_part_penn_ars/call_gene/cnv_annotation.txt", 
                 target_gene = consensus_gene$name2_1[i],
                 folder = "2_cnv_results/12_visual_gene_ars"),silent = T)
}
for(i in 1:nrow(consensus_gene)){
  try(cnv_visual(clean_cnv = "2_cnv_results/9_cnvr_combine_part_penn_umd/call_gene/cnv_annotation.txt", 
                 target_gene = consensus_gene$name2_1[i], 
                 folder = "2_cnv_results/13_visual_gene_umd"),silent = T)
}
```


# 10. Check unique sample size of CNV-Gene in four results

```{r}
consensus_gene <- read.table(file = "2_cnv_results/11_compare_gene/four_gene_comparison.txt", header = T, sep = "\t") %>%
                  filter(Common_check == "Common_high" | Common_check == "High_and_Low") 
gene_list <- as.vector(consensus_gene$name2_1)
for (i in 1:length(gene_list)) {
  penn_ars <- get_samples(annotated_cnvlist = "2_cnv_results/4_penn_ars/call_gene/cnv_annotation.txt", gene_name = gene_list[i])
  penn_umd <- get_samples(annotated_cnvlist = "2_cnv_results/3_penn_umd/call_gene/cnv_annotation.txt", gene_name = gene_list[i])
  part_ars <- get_samples(annotated_cnvlist = "2_cnv_results/2_part_ars/call_gene/cnv_annotation.txt", gene_name = gene_list[i])
  part_umd <- get_samples(annotated_cnvlist = "2_cnv_results/1_part_umd/call_gene/cnv_annotation.txt", gene_name = gene_list[i])
  all <- rbind(penn_ars, penn_umd, part_ars, part_umd)
  max_sample <- max(nrow(penn_ars), nrow(penn_umd), nrow(part_ars), nrow(part_umd))
  cat(paste0(gene_list[i], " gene with ", length(unique(all$Sample_ID)), " unqiue samples, maximal samples"), max_sample, " in list.\n")
  }
```

# 11. Check if any CNVRs located in known false positive regions?
There are 9 false positive CNVR regions caused by genome assembly errors, reported by following reference:
Zhou, Y., Utsunomiya, Y., Xu, L. et al. Comparative analyses across cattle genders and breeds reveal the pitfalls caused by false positive and lineage-differential copy number variations. Sci Rep 6, 29219 (2016). https://doi.org/10.1038/srep29219
```{r}
false_cnvr <- fread("2_cnv_results/false_positive_CNVR_by_mis_assembly.csv")
compare_interval(interval_1 = "2_cnv_results/9_cnvr_combine_part_penn_umd/consensus_cnvr_umd.txt", 
                 interval_2 = "2_cnv_results/false_positive_CNVR_by_mis_assembly.csv")
false_cnvr <- compare_interval(interval_1 = "2_cnv_results/9_cnvr_combine_part_penn_umd/cnvr.txt", 
                               interval_2 = "2_cnv_results/false_positive_CNVR_by_mis_assembly.csv")
umd_ars_map <- fread(input = "convert_map/def_tar_map.map") 
map_info <- data.frame()
for (i in 1:nrow(false_cnvr)){
  start <- umd_ars_map %>% filter(Chr_def == false_cnvr$Chr_2[i] & Position_def == false_cnvr$Start[i])
  end <- umd_ars_map %>% filter(Chr_def == false_cnvr$Chr_2[i] & Position_def == false_cnvr$End[i])
  map_info <- rbind(map_info, start, end)
}
```
