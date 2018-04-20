##########final version of Nilogen cytokine&nanostring analysis#############################
# Author: Jingping Qiao
# Reference: Alex Rolfe; Bartholomew Naughton
# Date: 10/26/2017
############################################################################################
#setwd("~/Desktop/Jieqing/jingping/Nilogen/final")

############################################################################################
# Multiplex data
# Manipulate cytokines data in excel (OORs, *)
# load in cytokine data, add 10 (for detect range) and then do logfc per sample
# plot heatmap
############################################################################################
raw_multiplex <- readODS::read.ods("~/Desktop/Jieqing/jingping/Nilogen/final/jp_final_data.ods",sheet = 4)
rownames(raw_multiplex) <- make.names(raw_multiplex[,1],unique = T)
raw_multiplex <- raw_multiplex[,-1]
colnames(raw_multiplex) <- raw_multiplex[1,]
raw_multiplex <- raw_multiplex[-1,]
#detect range 10
log_norm_raw <- log(data.matrix(raw_multiplex)+10)
#regulated normalize
cyto_logfc_vs_veh <- sapply(colnames(log_norm_raw), function (cn) {
  ctrl_col <- gsub("(A|A[+]B|B)$", "Veh", cn)
  log_norm_raw[,cn] - log_norm_raw[,ctrl_col]
})
cyto_logfc_vs_veh <- cyto_logfc_vs_veh[,-grep("Veh",colnames(cyto_logfc_vs_veh))]
################################
# mess with colors if you want #
################################
# makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
# {
#   stopifnot(length(colors) == 4)
#   ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
#   ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
#   return(c(ramp1, ramp2))
# }
# 
# cutoff.distance <- 3 
# cols <- makeColorRampPalette(c("blue", "white",    # distances 0 to 3 colored from white to red
#                                "white", "red"), # distances 3 to max(distmat) colored from green to black
#                              cutoff.distance / max(cyto_logfc_vs_veh),
#                              40)
# 
# 
# col1 = colorRampPalette(c("lightblue4","lightblue3","lightblue2","lightblue1","lightblue"))(4) #set the order of greys
# col2 <- rep("white", 1)
# col3 = colorRampPalette(c("red","red1","red2","red3", "red4"))(15)
#colors2 <- c(col1, col2, col3)
#tre <- rep("B",length(colnames(cyto_logfc_vs_veh)))
#tre[grep("A$",colnames(cyto_logfc_vs_veh))] <- "A"
#tre[grep("A[+]B",colnames(cyto_logfc_vs_veh))] <- "AB"
#tre <- data.frame(ID = tre)
#tre_colr <- list(ID = c(A = "red", B="blue",AB="green"))
pheatmap::pheatmap(
  cyto_logfc_vs_veh,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "logfc_vs_veh_plex",
  border_color = FALSE
)



############################################################################################
# PMAI_17 data
# scale
# see if correlated with others
############################################################################################
library(RColorBrewer)
PMAI <- readODS::read.ods("~/Desktop/Jieqing/jingping/Nilogen/final/jp_final_data.ods",sheet = 1)
colnames(PMAI) <- PMAI[1,]
rownames(PMAI) <- PMAI[,1]
PMAI <- PMAI[-1,-1]
PMAI <- t(scale(data.matrix(PMAI)))
pheatmap::pheatmap(
  PMAI,
  cluster_cols = T,
  cluster_rows = T,
 breaks = c(-2,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "PMAI_scaled",
  border_color = FALSE
)




############################################################################################
# FACS data
# add ldh cytotoxicty 48h normalized data
# remove live cell number
# log and substract with veh
############################################################################################
raw_facs <- readODS::read.ods("~/Desktop/Jieqing/jingping/Nilogen/final/jp_final_data.ods",sheet = 4)
ldh_toxi <- readODS::read.ods("~/Desktop/Jieqing/jingping/Nilogen/final/jp_final_data.ods",sheet = 1)
colnames(raw_facs) <- raw_facs[1,]
raw_facs <- raw_facs[-1,]
colnames(raw_facs)[1] <- "Variables"
colnames(ldh_toxi) <-colnames(raw_facs)
raw_facs <- rbind(raw_facs,ldh_toxi)

remove_raw_facs <- raw_facs[-c(1:4),]
variables <- remove_raw_facs[,1]
log_facs <- log2(data.matrix(remove_raw_facs[,-1])+1)
logfacs_vs_veh <- sapply(colnames(log_facs), function (cn) {
  ctrl_col <- gsub("(A|A[+]B|B)$", "Veh", cn)
  log_facs[,cn] - log_facs[,ctrl_col]
})
logfacs_vs_veh <- logfacs_vs_veh[,-grep("Veh",colnames(logfacs_vs_veh))]

pheatmap::pheatmap(
  logfacs_vs_veh,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-5,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,4),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "logfc_vs_veh_facs",
  labels_row = variables,
  border_color = FALSE
)



############################################################################################
# NANOSTRING PART from ALEX
# This part loads in rcc nanostring files and load signature score functions
#also normalized/logged raw data
############################################################################################

source(sprintf("%s/bmto-ngs/R/R/nanostring.R", Sys.getenv("GITCLONES")), chdir=T)
source(sprintf("%s/bmto-ngs/R/R/expression-outliers.R", Sys.getenv("GITCLONES")), chdir=T)

rcc_files <- load_rcc_files(list.files("~/GITCLONES/tim3-nilogen-tumor-explant-may-2017/data/NanoString/NanoString RCC data/", pattern="*.RCC", full.names=T))
expr_matrix <- rcc_list_to_matrix(rcc_files,
                                  gsub(".RCC", "", gsub(".*/", "", names(rcc_files))))

bg_normed_data <- bg_norm(expr_matrix)
pos_normed_data <- pos_norm(bg_normed_data)
nilogen_explant_hk_normed_data <- hk_norm(pos_normed_data,
                                          housekeeping_probes=intersect(nanostring_pan_cancer_housekeeping_genes,
                                                                        rownames(pos_normed_data)))

nilogen_explant_tumor_number <- gsub("EMD Tumor (\\d+) .*", "\\1", colnames(nilogen_explant_hk_normed_data))
nilogen_explant_treatment <- gsub("EMD Tumor \\d+ (.*)", "\\1", colnames(nilogen_explant_hk_normed_data))
nilogen_explant_tumor_treat <- sprintf("%s_%s", nilogen_explant_tumor_number, nilogen_explant_treatment)

nilogen_explant_log_counts <- log2(1 + nilogen_explant_hk_normed_data)
nilogen_explant_logfc <- nilogen_explant_log_counts - apply(nilogen_explant_log_counts, 1, median)
## spia and camera need log fold changes or similar.  let's compute against the
## tumor only samples.  Other option would be 0h.
nilogen_explant_logfc_vs_tu <- nilogen_explant_log_counts - apply(nilogen_explant_log_counts[,nilogen_explant_treatment == "Tu"], 1, median)

nilogen_explant_per_sample_logfc_vs_tu <- sapply(colnames(nilogen_explant_log_counts), function (cn) {
  ctrl_col <- gsub("\\s(A|AB|B|0h|PMAI)$", " Tu", cn)
  nilogen_explant_log_counts[,cn] - nilogen_explant_log_counts[,ctrl_col]
})

nilogen_explant_rtg <- rownames(nilogen_explant_logfc)
names(nilogen_explant_rtg) <- rownames(nilogen_explant_logfc)


discard <- (rowSums(expr_matrix) < 20) |
  (rowSums(expr_matrix > 0) < 10) |
  grepl("^(POS|NEG)", rownames(expr_matrix)) |
  rownames(expr_matrix) %in% nanostring_pan_cancer_housekeeping_genes

treatment_order <- c("0h"=0,
                     "PMAI"=1,
                     "Tu"=2,
                     "A"=3,
                     "B"=4,
                     "AB"=5)

nilogen_explant_treatment_order <- order(treatment_order[nilogen_explant_treatment],
                                         nilogen_explant_tumor_number)
nilogen_explant_tumor_order <- order(nilogen_explant_tumor_number,
                                     treatment_order[nilogen_explant_treatment])


########################################################################################################
# use nilogen_explant_per_sample_logfc_vs_tu
# take Amit's suggestion, remove PMAI
# plot heatmaps, do limma 

########################################################################################################
nano_sample_logfc <- nilogen_explant_per_sample_logfc_vs_tu
#remove PMAI and Tu, then we get same
nano_sample_logfc <- nano_sample_logfc[,-grep("0h$|Tu$|PMAI$",colnames(nano_sample_logfc))]
tu_number <- nilogen_explant_tumor_number[-grep("0h$|Tu$|PMAI$",colnames(nilogen_explant_per_sample_logfc_vs_tu))]

#not working
pheatmap::pheatmap(
  nano_sample_logfc,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-10,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1,1.5,2,3,10),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "nanostring_logfc_sample",
  border_color = FALSE
)
#get clusters
library(pheatmap)
res <- pheatmap(nano_sample_logfc)
res.clust <- cbind(nano_sample_logfc, 
                      cluster = cutree(res$tree_row, 
                                       k = 10))
head(res.clust)
anno <- data.frame(cluster = cluster)
pheatmap( nano_sample_logfc,
          annotation_row = anno)
##get gene in cluster 1
test_gene <- rownames(anno)[anno$cluster %in% c(1)]
test_nano <- nano_sample_logfc[test_gene,]
pheatmap::pheatmap(
  test_nano,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-10,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1,1.5,2,3,10),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "testnano_logfc_sample",
  border_color = FALSE
)



##############do it again
res <- pheatmap( test_nano)
res.clust <- cbind(test_nano, 
                   cluster = cutree(res$tree_row, 
                                    k = 8))
head(res.clust)
anno <- data.frame(cluster = cluster)
##get gene in cluster 1
test_gene2 <- rownames(anno)[anno$cluster %in% c(1)]

test_nano2 <- nano_sample_logfc[test_gene2,]
pheatmap::pheatmap(
  test_nano2,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-10,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1,1.5,2,3,10),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "testnano2_logfc_sample",
  border_color = FALSE
)

res <- pheatmap( test_nano2)
res.clust <- cbind(test_nano2, 
                   cluster = cutree(res$tree_row, 
                                    k = 8))
head(res.clust)
anno <- data.frame(cluster = cluster)
##get gene in cluster 1
test_gene3 <- rownames(anno)[anno$cluster ==1]
test_nano3 <- nano_sample_logfc[test_gene3,]
pheatmap::pheatmap(
  test_nano3,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-10,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1,1.5,2,3,10),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "testnano2_logfc_sample",
  border_color = FALSE
)

res <- pheatmap( test_nano3)
res.clust <- cbind(test_nano3, 
                   cluster = cutree(res$tree_row, 
                                    k = 8))
head(res.clust)
anno <- data.frame(cluster = cluster)

test_gene4 <- rownames(anno)[anno$cluster ==2]
test_nano4 <- nano_sample_logfc[test_gene4,]
pheatmap::pheatmap(
  test_nano4,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-10,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1,1.5,2,3,10),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "testnano2_logfc_sample",
  border_color = FALSE
)

res <- pheatmap( test_nano5)
res.clust <- cbind(test_nano5, 
                   cluster = cutree(res$tree_row, 
                                    k = 8))
head(res.clust)
anno <- data.frame(cluster = cluster)
##get gene in cluster 1
test_gene5 <- rownames(anno)[anno$cluster ==1]
test_nano5 <- nano_sample_logfc[test_gene5,]
pheatmap::pheatmap(
  test_nano5,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-10,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1,1.5,2,3,10),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "testnano2_logfc_sample",
  border_color = FALSE
)

res <- pheatmap( test_nano5)
res.clust <- cbind(test_nano5, 
                   cluster = cutree(res$tree_row, 
                                    k = 8))
head(res.clust)
anno <- data.frame(cluster = cluster)
##get gene in cluster 1
test_gene6 <- rownames(anno)[anno$cluster ==2]
test_nano6 <- nano_sample_logfc[test_gene6,]
pheatmap::pheatmap(
  test_nano6,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-10,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1,1.5,2,3,10),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "testnano2_logfc_sample",
  border_color = FALSE
)
#############################make a gene list
gene_list <- c("FOXJ1","FEZ1","C4BPA","BCL2","CFI","NOD1","PLA2G6","FLT3LG","INPP5D","MR1",
               "NOD2","IRF3","HRAS","SIGLEC1","CFI","CC2D1B","RRAD","ITGAE","CD3G","FLT3LG",
               "IL2RB","PECAM1","CD79B","GZMK","IL2RB","MR1","CD3G","SELPLG","IL8","ITGAM",
               "TIGIT","LY9","DUSP6","CASP1","KLRD1","CMKLR1","CLEC7A","TNF","CSF2RB","IL8",
               "CXCL2","IL32","TARP","EGR1","PTGS2","CD14","C3AR1","CCR3","PTGS2","THBS1",
               "MRC1","CCL26","CCL5","CD3E","ICAM3","CCL20","PTGS2","EGR1","IRAK2","NFKB2",
               "CD5","TNFSF12",
               "ITGA4","PSMB9","IL1RN","ATG7","IKBKE","CD3D","TNFSF15","PSMB9","TAP2","CSF1","CNOT10",
               "CTSH","IL15RA","RIPK2","IFI16","ICAM1",
               "TLR10","CX3CL1","MAPK11","CXCL6","CMA1","CSF2","CD1D","IL1B","IDO1","CD40","CD80","GZMB","PRF1","GZMA","COLEC12",
               "DDX43","LILRB3","FCER2","HRAS","CKLF","CD3G","EGR1","PTGS2","TNFSF12","CFP","NFKBIA","S100A8","REL","TANK",
               "IL22RA2","IL21R","F12","CCL21","SLAMF1","FOXP3","TNFRSF18","IL3","MS4A1","PPBP","GZMA","XCL2","PRF1",
               "DDX43","IL11RA","ARG1",
               "SLAMF6","CREBBP","CXCR2","GZMH","TIRAP","IL12RB2","CCR2","PPBP","CXCL3","CXCL5","CXCL6","DDX43",
               "COLEC12","CXCR5;","KLRF1")
#HEATMAP WITH SELECT GENE
select_gene <- nano_sample_logfc[rownames(nano_sample_logfc) %in% gene_list,]
pheatmap::pheatmap(
  select_gene,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-10,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1,1.5,2,3,10),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "select gene",
  border_color = FALSE
)


#########################################try perform logestic regression
#label first:
clust <- c(1,2,1,3,3,3,4,4,5,6,6,6,7,7,8,9,10,10,11,11,11,12,12,12,13,13,13,14,14,14)
library(MASS)
library(VGAM)
test <- data.frame(t(nano_sample_logfc))
test$clust <- as.factor(clust)
#remove constant
test <- test[,-c(60,219,628,653,737,740,742,744)]
fit <- lda(clust~.,test)
fit2 <- vglm(clust~.,family = multinomial, data=test)
plda <- predict(object = fit,
                newdata = test)
pca <- prcomp(test[,-777],
              center = T)
dataset = data.frame(clust = clust,
                     pca = pca$x, lda = plda$x)
require(ggplot2)
require(scales)
require(gridExtra)

 ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = clust), size = 1) 

#TRY CLUSTER ON GENE
 pheatmap::pheatmap(
   t(nano_sample_logfc),
   cluster_cols = T,
   cluster_rows = T,
   breaks = c(-10,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1,1.5,2,3,10),
   color = colorRampPalette(rev(brewer.pal(
     name = 
       "RdBu",
     n=9
   )))(15),
   #annotation_colors = tre_colr[1],
   #annotation_col = tre,
   #gaps_col = 1,
   main = "select gene",
   border_color = FALSE
 )
 

 
 



#######try sig scores
#######confirmed sigs by Alex
confirmed_nanostring_signatures <- c(
  "angelova2015.DC",
  "bindea2013.Neutrophils",
  "bindea2013.B.cells",
  "bindea2013.Cytotoxic",    
  "bindea2013.Macrophages",
  "bindea2013.T.cells",      
  "bindea2013.Treg",
  "S1_NK_cells",             
  "S10_B_cells",
  "angelova2015.Macrophages",
  "msd.asco.2016.TcellExhaustion",    
  "msd.asco.2016.TcellAndNK",
  "msd.asco.2016.Ifng",               
  "msd.asco.2016.AntigenPresentation",
  "msd.asco.2015.TCR",
  "msd.asco.2015.melanoma",
  "msd.asco.2015.Immune",             
  "msd.asco.2015.Ifng",
  "msd.asco.2015.DeNovo"                 
)

source(sprintf("%s/bioinformatics-signatures/load.R", Sys.getenv("GITCLONES")), chdir=T)
library(pheatmap)
rtg <- rownames(nano_sample_logfc)
names(rtg) <- rownames(nano_sample_logfc)
nanostring_sigscores <- scoreSignatures(nano_sample_logfc,
                                        rtg,
                                        signature.list[confirmed_nanostring_signatures])
library(RColorBrewer)
pheatmap::pheatmap(
  t(nanostring_sigscores),
  cluster_cols = T,
  cluster_rows = T,
 breaks = c(-5,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,4),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "sigscore_logfc_nanostring",
  border_color = FALSE
)

##########################################################IHC
raw_ihc <- readODS::read.ods("~/Desktop/Jieqing/jingping/Nilogen/cytokine/analysis/summary_parsed.ods",sheet = 7)
rownames(raw_ihc) <- raw_ihc[,1]
colnames(raw_ihc) <- raw_ihc[1,]
raw_ihc <- raw_ihc[-1,-1]

ihc_p <- raw_ihc[c(6,7),]
ihc_wo_p <- raw_ihc[c(1:5),]
library(RColorBrewer)
pheatmap::pheatmap(
 data.matrix(ihc_wo_p),
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-1.5,-1,-0.75,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.2,0.4,0.6,0.8,1,1.5,4),
 #color = colorRampPalette(rev(brewer.pal(
#    name = 
 #     "RdBu",
  #  n=9
  #)))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "IHC scores",
  border_color = FALSE
)

raw_ldhglu<- readODS::read.ods("~/Desktop/Jieqing/jingping/Nilogen/cytokine/analysis/summary_parsed.ods",sheet = 8)
rownames(raw_ldhglu) <- raw_ldhglu[,1]
colnames(raw_ldhglu) <- raw_ldhglu[1,]
raw_ldhglu <- data.matrix(raw_ldhglu[-1,-1])


ldhglu_change <- sapply(colnames(raw_ldhglu), function (cn) {
  ctrl_col <- gsub("(A|A[+]B|B|PMA/I)$", "Veh", cn)
  raw_ldhglu[,cn] - raw_ldhglu[,ctrl_col]
})
ldhglu_change_remove <-  ldhglu_change[,-grep("Veh$|PMA/I$",colnames(ldhglu_change))]




library(RColorBrewer)
pheatmap::pheatmap(
  data.matrix(ldhglu_change_remove[c(1:2),]),
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-1.5,-1,-0.75,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.2,0.4,0.6,0.8,1,1.5,4),
  #color = colorRampPalette(rev(brewer.pal(
  #    name = 
  #     "RdBu",
  #  n=9
  #)))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "LDH foldchange vs Veh normalized to media",
  border_color = FALSE
)

pheatmap::pheatmap(
  t(scale(t(data.matrix(ldhglu_change_remove[c(3,4),])))),
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-1.5,-1,-0.75,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.2,0.4,0.6,0.8,1,1.5,4),
  #color = colorRampPalette(rev(brewer.pal(
  #    name = 
  #     "RdBu",
  #  n=9
  #)))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Scaled foldchange Glucose",
  border_color = FALSE
)





