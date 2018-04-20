###############################Request from Biologist after presentation##########
#1. veh for all
#2. Treatment specific comparison vs veh, no fc
#3. signatured gene, heatmap plot on all treatment, 0h, pmai and veh
#4. all clust and unclust
##################################################################################


#################1.veh for all####################################################
#for multiplex
raw_multiplex <- readODS::read.ods("~/Desktop/Jieqing/jingping/Nilogen/final/jp_final_data.ods",sheet = 4)
rownames(raw_multiplex) <- make.names(raw_multiplex[,1],unique = T)
raw_multiplex <- raw_multiplex[,-1]
colnames(raw_multiplex) <- raw_multiplex[1,]
raw_multiplex <- raw_multiplex[-1,]
#detect range 10
log_norm_raw <- log(data.matrix(raw_multiplex)+10)
library(RColorBrewer)
###############veh logged multiplex, clust and unclust
setwd("~/Desktop/Jieqing/jingping/Nilogen/final/after")
png("Plex_Veh_logged_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  log_norm_raw[,grep("Veh",colnames(log_norm_raw))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Plex_Veh_logged",
  border_color = FALSE
)
dev.off()
png("Plex_unclust_Veh_logged_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  log_norm_raw[,grep("Veh",colnames(log_norm_raw))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Plex_unclust_Veh_logged",
  border_color = FALSE
)
dev.off()
################################A vs veh, clust and unclust
png("Plex_unclust_A_vs_Veh_logged_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  log_norm_raw[,grep("A$|Veh",colnames(log_norm_raw))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Plex_unclust_A_vs_Veh_logged",
  border_color = FALSE
)
dev.off()

png("Plex_A_vs_Veh_logged_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  log_norm_raw[,grep("A$|Veh",colnames(log_norm_raw))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Plex_A_vs_Veh_logged",
  border_color = FALSE
)
dev.off()

#######################################################B vs Veh multiplex, clust and unclust
png("Plex_unclust_B_vs_Veh_logged_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  log_norm_raw[,grep("[_]B$|Veh",colnames(log_norm_raw))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Plex_unclust_B_vs_Veh_logged",
  border_color = FALSE
)
dev.off()

png("Plex_B_vs_Veh_logged_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  log_norm_raw[,grep("[_]B$|Veh",colnames(log_norm_raw))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Plex_B_vs_Veh_logged",
  border_color = FALSE
)
dev.off()
#####################################################A+B vs Veh, clust and unclust
png("Plex_unclust_AB_vs_Veh_logged_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  log_norm_raw[,grep("A[+]B$|Veh",colnames(log_norm_raw))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Plex_unclust_AB_vs_Veh_logged",
  border_color = FALSE
)
dev.off()

png("Plex_AB_vs_Veh_logged_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  log_norm_raw[,grep("A[+]B$|Veh",colnames(log_norm_raw))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Plex_AB_vs_Veh_logged",
  border_color = FALSE
)
dev.off()
#########################all in one for plex
png("Plex_unclust_all_logged_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  log_norm_raw,
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Plex_unclust_all_logged",
  border_color = FALSE
)
dev.off()

png("Plex_all_logged_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  log_norm_raw,
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Plex_all_logged",
  border_color = FALSE
)
dev.off()


#################################Facs data veh , and veh vs treatment
raw_facs <- readODS::read.ods("~/Desktop/Jieqing/jingping/Nilogen/final/jp_final_data.ods",sheet = 5)
rownames(raw_facs) <- make.names(raw_facs[,1],unique = T)
raw_facs <- raw_facs[,-1]
colnames(raw_facs) <- raw_facs[1,]
raw_facs <- raw_facs[-1,]
log_norm_facs <- log2(data.matrix(raw_facs)+1)
#scale gives a smaller range,use scale
scale_facs <- t(scale(t(data.matrix(raw_facs))))


#veh
png("Facs_Veh_scaled_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  scale_facs[,grep("Veh",colnames(scale_facs))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Facs_Veh_scaled",
  border_color = FALSE
)
dev.off()
png("Facs_unclust_Veh_scaled_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  scale_facs[,grep("Veh",colnames(scale_facs))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Facs_unclust_Veh_scaled",
  border_color = FALSE
)
dev.off()

#a vs veh
png("Facs_A_vs_Veh_scaled_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  scale_facs[,grep("A$|Veh",colnames(scale_facs))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Facs_A_vs_Veh_scaled",
  border_color = FALSE
)
dev.off()

##A vs veh
png("Facs_unclust_A_vs_Veh_scaled_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  scale_facs[,grep("A$|Veh",colnames(scale_facs))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Facs_unclust_A_vs_Veh_scaled",
  border_color = FALSE
)
dev.off()

#B vs veh
png("Facs_B_vs_Veh_scaled_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  scale_facs[,grep("[_]B$|Veh",colnames(scale_facs))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Facs_B_vs_Veh_scaled",
  border_color = FALSE
)
dev.off()

png("Facs_unclust_B_vs_Veh_scaled_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  scale_facs[,grep("[_]B$|Veh",colnames(scale_facs))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Facs_unclust_B_vs_Veh_scaled",
  border_color = FALSE
)
dev.off()

#veh vs A+B

png("Facs_AB_vs_Veh_scaled_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  scale_facs[,grep("A[+]B$|Veh",colnames(scale_facs))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Facs_AB_vs_Veh_scaled",
  border_color = FALSE
)
dev.off()

png("Facs_unclust_AB_vs_Veh_scaled_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  scale_facs[,grep("A[+]B$|Veh",colnames(scale_facs))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Facs_unclust_AB_vs_Veh_scaled",
  border_color = FALSE
)
dev.off()

###########facs all
png("Facs_all_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  scale_facs,
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Facs_all_scaled",
  border_color = FALSE
)
dev.off()

png("Facs_unclust_all_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  scale_facs[,grep("[_]B$|Veh",colnames(scale_facs))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Facs_unclust_all_scaled",
  border_color = FALSE
)
dev.off()


################################NANOstring
#####################################load first

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

#use nilogen_explant_log_counts
#0h only

png("Nano_unclust_0h_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[,grep("0h$",colnames(nilogen_explant_log_counts))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_unclust_0h_logged",
  border_color = FALSE
)
dev.off()

png("Nano_0h_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[,grep("0h$",colnames(nilogen_explant_log_counts))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_0h_logged",
  border_color = FALSE
)
dev.off()


###################select signatures gene for 0h
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
select_confirmed_nanostring_signatures <- c(
  "angelova2015.DC",
  "angelova2015.Macrophages",
  "bindea2013.Treg",
  "msd.asco.2016.Ifng", 
  "bindea2013.T.cells", 
  "bindea2013.B.cells"
)
source(sprintf("%s/bioinformatics-signatures/load.R", Sys.getenv("GITCLONES")), chdir=T)
sig_genes <- signature.list[select_confirmed_nanostring_signatures]
sig_genes <- unlist(sig_genes)

png("Nano_unclust_0h_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("0h$",colnames(nilogen_explant_log_counts))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_unclust_0h_select_gene_logged",
  border_color = FALSE
)
dev.off()

png("Nano_0h_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("0h$",colnames(nilogen_explant_log_counts))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_0h_select_gene_logged",
  border_color = FALSE
)
dev.off()


###############################veh and 0h
png("Nano_unclust_veh_0h_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("Tu$|0h$",colnames(nilogen_explant_log_counts))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_unclust_veh$_0h_select_gene_logged",
  border_color = FALSE
)
dev.off()

png("Nano_Veh_0h_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("Tu$|0h$",colnames(nilogen_explant_log_counts))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_Veh_0h_select_gene_logged",
  border_color = FALSE
)
dev.off()


###############################A vs veh

png("Nano_unclust_veh_A_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("A$|Tu$",colnames(nilogen_explant_log_counts))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_unclust_veh_A_select_gene_logged",
  border_color = FALSE
)
dev.off()

png("Nano_Veh_A_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("A$|Tu$",colnames(nilogen_explant_log_counts))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_Veh_A_select_gene_logged",
  border_color = FALSE
)
dev.off()

###################################B vs Veh
png("Nano_unclust_veh_B_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("[ ]B$|Tu$",colnames(nilogen_explant_log_counts))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_unclust_veh_B_select_gene_logged",
  border_color = FALSE
)
dev.off()

png("Nano_Veh_B_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("[ ]B$|Tu$",colnames(nilogen_explant_log_counts))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_Veh_B_select_gene_logged",
  border_color = FALSE
)
dev.off()

#######################AB vs Veh
png("Nano_unclust_veh_AB_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("AB$|Tu$",colnames(nilogen_explant_log_counts))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_unclust_veh_AB_select_gene_logged",
  border_color = FALSE
)
dev.off()

png("Nano_Veh_AB_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("AB$|Tu$",colnames(nilogen_explant_log_counts))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_Veh_AB_select_gene_logged",
  border_color = FALSE
)
dev.off()

####################################PMAI only
png("Nano_unclust_PMAI_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("PMAI$",colnames(nilogen_explant_log_counts))],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_unclust_PMAI_select_gene_logged",
  border_color = FALSE
)
dev.off()

png("Nano_PMAI_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),grep("PMAI$",colnames(nilogen_explant_log_counts))],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_PMAI_select_gene_logged",
  border_color = FALSE
)
dev.off()

#################################################all in one for nano
png("Nano_unclust_all_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),],
  cluster_cols = F,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_unclust_all_select_gene_logged",
  border_color = FALSE
)
dev.off()

png("Nano_all_select_gene_heatmap.png", width=1200, height=1000)
pheatmap::pheatmap(
  nilogen_explant_log_counts[rownames(nilogen_explant_log_counts) %in% unname(sig_genes),],
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Nano_all_select_gene_logged",
  border_color = FALSE
)
dev.off()








##################################################try perform WGCNA
oh_logged <- nilogen_explant_log_counts[,grep("0h$",colnames(nilogen_explant_log_counts))]
pmai_logged <- nilogen_explant_log_counts[,grep("PMAI$",colnames(nilogen_explant_log_counts))]
TA_logfc <- nilogen_explant_per_sample_logfc_vs_tu[,grep("A$",colnames(nilogen_explant_per_sample_logfc_vs_tu))]
TB_logfc <-  nilogen_explant_per_sample_logfc_vs_tu[,grep("[ ]B$",colnames(nilogen_explant_per_sample_logfc_vs_tu))]
TAB_logfc <- nilogen_explant_per_sample_logfc_vs_tu[,grep("AB$",colnames(nilogen_explant_per_sample_logfc_vs_tu))]
allT_logfc <- nilogen_explant_per_sample_logfc_vs_tu[,-grep("0h$|PMAI$|Tu$",colnames(nilogen_explant_per_sample_logfc_vs_tu))]
tumor_number <- c(1,10,paste(seq(from=2,to=9)))
tumor_number_all <- c(rep(1,3),rep(10,3),rep(seq(from=2,to=9),each=3))

library(WGCNA)
library(flashClust)
library(doMC)
library("parallel")
Cl <- makeCluster(detectCores()-1)
load("/home/x189259/GITCLONES/bioinformatics-signatures/bioinformatics-signatures.rda")


soft <- NULL
gsg <- NULL
datExpr = NULL
firstWGCNA <- function(datExpr){
  registerDoMC(8)
  datExpr2 <-datExpr
  ##reverse dataframe so it fit feeds
  row.names(datExpr) = datExpr$X
  datExpr$X = NULL
  datExpr = as.data.frame(t(datExpr)) # now samples are rows and genes are columns
  dim(datExpr)
  
  ## Run this to check if there are gene outliers
  gsg = goodSamplesGenes(datExpr, verbose = 3)
  gsg$allOK 
  
  #If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
  if (!gsg$allOK)
  {if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
    datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  # Call the network topology analysis function
  soft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  yd <- -sign(soft$fitIndices[,3])*soft$fitIndices[,2]
  if(max(yd>0.9)){
    index_yd <- which(yd>0.9)
    # id <- diff(index_yd)
    # if (all(id == c(rep(1,length(id))))  ){
    sft <- soft$fitIndices$Power[index_yd[1]]}else{
      index_yd <- which(yd>0.8)
      sft <- soft$fitIndices$Power[index_yd[1]]
    }
  print(paste("The first WGCNA softPower we choose is :", sft))
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  # Mean connectivity as a function of the soft-thresholding power
  plot(soft$fitIndices[,1], soft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(soft$fitIndices[,1], soft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  
  ## Construct Networks------------------------------------------------------------------------------
  #seems that 4 would be a nice sft power. build a adjacency "correlation" matrix with power=4
  enableWGCNAThreads()
  softPower = sft
  adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type
  head(adjacency)
  
  # Construct Networks
  #translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
  TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
  dissTOM = 1-TOM
  
  # Generate Modules ------------------------------------------------------------------
  # Generate a clustered gene tree
  geneTree = flashClust(as.dist(dissTOM), method="average")
  plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
  #This sets the minimum number of genes to cluster into a module
  minModuleSize = 30
  dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=4, pamRespectsDendro= TRUE, minClusterSize = minModuleSize)
  dynamicColors= labels2colors(dynamicMods)
  MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = sft)
  MEs= MEList$eigengenes
  MEDiss= 1-cor(MEs)
  METree= flashClust(as.dist(MEDiss), method= "average")
  
  
  #plots tree showing how the eigengenes cluster together
  png(file="clusterw_gene.png")
  plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
  #set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
  MEDissThres = 0.0
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
  mergedColors = merge$colors
  mergedMEs = merge$newMEs
  dev.off()
  
  #plot dendrogram with module colors below it
  png(file="cluster_gene.png")
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
  moduleColors = mergedColors
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  dev.off()
  
  
  
  ###########################################################################################
  #use WGCNA again to cluster Jingping_cluster, mgsidb, and built in sig
  #sort(as.character(stopsle.rtg[cluster_list[["Jingping_7"]]]))
  
  #found that these are 11 clusters in CD16, so we label clusters with name, and change emsbl to symbles
  names_clusters <- dynamicMods
  genes <- as.vector(rownames(datExpr2)[gsg$goodGenes])
  names(genes) <- names_clusters
  cluster_list <- as.list(split(genes,names(genes)))
  cluster_list <- lapply(1:length(cluster_list),function(x) unname(cluster_list[[x]]))
  cluster_list <- cluster_list[order(sapply(cluster_list,length),decreasing=T)]
  
  #since the order was changed, manully change names
  #check order with summary(cluster_list)
  
  a <-  NULL
  for (i in 1: length(cluster_list)){
    b <- paste("Jingping_",i,sep = "")
    a <- c(a,b)
  }

  names(cluster_list) <- a
  for( i in 1:length(cluster_list)){
    cluster_list[[i]]  = as.character(nilogen_explant_rtg[cluster_list[[i]]])
  }
  y2 <- cluster_list
return(y2)
}
#test removed genes
#table(unname(unlist(test)) %in% rownames(TA_logfc)[gsg$goodGenes] )
#table(rownames(TA_logfc)[gsg$goodGenes] %in% unname(unlist(test)))
pheatmap::pheatmap(
  allT_logfc[gene_wgcna,],
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-9,-2,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.5,5.5),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "test for TA WGCNA",
  border_color = FALSE
)
gene_wgcna <- c("CXCL10","PPBP","CXCL1","IL1B","CCL7","CCL3L1","CXCL3","CCL2","CXCL5","SELL")
gene_test <- c("CCL2","CCL7","CCL3L1","CXCL5","CXCL1","CXCL3",
               "IL1B","CXCL10","CXCL15")
#"LILRA5","SELL","RORC",,"FCGR1A","PPBP","CCL8","CXCL10",,"HSD11B1","IL21""LILRA5","SELL",
#"HSD11B1","FCGR1A", "RORC","IL21","LILRA5","CXCL10","PPBP","CXCL1","IL1B","CCL7","CCL3L1","CXCL3","CCL2","CXCL5","SELL","CCL8"
#test this cluster on logfoldchange

wgcna <- list(gene_wgcna)
names(wgcna) <- "wgcna"
new.sig <- c(signature.list,wgcna)

confirmed_nanostring_signatures[20] <- "wgcna"
rtg <- rownames(allT_logfc)
names(rtg) <- rownames(allT_logfc)
nanostring_sigscores <- scoreSignatures(allT_logfc,
                                        rtg,
                                        new.sig[confirmed_nanostring_signatures])
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



#try sig genes
new_sig_genes <- sig_genes[unname(sig_genes) %in% rownames(allT_logfc)]
pheatmap::pheatmap(
  allT_logfc[new_sig_genes,],
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
  main = "select_genes_logfc_nanostring",
  border_color = FALSE
)
#####################write list into txt
library(marray)
write.list(signature.list[select_confirmed_nanostring_signatures],filename = "responding_nanostring_gene_from_signature_score.txt")
write.list(signature.list[confirmed_nanostring_signatures],filename = "all_signatures_used_in_heatmaps.txt")



















########################################heatmaps for each signatures####################
setwd("~/Desktop/Jieqing/jingping/Nilogen/final/after/sig_heatmaps")
#0h
baseline <- nilogen_explant_log_counts[,grep("0h",colnames(nilogen_explant_log_counts))]
#treatments and veh
treatments <- nilogen_explant_log_counts[,-grep("0h|PMAI",colnames(nilogen_explant_log_counts))]

plotHeatmap <- function(df){
for (i in 1:length(confirmed_nanostring_signatures)){
  if(length(unname(unlist(signature.list[confirmed_nanostring_signatures[i]])))>1){
  test <- df[rownames(df) %in% unname(unlist(signature.list[confirmed_nanostring_signatures[i]])),]
 png(filename = paste(confirmed_nanostring_signatures[i],".png"))
  pheatmap::pheatmap(
    test,
    cluster_cols = T,
    cluster_rows = T,
    breaks = c(0,1,3,4,4.5,5,5.5,6,6.5,7,7.5,8,9,15,21),
    color = colorRampPalette(rev(brewer.pal(
      name = 
        "RdBu",
      n=9
    )))(15),
    #annotation_colors = tre_colr[1],
    #annotation_col = tre,
    #gaps_col = 1,
    main = paste(confirmed_nanostring_signatures[i],"all Treatment, Gene Specific Log count"),
    border_color = FALSE
  )
  dev.off()
  
}else{
  print(paste(confirmed_nanostring_signatures[i],"cannot print"))
}
}
}
  
  
  
##########################################multiplex
#########################################try compare a vs b, a vs a+b
raw_multiplex <- readODS::read.ods("/home/x189259/Desktop/Jieqing/jingping/Nilogen/final/Nilogen_jp_push/jp_edited_final_data.ods",sheet = 4)
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
logfc_b_vs_a <- sapply(colnames(cyto_logfc_vs_veh), function (cn) {
  ctrl_col <- gsub("(A[+]B|B)$", "A", cn)
  cyto_logfc_vs_veh[,cn] - cyto_logfc_vs_veh[,ctrl_col]
})
logfc_b_vs_a <- logfc_b_vs_a[,-grep("A$",colnames(logfc_b_vs_a))]

  
pheatmap::pheatmap(
  logfc_b_vs_a,
  cluster_cols = F,
  cluster_rows = T,
  breaks = c(-4.5,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.1,0.2,0.3,0.4,0.5,0.7,1.1),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "logfc_vs_logfc_veh-A",
  border_color = FALSE
)
  
  
###################################################nanostring
##################################################logfc_veh_0h
logfc_0h <- sapply(colnames(nilogen_explant_log_counts), function (cn) {
  ctrl_col <- gsub("(PMAI|0h|Tu|A|AB|B)$", "0h", cn)
  nilogen_explant_log_counts[,cn] - nilogen_explant_log_counts[,ctrl_col]
})
logfc_Tu_0h <- logfc_0h[,grep("Tu$",colnames(logfc_0h))]
select_gene <- c("CCL3","CCL4","IL6","CSF3")
pheatmap::pheatmap(
  logfc_Tu_0h,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-10,-5,-1.5,-1,-0.75,-0.5,-0.25,0,0.1,0.2,0.3,0.4,0.5,1,6,15),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "Tu_logfc_vs_0h",
  border_color = FALSE
)
###################################################try logfc vs veh vs a

logfc_a <- nilogen_explant_log_counts[,-grep("(PMAI|Tu|0h)$",colnames(nilogen_explant_log_counts))]
logfc_a <- sapply(colnames(logfc_a), function (cn) {
  ctrl_col <- gsub("(AB|B)$", "A", cn)
  logfc_a[,cn] - logfc_a[,ctrl_col]
})
logfc_a <- logfc_a[,-grep("A$",colnames(logfc_a))]
select_gene <- c("CCL3","CCL4","IL6","CSF3")

pheatmap::pheatmap(
  logfc_a,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-5.1,-1.5,-1,-0.75,-0.5,-0.4,-0.2,0,0.1,0.2,0.3,0.4,0.5,0.8,1,5.1),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "double logfc vs a",
  border_color = FALSE
)
rtg <- rownames(logfc_a)
names(rtg) <- rownames(logfc_a)
test_sig <- scoreSignatures(logfc_a,rtg,signatures = signature.list[confirmed_nanostring_signatures])
pheatmap::pheatmap(
  t(test_sig),
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-1.5,-1,-0.75,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.8,1.4,2.0),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = "double logfc vs a",
  border_color = FALSE
)








#######################try wcgna

logfc <- nilogen_explant_per_sample_logfc_vs_tu[,-grep("0PMAI|h|Tu$",colnames(nilogen_explant_per_sample_logfc_vs_tu))]
test_logfc <- firstWGCNA(data.frame(logfc))



for (i in 1:length(test_logfc)){
pheatmap::pheatmap(
  logfc_a[unname(unlist(test_logfc[i])),],
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-5.1,-1.5,-1,-0.75,-0.5,-0.4,-0.2,0,0.1,0.2,0.3,0.4,0.5,0.8,1,5.1),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  #annotation_col = tre,
  #gaps_col = 1,
  main = paste("cluster",i,"double logfc vs a"),
  border_color = FALSE
)
}

#################################################perform paired test on logfc data
nilogen_logfc <- cyto_logfc_vs_veh
a_samples <- grep("A$",colnames(nilogen_logfc))
b_samples <- grep("[_]B$",colnames(nilogen_logfc))
ab_samples <- grep("A[+]B$",colnames(nilogen_logfc))
pval_a <- sapply(rownames(nilogen_logfc), function (rn) {
  wilcox.test(nilogen_logfc[rn, a_samples], nilogen_logfc[rn, -a_samples],alternative = "greater")$p.value
})
#gene list affect m score
names(pval_a[pval_a <0.05])
pval_less_a <- sapply(rownames(nilogen_logfc), function (rn) {
  wilcox.test(nilogen_logfc[rn, a_samples], nilogen_logfc[rn, -a_samples],alternative = "less")$p.value
})
names(pval_less_a[pval_less_a <0.05])

pval_b <- sapply(rownames(nilogen_logfc), function (rn) {
  wilcox.test(nilogen_logfc[rn, b_samples], nilogen_logfc[rn, -b_samples],alternative = "greater")$p.value
})
#gene list affect m score
names(pval_b[pval_b <0.05])
pval_less_b <- sapply(rownames(nilogen_logfc), function (rn) {
  wilcox.test(nilogen_logfc[rn, b_samples], nilogen_logfc[rn, -b_samples],alternative = "less")$p.value
})
names(pval_less_b[pval_less_b <0.05])

pval_ab <- sapply(rownames(nilogen_logfc), function (rn) {
  wilcox.test(nilogen_logfc[rn, ab_samples], nilogen_logfc[rn, -ab_samples],alternative = "greater")$p.value
})
#gene list affect m score
names(pval_ab[pval_ab <0.05])
pval_less_ab <- sapply(rownames(nilogen_logfc), function (rn) {
  wilcox.test(nilogen_logfc[rn, ab_samples], nilogen_logfc[rn, -ab_samples],alternative = "less")$p.value
})
names(pval_less_ab[pval_less_ab <0.05])


#try to label 3b,5ab,9ab,4ab,7ab
respond <- rep(0,dim(nilogen_logfc)[2])
respond[grep("(5[_]A[+]B|9[_]A[+]B|7[_]A[+]B|4[_]A[+]B)$",colnames(nilogen_logfc))] <- 1
pval_g <- sapply(rownames(nilogen_logfc), function (rn) {
  wilcox.test(nilogen_logfc[rn, grep("(5[_]A[+]B|9[_]A[+]B|7[_]A[+]B|4[_]A[+]B)$",colnames(nilogen_logfc))], nilogen_logfc[rn, -grep("(3[_]B|5[_]A[+]B|9[_]A[+]B|7[_]A[+]B|4[_]A[+]B)$",colnames(nilogen_logfc))],alternative = "greater")$p.value
})
names(pval_g[pval_g <0.05])
pval_g_less <- sapply(rownames(nilogen_logfc), function (rn) {
  wilcox.test(nilogen_logfc[rn, grep("(5[_]A[+]B|9[_]A[+]B|7[_]A[+]B|4[_]A[+]B)$",colnames(nilogen_logfc))], nilogen_logfc[rn, -grep("(5[_]A[+]B|9[_]A[+]B|7[_]A[+]B|4[_]A[+]B)$",colnames(nilogen_logfc))],alternative = "less")$p.value
})
names(pval_g_less[pval_g_less <0.05])


########################plot treatment
ta <- rep(0,dim(nilogen_logfc)[2])
ta[grep("A$",colnames(nilogen_logfc))] <- 1
tb <- rep(0,dim(nilogen_logfc)[2])
tb[grep("[_]B$",colnames(nilogen_logfc))] <- 1
tab <- rep(0,dim(nilogen_logfc)[2])
tab[grep("[_]A[+]B)$",colnames(nilogen_logfc))] <- 1

A <- data.frame(A=as.factor(ta))
rownames(A) <- colnames(nilogen_logfc)
B <- data.frame(B=as.factor(tb))
rownames(B) <- colnames(nilogen_logfc)
AB <- data.frame(AB=as.factor(tab))
rownames(AB) <- colnames(nilogen_logfc)
pheatmap::pheatmap(
  nilogen_logfc[names(pval_a[pval_a <0.05]),],
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-16,-6,-3,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,3,14),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment A tested cytokines",
  annotation_col = A
)

pheatmap::pheatmap(
  nilogen_logfc[names(pval_less_b[pval_less_b <0.05]),],
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-5,-3,-1,-0.75,-0.5,-0.25,-0.1,0,0.1,0.2,0.4,0.6,0.8,1,3),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment B tested cytokines",
  annotation_col = B
)


########################treatment only
Nilogen_TA <- nilogen_logfc[,grep("[_]A$",colnames(nilogen_logfc))]
Nilogen_TB <- nilogen_logfc[,grep("[_]B$",colnames(nilogen_logfc))]
Nilogen_TAB <- nilogen_logfc[,grep("A[+]B$",colnames(nilogen_logfc))]
pheatmap::pheatmap(
  Nilogen_TA,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-4,-1,-0.75,-0.5,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.6,0.8,1,3),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment A"
)
pheatmap::pheatmap(
  Nilogen_TB,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-4,-1,-0.75,-0.5,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.6,0.8,1,3),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment B"
)
pheatmap::pheatmap(
  Nilogen_TAB,
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-4,-1,-0.75,-0.5,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.6,0.8,1,3),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment AB"
)


##################################################try same thing on facs

raw_facs <- readODS::read.ods("/home/x189259/Desktop/Jieqing/jingping/Nilogen/final/Nilogen_jp_push/jp_edited_final_data.ods",sheet = 5)
colnames(raw_facs) <- raw_facs[1,]
rownames(raw_facs) <- raw_facs[,1]
raw_facs <- raw_facs[-1,-1]
log_facs <- log2(data.matrix(raw_facs)+1)
facs_logfc_vs_veh <- sapply(colnames(log_facs), function (cn) {
  ctrl_col <- gsub("(A|A[+]B|B)$", "Veh", cn)
  log_facs[,cn] - log_facs[,ctrl_col]
})
facs_logfc_vs_veh <- facs_logfc_vs_veh[,-grep("Veh$",colnames(facs_logfc_vs_veh))]
ta <- rep(0,dim(facs_logfc_vs_veh)[2])
ta[grep("A$",colnames(facs_logfc_vs_veh))] <- 1
tb <- rep(0,dim(facs_logfc_vs_veh)[2])
tb[grep("[_]B$",colnames(facs_logfc_vs_veh))] <- 1
tab <- rep(0,dim(facs_logfc_vs_veh)[2])
tab[grep("A[+]B$",colnames(facs_logfc_vs_veh))] <- 1



pval_a <- sapply(rownames(facs_logfc_vs_veh), function (rn) {
  wilcox.test(facs_logfc_vs_veh[rn, grep("A$",colnames(facs_logfc_vs_veh))], facs_logfc_vs_veh[rn, - grep("A$",colnames(facs_logfc_vs_veh))],alternative = "greater")$p.value
})
names(pval_a[pval_a <0.05])
pval_a_less <- sapply(rownames(facs_logfc_vs_veh), function (rn) {
  wilcox.test(facs_logfc_vs_veh[rn, grep("A$",colnames(facs_logfc_vs_veh))], facs_logfc_vs_veh[rn, - grep("A$",colnames(facs_logfc_vs_veh))],alternative = "less")$p.value
})
names(pval_a_less[pval_a_less <0.05])


pval_b <- sapply(rownames(facs_logfc_vs_veh), function (rn) {
  wilcox.test(facs_logfc_vs_veh[rn, grep("[_]B$",colnames(facs_logfc_vs_veh))], facs_logfc_vs_veh[rn, - grep("[_]B$",colnames(facs_logfc_vs_veh))],alternative = "greater")$p.value
})
names(pval_b[pval_b <0.05])
pval_b_less <- sapply(rownames(facs_logfc_vs_veh), function (rn) {
  wilcox.test(facs_logfc_vs_veh[rn, grep("[_]B$",colnames(facs_logfc_vs_veh))], facs_logfc_vs_veh[rn, - grep("[_]B$",colnames(facs_logfc_vs_veh))],alternative = "less")$p.value
})
names(pval_b_less[pval_b_less <0.05])



pval_ab <- sapply(rownames(facs_logfc_vs_veh), function (rn) {
  wilcox.test(facs_logfc_vs_veh[rn, grep("A[+]B$",colnames(facs_logfc_vs_veh))], facs_logfc_vs_veh[rn, - grep("A[+]B$",colnames(facs_logfc_vs_veh))],alternative = "greater")$p.value
})
names(pval_ab[pval_ab <0.05])
pval_ab_less <- sapply(rownames(facs_logfc_vs_veh), function (rn) {
  wilcox.test(facs_logfc_vs_veh[rn, grep("A[+]B$",colnames(facs_logfc_vs_veh))], facs_logfc_vs_veh[rn, - grep("A[+]B$",colnames(facs_logfc_vs_veh))],alternative = "less")$p.value
})
names(pval_ab_less[pval_ab_less <0.05])






#conclusion:
names(pval_a[pval_a <0.05])
names(pval_b_less[pval_b_less <0.05])



################nano
nanologfc <- nilogen_explant_per_sample_logfc_vs_tu[,-grep("(0h|PMAI|Tu)$",colnames(nilogen_explant_per_sample_logfc_vs_tu))]
pval_a <- sapply(rownames(nanologfc), function (rn) {
  wilcox.test(nanologfc[rn, grep("A$",colnames(nanologfc))], nanologfc[rn, - grep("A$",colnames(nanologfc))],alternative = "greater")$p.value
})
names(pval_a[pval_a <0.05])
pval_a_less <- sapply(rownames(nanologfc), function (rn) {
  wilcox.test(nanologfc[rn, grep("A$",colnames(nanologfc))], nanologfc[rn, - grep("A$",colnames(nanologfc))],alternative = "less")$p.value
})
names(pval_a_less[pval_a_less <0.05])

pval_b <- sapply(rownames(nanologfc), function (rn) {
  wilcox.test(nanologfc[rn, grep("[ ]B$",colnames(nanologfc))], nanologfc[rn, - grep("[ ]B$",colnames(nanologfc))],alternative = "greater")$p.value
})
names(pval_b[pval_b <0.05])
pval_b_less <- sapply(rownames(nanologfc), function (rn) {
  wilcox.test(nanologfc[rn, grep("[ ]B$",colnames(nanologfc))], nanologfc[rn, - grep("[ ]B$",colnames(nanologfc))],alternative = "less")$p.value
})
names(pval_b_less[pval_b_less <0.05])


pval_ab <- sapply(rownames(nanologfc), function (rn) {
  wilcox.test(nanologfc[rn, grep("AB$",colnames(nanologfc))], nanologfc[rn, - grep("AB$",colnames(nanologfc))],alternative = "greater")$p.value
})
names(pval_ab[pval_ab <0.05])
pval_ab_less <- sapply(rownames(nanologfc), function (rn) {
  wilcox.test(nanologfc[rn, grep("AB$",colnames(nanologfc))], nanologfc[rn, - grep("AB$",colnames(nanologfc))],alternative = "less")$p.value
})
names(pval_ab_less[pval_ab_less <0.05])

A <- data.frame(A=as.factor(ta))
rownames(A) <- colnames(nanologfc)
B <- data.frame(B=as.factor(tab))
rownames(B) <- colnames(nanologfc)
AB <- data.frame(AB=as.factor(tb))
rownames(AB) <- colnames(nanologfc)
library(RColorBrewer)
pheatmap::pheatmap(
  nanologfc[names(pval_a[pval_a <0.05]),],
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-8,-3,-1,-0.75,-0.5,-0.3,-0.15,0,0.2,0.4,0.6,0.8,1,2,6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment A tested nanostring (greater)",
  annotation_col = A
)
pheatmap::pheatmap(
  nanologfc[names(pval_a_less[pval_a_less <0.05]),],
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-8,-3,-1,-0.75,-0.5,-0.3,-0.15,0,0.2,0.4,0.6,0.8,1,2,6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment A tested nanostring (less)",
  annotation_col = A
)
pheatmap::pheatmap(
  nanologfc[names(pval_b[pval_b <0.05]),],
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-8,-3,-1,-0.75,-0.5,-0.3,-0.15,0,0.2,0.4,0.6,0.8,1,2,6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment B tested nanostring (greater)",
  annotation_col = B
)
pheatmap::pheatmap(
  nanologfc[names(pval_b_less[pval_b_less <0.05]),],
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-8,-3,-1,-0.75,-0.5,-0.3,-0.15,0,0.2,0.4,0.6,0.8,1,2,6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment B tested nanostring (less)",
  annotation_col = B
)

pheatmap::pheatmap(
  nanologfc[names(pval_ab[pval_ab <0.05]),],
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-8,-3,-1,-0.75,-0.5,-0.3,-0.15,0,0.2,0.4,0.6,0.8,1,2,6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment AB tested nanostring (greater)",
  annotation_col = AB
)
pheatmap::pheatmap(
  nanologfc[names(pval_ab_less[pval_ab_less <0.05]),],
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-8,-3,-1,-0.75,-0.5,-0.3,-0.15,0,0.2,0.4,0.6,0.8,1,2,6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  main = "Treatment AB tested nanostring (less)",
  annotation_col = AB
)
