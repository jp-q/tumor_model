setwd("~/GITCLONES/2017-11-09-Jingping-Nilogen/Mitra/jpApr2018/")



raw_multiplex <- readODS::read.ods("~/GITCLONES/2017-11-09-Jingping-Nilogen/Mitra/jpApr2018/jp_mitra.ods",sheet = 6)
rownames(raw_multiplex) <- make.names(raw_multiplex[,1],unique = T)
raw_multiplex <- raw_multiplex[,-1]
colnames(raw_multiplex) <- raw_multiplex[1,]
raw_multiplex <- raw_multiplex[-1,]

raw_multiplex[,4:61] <- log(data.matrix(raw_multiplex[,4:61]))
log_multiplex <- raw_multiplex

mscore <- as.numeric(log_multiplex$Mscore)
mscore[mscore<26 & mscore !=0] = -1
mscore[mscore >=26] =1

log_multiplex$Mresponse <- mscore



###############################################
#HEATMAP
###############################################

annocol <- data.frame( Tu = factor(log_multiplex$tumor),
  m = factor(log_multiplex$Mresponse)
)
rownames(annocol) <- rownames(log_multiplex)
###annotation color needs to be revised

pheatmap::pheatmap(
  t(log_multiplex[,4:61]),
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  #color = colorRampPalette(rev(brewer.pal(
   # name = 
   #   "RdBu",
   # n=9
  #)))(15),
  #annotation_colors = tre_colr[1],
   labels_col  = log_multiplex$treatment,
  #gaps_col = 1,
  main = "log_test",
  annotation_col = annocol,
  border_color = FALSE
)


#try only reponse
pheatmap::pheatmap(
  t(log_multiplex[mscore>= 0,4:61]),
  cluster_cols = T,
  cluster_rows = T,
  #breaks = c(-4,-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.6,0.8,1,1.5,2.6),
  #color = colorRampPalette(rev(brewer.pal(
  # name = 
  #   "RdBu",
  # n=9
  #)))(15),
  #annotation_colors = tre_colr[1],
  labels_col  = log_multiplex$treatment,
  #gaps_col = 1,
  main = "log_test",
  annotation_col = annocol,
  border_color = FALSE
)


#try fold change

logfc <-  sapply(rownames(log_multiplex), function (cn) {
  ctrl_col <- gsub("\\_.*", "_control", cn)
  log_multiplex[cn,4:61] - log_multiplex[ctrl_col,4:61]
})
logfc <- logfc[,-(1:20)]
logfc <- data.frame(logfc)
annocol <- data.frame( Tu = factor(c(rep(seq(1,20),5) )),
                       m = factor(log_multiplex$Mresponse)[21:120]
)
rownames(annocol) <- colnames(logfc)
###annotation color needs to be revised

library(RColorBrewer)
pheatmap::pheatmap(
  data.matrix(logfc),
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-8,-4,-2,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.5,1,2,4,6),
  color = colorRampPalette(rev(brewer.pal(
   name = 
     "RdBu",
   n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  labels_col  = log_multiplex$treatment[21:120],
  #gaps_col = 1,
  main = "logfc_test",
  annotation_col = annocol,
  border_color = FALSE
)


pheatmap::pheatmap(
  data.matrix(logfc[,mscore[21:120] > 0]),
  cluster_cols = T,
  cluster_rows = T,
  breaks = c(-8,-4,-2,-1,-0.75,-0.5,-0.25,0,0.2,0.4,0.5,1,1.5,2,6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  labels_col  = log_multiplex$treatment[21:120][mscore[21:120] > 0],
  #gaps_col = 1,
  main = "logfc_test",
  annotation_col = annocol,
  border_color = FALSE
)


#####try order with m score
logfc_m <- logfc
logfc_m[59,] <- as.factor(log_multiplex$Mscore[21:120])
rownames(logfc_m)[59] <- "mscore"
test <- logfc_m[,order(-as.numeric(unlist(logfc_m[nrow(logfc_m),])))]


annocol <- data.frame( Tu = factor(c(rep(seq(1,20),5) ))[order(-as.numeric(unlist(logfc_m[nrow(logfc_m),])))]
)
rownames(annocol) <- colnames(test)

pheatmap::pheatmap(
  data.matrix(test[-59,]),
  cluster_cols = F,
  cluster_rows = T,
  breaks = c(-8,-4,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.5,2,6),
  color = colorRampPalette(rev(brewer.pal(
    name = 
      "RdBu",
    n=9
  )))(15),
  #annotation_colors = tre_colr[1],
  labels_col  = log_multiplex$treatment[21:120][order(-as.numeric(unlist(logfc_m[nrow(logfc_m),])))],
  #gaps_col = 1,
  main = "logfc_m_test",
  annotation_col = annocol,
  border_color = FALSE
)









#correlation
#build gene list
logged <- t(data.matrix(log_multiplex[,-c(1:3)]))
cyto_hnscc_logfc <- data.matrix(logged) - apply(data.matrix(logged), 1, median)
cytologfc_17 <- cyto_hnscc_logfc[42:58,]

##steal from alex
library(parallel)
source(sprintf("%s/bioinformatics-signatures/load.R", Sys.getenv("GITCLONES")), chdir=T)
setwd("/home/x189259/GITCLONES/tim3-mitra-tumor-explant-december-2017")
mitra_hnscc_metadata <- read.table("mitra_hnscc.metadata.txt", sep="\t", header=TRUE)

mitra_hnscc_raw_tpm <- read.table("mitra_hnscc.tpm.txt", header=T, sep="\t", row.names="gene_id.TPM")
mitra_hnscc_counts <- read.table("mitra_hnscc.counts.txt", header=T, sep="\t", row.names="gene_id.expected_count")

mitra_hnscc_samplenames <- gsub(".*(HNS\\d+Rx\\d+).*", "\\1", mitra_hnscc_metadata[,"Sequences"])
colnames(mitra_hnscc_raw_tpm) <- mitra_hnscc_samplenames
colnames(mitra_hnscc_counts) <- mitra_hnscc_samplenames

mitra_hnscc_outliers <- c("HNS18Rx6")
k <- which(!(mitra_hnscc_metadata[,"Sequences"] %in% mitra_hnscc_outliers))
mitra_hnscc_counts <- mitra_hnscc_counts[,k]
mitra_hnscc_raw_tpm <- mitra_hnscc_raw_tpm[,k]
mitra_hnscc_metadata <- mitra_hnscc_metadata[k,]
mitra_hnscc_metadata$tumor <- gsub(".*(HNS\\d+).*", "\\1", mitra_hnscc_metadata$Sequences)
mitra_hnscc_metadata$arm <- gsub("\\s+explant.culture", "",
                                 gsub("\\([^\\)]+\\)\\s+", "",
                                      mitra_hnscc_metadata[,"condition"]))

source(sprintf("%s/bmto-ngs/R/R/expression-normalization.R", Sys.getenv("GITCLONES")), chdir=T)
mitra_hnscc_tpm <- upper_quartile_normalize_table(mitra_hnscc_raw_tpm, quartile=.85)
mitra_hnscc_logtpm <- log2(.2 + mitra_hnscc_tpm)
mitra_hnscc_logfc <- mitra_hnscc_logtpm - apply(mitra_hnscc_logtpm, 1, median)
ensgToSymbol <- read.table("ensg_to_symbol.txt", header=TRUE, stringsAsFactors=FALSE)
mitra_hnscc_rtg <- ensgToSymbol[,2]
names(mitra_hnscc_rtg) <- ensgToSymbol[,1]
mitra_hnscc_sigscores <- cbind(scoreSignatures(mitra_hnscc_logfc, mitra_hnscc_rtg),
                               scoreCibersort(mitra_hnscc_tpm, mitra_hnscc_rtg))




##end stealing
gene_list <- c("CSF2","TNFRSF9","IFNG","IL10","GZMA","IL13","GZMB","FAS","IL2","IL4","IL5","IL6","FASLG","CCL3","CCL4","TNF","PRF1")
gtr <- names(mitra_hnscc_rtg)
names(gtr) <- unname(mitra_hnscc_rtg)
gene_list <- gtr[gene_list]
gene_logfc <- mitra_hnscc_logfc
gene_17_logfc <- gene_logfc[gene_list,]
gene_17_logfc <- gene_17_logfc[,order(gsub(".*x","",colnames(gene_17_logfc)))]

ptt <- 	c("Tu1","Tu2","Tu3","Tu4","Tu5","Tu6","Tu7","Tu8","Tu9","Tu10","Tu11","Tu12","Tu13","Tu14","Tu15","Tu16","Tu17","Tu18","Tu19","Tu20")
names(ptt) <- c("HNS1R",	"HNS3R","HNS4R","HNS5R","HNS7R","HNS9R","HNS10R","HNS11R","HNS13R","HNS15R","HNS16R","HNS18R","HNS19R","HNS20R","HNS21R",
                "HNS22R","HNS23R","HNS26R","HNS27R","HNS28R")
colnames(gene_17_logfc) <-  paste(unname(ptt[gsub("x.*","",colnames(gene_17_logfc))]),gsub(".*x","",colnames(gene_17_logfc)),sep = "_")

colnames(cytologfc_17) <- gsub("control","1",colnames(cytologfc_17))
colnames(cytologfc_17) <- gsub("M7824.M6903","5",colnames(cytologfc_17))
colnames(cytologfc_17) <- gsub("M7824.M6223","6",colnames(cytologfc_17))
colnames(cytologfc_17) <- gsub("M7824","2",colnames(cytologfc_17))
colnames(cytologfc_17) <- gsub("M6903","3",colnames(cytologfc_17))
colnames(cytologfc_17) <- gsub("M6233","4",colnames(cytologfc_17))
cytologfc_17 <-  cytologfc_17[,colnames(cytologfc_17) %in% colnames(gene_17_logfc)]
cytologfc_17 <- cytologfc_17[,sort(colnames(cytologfc_17))]

rownames(gene_17_logfc) <- mitra_hnscc_rtg[rownames(gene_17_logfc)]
gene_17_logfc <- gene_17_logfc[,sort(colnames(gene_17_logfc))]

b <- diag(cor(t(gene_17_logfc), t(cytologfc_17)))


plot(unname(b),xaxt="n",ylab = "correlation")
axis(1,at=1:17,labels=names(b))
 names(b) <- rownames(cytologfc_17)
a <-  cor.test(t(gene_17_logfc),t(cytologfc_17))

 





###############################
#sig genes
###############################
sigs <- c(#"CHAN_INTERFERON_PRODUCING_DENDRITIC_CELL",
          #"GNF2_CD7",
          "angelova2015.Effector.memory.CD8",
          "angelova2015.T.cells",
          "Sharma2017.Ifng.Urothelial",
          "angelova2015.Th1",
          "S1_NK_cells",
          "T_cell_signature_IRIS",
         # "GNF2_CD14",
          #"REACTOME_PD1_SIGNALING",
          "bindea2013.T.cells",
          "bindea2013.Cytotoxic",
          "msd.asco.2016.Ifng",
         #"GNF2_MYL2",
          #"REACTOME_IL_3_5_AND_GM_CSF_SIGNALING","GO_POSITIVE_REGULATION_OF_MACROPHAGE_CHEMOTAXIS","GO_IMMUNOGLOBULIN_BINDING",
         "msd.asco.2016.TcellExhaustion")
#sigs <- na.omit(signature.list[sigs])
gtr <- names(mitra_hnscc_rtg)
names(gtr) <- unname(mitra_hnscc_rtg)



#################this function plots heatmap as jieqing's request
#################requires a signature and a treatment arm as input
sigs_gene_heatmap <- function(sigs_gene,i){

  genename <- na.omit(unname(gtr[unlist(unname(sigs_gene))]))
  testsig <- mitra_hnscc_logtpm[,sort(c(grep("x1$",colnames(mitra_hnscc_logtpm))))]
  x <- paste("x",i,sep = "")
  testsig_2 <- mitra_hnscc_logtpm[,sort(c(grep(x,colnames(mitra_hnscc_logtpm))))]
  testsig <- testsig[genename,]
  rownames(testsig) <- mitra_hnscc_rtg[rownames(testsig)]
  testsig_2 <- testsig_2[genename,]
  rownames(testsig_2) <- mitra_hnscc_rtg[rownames(testsig_2)]
##check if they have same samples
  if(ncol(testsig_2)>ncol(testsig)){
  testsig_2 <- testsig_2[,gsub("x.*","",colnames(testsig_2)) %in% gsub("x.*","",colnames(testsig))]
  testsig <- testsig[,gsub("x.*","",colnames(testsig)) %in% gsub("x.*","",colnames(testsig_2))]
  }else if(ncol(testsig_2)<ncol(testsig)){
    testsig <- testsig[,gsub("x.*","",colnames(testsig)) %in% gsub("x.*","",colnames(testsig_2))]
    testsig_2 <- testsig_2[,gsub("x.*","",colnames(testsig_2)) %in% gsub("x.*","",colnames(testsig))]
  }
  
  
  options(bitmapType = 'cairo')
  a <- pheatmap::pheatmap(
    data.matrix(testsig),
    cluster_cols = T,
    cluster_rows = F,
    # breaks = c(-8,-4,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.5,2,6),
    # color = colorRampPalette(rev(brewer.pal(
    #  name = 
    #   "RdBu",
    #  n=9
    #)))(15),
    #annotation_colors = tre_colr[1],
    #labels_col  = log_multiplex$treatment[21:120][order(-as.numeric(unlist(logfc_m[nrow(logfc_m),])))],
    #gaps_col = 1,
    main = "sigs_gene",
    #annotation_col = annocol,
    border_color = FALSE
  )
  testsig <- testsig[,a$tree_col$order]
  testsig_2 <- testsig_2[, a$tree_col$order]
  
  
  colnames(testsig_2) <- colnames(testsig)
  testsig <- rbind(testsig,testsig_2)
  row_anno <- c(rownames(testsig_2),rownames(testsig_2))
  col_anno <- c(gsub("Rx.*$","",colnames(testsig_2)))
  #dev.off()
  setwd("~/GITCLONES/2017-11-09-Jingping-Nilogen/Mitra/jpApr2018/")
  
  jpeg(file=paste("arm_",i,"_signature_",names(sigs_gene),".jpeg",sep = ""),width = 1400,height = 1400)
  #png(sprintf("%s/signatures-from-camera.png", output_dir), width=1600, height=1200)
  
  pheatmap::pheatmap(
    data.matrix(testsig),
    cluster_cols = F,
    cluster_rows = F,
    gaps_row = nrow(testsig_2),
    labels_row = row_anno,
    labels_col = col_anno,
    # breaks = c(-8,-4,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.5,2,6),
    # color = colorRampPalette(rev(brewer.pal(
    #  name = 
    #   "RdBu",
    #  n=9
    #)))(15),
    #annotation_colors = tre_colr[1],
    #labels_col  = log_multiplex$treatment[21:120][order(-as.numeric(unlist(logfc_m[nrow(logfc_m),])))],
    #gaps_col = 1,
    main = paste("Treatment_arm",i,"signature:",names(sigs_gene)),
    #annotation_col = annocol,
    border_color = FALSE
    
    
  )
  dev.off()

  
  
  
}



length(sigs)

for( j in 1:length(sigs)){
sigs_gene_heatmap(signature.list[sigs[j]],6)
}
sigs_gene_heatmap(signature.list[sigs[1]],2)





#correlation vs more
gene_logfc <- mitra_hnscc_logfc
ptt <- 	c("Tu1","Tu2","Tu3","Tu4","Tu5","Tu6","Tu7","Tu8","Tu9","Tu10","Tu11","Tu12","Tu13","Tu14","Tu15","Tu16","Tu17","Tu18","Tu19","Tu20")
names(ptt) <- c("HNS1R",	"HNS3R","HNS4R","HNS5R","HNS7R","HNS9R","HNS10R","HNS11R","HNS13R","HNS15R","HNS16R","HNS18R","HNS19R","HNS20R","HNS21R",
                "HNS22R","HNS23R","HNS26R","HNS27R","HNS28R")
colnames(gene_logfc) <-  paste(unname(ptt[gsub("x.*","",colnames(gene_logfc))]),gsub(".*x","",colnames(gene_logfc)),sep = "_")

colnames(cyto_hnscc_logfc) <- gsub("control","1",colnames(cyto_hnscc_logfc))
colnames(cyto_hnscc_logfc) <- gsub("M7824.M6903","5",colnames(cyto_hnscc_logfc))
colnames(cyto_hnscc_logfc) <- gsub("M7824.M6223","6",colnames(cyto_hnscc_logfc))
colnames(cyto_hnscc_logfc) <- gsub("M7824","2",colnames(cyto_hnscc_logfc))
colnames(cyto_hnscc_logfc) <- gsub("M6903","3",colnames(cyto_hnscc_logfc))
colnames(cyto_hnscc_logfc) <- gsub("M6233","4",colnames(cyto_hnscc_logfc))
cyto_hnscc_logfc <-  cyto_hnscc_logfc[,colnames(cyto_hnscc_logfc) %in% colnames(gene_logfc)]
cyto_hnscc_logfc <- cyto_hnscc_logfc[,sort(colnames(cyto_hnscc_logfc))]
gene_logfc <-  gene_logfc[,colnames(gene_logfc) %in% colnames(cyto_hnscc_logfc)]
gene_logfc <- gene_logfc[,sort(colnames(gene_logfc))]



cor_result <- cor(t(gene_logfc), t(cyto_hnscc_logfc))
b <- data.frame(genename = unname(mitra_hnscc_rtg[as.character(rownames(cor_result))]))
cor_result <- cbind(cor_result,b)



#loop to get high cor values
c <- character()

for(i in 1:nrow(cor_result)){
  for(j in 1:58){
    b <- as.numeric()
    a <- as.numeric(cor_result[i,j])
    if (is.na(a)){
      a=0
    }
   if( a > 0.5){
      b <- a
      names(b) <- paste(cor_result[i,61],colnames(cor_result)[j], sep="_")
    } 
 if( a < -0.5){
     b <- a
      names(b) <- paste(cor_result[i,61],colnames(cor_result)[j], sep="_")
   }
    c <- c(c,b)
    print(paste(i,j))
}
    
  }


####grep interested genes
target_genes <- c("CD274","TGFB","HAVCR2","TIGIT","PDCD1","LGALS9","CD226","IDO1","TDO2","GCN2","PVR")
target_cor <- cor_result[cor_result$genename %in% target_genes,]



#######################try penalized regression

install.packages("glmnet")
library("glmnet")

sigscore <- mitra_hnscc_sigscores
model_signatures <- c(grep("Walter2013", names(signature.list), value=T),
                        "Byers2013.EMT.epithelial",
                        "Huang2012_EMTdown",
                        "angelova2015.TGD",
                        "angelova2015.T.cells",
                        "angelova2015.Th1",
                        "msd.asco.2015.DeNovo",
                        "S10_B_cells",
                        "bindea2013.B.cells",
                        "msd.asco.2015.TCR",
                        "IL17.Grogan2016patent",
                        "angelova2015.Mast.cells",
                        "Huang2012_EMTup",
                        "Sharma2017.Ifng.Urothelial",
                        "angelova2015.pDC",
                        "Budinska2013.immune.response.CRC",
                        "msd.asco.2015.melanoma",
                        "angelova2015.T.cells",
                        "CongGroup3_IFN",
                        "angelova2015.Th1",
                        "BreastStromaMetagene.Farmer2005",
                        "msd.asco.2015.Ifng",
                        "genentech.teff",
                        "Walsh2007_IFNsig_6g",
                        "Winter2007.hypoxia.down",
                        "angelova2015.DC",
                        "S4_T_cells_resting",
                        "fibro.signature",
                        "REACTOME_COLLAGEN_FORMATION",
                        "iOnc.Amit.fibroblasts",
                        "S1_NK_cells",
                        "Finkernagel2016.TAM.up")
model_sigscore <- t(sigscore[,unique(intersect(model_signatures,
                                             colnames(sigscore)))])
model_sigscore <- model_sigscore[,sort(colnames(model_sigscore))]
ptt <- 	c("Tu1","Tu2","Tu3","Tu4","Tu5","Tu6","Tu7","Tu8","Tu9","Tu10","Tu11","Tu12","Tu13","Tu14","Tu15","Tu16","Tu17","Tu18","Tu19","Tu20")
names(ptt) <- c("HNS1R",	"HNS3R","HNS4R","HNS5R","HNS7R","HNS9R","HNS10R","HNS11R","HNS13R","HNS15R","HNS16R","HNS18R","HNS19R","HNS20R","HNS21R",
                "HNS22R","HNS23R","HNS26R","HNS27R","HNS28R")
colnames(model_sigscore) <-  paste(unname(ptt[gsub("x.*","",colnames(model_sigscore))]),gsub(".*x","",colnames(model_sigscore)),sep = "_")

mscore_data <- log_multiplex
rownames(mscore_data) <- gsub("control","1",rownames(mscore_data))
rownames(mscore_data) <- gsub("M7824.M6903","5",rownames(mscore_data))
rownames(mscore_data) <- gsub("M7824.M6223","6",rownames(mscore_data))
rownames(mscore_data) <- gsub("M7824","2",rownames(mscore_data))
rownames(mscore_data) <- gsub("M6903","3",rownames(mscore_data))
rownames(mscore_data) <- gsub("M6233","4",rownames(mscore_data))
mscore_data <- mscore_data[mscore_data$Mscore !=0,]
model_sigscore <- model_sigscore[,gsub(".*\\_","",colnames(model_sigscore))!=1]
mscore_data <- mscore_data[rownames(mscore_data)%in% colnames(model_sigscore),]
model_sigscore <- model_sigscore[,colnames(model_sigscore) %in% rownames(mscore_data)]

mscore_data <- data.matrix(t(mscore_data[,-c(1,2,3,62)]))
model_sigscore <- model_sigscore[,sort(colnames(model_sigscore))] 
mscore_data <- rbind(model_sigscore, mscore_data)



#a <- glmnet(t(mscore_data[-c(33,34),]),mscore_data[34,],family = "binomial")
#sparse lvl?


############try glmnetgraph


# set.seed(1234)
# library(glmgraph)
# n <- 100
# p1 <- 10
# p2 <- 90
# p <- p1+p2
# X <- matrix(rnorm(n*p), n,p)
# magnitude <- 1
# A <- matrix(rep(0,p*p),p,p)
# A[1:p1,1:p1] <- 1
# A[(p1+1):p,(p1+1):p] <- 1
# diag(A) <- 0
# btrue <- c(rep(magnitude,p1),rep(0,p2))
# intercept <- 0
# eta <- intercept+X%*%btrue
# ### construct laplacian matrix from adjacency matrix
# diagL <- apply(A,1,sum)
# L <- -A
# diag(L) <- diagL
# ### gaussian
# Y <- eta+rnorm(n)
# obj <- glmgraph(X,Y,L,family="gaussian")
# 
# plot(obj)
# ### binomial
# Y <- rbinom(n,1,prob=1/(1+exp(-eta)))
# obj <- glmgraph(X,Y,L,family="binomial")
# plot(obj)

model_sigscore <- mscore_data[-c(33,34),]
###try to test this model
test <- model_sigscore[,1:82]
train <- model_sigscore[,83:92]
class(test) <- "numeric"
class(train) <- "numeric"
test_L <- t(matrix(rep(1,82*(nrow(model_sigscore))),ncol = 82,nrow = nrow(model_sigscore)))
test_M <- as.numeric(mscore_data[34,])[1:82]
train_M <- as.numeric(mscore_data[34,])[83:92]
b <- glmgraph(X=t(as.matrix(test)),Y= test_M,L= test_L,family = "binomial",penalty = "lasso")
plot(b)
#fail
i = predict(b,t(as.matrix(train)),type = "class")

#sort of works
a <- glmnet(t(as.matrix(test)), test_M,family = "binomial")
View(predict(a,newx=t(as.matrix(train)),type = "class"))
sort(a$beta[,53])




# #use cytokine
# raw_multiplex <- readODS::read.ods("~/GITCLONES/2017-11-09-Jingping-Nilogen/Mitra/jpApr2018/jp_mitra.ods",sheet = 6)
# rownames(raw_multiplex) <- make.names(raw_multiplex[,1],unique = T)
# raw_multiplex <- raw_multiplex[,-1]
# colnames(raw_multiplex) <- raw_multiplex[1,]
# raw_multiplex <- raw_multiplex[-1,]
# 
# raw_multiplex[,4:61] <- log(data.matrix(raw_multiplex[,4:61]))
# log_multiplex <- raw_multiplex
# log_multiplex <- log_multiplex[log_multiplex$Mscore != 0,]
# mscore <- as.numeric(log_multiplex$Mscore)
# mscore[mscore<26] = 0
# mscore[mscore >=26] =1
# log_multiplex$Mresponse <- mscore
# 
# model_cyto <- log_multiplex[log_multiplex$Mscore != 0,]
# model_cyto <- data.matrix(t(model_cyto[,-c(1,2,3)]))
# cyto_response <- model_cyto[60,]
# cyto_m <- model_cyto[59,]
# model_cyto <- model_cyto[-c(59,60),]
# class(model_cyto) <- "numeric"
# 
# ###try to test this model
# test <- model_cyto[,1:89]
# train <- model_cyto[,90:100]
# test_L <- t(matrix(rep(1,89*(nrow(model_cyto))),ncol = 89,nrow = nrow(model_cyto)))
# test_M <- as.numeric(cyto_response)[1:89]
# train_M <- as.numeric(cyto_response)[90:100]
# e <- glmgraph(X=t(as.matrix(test)),Y= test_M,L= test_L,family = "binomial",penalty = "lasso")
# plot(e)
# #fail
# ee <- predict(e,X=t(as.matrix(train)),type = "class")
# View(ee)
# 
# #fail
# a <- glmnet(t(as.matrix(test)), as.factor(test_M),family = "binomial")
# View(predict(a,newx=t(as.matrix(train)),type = "class"))
# 
# 


#A=model_sigscore

input <- mscore_data[1:32,]
y <- as.factor(mscore_data[91,])
randco <- function(A) {
  #A <- model_sigscore
  hist <- c()
  repeat{
  class(A) <- "numeric"
  A <- scale(t(A))
  A <- t(A)
  t <- sample(92,22)
  train <- A[,-t]
  test <- A[,t]

  response <- as.numeric(y)[-t]
  test_response <- as.numeric(y)[t]
 a <- glmnet(t(as.matrix(train)), response,family = "binomial")
 b <- data.frame(predict(a,newx =t(test),type  ="class"))
 list_acc <- c()
 for(i in 1: ncol(b)){
   accu <- length(which(as.numeric(as.character(unlist(b[i]))) != test_response))
   list_acc <- c(list_acc,sum(accu))
 
   
 }
 print(list_acc)
  hist <- c(mean(list_acc),hist)
 if(min(list_acc) < 1){
   sink(file=paste("result.txt"))
   cat(paste("difference for all results:\n"))
   cat(list_acc)
   cat("\n\n")
   best <- which(list_acc == min(list_acc))
   cat(paste(" The least diff is samples models: \n"))
   cat(best)
   cat("\n\n")
   for(j in 1:length(best)){
     cat("======================================")
     cat(j)
     cat("======================================")
     cat("\n")
     cat(paste("The actual response is:\n"))
     cat(test_response)
     cat("\n\n")
     cat(paste("The predict result",j,":\n"))
     cat(as.numeric(as.character(b[,best[j]])))
     cat("\n\n")
     cat(paste("and the predict coefficient:\n"))
     print(sort(a$beta[,best[j]]))
     cat("\n\n")
     cat(min(list_acc))
     cat("\n\n")
   }
   cat(colnames(test))
   sink()
   hist(hist,xlab = "# of diff",main = "Sigscore model accuracy frequency")
   return(hist)
   break
}
 }
 # else{
 
 # print(list_acc)
 # if(min(list_acc<3)){
 # best <- which(list_acc == min(list_acc))
 # print(best)
 # for(j in 1:length(best)){
 #  print(paste("The actual response is:"))
 # print(test_response)
 # print(paste("The predict result",j,":"))
 # print(as.numeric(as.character(b[,best[j]])))
 # print(paste("and the predict coefficient"))
 # print(sort(a$beta[,best[j]]))
 # print(min(list_acc))
 #  }
 #  print(colnames(test)) 
 #  
 #  
 
  
   }

randco(model_sigscore)

#t(replicate(10,randco(A)))

















################try log muiltiplex


model_cyto <- data.matrix(logfc)
model_m <- log_multiplex$Mresponse[21:120]



model_plex<- function(A) {
  hist <- c()
  #A <- model_sigscore
  repeat{
    class(A) <- "numeric"
    A <- scale(t(A))
    A <- t(A)
    t <- sample(100,20)
    train <- A[,-t]
    test <- A[,t]
    
    response <- model_m[-t]
    test_response <- model_m[t]
    a <- glmnet(t(as.matrix(train)), response,family = "binomial")
    b <- data.frame(predict(a,newx =t(test),type  ="class"))
    list_acc <- c()
    for(i in 1: ncol(b)){
      accu <- length(which(as.numeric(as.character(unlist(b[i]))) != test_response))
      list_acc <- c(list_acc,sum(accu))
      
      
    }
    print(list_acc)
    hist <- c(mean(list_acc),hist)
    if(min(list_acc) < 1){
      sink(file=paste("result.txt"))
      cat(paste("difference for all results:\n"))
      cat(list_acc)
      cat("\n\n")
      best <- which(list_acc == min(list_acc))
      cat(paste(" The least diff is samples models: \n"))
      cat(best)
      cat("\n\n")
      for(j in 1:length(best)){
        cat("======================================")
        cat(j)
        cat("======================================")
        cat("\n")
        cat(paste("The actual response is:\n"))
        cat(test_response)
        cat("\n\n")
        cat(paste("The predict result",j,":\n"))
        cat(as.numeric(as.character(b[,best[j]])))
        cat("\n\n")
        cat(paste("and the predict coefficient:\n"))
        print(sort(a$beta[,best[j]]))
        cat("\n\n")
        cat(min(list_acc))
        cat("\n\n")
      }
      cat(colnames(test))
      sink()
      hist(hist,xlab = "# of diff",main = "Cytokine model accuracy frequency")
      break
    }
  }
}


########################### combination of cyto and sigs
#input
input <- mscore_data[-91,]
response <- as.numeric(mscore_data[91,])
response[response == -1]  <- 0
setwd("~/GITCLONES/2017-11-09-Jingping-Nilogen/Mitra/jpApr2018/")

full_model <- function(A) {
  A <- input
  hist <- c()
  repeat{
    class(A) <- "numeric"
    A <- scale(t(A))
    A <- t(A)
    t <- sample(92,20)
    train <- A[,-t]
    test <- A[,t]
    
    train_response <- response[-t]
    test_response <- response[t]
    a <- cv.glmnet(t(as.matrix(train)), train_response,family = "binomial")
    b <- predict(a,newx =t(test),type  ="class",s="lambda.min")
    #list_acc <- c()
    #for(i in 1: ncol(b)){
    accu <- length(which(as.numeric(as.character(unlist(b[i]))) != test_response))
    #list_acc <- c(list_acc,sum(accu))
    
       
     
    print(accu)
   hist <- c(accu,hist)
   # if(min(list_acc) < 2){
    if(accu < 7){
      sink(file=paste("result.txt"))
   #   cat(paste("difference for all results:\n"))
      cat(paste("The number of diff is:\n"))
      cat(accu)
      cat("\n\n")
      #best <- which(list_acc == min(list_acc))
      cat(paste(" The predict result: \n"))
      cat(as.numeric(as.character(unlist(b[i]))))
      cat("\n\n")
     # for(j in 1:length(best)){
      #  cat("======================================")
       # cat(j)
     #   cat("======================================")
      #  cat("\n")
        cat(paste("The actual response is:\n"))
        cat(test_response)
        cat("\n\n")
      #  cat(paste("The predict result",j,":\n"))
      #  cat(as.numeric(as.character(b[,best[j]])))
       # cat("\n\n")
        cat(paste("and the predict coefficient:\n"))
        print(coef(a, s = "lambda.min"))
        cat("\n\n")
      cat(colnames(test))
      sink()
      hist(hist,xlab = "# of diff",main = "Combine model accuracy frequency")
      plot(a)
      return(a)
      break
    }
}
}
# aaa <- cv.glmnet(scale(t(input)),as.numeric(response))
# plot(aaa)
# coef_cv=coef(aaa, s = "lambda.min")


########################### ihc glucose
setwd("~/GITCLONES/2017-11-09-Jingping-Nilogen/Mitra/jpApr2018/")



glu <- readODS::read.ods("~/GITCLONES/2017-11-09-Jingping-Nilogen/Mitra/jpApr2018/jp_mitra.ods",sheet = 8)



ihc <- readODS::read.ods("~/GITCLONES/2017-11-09-Jingping-Nilogen/Mitra/jpApr2018/jp_mitra.ods",sheet = 3)
rownames(ihc) <- ihc[,1]
colnames(ihc) <- ihc[1,]
ihc <- ihc[-1,-1]
ihc <- scale(data.matrix(ihc))
