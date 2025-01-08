# Spatial DE_ meta

# load Libraries --------------------------
library(EnhancedVolcano)
library(gp.sa.diff)
library(SingleCellExperiment)
library(scater)
library(scran.chan)
library(Seurat)
library(dplyr)
library(tidyverse)
library(metap)
library(metafor)
library(scater)
sc <- import("scanpy")

H5toSeurat <- function(adata=adata){
  counts <- t(as.matrix(adata$layers["counts"]))
  colnames(counts) <- adata$obs_names$to_list()
  rownames(counts) <- adata$var_names$to_list()
  counts <- Matrix::Matrix(as.matrix(counts), sparse = T)
  
  data <- t(as.matrix(adata$layers["logcounts"]))
  colnames(data) <- adata$obs_names$to_list()
  rownames(data) <- adata$var_names$to_list()
  data <- Matrix::Matrix(as.matrix(data), sparse = T)
  
  print( "Creating Seurat")
  seurat <- CreateSeuratObject(counts,assay = "Spatial")
  seurat <- SetAssayData(seurat, layer = "data", data,assay = "Spatial")
  seurat <- AddMetaData(seurat, data.frame(adata$obs))
  
  print( "Adding Embeddings")
  embedding <- adata$obsm["X_umap"]
  rownames(embedding) <- adata$obs_names$to_list()
  colnames(embedding) <- c("umap_1", "umap_2")
  seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_",assay = "Spatial")
  
  embedding <- adata$obsm["X_pca"]
  rownames(embedding) <- adata$obs_names$to_list()
  colnames(embedding) <- paste0("PC_",1:dim(embedding)[2])
  seurat[["pca"]] <- CreateDimReducObject(embedding, key = "pca_",assay = "Spatial")
  
  features <- data.frame(gene_ids=adata$var$gene_ids,
                         symbols=adata$var_names$to_list(),
                         row.names = adata$var_names$to_list())
  seurat@misc <- data.frame(adata$var)
  return(seurat)
  
}
# Set paths --------
RESULTS_PATH = "/Users/lucast3/Documents/VisiumAnalysis/results"
data.dir <- "/Users/lucast3/Documents/VisiumAnalysis/data"

# Load Data -----------------------------------------------------
sobj_merge <- readRDS(file.path(data.dir,"MergedSpatialObject.rds"))

# Puesdobulk DE AD vs Control Analysis =================================
library(SingleCellExperiment)
library(wesanderson)
library(scater)


counts <- GetAssayData(object = sobj_merge,
                       layer = "data",
                       assay="SCT") # get counts table
se <- SingleCellExperiment(assays = list(counts = counts),
                           colData = sobj_merge@meta.data)

set.seed(824)
sum_by <- c("Abeta","BrainRegion","orig.ident")
summed <- aggregateAcrossCells(se, id=colData(se)[,sum_by])

condition = c('TauPS2APP','TauPS2APP','TauPS2APP','NonTG','TauPS2APP','NonTG')

# DE for all Pathology vs no ===========
AllPath_res <- list()
for (i in 1:length(condition)) {
  study <- datasets[i]
  out <- runVoom(summed, groups="Abeta", comparisons=list(c("Abeta", "NoPathology")),
                 subset.factor="orig.ident",subset.levels=c(study), commit="never")
  
  out_final <- as.data.frame(out[[1]]) %>% mutate(symbol = rowData(summed)$symbol)
  
  AllPath_res[[study]] <- out_final
}

# DE for all TauPS2APP Pathology vs no ===========
subset_se <- se[,se$Genotype=="TauPS2APP"]
sum_by <- c("Abeta","BrainRegion","orig.ident")
summed <- aggregateAcrossCells(subset_se, id=colData(subset_se)[,sum_by])


TauPS2APP_res <- list()
for (i in levels(factor(summed$orig.ident))) {
  study <- i
  out <- runVoom(summed, groups="Abeta", comparisons=list(c("Abeta", "NoPathology")),
                 subset.factor="orig.ident",subset.levels=c(study), commit="never")
  
  out_final <- as.data.frame(out[[1]]) %>% mutate(symbol = rowData(summed)$symbol)
  
  TauPS2APP_res[[study]] <- out_final
}

# DE for all Cortex: TauPS2APP Pathology vs no ===========
subset_cortex <- subset_se[,subset_se$BrainRegion=="Cortex"]
sum_by <- c("Abeta","seurat_clusters","orig.ident")
summed <- aggregateAcrossCells(subset_cortex, id=colData(subset_cortex)[,sum_by])
Cortex_res <- list()
for (i in levels(factor(summed$orig.ident))) {
  study <- i
  out <- runVoom(summed, groups="Abeta", comparisons=list(c("Abeta", "NoPathology")),
                 subset.factor="orig.ident",subset.levels=c(study), commit="never")
  
  out_final <- as.data.frame(out[[1]]) %>% mutate(symbol = rowData(summed)$symbol)
  
  Cortex_res[[study]] <- out_final
}

# DE for all Hippocampus: TauPS2APP Pathology vs no ===========
subset_hipp <- subset_se[,subset_se$BrainRegion=="Hipp"]
sum_by <- c("Abeta","seurat_clusters","orig.ident")
summed <- aggregateAcrossCells(subset_hipp, id=colData(subset_hipp)[,sum_by])
Hipp_res <- list()
for (i in levels(factor(summed$orig.ident))) {
  study <- i
  out <- runVoom(summed, groups="Abeta", comparisons=list(c("Abeta", "NoPathology")),
                 subset.factor="orig.ident",subset.levels=c(study), commit="never")
  
  out_final <- as.data.frame(out[[1]]) %>% mutate(symbol = rowData(summed)$symbol)
  
  Hipp_res[[study]] <- out_final
}

# P-value meta----------
res_list <- Hipp_res
feat <- as.data.frame(genomitory::getFeatures("GMTY17:GRCm38/GRCm38.IGIS4.0.genes.rds@REVISION-3"))

fc <- matrix(0, ncol=length(res_list), nrow=dim(feat)[1])
rownames(fc) <- feat$symbol

pval <- matrix(NA, ncol=length(res_list), nrow=dim(feat)[1])
rownames(pval) <- feat$symbol

mexp <- matrix(NA, ncol=length(res_list), nrow=dim(feat)[1])
rownames(mexp) <- feat$symbol

sefc <- matrix(NA, ncol=length(res_list), nrow=dim(feat)[1])
rownames(sefc) <- feat$symbol

for(i in 1:length(res_list)){
  m <- match(rownames(fc), rownames(res_list[[i]]))
  f.a = !is.na(m)
  f.b = m[f.a]
  fc[f.a,i] <- res_list[[i]][f.b,"LogFC"]
  pval[f.a,i] <- res_list[[i]][f.b,"PValue"]
  mexp[f.a,i] <- res_list[[i]][f.b,"AveExpr"]
  sefc[f.a,i] <- res_list[[i]][f.b,"LogFC"]/res_list[[i]][f.b,"t"]
}

x<-strsplit(names(res_list),".", fixed=T)
x2 <- sapply(x, function(x) purrr::pluck(x,1))
colnames(pval) = x2
colnames(fc) = x2
colnames(mexp) = x2
colnames(sefc) = x2


deg_master <- list(pval=pval, fc=fc, se_fc=sefc, mean_exp=mexp)

# Set up for FC meta ---------------------------------
f <- colnames(deg_master$pval)
pval <- deg_master$pval[, f]
fc <- deg_master$fc[, f]
mexp <- deg_master$mean_exp[, f]
sefc <- cbind(fc,deg_master$se_fc[,f])

pval_up <- pval * 0
pval_down <- pval * 0

for (i in 1:dim(fc)[2]) {
  f1 <- which(fc[, i] > 0)
  f2 <- which(fc[, i] < 0)
  pval_up[f1, i] <- pval[f1, i] * 0.5
  pval_up[f2, i] <- 1 - (pval[f2, i] * 0.5)
  pval_down[f1, i] <- 1 - (pval[f1, i] * 0.5)
  pval_down[f2, i] <- pval[f2, i] * 0.5
  
}


# Syntax
x <- pval_down[rowSums(is.na(pval_down)) != ncol(pval_down), ]

#combine pvalues by the sum of logs for each gene (Fishers meta pvalue)
meta_up <- suppressWarnings(apply(pval_up, 1, function(y) sumlog(y[!is.na(y)])$p))
meta_down <- suppressWarnings(apply(pval_down, 1, function(y) sumlog(y[!is.na(y)])$p))
#Average fc an expression
mean_fc <- rowMeans(fc, na.rm = T)
mean_exp <- rowMeans(mexp, na.rm = T)

meta_down <- meta_down[-which(duplicated(names(meta_down))==TRUE)]
meta_up <- meta_up[-which(duplicated(names(meta_up))==TRUE)]
mean_fc <- mean_fc[-which(duplicated(names(mean_fc))==TRUE)]
mean_exp <- mean_exp[-which(duplicated(names(mean_exp))==TRUE)]





f_se_fc <- apply(sefc[,1:length(res_list)],1,function(x) sum(!is.na(x))>=0.1)
f_se_se <- apply(sefc[,c(length(res_list)+1):(length(res_list)*2)],1,function(x) sum(!is.na(x))>=0.1)
f_se <- f_se_fc & f_se_se

# Calulate meta using different approaces -------------------------------------
dsl_res <- apply(sefc[f_se,], 1, function(row) rma(yi=as.vector(row[1:length(res_list)]),vi=as.vector(row[(length(res_list)+1):(length(res_list)*2)]),method="DL"))
hsk_res <- apply(sefc[f_se,], 1, function(row) rma(yi=as.vector(row[1:length(res_list)]),vi=as.vector(row[(length(res_list)+1):(length(res_list)*2)]),method="HSk"))
reml_res <- apply(sefc[f_se,], 1, function(row) rma(yi=as.vector(row[1:length(res_list)]),vi=as.vector(row[(length(res_list)+1):(length(res_list)*2)]),method="REML"))
sj_res <- apply(sefc[f_se,], 1, function(row) rma(yi=as.vector(row[1:length(res_list)]),vi=as.vector(row[(length(res_list)+1):(length(res_list)*2)]),method="SJ"))
hksj_res <- apply(sefc[f_se,], 1, function(row) rma(yi=as.vector(row[1:length(res_list)]),vi=as.vector(row[(length(res_list)+1):(length(res_list)*2)]),method="SJ",test = "knha"))

meta_list <- list(DL=dsl_res,HSk=hsk_res,
                 REML=reml_res,SJ=sj_res,
                 HKSJ=hksj_res)
generateDataFrame <- function(res,sefc,f_se,fc,pval_up,meta_up,meta_down){
  betas <- sapply(res,function(x) purrr::pluck(x,"beta"))
  ses <- sapply(res,function(x) purrr::pluck(x,"se"))
  dl_ps <- sapply(res,function(x) purrr::pluck(x,"pval"))
  dl_ps <- sapply(res,function(x) purrr::pluck(x,"pval"))
  
  dl_sub <- data.frame(ID=rownames(sefc)[f_se],beta=betas,se=ses,pval=dl_ps)
  
  dl_all <- data.frame(ID=rownames(sefc))
  dl_all <- left_join(dl_all,dl_sub,by=c("ID"="ID"))
  
  dl_all <- dl_all[-which(duplicated(dl_all$ID)==TRUE),]
  n_tested = apply(pval_up, 1, function(y)
    sum(!is.na(y)))
  n_tested <- n_tested[-which(duplicated(names(n_tested))==TRUE)]
  n_up = apply(fc, 1, function(y)
    sum(y > 0, na.rm = T))
  n_up <- n_up[-which(duplicated(names(n_up))==TRUE)]
  n_down = apply(fc, 1, function(y)
    sum(y < 0, na.rm = T))
  n_down <- n_down[-which(duplicated(names(n_down))==TRUE)]
  
  temp <-
    data.frame(
      ID = dl_all$ID,
      AveExpr = as.numeric(mean_exp),
      PValue = dl_all$pval,
      FDR = p.adjust(dl_all$pval,method="BH"),
      LogFC = as.numeric(mean_fc),
      dl_mu = dl_all$beta,
      dl_s = dl_all$se,
      metap_up = meta_up,
      metap_down = meta_down,
      adj_metap_up = p.adjust(meta_up, method = "BH"),
      adj_metap_down = p.adjust(meta_down, method = "BH"),
      n_tested = n_tested,
      n_up = n_up,
      n_down = n_down
    )
  
  x <- match(temp$ID,feat$symbol)
  features<- feat[x,]
  temp <- left_join(temp, features,by=c("ID"="symbol"))
  temp<- temp[which(is.nan(temp$AveExpr)==FALSE),]
  ord <- order(temp$metap_up)
  
  temp <- temp[ord, ]
  x <-
    c(
      "ID.y",
      "ID",
      "AveExpr",
      "PValue",
      "FDR",
      "LogFC",
      "dl_mu",
      "dl_s",
      "metap_up",
      "metap_down",
      "adj_metap_up",
      "adj_metap_down",
      "n_tested",
      "n_up",
      "n_down","desc"
    )
  
  res_dsl <- temp %>% dplyr::select(x)
  colnames(res_dsl)[1:2] <- c("ID","symbol")
  return(res_dsl)
}
generatPlotData <- function(res_dsl){
  plotdata_up <- res_dsl %>% filter(LogFC > 0)
  res_dsl_plot <- data.frame(ID=plotdata_up$ID,
                             symbol=plotdata_up$symbol,
                             dl_mu=plotdata_up$dl_mu,
                             LogFC = plotdata_up$LogFC,
                             PValue=plotdata_up$metap_up,
                             FDR=plotdata_up$adj_metap_up,
                             n_tested=plotdata_up$n_tested,
                             n_up=plotdata_up$n_up,
                             n_down=plotdata_up$n_down)
  plotdata_down <- res_dsl %>% filter(LogFC < 0)
  res_dsl_plot <- rbind(res_dsl_plot,data.frame(ID=plotdata_down$ID,
                                                symbol=plotdata_down$symbol,
                                                dl_mu=plotdata_down$dl_mu,
                                                LogFC = plotdata_down$LogFC,
                                                PValue=plotdata_down$metap_down,
                                                FDR=plotdata_down$adj_metap_down,
                                                n_tested=plotdata_down$n_tested,
                                                n_up=plotdata_down$n_up,
                                                n_down=plotdata_down$n_down))
  res_dsl_plot <- res_dsl_plot %>% mutate(sig = ifelse((FDR<=0.05 | FDR <= 0.05) & n_up | n_down >= n_tested/2 & abs(dl_mu) >= 0.5, "yes", "no"))
  
  
  res_dsl_plot$diffexpressed <- "unchanged"
  
  res_dsl_plot$diffexpressed[res_dsl_plot$dl_mu >= 0.5 &
                               res_dsl_plot$FDR<=0.05 &
                               res_dsl_plot$sig=="yes" ] <- "up"
  
  res_dsl_plot$diffexpressed[res_dsl_plot$dl_mu <= -0.5 &
                               res_dsl_plot$FDR<=0.05 &
                               res_dsl_plot$sig=="yes" ] <- "down"
  return(res_dsl_plot)
}

df <- lapply(meta_list, function(method){
  generateDataFrame(method,sefc,f_se,fc,pval_up,meta_up,meta_down)
})

plotdata <- lapply(df, function(method){
  generatPlotData(method)
})
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("down", "up", "unchanged")


# Method compare ------------
m.plots <- list()
for(study in 1:length(plotdata)){
  keyvals <- ifelse(
    plotdata[[study]]$LogFC < -0.5 & plotdata[[study]]$PValue < 0.05 , 'royalblue',
    ifelse(plotdata[[study]]$LogFC  > 0.5 & plotdata[[study]]$PValue < 0.05, 'red',
           'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'red'] <- 'up'
  names(keyvals)[keyvals == 'black'] <- 'NC'
  names(keyvals)[keyvals == 'royalblue'] <- 'down'
  p <- EnhancedVolcano(plotdata[[study]],
                       lab = plotdata[[study]]$symbol,
                       x = 'LogFC',
                       y = 'PValue',
                       colAlpha = 1,pCutoff = 0.05,
                       FCcutoff = 0.5,colCustom = keyvals,
                       title= names(plotdata[study]),
                       subtitle = "", drawConnectors = T) + theme_classic() +
    theme(legend.position = "none")
  m.plots[[names(plotdata[study])]] <- p
}

do.call("grid.arrange", c(m.plots, ncol=3))

m.plots[[1]]
pd$diffexpressed[pd$logFC >= 0.5 &
                           pd$FDR<=0.05 ] <- "Near Plaque"

pd$diffexpressed[pd$logFC <= -0.5 &
                           pd$FDR<=0.05] <- "No Pathology"

mycolors <- c("red", "blue", "black")
names(mycolors) <- c("Near Plaque", "No Pathology", "unchanged")

png(file.path(RESULTS_PATH,"Volcano_TauPS2APPNoPathologyvsNoPathology_Hippocampus.png"),
    width =1000,height = 800)
ggplot(data=pd, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=rownames(pd))) + 
  geom_point(size=2) + 
  theme_classic()+ geom_text_repel(aes(size=8),max.overlaps=2) +
  scale_colour_manual(values = mycolors) + ggtitle("PseudoBulk: TauPS2APP Pathology vs No Pathology Hippocampus")+
  theme(legend.text = element_text(size=16,face = "bold"),
        legend.title = element_text(size=18,face='bold'),
        legend.key.size = unit(1 ,'cm'),
        axis.title = element_text(size=20,face='bold'),
        plot.title = element_text(size=20,face='bold',hjust = 0.65)) +
  xlab("LogFoldChange") + ylab(bquote(bold(-log[10](PValue)))) +
  guides(fill = guide_legend(override.aes = aes(label = "")))
dev.off()
