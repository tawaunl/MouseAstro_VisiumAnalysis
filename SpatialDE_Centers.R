# Spatial DE_centroids

# Package Load ==================================================
suppressPackageStartupMessages({
  library("reticulate")
  library("ggplot2")
  library("scater")
  library("Seurat")
  library("ggrepel")
  library("tidyverse")
})
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


# Read in AnnData and Convert to SEurat objects ===============================
datasets = c('LIB5472833_SAM24434294','LIB5472834_SAM24434295', 'LIB5472835_SAM24434296',
             'LIB5472836_SAM24434297','LIB5472837_SAM24434298','LIB5472838_SAM24434299')

data_list <- list()
for(dataset in datasets){
  file <- file.path(data.dir,paste0(dataset,"_Spotdata.h5ad"))
  adata <- sc$read_h5ad(file)
  
  sdata <- H5toSeurat(adata = adata)
  data_list[[dataset]] <- sdata
}

# Process data =========================

for (n in names(data_list)) {
  data_list[[n]] <- subset(data_list[[n]], nCount_Spatial>0) #filter object that has spots with no counts for SCtransform
  data_list[[n]] <- SCTransform(data_list[[n]], assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
  # also run standard log normalization for comparison
  data_list[[n]] <- NormalizeData(data_list[[n]], verbose = FALSE, assay = "Spatial")
  
  # HVGs
  ignore_regex <- "(^Rp[sl]\\d+)|(^mt-)"
  sobj <- data_list[[n]]
  g <- sobj@assays$SCT@var.features
  sobj@assays$SCT@var.features <- g[!grepl(ignore_regex, g, ignore.case=T)]
  
  # Dimensional reduction and clustering
  sobj <- RunPCA(sobj, npcs=30,assay = "SCT",verbose = FALSE)
  sobj <- FindNeighbors(sobj, reduction="pca", dims=1:30)
  sobj <- FindClusters(sobj, resolution=c(0.8, 1.2, 1.6))
  sobj <- RunUMAP(sobj, assay="SCT", reduction="pca",dims = 1:30)
  sobj <- RunTSNE(sobj, assay="SCT", reduction="pca", dims=1:30)
  
  # Save
  data_list[[n]] <- sobj
  
}

nontg <- c("LIB5472836_SAM24434297","LIB5472838_SAM24434299")
# Add Sample information to cells
s <- names(data_list)
samples <- sapply(s, function(x){
  str_split(x,"_")[[1]][2]
})

for(dataset in 1:length(data_list)){ #create unique colnames 
  colnames(data_list[[dataset]]) <- paste0(samples[dataset],"_",colnames(data_list[[dataset]]))
  data_list[[dataset]]$orig.ident <- data_list[[dataset]]$region
  if(data_list[[dataset]]$orig.ident[1] %in% nontg){
    data_list[[dataset]]$Genotype <- "NonTG"
  } else{
    data_list[[dataset]]$Genotype <- "TauPS2APP"
  }
}

# Save indivdual Seurats ================
saveRDS(data_list,file.path(data.dir,"ProcessedSpatialObjects.rds"))
# Merge Objects ==========================
sobj_merge <- merge(x=data_list[[1]], y=c(data_list[[2]],data_list[[3]],
                                          data_list[[4]],data_list[[5]],data_list[[6]]))

DefaultAssay(sobj_merge) <- "SCT"
VariableFeatures(sobj_merge) <- c(VariableFeatures(data_list[[1]]), VariableFeatures(data_list[[2]]),
                                  VariableFeatures(data_list[[3]]), VariableFeatures(data_list[[4]]),
                                  VariableFeatures(data_list[[5]]), VariableFeatures(data_list[[6]]))
sobj_merge <- RunPCA(sobj_merge,assay = "SCT", verbose = FALSE)
sobj_merge <- FindNeighbors(sobj_merge, dims = 1:30)
sobj_merge <- FindClusters(sobj_merge, verbose = FALSE)
sobj_merge <- RunUMAP(sobj_merge, assay="SCT", reduction="pca",dims = 1:30)                   

DimPlot(sobj_merge, reduction = "umap", group.by = c("ident", "orig.ident"))
# Save Merged Seurat ================
saveRDS(sobj_merge,file.path(data.dir,"MergedSpatialObject.rds"))

# Puesdobulk DE AD vs Control Analysis =================================
library(SingleCellExperiment)
library(wesanderson)
library(scater)
library(scales)
library(edgeR)


counts <- GetAssayData(object = sobj_merge,
                       layer = "data",
                       assay="SCT") # get counts table
se <- SingleCellExperiment(assays = list(counts = counts),
                           colData = sobj_merge@meta.data)

set.seed(824)

# Puesdobulk DE Pathology Centers vs No ALL Analysis =================================
sum_by <- c("Abeta_centers","orig.ident")
summed <- aggregateAcrossCells(se, id=colData(se)[,sum_by])

# Creating up a DGEList object for use in edgeR:

y <- DGEList(counts(summed), samples=colData(summed))
discarded <- summed$ncells < 10
y <- y[,!discarded]

keep <- filterByExpr(y, group=summed$Abeta_centers)
y <- y[keep,]
y <- calcNormFactors(y)

design <- model.matrix(~0+factor(Abeta_centers), y$samples)
colnames(design) <- c("Abeta","NoPathology")
contr <- makeContrasts(Abeta-NoPathology,levels=design ) # NoPathology as reference
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust=TRUE)
AllPath_de <- glmQLFTest(fit, coef=ncol(design), contrast = contr)
AllPath_de <- data.frame(topTags(AllPath_de, n=50000,p.value = 10))


## Volcano Plot =============================== 
AllPath_de$diffexpressed <- "unchanged"

AllPath_de$diffexpressed[AllPath_de$logFC >= 0.5 &
                           AllPath_de$FDR<=0.05 ] <- "Near Plaque"

AllPath_de$diffexpressed[AllPath_de$logFC <= -0.5 &
                           AllPath_de$FDR<=0.05] <- "No Pathology"

mycolors <- c("red", "blue", "black")
names(mycolors) <- c("Near Plaque", "No Pathology", "unchanged")
p <- ggplot(data=AllPath_de, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=rownames(AllPath_de))) + 
  geom_point(size=2) + 
  theme_classic() +
  scale_colour_manual(values = mycolors) + ggtitle("PseudoBulk: Pathology vs No Pathology All Spots Centers")+
  theme(legend.text = element_text(size=16,face = "bold"),
        legend.title = element_text(size=18,face='bold'),
        legend.key.size = unit(1 ,'cm'),
        axis.title = element_text(size=20,face='bold'),
        plot.title = element_text(size=20,face='bold',hjust = 0.5)) +
  xlab("LogFoldChange") + ylab(bquote(bold(-log[10](PValue)))) + 
  geom_text_repel(max.overlaps=25,box.padding = 4,point.padding = 2)

png(file.path(RESULTS_PATH,"Volcano_AllNoPathologyvsNoPathology_Centers.png"),
    width =1000,height = 800)
p + guides(fill = guide_legend(override.aes = aes(color = NA)))
dev.off()

write.csv(AllPath_de,file.path(RESULTS_PATH,"AllSpots_DE_PathologyvsNoPathology.csv"))

# Puesdobulk DE: TauPS2APP Pathology vs No Analysis =================================
subset_se <- se[,se$Genotype=="TauPS2APP"]
sum_by <- c("Abeta_centers","orig.ident")
summed <- aggregateAcrossCells(subset_se, id=colData(subset_se)[,sum_by])

# Creating up a DGEList object for use in edgeR:

y <- DGEList(counts(summed), samples=colData(summed))
discarded <- summed$ncells < 10
y <- y[,!discarded]

keep <- filterByExpr(y, group=summed$Abeta_centers)
y <- y[keep,]
y <- calcNormFactors(y)

design <- model.matrix(~0+factor(Abeta_centers), y$samples)
colnames(design) <- c("Abeta","NoPathology")
contr <- makeContrasts(Abeta-NoPathology,levels=design ) # NoPathology as reference
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust=TRUE)
AllPath_de <- glmQLFTest(fit, coef=ncol(design), contrast = contr)
AllPath_de <- data.frame(topTags(AllPath_de, n=50000,p.value = 10))

## Volcano Plot =============================== 
AllPath_de$diffexpressed <- "unchanged"

AllPath_de$diffexpressed[AllPath_de$logFC >= 0.5 &
                           AllPath_de$FDR<=0.05 ] <- "Near Plaque"

AllPath_de$diffexpressed[AllPath_de$logFC <= -0.5 &
                           AllPath_de$FDR<=0.05] <- "No Pathology"

mycolors <- c("red", "blue", "black")
names(mycolors) <- c("Near Plaque", "No Pathology", "unchanged")

png(file.path(RESULTS_PATH,"Volcano_TauPS2APPNoPathologyvsNoPathology_AllRegions_Centers.png"),
    width =1000,height = 800)
ggplot(data=AllPath_de, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=rownames(AllPath_de))) + 
  geom_point(size=2) + 
  theme_classic()+ geom_text_repel(aes(size=8),max.overlaps=2) +
  scale_colour_manual(values = mycolors) + ggtitle("PseudoBulk: TauPS2APP Pathology vs No Pathology All Spots Centers")+
  theme(legend.text = element_text(size=16,face = "bold"),
        legend.title = element_text(size=18,face='bold'),
        legend.key.size = unit(1 ,'cm'),
        axis.title = element_text(size=20,face='bold'),
        plot.title = element_text(size=20,face='bold',hjust = 0.65)) +
  xlab("LogFoldChange") + ylab(bquote(bold(-log[10](PValue)))) +
  guides(fill = guide_legend(override.aes = aes(label = "")))
dev.off()

write.csv(AllPath_de,file.path(RESULTS_PATH,"TauPS2APP_DE_PathologyvsNoPathology_AllRegions_Centers.csv"))





# Puesdobulk DE: Cortex TauPS2APP Pathology vs No Analysis =================================
cortex <- subset_se[,subset_se$BrainRegion=="Cortex"]
sum_by <- c("Abeta_centers","orig.ident")
summed <- aggregateAcrossCells(cortex, id=colData(cortex)[,sum_by])

# Creating up a DGEList object for use in edgeR:

y <- DGEList(counts(summed), samples=colData(summed))
discarded <- summed$ncells < 10
y <- y[,!discarded]

keep <- filterByExpr(y, group=summed$Abeta)
y <- y[keep,]
y <- calcNormFactors(y)

design <- model.matrix(~0+factor(Abeta_centers), y$samples)
colnames(design) <- c("Abeta","NoPathology")
contr <- makeContrasts(Abeta-NoPathology,levels=design ) # NoPathology as reference
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust=TRUE)
AllPath_de <- glmQLFTest(fit, coef=ncol(design), contrast = contr)
AllPath_de <- data.frame(topTags(AllPath_de, n=50000,p.value = 10))

## Volcano Plot =============================== 
AllPath_de$diffexpressed <- "unchanged"

AllPath_de$diffexpressed[AllPath_de$logFC >= 0.5 &
                           AllPath_de$FDR<=0.05 ] <- "Near Plaque"

AllPath_de$diffexpressed[AllPath_de$logFC <= -0.5 &
                           AllPath_de$FDR<=0.05] <- "No Pathology"

mycolors <- c("red", "blue", "black")
names(mycolors) <- c("Near Plaque", "No Pathology", "unchanged")

png(file.path(RESULTS_PATH,"Volcano_TauPS2APPNoPathologyvsNoPathology_Cortex_centers.png"),
    width =1000,height = 800)
ggplot(data=AllPath_de, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=rownames(AllPath_de))) + 
  geom_point(size=2) + 
  theme_classic()+ geom_text_repel(aes(size=8),max.overlaps=2) +
  scale_colour_manual(values = mycolors) + ggtitle("PseudoBulk: TauPS2APP Pathology vs No Pathology Cortex Centers")+
  theme(legend.text = element_text(size=16,face = "bold"),
        legend.title = element_text(size=18,face='bold'),
        legend.key.size = unit(1 ,'cm'),
        axis.title = element_text(size=20,face='bold'),
        plot.title = element_text(size=20,face='bold',hjust = 0.65)) +
  xlab("LogFoldChange") + ylab(bquote(bold(-log[10](PValue)))) +
  guides(fill = guide_legend(override.aes = aes(label = "")))
dev.off()

write.csv(AllPath_de,file.path(RESULTS_PATH,"TauPS2APP_DE_PathologyvsNoPathology_Cortex_Centers.csv"))

# Puesdobulk DE: Hippocampus TauPS2APP Pathology vs No Analysis =================================
Hippo <- subset_se[,subset_se$BrainRegion=="Hipp"]
sum_by <- c("Abeta_centers","orig.ident")
summed <- aggregateAcrossCells(Hippo, id=colData(Hippo)[,sum_by])

# Creating up a DGEList object for use in edgeR:

y <- DGEList(counts(summed), samples=colData(summed))
discarded <- summed$ncells < 10
y <- y[,!discarded]

keep <- filterByExpr(y, group=summed$Abeta_centers)
y <- y[keep,]
y <- calcNormFactors(y)

design <- model.matrix(~0+factor(Abeta_centers), y$samples)
colnames(design) <- c("Abeta","NoPathology")
contr <- makeContrasts(Abeta-NoPathology,levels=design ) # NoPathology as reference
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust=TRUE)
AllPath_de <- glmQLFTest(fit, coef=ncol(design), contrast = contr)
AllPath_de <- data.frame(topTags(AllPath_de, n=50000,p.value = 10))

## Volcano Plot =============================== 
AllPath_de$diffexpressed <- "unchanged"

AllPath_de$diffexpressed[AllPath_de$logFC >= 0.5 &
                           AllPath_de$FDR<=0.05 ] <- "Near Plaque"

AllPath_de$diffexpressed[AllPath_de$logFC <= -0.5 &
                           AllPath_de$FDR<=0.05] <- "No Pathology"

mycolors <- c("red", "blue", "black")
names(mycolors) <- c("Near Plaque", "No Pathology", "unchanged")

png(file.path(RESULTS_PATH,"Volcano_TauPS2APPNoPathologyvsNoPathology_Hippocampus_Centers.png"),
    width =1000,height = 800)
ggplot(data=AllPath_de, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=rownames(AllPath_de))) + 
  geom_point(size=2) + 
  theme_classic()+ geom_text_repel(aes(size=8),max.overlaps=2) +
  scale_colour_manual(values = mycolors) + ggtitle("PseudoBulk: TauPS2APP Pathology vs No Pathology Hippocampus Centers")+
  theme(legend.text = element_text(size=16,face = "bold"),
        legend.title = element_text(size=18,face='bold'),
        legend.key.size = unit(1 ,'cm'),
        axis.title = element_text(size=20,face='bold'),
        plot.title = element_text(size=20,face='bold',hjust = 0.65)) +
  xlab("LogFoldChange") + ylab(bquote(bold(-log[10](PValue)))) +
  guides(fill = guide_legend(override.aes = aes(label = "")))
dev.off()

write.csv(AllPath_de,file.path(RESULTS_PATH,"TauPS2APP_DE_PathologyvsNoPathology_Hippocampus_Centers.csv"))

# 4 way plot Hipp vs Cortex ================================
library(gg4way)
hipp <- read.csv(file.path(RESULTS_PATH,"TauPS2APP_DE_PathologyvsNoPathology_Hippocampus_Centers.csv"))
cortex <- read.csv(file.path(RESULTS_PATH,"TauPS2APP_DE_PathologyvsNoPathology_Cortex_Centers.csv"))
colnames(hipp)[1] <- "symbol"
colnames(cortex)[1] <- "symbol"

feat <- as.data.frame(genomitory::getFeatures("GMTY17:GRCm38/GRCm38.IGIS4.0.genes.rds@REVISION-2"))
hipp$ID   <-  feat$ID[match(hipp$symbol,feat$symbol)]
cortex$ID <-  feat$ID[match(cortex$symbol,feat$symbol)]

hipp <- hipp %>% dplyr::rename(adj.P.Val = FDR)
cortex <- cortex %>% dplyr::rename(adj.P.Val = FDR)

x <- list("Cortex:Pathology vs No Pathology" = cortex,
          "Hippocampus:Pathology vs No Pathology" = hipp)
DAM <- c("Apoe","Axl","Bhlhe40","Clec7a","Csf1","Cst7","Ctsb","Ctsd",
         "Ctsl","Cybb","Fabp5","Fth1","Itgax","Gnas","Gpnmb",
         "Grn","Il1b","Lgals3","Lilrb4","Lpl","Lyz2","Mir155",
         "Msr1","Nos2","Spp1","Tfec","Trem2","Tyrobp","Vegfa")
DAA1 <- c("Gfap", "Id3", "Mt2", "Id4", "Cd81","Aqp4", "Myoc", "Igfbp5", "Prdx6", "Gja1")
DAA2 <- c("Cst3", "Apoe", "Mt1", 'Aldoc', "Clu", "Ckb", "Sparcl1", "Mt3", "Glu1" , "Slc1a2")

labels <- c(DAM,DAA1,DAA2)

gg4way(x,
       x = "Cortex:Pathology vs No Pathology", y= "Hippocampus:Pathology vs No Pathology",
       sep = " vs ",
       FDR = "adj.P.Val", logFCcutoff=0.25,FDRcutoff=.1,
       label=labels, colorVector=c("darkgrey", "firebrick", "forestgreen", "mediumblue")) +
  xlab(expression(atop(
    paste("Higher in No Pathology" %<->% "Higher in Pathology Centers"),
    paste("Cortex LogFC")))) +
  ylab(expression(atop(
    paste("Hippocampus LogFC"),
    paste("Higher in No Pathology" %<->% "Higher in Pathology Centers"))))


# Module scoring by Region ========================================
markers <- readRDS("~/Documents/AstrocytePaper/IntegratedMarkers.rds")
ignore_regex <- "(^Rp[sl]\\d+)|(^mt-)"

# Homeostatic
genes_homeo <- data.frame(markers$statistics$`1`)
genes_homeo <- genes_homeo[order(genes_homeo$logFC,decreasing = TRUE),]
genes_homeo <- rownames(genes_homeo)[1:20]
genes_homeo <- genes_homeo[!grepl(ignore_regex,genes_homeo,ignore.case = T)]

# DAA1
genes_DAA1 <- data.frame(markers$statistics$`2`)
genes_DAA1 <- genes_DAA1[order(genes_DAA1$logFC,decreasing = TRUE),]
genes_DAA1 <- rownames(genes_DAA1)[1:20]
genes_DAA1 <- genes_DAA1[!grepl(ignore_regex,genes_DAA1,ignore.case = T)]


# DAA2
genes_DAA2 <- data.frame(markers$statistics$`3`)
genes_DAA2 <- genes_DAA2[order(genes_DAA2$logFC,decreasing = TRUE),]
genes_DAA2 <- rownames(genes_DAA2)[1:20]
genes_DAA2 <- genes_DAA2[!grepl(ignore_regex,genes_DAA2,ignore.case = T)]


# Synapse-related
genes_synapse <- data.frame(markers$statistics$`4`)
genes_synapse <- genes_synapse[order(genes_synapse$logFC,decreasing = TRUE),]
genes_synapse <- rownames(genes_synapse)[1:25]
genes_synapse <- genes_synapse[!grepl(ignore_regex,genes_synapse,ignore.case = T)]


# Astro Substates
feature_list <- list(HomeostaticModule=genes_homeo, DAA1Module=DAA1,
                     DAA2Module=DAA2, SynapseModule=genes_synapse)
for (n in names(data_list)) {
  data_list[[n]] <- AddModuleScore(data_list[[n]], features=feature_list,  name=names(feature_list),
                                   assay = "SCT", slot = "data")
  names(data_list[[n]]@meta.data) <-  str_replace(names(data_list[[n]]@meta.data), "Module\\d", "")
}

sobj_merge <- AddModuleScore(sobj_merge, features=feature_list,  name=names(feature_list),
                             assay = "SCT", slot = "data")
names(sobj_merge@meta.data) <- str_replace(names(feature_list), "Module", "")

VlnPlot(sobj_merge,features = str_replace(names(feature_list), "Module", ""),
        group.by = "BrainRegion",ncol=2,split.by = "Genotype",assay = "SCT",
        layer="data") + 
  stat_summary(fun= median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5)


library(RColorBrewer)
library(viridis)
coldata <- as.data.frame(sobj_merge@meta.data)

abundances <- table(sobj_merge$BrainRegion,sobj_merge$Sample) 
abundances <- unclass(abundances) 

extra.info <- coldata[match(colnames(abundances), sobj_merge$Sample),]
d <- DGEList(abundances, samples=extra.info)
d = calcNormFactors(d)
d= estimateCommonDisp(d, verbose=TRUE)

norm_counts <- as.data.frame(t(d$counts)) 
norm_counts <- norm_counts %>% mutate(Sample = rownames(.))

coldata_short <- coldata %>% dplyr::select(Sample,Abeta_centers,Condition,Homeostatic,DAA1,DAA2,Synapse,BrainRegion) %>% unique


df_long_final <- left_join(norm_counts,coldata_short, by="Sample") 

plotdata <- aggregate(cbind(Homeostatic,DAA1,DAA2,Synapse) ~Sample + BrainRegion +Condition+Abeta_centers, data = df_long_final, median, na.rm = TRUE)
data <- pivot_longer(plotdata,cols = c(Homeostatic,DAA1,DAA2,Synapse))
data_tau <- data[data$Condition=="TauPS2APP",]


ggplot(data, aes(x=BrainRegion, y=value, fill=Abeta_centers)) + theme_classic() +
  geom_boxplot(outlier.shape=NA) + scale_fill_viridis(discrete = TRUE, alpha=1) +
  geom_jitter(width = 0.1,size=1) + 
  facet_wrap(~factor(name),scales = "free", ncol=2) + xlab("Brain region") +
  theme(axis.text.x = element_text(angle=75,vjust = 0.5),axis.text.y = element_text(size=20)) +
  ylab("Module Score") +  
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),legend.text = element_text(size=14),
        legend.title = element_text(size=16)) 

