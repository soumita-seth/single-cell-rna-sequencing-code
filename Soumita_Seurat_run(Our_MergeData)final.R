setwd("D:/PhD Work/3rd Work/Restart Work")

#####https://github.com/debsin/dropClust####
#source("https://bioconductor.org/biocLite.R") 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install()
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocVersion")
BiocManager::install(version ='3.12')

#install.packages("devtools")
BiocManager::install(c("devtools"))
BiocManager::install(c("usethis"))
BiocManager::install(c("fs"))
BiocManager::install(c("RANN"))
BiocManager::install(c("cli"))
BiocManager::install(c("memoise"))
BiocManager::install(c("cachem"))
BiocManager::install(c("pkgbuild"))
BiocManager::install(c("callr"))
BiocManager::install(c("processx"))
BiocManager::install(c("ps"))
BiocManager::install(c("pkgload"))
BiocManager::install(c("withr"))
BiocManager::install(c("desc"))
BiocManager::install(c("testthat"))
BiocManager::install(c("remotes"))
BiocManager::install(c("sessioninfo"))

library(usethis)

library(devtools)

#install_github("debsin/dropClust", dependencies = T,force = TRUE)

#BiocManager::install(c("dropClust"))

#library(dropClust)

library(ggplot2)

####https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html#####
BiocManager::install(c("dplyr"))
library(dplyr)
BiocManager::install(c("Seurat"))

BiocManager::install(c("RColorBrewer"))
BiocManager::install(c("Rtsne"))
BiocManager::install(c("gtable"))
BiocManager::install(c("lifecycle"))
BiocManager::install(c("munsell"))
BiocManager::install(c("ellipsis"))
BiocManager::install(c("stringr"))
BiocManager::install(c("gridExtra"))
BiocManager::install(c("httr"))
BiocManager::install(c("ica"))
BiocManager::install(c("irlba"))
BiocManager::install(c("xtable"))
BiocManager::install(c("fastmap"))
BiocManager::install(c("lazyeval"))
BiocManager::install(c("png"))
BiocManager::install(c("abind"))

library(Seurat)

###download data from "https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html" by clicking on the 'here' word of "The raw data can be found here."##
###then unzip and save the folder in "R Script files" folder and set the path into the following command to read###
##pbmc.data <- Read10X(data.dir = "D:/PhD Work/3rd Work/pbmc3k_filtered_gene_bc_matrices.tar/filtered_gene_bc_matrices/hg19/")
#pbmc.data <-Read10X(data.dir = "D:/PhD Work/3rd Work/Final_Input_Matrix.csv")
scRNAseq.data <- read.csv("D:/PhD Work/3rd Work/Restart Work/merged_data.csv",row.names = 1)
#Initialize the Seurat object with the raw (non-normalized data).
class(scRNAseq.data)
scRNAseq.mat<-as.matrix(scRNAseq.data)
dim(scRNAseq.mat)
scRNAseq <- CreateSeuratObject(counts = scRNAseq.mat, project = "scRNA-seq", min.cells = 3, min.features = 200)
#scRNAseq
#class(pbmc.mat)
#dim(pbmc.mat)
dim(scRNAseq)
###pre-processing########

scRNAseq[["percent.mt"]] <- PercentageFeatureSet(scRNAseq, pattern = "^Mt")
class(scRNAseq[["percent.mt"]])
dim(scRNAseq[["percent.mt"]])
View(scRNAseq[["percent.mt"]])
View(scRNAseq[["nFeature_RNA"]])
scRNAseq[["nCount_RNA"]]<-as.integer(unlist(scRNAseq[["nCount_RNA"]]))
View(scRNAseq[["nCount_RNA"]])
#scRNAseq[["nCount_RNA"]]<-as.numeric(scRNAseq[["nCount_RNA"]])
# Visualize QC metrics as a violin plot
VlnPlot(scRNAseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
View(scRNAseq@meta.data)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(scRNAseq, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNAseq, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

scRNAseq <- subset(scRNAseq, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
class(scRNAseq)
dim(scRNAseq)

####Normalizing the data####
scRNAseq<- NormalizeData(scRNAseq, normalization.method = "LogNormalize", scale.factor = 10000)
#scRNAseq<- NormalizeData(scRNAseq)
scRNAseq1<-scRNAseq[["RNA"]]@data
class(scRNAseq1)
dim(scRNAseq1)
View(scRNAseq1)
scRNAseq2<-as.matrix(scRNAseq1)
class(scRNAseq2)
dim(scRNAseq2)
View(scRNAseq2)
###Identification of highly variable features (feature selection)###
scRNAseq <- FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNAseq), 10)
View(top10)
BiocManager::install(c("ggrepel"))

library(ggrepel)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scRNAseq)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge=0,ynudge=0)
plot2
#CombinePlots(plots = list(plot1, plot2))
#BiocManager::install(c("patchwork"))
#library(patchwork)
#plot1+plot2
#####Scaling the data###

all.genes <- rownames(scRNAseq)
length(all.genes)
class(all.genes)
View(all.genes)
scRNAseq <- ScaleData(scRNAseq, features = all.genes)
length(VariableFeatures(object = scRNAseq))
####Perform linear dimensional reduction###
scRNAseq <- RunPCA(scRNAseq, features = VariableFeatures(object = scRNAseq))
class(scRNAseq)
View(scRNAseq@meta.data)
scRNAseq[["pca"]]
View(scRNAseq[["pca"]])
# Examine and visualize PCA results a few different ways
#print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
print(scRNAseq[["pca"]], dims = 1:50, nfeatures = 5)
VizDimLoadings(scRNAseq, dims = 1:2, reduction = "pca")

DimPlot(scRNAseq, reduction = "pca")

DimHeatmap(scRNAseq, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(scRNAseq, dims = 1:15, cells = 500, balanced = TRUE)


####Determine the 'dimensionality' of the dataset###
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
scRNAseq <- JackStraw(scRNAseq, num.replicate = 100)
scRNAseq <- ScoreJackStraw(scRNAseq, dims = 1:20)
JackStrawPlot(scRNAseq, dims = 1:15)
dim(scRNAseq)

scRNAseq <- FindNeighbors(scRNAseq, dims = 1:10)
View(scRNAseq@meta.data)
scRNAseq_cluster <- FindClusters(scRNAseq, resolution = 0.5) ####need to add new multiobjective clustering algo (moga-svm)

View(scRNAseq_cluster@meta.data)
#print(scRNAseq_cluster@graph)
# Look at cluster IDs of the first 5 cells
head(Idents(scRNAseq_cluster), 5)
View(Idents(scRNAseq_cluster))
#View(scRNAseq)
#View(scRNAseq_cluster@meta.data)
#get.centroids(scRNAseq$seurat_clusters)
#View(scRNAseq_cluster$seurat_clusters)
#cluster1<-subset(scRNAseq_cluster,seurat_clusters==1)

#class(cluster1@meta.data)
#cluster1.points<-cluster1@meta.data[,2:4]
#View(cluster1.points)
#BiocManager::install(c("LICORS"))
#BiocManager::install(c("fields"))
#BiocManager::install(c("maps"))
#BiocManager::install(c("viridis"))
#BiocManager::install(c("FNN"))
#BiocManager::install(c("mvtnorm"))
#library(LICORS)
#BiocManager::install(c("pracma"))
#library(pracma)
#class(cluster1.points)
#cluster1.points.numeric <-data.matrix(cluster1.points, rownames.force = NA)
#View(cluster1.points.numeric)
#cluster1.Result<-kmeanspp(cluster1.points.numeric, k=1)
#View(cluster1.Result)
#View(cluster1.Result$cluster)
#View(cluster1.Result$centers)
#aggregate(cluster1.points.numeric, by=list(cluster=cluster1.Result$cluster), mean)
#dd <- cbind(cluster1.points.numeric, cluster = cluster1.Result$cluster )
#View(dd)
#BiocManager::install(c("CountClust"))
#library(CountClust)
#install.packages("remotes")
#remotes::install_github("paodan/studySeu")
 #doKmeans(cluster1,k.cells=0)
#install.packages("remotes")
#BiocManager::install(c("ncchung/SeuratAddon"))
#seurat_clusters <- as.data.frame(Idents(scRNAseq_cluster))
#dim(seurat_clusters)
#View(seurat_clusters)
#scRNAseq2.df<-as.data.frame(t(scRNAseq2))
#cluster_center <- aggregate(scRNAseq2.df,list(cluster=seurat_clusters),mean)
#cluster_center[,1:4]


#clust.centroid = function(i, dat, seurat_clusters) {
 # ind = (seurat_clusters == i)
 # colMeans(dat[ind,])
#}

#sapply(unique(seurat_clusters), clust.centroid,scRNAseq2.df , seurat_clusters)

#remotes::install_github("ncchung/SeuratAddon")
#remotes::install_github()
#Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
write.csv(Idents(scRNAseq_cluster),'D:\\PhD Work\\3rd Work\\Restart Work\\Cell_Cluster_Set_Latest.csv', row.names = TRUE)
####Run non-linear dimensional reduction (UMAP/tSNE)###
###download and install anaconda python software 3.7 version###
###reticulate::py_install(packages ='umap-learn')
#reticulate::py_install(packages ='umap-learn')
scRNAseq_new <- RunUMAP(scRNAseq_cluster, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
#DimPlot(pbmc, reduction = "umap")
DimPlot(scRNAseq_new, reduction = "umap", label=TRUE)

saveRDS(scRNAseq_new, file = "D:/PhD Work/3rd Work/Restart Work/scRNASeq_tutorial_latest.rds")
#cluster_pbmc<-readRDS(file="D:/PhD Work/3rd Work/pbmc_tutorial2.rds")
#View(cluster_pbmc)

####Finding differentially expressed features (cluster biomarkers)###
# find all markers of cluster 1
#BiocManager::install("DESeq2")
#library(DESeq2)
#class(scRNAseq_new[["RNA"]]@counts)
#as.integer(scRNAseq_new[["RNA"]]@counts)
#View(scRNAseq_new@meta.data)
#scRNAseq_new@active.assay
#expr_raw <- as.matrix(GetAssayData(object = scRNAseq_new, slot = "counts"))
#View(expr_raw)
BiocManager::install("MAST")
BiocManager::install("progress")
library(MAST)
cluster1.markers <- FindMarkers(scRNAseq_new, ident.1 = 1, test.use="MAST", min.pct = 0.25)
head(cluster1.markers, n = 5)
class(cluster1.markers)
dim(cluster1.markers)
View(cluster1.markers)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(scRNAseq, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
class(cluster5.markers)
dim(cluster5.markers)
View(cluster5.markers)
# find markers for every cluster compared to all remaining cells, report only the positive ones
scRNAseq.markers <- FindAllMarkers(scRNAseq_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST")
#pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
View(scRNAseq.markers)
dim(scRNAseq.markers)
class(scRNAseq.markers)
scRNAseq.markers.name<- scRNAseq.markers[,7]
View(scRNAseq.markers.name)
class(scRNAseq.markers.name)
scRNAseq.markers.name.df<-as.data.frame(scRNAseq.markers.name)
View(scRNAseq.markers.name.df)
dim(scRNAseq.markers.name.df)
#Occurences<-table(unlist(scRNAseq.markers.name.df)) // It also works correctly
#View(Occurences)
#dim(Occurences)
#scRNAseq.markers.name.df.temp<-as.data.frame(scRNAseq.markers.name.df[1:1000,])
#scRNAseq.markers.temp<-as.data.frame(scRNAseq.markers[1:1000,])
#View(scRNAseq.markers.temp)
#View(scRNAseq.markers.name.df.temp)
library(stringr)
#loop over the strings against the pattern from df2
scRNAseq.markers.name.df$Counts<-sapply(scRNAseq.markers$gene, function(x){
  sum(str_count(x, scRNAseq.markers.name.df$scRNAseq.markers.name))
})
View(scRNAseq.markers.name.df)
scRNAseq.markers.name.omitted_duplicated<-scRNAseq.markers.name.df[!duplicated(scRNAseq.markers.name.df$scRNAseq.markers.name),]
View(scRNAseq.markers.name.omitted_duplicated)
dim(scRNAseq.markers.name.omitted_duplicated)
scRNAseq.markers.name.omitted_duplicated_subset<-scRNAseq.markers.name.omitted_duplicated[scRNAseq.markers.name.omitted_duplicated$Counts!=1,]
View(scRNAseq.markers.name.omitted_duplicated_subset)
dim(scRNAseq.markers.name.omitted_duplicated_subset)
frequent_marker<-scRNAseq.markers.name.omitted_duplicated_subset$scRNAseq.markers.name
scRNAseq.markers.subset<-scRNAseq.markers[scRNAseq.markers$gene%in%frequent_marker,]
dim(scRNAseq.markers.subset)
View(scRNAseq.markers.subset)
scRNAseq.markers.name.df.subset<-scRNAseq.markers.name.df[scRNAseq.markers.name.df$Counts!=1,]
dim(scRNAseq.markers.name.df.subset)
scRNAseq.markers.subset.final<-cbind(scRNAseq.markers.subset,scRNAseq.markers.name.df.subset$Counts)
View(scRNAseq.markers.subset.final)
dim(scRNAseq.markers.subset.final)
write.csv(scRNAseq.markers.name.omitted_duplicated_subset,'D:\\PhD Work\\3rd Work\\Restart Work\\Cluster_specified_frequent_DEGs_Latest.csv', row.names = FALSE)
#count(scRNAseq.markers.name.df, vars = "scRNAseq.markers.name")
write.csv(scRNAseq.markers,'D:\\PhD Work\\3rd Work\\Restart Work\\Cluster_specified_DEGList_Latest.csv', row.names = FALSE)
write.csv(scRNAseq.markers.subset.final,'D:\\PhD Work\\3rd Work\\Restart Work\\Cluster_specified_frequent_DEGList_fullresult.csv', row.names = FALSE)
##ROC test returns the 'classification power' for any individual marker (ranging from 0 - random, to 1 - perfect)##
#cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#View(cluster1.markers)
#res<-scRNAseq.markers.name.omitted_duplicated_subset.merge(scRNAseq.markers, how="inner", left_on='scRNAseq.markers.name', right_on='gene')
#res<-merge(scRNAseq.markers,scRNAseq.markers.name.omitted_duplicated_subset,by=intersect(names(scRNAseq.markers),names(scRNAseq.markers.name.omitted_duplicated_subset)))
#dim(res)
#scRNAseq.markers.name.omitted_duplicated_subset1<-scRNAseq.markers.name.omitted_duplicated[scRNAseq.markers.name.omitted_duplicated$Counts>3,]
#dim(scRNAseq.markers.name.omitted_duplicated_subset1)
scRNAseq.markers.name.omitted_duplicated_subset2<-scRNAseq.markers.name.omitted_duplicated[scRNAseq.markers.name.omitted_duplicated$Counts>2,]
dim(scRNAseq.markers.name.omitted_duplicated_subset2)
#most_frequent_marker<-scRNAseq.markers.name.omitted_duplicated_subset1$scRNAseq.markers.name
most_frequent_marker<-scRNAseq.markers.name.omitted_duplicated_subset2$scRNAseq.markers.name
#scRNAseq_Frequency_filtered_Mat<-scRNAseq2[frequent_marker,]
scRNAseq_Frequency_filtered_Mat<-scRNAseq2[most_frequent_marker,]
dim(scRNAseq_Frequency_filtered_Mat)
View(scRNAseq_Frequency_filtered_Mat)
scRNAseq_Frequency_filtered_Mat_Transpose<-t(scRNAseq_Frequency_filtered_Mat)
dim(scRNAseq_Frequency_filtered_Mat_Transpose)
scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames<-scRNAseq_Frequency_filtered_Mat_Transpose
rownames(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames)<-NULL
View(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames)
dim(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames)
scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor<-cor(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames,method="spearman")
dim(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor)##[1892, 1892]|[85,85]|[543,543]
class(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor)
save(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor, file = "correlation_matrix_MoreFrequent_Markers_spearman.rda")
write.csv(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor,'D:\\PhD Work\\3rd Work\\Restart Work\\correlation_matrix_MoreFrequent_Markers_spearman.csv', row.names = TRUE)
scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1<-scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor
dim(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1)##[1892, 1892]|[85,85]|[543,543]
for (i in 1:nrow(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1)){
  #for(i in 1:5 ){
  print(paste0("i=",i))
  ##for(j in 1:i){###???wrong
  for(j in 1:nrow(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1)){
    if(i!=j){
      if(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1[i,j]>=0.5 | scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1[i,j]<=-0.5){
        scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1[i,j]=1
      }
      else{
        scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1[i,j]=0
      }
    }
    else{
      scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1[i,j]=-99
    }
  }
}
#View(myExpr_FDR_filtered_transpose_no_rownames.cor1[1:5,1:5])
#n1<-as.double(myExpr_FDR_filtered_transpose_no_rownames.cor1[1:5,1:5])
#View(n1)
##check how many 1s in matrix####newly added
howmany1<-sum(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1==1)
howmany1##6740|54|1198
totalelement_corrmat<-nrow(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1)*nrow(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1)
totalelement_corrmat##3579664|7225|294849
fraction_howmany1<-howmany1*100/totalelement_corrmat
fraction_howmany1##0.188%|0.747%|0.406%
################
View(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1)
save(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1, file = "Boolean_converted_corr_matrix_Morefrequent_Markers_spearman.rda")##newly added
write.csv(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1,'D:\\PhD Work\\3rd Work\\Restart Work\\Boolean_converted_corr_matrix_Morefrequent_Markers_spearman.csv', row.names = TRUE)
#write.table(myExpr_FDR_filtered_transpose_no_rownames.cor1,file="Boolean_converted_corr_matrix_DEG_pearson.txt",quote=F,sep="\t",row.names=TRUE,col.name=T)
#myExpr_FDR_filtered_transpose_no_rownames.cor1[1:5,1:5]
#colnames(myExpr_FDR_filtered_transpose_no_rownames.cor1)[2]
#rownames(myExpr_FDR_filtered_transpose_no_rownames.cor1)[3]
#strongly_correlated_genes <- matrix(0,nrow=nrow(myExpr_FDR_filtered_transpose_no_rownames.cor1),ncol=2)
####strongly_correlated_genes <- matrix(0, nrow = 2000000, ncol = 2, byrow = FALSE, dimnames = NULL)
###newly added####
strongly_correlated_genes <- matrix(NA, nrow = (howmany1/2), ncol = 2, byrow = FALSE, dimnames = NULL)
#View(strongly_correlated_genes)
#dim(strongly_correlated_genes)
row<-1
###col<-1##no need
for(i in 2:nrow(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1)){###corrected
  ##for(i in 2:5){###corrected
  print(paste0("i=",i))
  ###for(j in 1:i){###wrong
  for(j in 1:(i-1)){
    ##if(i!=j){
    print(paste0("j=",j))
    #print(paste(myExpr_FDR_filtered_transpose_no_rownames.cor1[i,j]))  
    if(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1[i,j]==1){
      s<-rownames(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1)[i]
      #print(paste("s=",s))
      ##strongly_correlated_genes[row,col]<-s###wrong
      strongly_correlated_genes[row,1]<-s##newly added
      #print(paste("row=",row))
      #print(paste("col=",col))
      ##col<-col+1##wrong
      d<-colnames(scRNAseq_Frequency_filtered_Mat_Transpose_no_rownames.cor1)[j]
      #print(paste("d=",d))
      ##strongly_correlated_genes[row,col]<-d###wrong
      strongly_correlated_genes[row,2]<-d
      #print(paste("row=",row))
      #print(paste("col=",col))
      ##col<-1##no need
      row<-row+1
    }
    ##}  
    
  }
}
dim(strongly_correlated_genes)##[3370, 2]|[27,2]|[599,2]
write.csv(strongly_correlated_genes,'D:\\PhD Work\\3rd Work\\Restart Work\\Strongly_Correlated_Markers(More Frequent).csv', row.names = FALSE)




















