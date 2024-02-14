#human1 scRNA data set
library(readr)
library(Seurat)
library(Matrix)





#READING THE DATA SET and Removing the X and assigned_cluster columns.

human1 <- read.csv("C:/ASHIM/scRNA/DATA_SET_HUMAN1/Human1.csv")
human1_data = human1
human1_data $ X <- NULL
human1_data $ assigned_cluster <- NULL

# just removing the "-" in barcode.
human1_data $barcode <- gsub("-", "", human1_data $barcode)
human1_data $ barcode <- row.names(human1_data)



#PRE-PROCESSING WORKFLOW FOR FURTHER ANALYSIS

# 1.extracting gene expression
gene_expression_data<-human1_data[,2:ncol(human1_data)]
gene_expression_matrix<- as.matrix(gene_expression_data)


# 2.creating a sparse matrix to count non zero cells
sparse_matrix <- as(gene_expression_matrix, "sparseMatrix")
sparse_matrix <- t(sparse_matrix)
sparse_matrix_dgc <- as(sparse_matrix, "dgCMatrix")



# 3.quality control(QC Metrics), REMOVING THE LOW QUALITY scRNA DATA
#creating a Seurat object
# min.cells = 3, this is a filtering parameter that filters out the Gene that is expressed in less than 3 cells.
# min.features = 200, this filters the cells that is not describing atleast 200 genes.

objHuman1<-CreateSeuratObject(counts = sparse_matrix_dgc, project = "HUMAN1", min.cells = 3, min.features = 200 )
objHuman1[c("A1BG", "A1CF", "A2M"), 1:30]
VlnPlot(objHuman1,features = c("nCount_RNA", "nFeature_RNA"),ncol = 2)






# NORMALIZING THE DATA SET

objHuman1 <- NormalizeData(objHuman1, normalization.method = "LogNormalize", scale.factor = 10000)




# Identifying the higg vAriable Genes

objHuman1 <- FindVariableFeatures(objHuman1, selection.method = "vst", nfeatures = 2000)
#identifying top 4 genes
top4 <- head(VariableFeatures(objHuman1), 4)
top4
#visualizing 
plot_variable <-VariableFeaturePlot(objHuman1)
plot_variable2 <- LabelPoints(plot = plot_variable, points = top4, repel = TRUE, xnudge = 0, ynudge = 0)
plot_variable2





# Scaling the data

objHuman1 <- ScaleData(objHuman1, features = rownames(objHuman1))





#perfirm linear dimes reduction

objHuman1 <- RunPCA(objHuman1, features = VariableFeatures(object = objHuman1))
#visualizing 
VizDimLoadings(objHuman1, dims = 1:2, reduction = "pca")
DimPlot(objHuman1, reduction = "pca") + NoLegend()






#Deetminig the dimension of the datasert

plotElbow <- ElbowPlot(objHuman1)
plotElbow
noPC <-7







#cluster the cells Using Seurat

objHuman1 <- FindNeighbors(objHuman1, dims = 1: noPC)
objHuman1 <- FindClusters(objHuman1, resolution = 1.0)
head(Idents(objHuman1), 5)





#Non-lineaer dimensional refuction

objHuman1 <- RunUMAP(objHuman1, dims = 1:noPC)
DimPlot(objHuman1, reduction = "umap")






