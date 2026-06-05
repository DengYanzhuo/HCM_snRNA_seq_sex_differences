#### Fig.S2A ----
options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) 
library(data.table) 
library(Seurat)
setwd("~/figure/大群/16_WGCNA/")

seurat.obj <- readRDS("~/data/seurat_obj_group1.rds")
seurat.obj <- NormalizeData(seurat.obj)
expr <- AverageExpression(seurat.obj, group.by = "biosample_id", assays = "RNA")[["RNA"]]
expr <- as.data.frame(expr)
dim(expr)
expr2 <- as.data.frame(AverageExpression(seurat.obj, group.by = c("group1","biosample_id"), assays = "RNA")[["RNA"]])

name <- as.data.frame(colnames(expr2))
name <- str_split_fixed(name$`colnames(expr2)`,"_",2)  # 以_为分割，分成两列
name <- as.data.frame(name)
table(name$V1)
expr <- t(expr)
expr <- as.data.frame(expr)
rownames(expr) <- paste0(rownames(expr),"#",name$V1)  # 添加分组信息
bioid <- t(as.data.frame(strsplit(rownames(expr),"#")))
bioid <- as.data.frame(bioid)
expr$biosampleid <- bioid$V1
colnames(name) <- c("group","biosampleid")
expr <- merge(expr,name,by="biosampleid")
rownames(expr) <- paste0(expr$biosampleid,"#",expr$group)
expr <- subset(expr, select = -c(biosampleid, group))
expr[1:4,1:4]

expr <- t(expr)
keep_data <- expr[order(apply(expr, 1, mad),decreasing = T)[1:7000],] # mad 中位数绝对偏差
dim(keep_data)
keep_data[1:4,1:4]
keep_data <- t(keep_data)
keep_data <- as.data.frame(keep_data)
dim(keep_data)
keep_data[1:4,1:4]

df <- data.frame(rownames(keep_data))
a <- as.data.frame(str_split_fixed(df$rownames.keep_data.,"#",2))
df <- cbind(df,a)
colnames(df) <- c("AverageExpression_group","Sample","group")

datTraits <- data.frame(row.names = rownames(keep_data),
                        group = df$group)
table(datTraits$group)

grouptype <- data.frame(group=sort(unique(datTraits$group)),
                        groupNo=1:length(unique(datTraits$group)))
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]}
head(datTraits)
table(datTraits$group)
table(datTraits$groupNo)

datExpr0 <- as.data.frame(keep_data)
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){ 
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

sampleTree <- hclust(dist(datExpr0), method = "average")
par(mar = c(0,5,2,0))
pdf("step1_Sample dendrogram.pdf",width = 12,height = 6)
p <- plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2,
          cex.axis = 1, cex.main = 1,cex.lab=1)
print(p)
dev.off()
sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)),
                                colors = rainbow(length(table(datTraits$group))),
                                signed = FALSE)
par(mar = c(1,4,3,1),cex=0.8)
pdf("step1_Sample dendrogram and trait.pdf",width = 12,height = 6)
p2 <- plotDendroAndColors(sampleTree, sample_colors,
                          groupLabels = "trait",
                          cex.dendroLabels = 0.8,
                          marAll = c(1, 4, 3, 1),
                          cex.rowText = 0.01,
                          main = "Sample dendrogram and trait")
print(p2)
dev.off()

group_list <- datTraits$group
dat.pca <- PCA(datExpr0, graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"), #"point","text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE,
                    col.ind = group_list,
                    axes.linetype=NA,
                    mean.point=F
) +
  theme(legend.position = "none")+
  coord_fixed(ratio = 1)
pca
ggsave(pca,filename= "step1_Sample PCA analysis2.pdf", width = 15, height = 15)

datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="~/data/WGCNA/step1_input2.Rdata")



#### Fig.S2B ----
load("~/data/WGCNA/step1_input2.Rdata")
R.sq_cutoff = 0.8
powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
sft <- pickSoftThreshold(datExpr, 
                         networkType = "unsigned",
                         powerVector = powers, 
                         RsquaredCut = R.sq_cutoff,  
                         verbose = 5)
pdf("step2_power-value.pdf",width = 16,height = 12)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=R.sq_cutoff ,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=100,col="red")
dev.off()

sft$powerEstimate
# 15
power = sft$powerEstimate
power 


save(sft, power, file = "~/data/WGCNA/step2_power_value2.Rdata")

rm(list = ls())  
load(file = "~/data/WGCNA/step1_input2.Rdata")
load(file = "~/data/WGCNA/step2_power_value2.Rdata")

net <- blockwiseModules(
  datExpr,
  power = power,
  maxBlockSize = ncol(datExpr),
  corType = "pearson",
  networkType = "unsigned",
  TOMType = "unsigned", 
  minModuleSize = 35,
  mergeCutHeight = 0.15,
  numericLabels = TRUE, 
  saveTOMs = T,
  verbose = 3 
)
table(net$colors) 

# R^2>0.8 power=15
#  0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
#  3954  806  614  354  255  166  133  107   77   70   66   65   64   62   58   55   55   39

moduleColors <- labels2colors(net$colors)
table(moduleColors)
pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
save(net, moduleColors, file = "~/data/WGCNA/step3_genes_modules.Rdata")


#### Fig.S2C ----
rm(list = ls())
load(file = '~/data/WGCNA/step1_input.Rdata')
load(file = "~/data/WGCNA/step2_power_value.Rdata")
load(file = "~/data/WGCNA/step3_genes_modules.Rdata")
load(file = "~/data/WGCNA/step4_design.Rdata")

MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
MET = orderMEs(MEs)

pdf("step13_module_cor_Eigengene-dendrogram.pdf",width = 8,height = 10)
plotEigengeneNetworks(MET, setLabels="",
                      marDendro = c(0,4,1,4),  # 留白：下右上左
                      marHeatmap = c(5,5,1,2), # 留白：下右上左
                      cex.lab = 0.8,
                      xLabelsAngle = 90)
dev.off()


#### Fig.S2D ----

load(file = '~/data/WGCNA/step1_input.Rdata')
load(file = "~/data/WGCNA/step2_power_value.Rdata")
load(file = "~/data/WGCNA/step3_genes_modules.Rdata")
load(file = "~/data/WGCNA/step4_design.Rdata")
load(file = "~/data/WGCNA/step6_module_GO_term.Rdata")
module_term <- formula_res@compareClusterResult

brown <- module_term[module_term$Description %in% c(  "muscle cell development",
                                                      "muscle system process",
                                                      "mitochondrial transmembrane transport",
                                                      "cardiac muscle cell action potential",
                                                      "energy derivation by oxidation of organic compounds") & module_term$Cluster == "brown",]


midnightblue <- module_term[module_term$Description %in% c( "extracellular matrix organization",
                                                            "regulation of integrin-mediated signaling pathway",
                                                            "regulation of plasma membrane organization",
                                                            "external encapsulating structure organization",
                                                            "extracellular structure organization") & module_term$Cluster == "midnightblue",]

salmon <- module_term[module_term$Description %in% c("extracellular matrix organization",
                                                     "morphogenesis of a branching epithelium",
                                                     "smoothened signaling pathway",
                                                     "positive regulation of cytosolic calcium ion concentration",
                                                     "positive regulation of chondrocyte differentiation") & module_term$Cluster == "salmon",]

tan <- module_term[module_term$Description %in% c("small GTPase-mediated signal transduction",
                                                  "cell-matrix adhesion",
                                                  "cellular response to vascular endothelial growth factor stimulus",
                                                  "calcium ion homeostasis",
                                                  "positive regulation of Ras protein signal transduction") & module_term$Cluster == "tan",]

## brown Cardiomyocyte
module <- "brown"
deg <- read.csv("~/data/gsea.data.H/Cardiomyocyte_input.csv")
brown_genes <- unlist(strsplit(brown$geneID, "/"))
brown_genes <- bitr(brown_genes,
                    fromType = "ENTREZID",
                    toType = "SYMBOL",
                    OrgDb = "org.Hs.eg.db")
brown_genes <- brown_genes$SYMBOL

target_genes <- brown_genes

gene <- colnames(datExpr)
inTarget <- gene %in% target_genes
sub_deg <- subset(deg,deg$gene %in% target_genes)

modgene <- gene[inTarget]

TOM <- TOMsimilarityFromExpr(datExpr,power=power)
modTOM <- TOM[inTarget,inTarget]
dimnames(modTOM) <- list(modgene,modgene)

modgene_log2FC <- data.frame(
  gene = modgene,
  log2FC = sub_deg$avg_log2FC[match(modgene, sub_deg$gene)],
  stringsAsFactors = FALSE
)

nTop = 30
IMConn = softConnectivity(datExpr[, modgene])
top = (rank(-IMConn) <= nTop)
filter_modTOM <- modTOM[top, top]

top_genes <- unlist(dimnames(filter_modTOM)[1])
top_log2FC <- modgene_log2FC$log2FC[match(top_genes, modgene_log2FC$gene)]

node_attr <- data.frame(
  moduleColor = moduleColors[inTarget][top],
  log2FC = top_log2FC,
  direction = ifelse(top_log2FC > 0, "Up", "Down"),
  abs_log2FC = abs(top_log2FC),
  stringsAsFactors = FALSE
)

cyt <- exportNetworkToCytoscape(filter_modTOM,
                                edgeFile = paste("~/data/WGCNA/step8_CytoscapeInput-edges-", paste(module, collapse="-"), "2.txt", sep=""),
                                nodeFile = paste("~/data/WGCNA/step8_CytoscapeInput-nodes-", paste(module, collapse="-"), "2.txt", sep=""),
                                weighted = TRUE,
                                threshold = 0.05,
                                nodeNames = modgene[top],
                                nodeAttr = node_attr)



## midnightblue Fibroblast

module <- "midnightblue"
deg <- read.csv("~/data/gsea.data.H/Fibroblast_input.csv")
midnightblue_genes <- unlist(strsplit(midnightblue$geneID, "/"))
midnightblue_genes <- bitr(midnightblue_genes,
                           fromType = "ENTREZID",
                           toType = "SYMBOL",
                           OrgDb = "org.Hs.eg.db")
midnightblue_genes <- midnightblue_genes$SYMBOL

target_genes <- midnightblue_genes

gene <- colnames(datExpr)
inTarget <- gene %in% target_genes
sub_deg <- subset(deg,deg$gene %in% target_genes)

modgene <- gene[inTarget]

modTOM <- TOM[inTarget,inTarget]
dimnames(modTOM) <- list(modgene,modgene)

modgene_log2FC <- data.frame(
  gene = modgene,
  log2FC = sub_deg$avg_log2FC[match(modgene, sub_deg$gene)],
  stringsAsFactors = FALSE
)
nTop = 30
IMConn = softConnectivity(datExpr[, modgene]) 
top = (rank(-IMConn) <= nTop)
filter_modTOM <- modTOM[top, top]

top_genes <- unlist(dimnames(filter_modTOM)[1])
top_log2FC <- modgene_log2FC$log2FC[match(top_genes, modgene_log2FC$gene)]


node_attr <- data.frame(
  moduleColor = moduleColors[inTarget][top],
  log2FC = top_log2FC,
  direction = ifelse(top_log2FC > 0, "Up", "Down"),
  abs_log2FC = abs(top_log2FC),
  stringsAsFactors = FALSE
)

cyt <- exportNetworkToCytoscape(filter_modTOM,
                                edgeFile = paste("~/data/WGCNA/step8_CytoscapeInput-edges-", paste(module, collapse="-"), "2.txt", sep=""),
                                nodeFile = paste("~/data/WGCNA/step8_CytoscapeInput-nodes-", paste(module, collapse="-"), "2.txt", sep=""),
                                weighted = TRUE,
                                threshold = 0.01,
                                nodeNames = modgene[top],
                                nodeAttr = node_attr)


## salmon Fibroblast
module <- "salmon"
deg <- read.csv("~/data/gsea.data.H/Fibroblast_input.csv")
salmon_genes <- unlist(strsplit(salmon$geneID, "/"))
salmon_genes <- bitr(salmon_genes,
                     fromType = "ENTREZID",
                     toType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db")
salmon_genes <- salmon_genes$SYMBOL

target_genes <- salmon_genes

gene <- colnames(datExpr)
inTarget <- gene %in% target_genes
sub_deg <- subset(deg,deg$gene %in% target_genes)

modgene <- gene[inTarget]

modTOM <- TOM[inTarget,inTarget]
dimnames(modTOM) <- list(modgene,modgene)

modgene_log2FC <- data.frame(
  gene = modgene,
  log2FC = sub_deg$avg_log2FC[match(modgene, sub_deg$gene)],
  stringsAsFactors = FALSE
)

nTop = 30
IMConn = softConnectivity(datExpr[, modgene])
top = (rank(-IMConn) <= nTop) 
filter_modTOM <- modTOM[top, top]

manual_top_genes <- c("GLI2", "FAP", "ABI3BP", "GLI3", "COL14A1",
                      "DZIP1","AEBP1","COL12A1","PRRX1","POSTN")

manual_top_genes <- manual_top_genes[manual_top_genes %in% modgene]

top_idx <- match(manual_top_genes, modgene)
filter_modTOM <- modTOM[top_idx, top_idx]

top_genes <- manual_top_genes
top_log2FC <- modgene_log2FC$log2FC[match(top_genes, modgene_log2FC$gene)]

node_attr <- data.frame(
  moduleColor = moduleColors[inTarget][top_idx],
  log2FC = top_log2FC,
  direction = ifelse(top_log2FC > 0, "Up", "Down"),
  abs_log2FC = abs(top_log2FC),
  stringsAsFactors = FALSE
)


cyt <- exportNetworkToCytoscape(filter_modTOM,
                                edgeFile = paste("~/data/WGCNA/step8_CytoscapeInput-edges-", paste(module, collapse="-"), "2.txt", sep=""),
                                nodeFile = paste("~/data/WGCNA/step8_CytoscapeInput-nodes-", paste(module, collapse="-"), "2.txt", sep=""),
                                weighted = TRUE,
                                threshold = 0.05,
                                nodeNames = top_genes,
                                nodeAttr = node_attr)


## tan Adipocyte

module <- "tan"
deg <- read.csv("~/data/gsea.data.H/Adipocyte_input.csv")
tan_genes <- unlist(strsplit(tan$geneID, "/"))
tan_genes <- bitr(tan_genes,
                  fromType = "ENTREZID",
                  toType = "SYMBOL",
                  OrgDb = "org.Hs.eg.db")
tan_genes <- tan_genes$SYMBOL

target_genes <- tan_genes

gene <- colnames(datExpr)
inTarget <- gene %in% target_genes
sub_deg <- subset(deg,deg$gene %in% target_genes)

modgene <- gene[inTarget]

modTOM <- TOM[inTarget,inTarget]
dimnames(modTOM) <- list(modgene,modgene)

modgene_log2FC <- data.frame(
  gene = modgene,
  log2FC = sub_deg$avg_log2FC[match(modgene, sub_deg$gene)],
  stringsAsFactors = FALSE
)

nTop = 30
IMConn = softConnectivity(datExpr[, modgene])
top = (rank(-IMConn) <= nTop)
filter_modTOM <- modTOM[top, top]

top_genes <- unlist(dimnames(filter_modTOM)[1])
top_log2FC <- modgene_log2FC$log2FC[match(top_genes, modgene_log2FC$gene)]


node_attr <- data.frame(
  moduleColor = moduleColors[inTarget][top],
  log2FC = top_log2FC,
  direction = ifelse(top_log2FC > 0, "Up", "Down"),
  abs_log2FC = abs(top_log2FC),
  stringsAsFactors = FALSE
)


cyt <- exportNetworkToCytoscape(filter_modTOM,
                                edgeFile = paste("~/data/WGCNA/step8_CytoscapeInput-edges-", paste(module, collapse="-"), "2.txt", sep=""),
                                nodeFile = paste("~/data/WGCNA/step8_CytoscapeInput-nodes-", paste(module, collapse="-"), "2.txt", sep=""),
                                weighted = TRUE,
                                threshold = 0.05,
                                nodeNames = modgene[top],
                                nodeAttr = node_attr)
