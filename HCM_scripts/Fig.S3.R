#### Fig.S3A ----
remotes::install_version("igraph", version = "1.2.11")
library(igraph, lib.loc = "/home/dengyz/R/x86_64-pc-linux-gnu-library/4.3")
library(monocle)
library(Seurat)
library(dplyr)
library(ggplot2)
library(CLL)
library(ggridges)
library(RColorBrewer)
library(tibble)
set.seed(717)


seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Cardiomyocyte")


n <- 30000

joint_ratio <- prop.table(table(
  seurat_obj$group1, 
  seurat_obj$cell_type_leiden0.6
))
joint_ratio<-as.data.frame(joint_ratio)
joint_ratio <- subset(joint_ratio,Freq != 0)
colnames(joint_ratio) <- c("group","cell_type","ratio")

joint_ratio$sample_sizes <- round(joint_ratio$ratio * n)


meta_data <- seurat_obj@meta.data
meta_data$cell_id <- rownames(meta_data)

sampled_cells <- c()
# i <- 1
for(i in 1:nrow(joint_ratio)) {
  group <- as.character(joint_ratio$group[i])
  cell_type <- as.character(joint_ratio$cell_type[i])
  n_sample <- joint_ratio$sample_sizes[i]
  
  # 获取符合条件的细胞
  target_cells <- meta_data$cell_id[
    meta_data$group1 == group & 
      meta_data$cell_type_leiden0.6 == cell_type
  ]
  
  n_actual <- min(n_sample, length(target_cells))
  
  sampled <- sample(target_cells, n_actual)
  sampled_cells <- c(sampled_cells, sampled)
}

seurat_sampled <- subset(seurat_obj, cells = sampled_cells)

saveRDS(seurat_sampled,"~/02_work_data/HCM_NF/raw_data/sampled_Car_seurat_30000.rds")


seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/sampled_Car_seurat_30000.rds")
seurat_obj <- NormalizeData(seurat_obj)

data <- as(as.matrix(seurat_obj@assays$RNA@counts),'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data) 
head(row.names(data)) 
fData <- 
  data.frame(gene_short_name = row.names(data),
             row.names = row.names(data))
fd <- new('AnnotatedDataFrame',data=fData)

cds <- newCellDataSet(
  data,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

save(cds, file = "temp_cds_1.RData")
load("temp_cds_1.RData")

cds <-
  detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed>=10)) 

diff <- differentialGeneTest(cds[expressed_genes,],
                             fullModelFormulaStr = "~cell_type_leiden0.6",
                             
)

deg <- read.csv("~/02_work_data/HCM_NF/monocle/deg_30000.csv")
deg <- subset(deg, qval < 0.01)
deg <- deg[order(deg$qval), ][1:2000, ]
head(deg)

ordergene <- deg$X
cds <- setOrderingFilter(cds, ordergene)

plot_ordering_genes(cds)
ggsave("~/03_figure/HCM/Car/06_monocle2/step1_order_deg.pdf",plot = last_plot(),width = 5,height = 4)

save(cds, file = "temp_cds_2.RData")
load("temp_cds_2.RData")
cds <- reduceDimension(cds, max_components = 2, 
                       reduction_method = "DDRTree",
                       verbose = T) 

load("temp_cds_4.RData")
cds <- orderCells(cds,root_state = 1)

cds$group1 <- ""
cds$group1[cds$sex == 'female' & cds$disease == 'HCM'] <- 'HF'  
cds$group1[cds$sex == 'male' & cds$disease == 'HCM'] <- 'HM'  
cds$group1[cds$sex == 'female' & cds$disease == 'NF'] <- 'NF'  
cds$group1[cds$sex == 'male' & cds$disease == 'NF'] <- 'NM' 
table(cds$group1)

current_state <- as.character(cds$State)

new_state <- current_state
new_state[current_state == "2"] <- "temp"
new_state[current_state == "3"] <- "2"
new_state[new_state == "temp"] <- "3"

cds$State <- factor(new_state)


plot_cell_trajectory(cds,
                     color_by = "Pseudotime",
                     size = 1,
                     show_backbone = T)
ggsave("~/03_figure/HCM/Car/06_monocle2/step2_orderCell_time.pdf",plot = last_plot(),width = 5,height = 4)



#### Fig.S3B ----
plot_cell_trajectory(cds,
                     color_by = "State",
                     size = 1,
                     show_backbone = TRUE)
ggsave("~/03_figure/HCM/Car/06_monocle2/step2_orderCell_state.pdf",plot = last_plot(),width = 5,height = 4)



#### Fig.S3C ----
plot_cell_trajectory(cds,
                     color_by = "cell_type_leiden0.6",
                     size = 1,
                     show_backbone = TRUE)
ggsave("~/03_figure/HCM/Car/06_monocle2/step2_orderCell_state.pdf",plot = last_plot(),width = 5,height = 4)



#### Fig.S3D ----
plot_cell_trajectory(cds, color_by = "State", cell_size = 0.5) +facet_wrap(~group1,nrow = 1)+
  theme(aspect.ratio = 1,
        legend.title = element_blank())
ggsave("~/03_figure/HCM/Car/06_monocle2/step2_orderCell_group.pdf",plot = last_plot(),width = 10,height = 4)


#### Fig.S3E ----
plotdf=Biobase::pData(cds)
colnames(plotdf)
plotdf <- plotdf[,c(7,18,22)]

ggplot(plotdf, aes(x=Pseudotime,y=group1,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(0.5,11.5,17.3),linetype=3)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=30,face = "bold"))

ggsave("~/03_figure/HCM/Car/06_monocle2/step4_expr_time_group.pdf",plot = last_plot(),width = 6,height = 4)



#### Fig.S3F ----
plotdf=Biobase::pData(cds)
colnames(plotdf)
plotdf <- plotdf[,c(7,18,22)]

ggplot(plotdf, aes(x=Pseudotime,y=cell_type_leiden0.6,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(0.5,12.2,18.3),linetype=3)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  ) + theme(axis.title.y=element_text(vjust=2, size=30,face = "bold"))
ggsave("~/03_figure/HCM/Car/06_monocle2/step4_expr_time_type.pdf",plot = last_plot(),width = 6,height = 4)



#### Fig.S3G ----
deg <- read.csv("~/data/monocle_Car/deg_30000.csv")
deg <- subset(deg, qval < 0.01)
deg <- deg[order(deg$qval), ][1:2000, ]
head(deg)
ordergene <- deg$X

Time_diff <- differentialGeneTest(cds[ordergene,], cores =1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
state_diff <- differentialGeneTest(cds[ordergene,],cores = 1,
                                   fullModelFormulaStr = "~State")

state_diff <- state_diff[order(state_diff$qval),]

BEAM_res=BEAM(cds,branch_point = 1,cores = 1,progenitor_method = "duplicate")

save(BEAM_res, file = "~/data/monocle_Car/BEAM_res.RData")
saveRDS(BEAM_res,"~/data/monocle_Car/BEAM_res.rds")


load("~/data/monocle_Car/BEAM_res.RData")
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]

head(BEAM_res)
table(BEAM_res$qval < 1e-4)
tmp1=plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, qval < 1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 3,
                                 cores = 1,
                                 use_gene_short_name = TRUE,
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "RdBu")))(62),
                                 show_rownames = F,
                                 return_heatmap = T)

pdf("branched_heatmap2.pdf",width = 5,height = 6)
print(tmp1$ph_res)
dev.off()

#### Fig.S3H ----

ph <- if (inherits(tmp1, "pheatmap"))tmp1 else tmp1$ph
k <- 3
gene_clusters <- cutree(ph$tree_row, k = k)
beam_modules <- split(names(gene_clusters), gene_clusters)
for (i in seq_along(beam_modules)) {
  writeLines(beam_modules[[i]],
             con = sprintf("BEAM_branch%s_module%s_genes.txt", 1, i))}

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(xlsx)

diff <- read.delim("BEAM_branch1_module3_genes.txt",header=F) # CarI
bp1 <- enrichGO(diff$V1,
                OrgDb = org.Hs.eg.db,
                keyType = 'SYMBOL',
                ont = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)
term <- bp1@result
saveRDS(term, "GO/term3.rds")

term <- readRDS("GO/term3.rds")

go <- c(
  "regulation of axonogenesis",
  "muscle cell development",
  "developmental growth involved in morphogenesis",
  "intracellular calcium ion homeostasis",
  "positive regulation of MAPK cascade",
  "regulation of cell junction assembly",
  "collagen-activated signaling pathway",
  "calcium-mediated signaling")
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
pdf("GO/GO-C3.pdf",height = 8,width = 6)
p <- ggplot(data = term,
            aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "#67a4cc",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=5.5,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
print(p)
dev.off()


#### Fig.S3I ----
diff <- read.delim("BEAM_branch1_module2_genes.txt",header=F) # CarI
bp1 <- enrichGO(diff$V1,
                OrgDb = org.Hs.eg.db,
                keyType = 'SYMBOL',
                ont = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)
term <- bp1@result
saveRDS(term, "GO/term2.rds")

term <- readRDS("GO/term2.rds")
go <- c(
  "cardiac muscle contraction",
  "heart process",
  "myofibril assembly",
  "cardiac muscle tissue development",
  "sarcomere organization"
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
pdf("GO/GO-C2.pdf",height = 5,width = 6)
p <- ggplot(data = term,
            aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "#af93c4",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=5.5,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
print(p)
dev.off()



#### Fig.S3J ----
diff <- read.delim("BEAM_branch1_module1_genes.txt",header=F) # CarII
bp1 <- enrichGO(diff$V1,
                OrgDb = org.Hs.eg.db,
                keyType = 'SYMBOL',
                ont = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)
term <- bp1@result
saveRDS(term, "GO/term1.rds")

term <- readRDS("GO/term1.rds")
go <- c(
  "aerobic respiration",
  "ATP synthesis coupled electron transport",
  "oxidative phosphorylation",
  "energy derivation by oxidation of organic compounds",
  "cell-substrate adhesion"
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
pdf("GO/GO-C1.pdf",height = 5,width = 6)
p <- ggplot(data = term,
            aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "#f38989",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=5.5,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
print(p)
dev.off()



#### Fig.S3K ----
a <- read.delim("BEAM_branch1_module3_genes.txt",header=F) # CarIII
b <- ordergene
c <- as.data.frame(intersect(a$V1,b))
test_genes=c("JAK2","NRXN3","EDA2R")
p <- plot_genes_branched_pseudotime(cds[test_genes,],
                                    branch_point = 1,
                                    color_by = "cell_type_leiden0.6")
p
ggsave("genes_branched_pseudotime.pdf",width = 5,height = 6)


#### Fig.S3L ----

a <- read.delim("BEAM_branch1_module3_genes.txt",header=F) # CarIII
b <- ordergene
c <- as.data.frame(intersect(a$V1,b))
test_genes=c("JAK2","NRXN3","EDA2R")
p <- plot_genes_branched_pseudotime(cds[test_genes,],
                                    branch_point = 1,
                                    color_by = "cell_type_leiden0.6")
p
df <- p$data
group_colors <- c("#E32d32","#7FC97F","#4dbbd5ff","darkgoldenrod1")

gene <- "EDA2R"
df1 <- df[df$f_id == gene,]
df12 <- df1[df1$Branch =="Y_16",]
ggplot(df12, aes(Pseudotime, adjusted_expression)) +
  geom_smooth(aes(colour = group1),size = 1,fill = "grey",
              alpha = 0.2,method = "loess", se = T)+ # se = FALSE
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 6),
        axis.title.x = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))+
  scale_color_manual(values = group_colors)+
  ggtitle(gene)
ggsave(paste0("gene/",gene,".pdf"),height = 3,width=4.5)


gene <- "JAK2"
df1 <- df[df$f_id == gene,]
df12 <- df1[df1$Branch =="Y_16",]
ggplot(df12, aes(Pseudotime, adjusted_expression)) +
  geom_smooth(aes(colour = group1),size = 1,fill = "grey",
              alpha = 0.2,method = "loess", se = T)+ # se = FALSE
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 6),
        axis.title.x = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))+
  scale_color_manual(values = group_colors)+
  ggtitle(gene)
ggsave(paste0("gene/",gene,".pdf"),height = 3,width=4.5)


gene <- "NRXN3"
df1 <- df[df$f_id == gene,]
df12 <- df1[df1$Branch =="Y_16",]
ggplot(df12, aes(Pseudotime, adjusted_expression)) +
  geom_smooth(aes(colour = group1),size = 1,fill = "grey",
              alpha = 0.2,method = "loess", se = T)+ # se = FALSE
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 6),
        axis.title.x = element_text(color = "black", size = 8),
        plot.title = element_text(color = "black", size = 8))+
  scale_color_manual(values = group_colors)+
  ggtitle(gene)
ggsave(paste0("gene/",gene,".pdf"),height = 3,width=4.5)


