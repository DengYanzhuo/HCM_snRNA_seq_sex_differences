#### Fig.5A ----
library(ggalluvial)
library(tidyverse)
library(cols4all)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
subcell_color <- c(
  'Macrophage' = "#f38989",
  'Proliferating_macrophage' = "#815e99"
)

seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Macrophage")
pdf("~/03_figure/HCM/Mac/01_umap1.pdf",width = 6, height = 4)
DimPlot(seurat_obj, 
        reduction = "umap", 
        group.by='cell_type_leiden0.6',
        label = T,
        cols = subcell_color,
        repel = TRUE,raster = FALSE)
dev.off()


#### Fig.5B ----
library(Seurat)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(dplyr)
setwd("~/03_figure/HCM/Mac/02_亚群鉴定/")


seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
Idents(seurat_obj) <- "cell_type"
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Macrophage")
seurat_obj <- NormalizeData(seurat_obj)

Idents(seurat_obj) <- "cell_type_leiden0.6"
all_mk <-FindAllMarkers(seurat_obj)
write.csv(all_mk,"~/02_work_data/HCM_NF/ScRNA.marker/Mac_sub.csv")
all_mk <- read.csv("~/02_work_data/HCM_NF/ScRNA.marker/Mac_sub.csv")

top_mk <- all_mk %>%
  group_by(cluster) %>%
  top_n(wt = avg_log2FC, n =5)%>%
  pull(gene)
top_mk <- c("CD86","CD74","HLA-DRB1","F13A1","LGMN",
            # Macrophage
            "TUBB","H2AFZ","STMN1","TOP2A","MKI67"
            # Proliferating_macrophage
)

desired_order <- c("Macrophage", "Proliferating_macrophage")

p <- DotPlot(seurat_obj,
             features = top_mk,
             scale.by ='size')+
  scale_size(range = c(1, 6))+
  theme( axis.text.x.bottom =element_text(angle =45,hjust =1,vjust =1))+
  scale_y_discrete(limits = desired_order)
p
ggsave("~/03_figure/HCM/Mac/02_亚群鉴定/marker.pdf",plot = p, width = 8,height = 3)




#### Fig.5C ----
library(ggalluvial)
library(tidyverse)
library(cols4all)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
subcell_color <- c(
  'Macrophage' = "#f38989",
  'Proliferating_macrophage' = "#815e99"
)

seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Macrophage")
pdf("~/03_figure/HCM/Mac/01_umap1.pdf",width = 6, height = 4)
DimPlot(seurat_obj, 
        reduction = "umap", 
        group.by='cell_type_leiden0.6',
        label = T,
        cols = subcell_color,
        repel = TRUE,raster = FALSE)
dev.off()
cellratio<-prop.table(table(seurat_obj$cell_type_leiden0.6,seurat_obj$group1),margin=2)
cellratio<-as.data.frame(cellratio)
colnames(cellratio) <- c("cell_type_leiden0.6","group","ratio")
cellratio <- subset(cellratio,ratio != 0)

p <- ggplot(cellratio, aes(x = group, y = ratio, fill = cell_type_leiden0.6)) +
  geom_bar(position = "fill", stat="identity", color = 'white', alpha = 5, width = 0.95) +
  scale_fill_manual(values = subcell_color) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  coord_flip() 
p

pp <- ggplot(cellratio, aes(x = group, y = ratio, fill = cell_type_leiden0.6,
                            stratum = cell_type_leiden0.6, alluvium = cell_type_leiden0.6)) +
  scale_fill_manual(values = subcell_color) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  coord_flip() 
p1 <- pp +
  geom_col(width = 0.6, color = NA) +
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0) 
p1
p2 <- pp +
  geom_col(width = 0.6,
           color = 'white', size = 0.5) +
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0,
            color = 'white', size = 0.5) 
p2
ggsave(file.path("~/03_figure/HCM/Mac/01_umap_ratio/01_ratio1.pdf"),
       ggplot2::last_plot(),
       height= 2,
       width= 10)



#### Fig.5D ----
library(Seurat)
library(tidyverse)
library(dplyr)

seurat.obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat.obj <- subset(seurat.obj,subset = cell_type %in% "Macrophage")
Idents(seurat.obj) <- "cell_type_leiden0.6"

seurat.obj <- subset(seurat.obj, group == "HCM")
Idents(seurat.obj) <- "cell_type_leiden0.6"

cell.types <- c("Macrophage","Proliferating_macrophage")
for (cell.type in cell.types) {
  seurat <- subset(seurat.obj, cell_type_leiden0.6 == cell.type)
  Idents(seurat) <- "group1"
  seurat <- NormalizeData(seurat)
  
  gsea.input <- FindMarkers(
    seurat,
    ident.1 = "HF",
    ident.2 = "HM",
    min.pct = 0,
    logfc.threshold = 0
  )
  gsea.input$celltype <- cell.type
  gsea.input$gene <- rownames(gsea.input)
  
  write.table(
    gsea.input,
    file = paste0("~/02_work_data/HCM_NF/gsea.input/gsea_H/Macrophage/",cell.type,"_input.csv"),
    quote = F,
    sep = ",",
    row.names = F
  )
}

category <- c("C5") 
genesets <- msigdbr(species = "Homo sapiens", category = category)
genesets <- subset(genesets, select = c("gs_name", "gene_symbol"))

i <- "Macrophage"

gsea.input <-read.csv(paste0("~/02_work_data/HCM_NF/gsea.input/gsea_H/Macrophage/",i,"_input.csv"), sep = ",", header = T)
rownames(gsea.input) <- gsea.input$gene
gsea.input <- gsea.input[order(gsea.input$avg_log2FC, decreasing = T),]
genelist <- structure(gsea.input$avg_log2FC, names = rownames(gsea.input))

write.table(genelist, quote = F, sep = "\t", col.names = F, file = paste0("~/02_work_data/HCM_NF/gsea_data/gsea_H/Macrophage/",i,"_H.rnk"))

genelist <- genelist[which(genelist!=0)]
res <- GSEA(genelist, TERM2GENE = genesets, eps = 0 ,pvalue= 1 ,minGSSize = 1)
saveRDS(res, paste0("~/02_work_data/HCM_NF/gsea_data/gsea_H/Macrophage/",i,"_H.rds"))

write.csv(res, paste0("~/02_work_data/HCM_NF/gsea_data/gsea_H/Macrophage/",i,"_H.csv"), row.names = F)

library(GseaVis)
a = readRDS("~/02_work_data/HCM_NF/gsea_data/gsea_H/Macrophage/Macrophage_H.rds")
b <- a@result
Mac = b[b$pvalue<0.05, ]

gseaNb(
  object = a,
  geneSetID = i,
  addPval = TRUE,
  legend.position = "right",
  newGsea = TRUE,
  newCurveCol = c("#815e99", "grey", "#f38989"),
  newHtCol = c("#f38989", "white", "#815e99")
)

ggsave(paste0("~/03_figure/HCM/Mac/08_GSEA/",i,".pdf"), width = 5, height = 4) 



#### Fig.5E ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(psych)
library(reshape2)
library(tidyverse)
library(SeuratDisk)
library(msigdbr)
library(AUCell)
library(pheatmap)
library(plyr)
library(devtools)
library(ggunchained)
library(ggpubr)

seurat_obj <- readRDS("~/data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Macrophage")
seurat_obj <- subset(seurat_obj,subset = cell_type_leiden0.6 %in% "Macrophage")
seurat_obj <- NormalizeData(seurat_obj)

cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data) 
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
Dataset <- Dataset %>% filter(str_detect(Dataset$gs_name,"GOBP_"))
Dataset$gs_name <- Dataset$gs_name|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
a <- as.data.frame(unique(Dataset$gs_name))

feature_list <- c("Adaptive immune response",
                  "Cell activation involved in immune response",
                  "Leukocyte mediated immunity",
                  "Lymphocyte activation involved in immune response",
                  "Mononuclear cell differentiation",
                  "T cell activation",
                  "Lymphocyte mediated immunity",
                  "Leukocyte activation involved in inflammatory response",
                  "Positive regulation of immune effector process",
                  "Leukocyte cell cell adhesion"
)

feature_list[!(feature_list %in% Dataset$gs_name)]

dataset <- Dataset %>% filter(Dataset$gs_name %in% feature_list)
length(unique(dataset$gs_name))

geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name)

cells_AUC <- AUCell_calcAUC(geneSets,
                            cells_rankings,
                            aucMaxRank=nrow(cells_rankings)*0.1)

AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
colnames(AUCell_socre) <- colnames(AUCell_socre)|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

sce <- seurat_obj
sce$ID <- rownames(sce@meta.data)

meta <- AUCell_socre
meta <- as.data.frame(meta)
meta$sub_cell_type <- sce$cell_type_leiden0.6
meta$group1 <- sce$group1

object1 <- meta

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

data_melt <- melt(object1,
                  id.vars = c("sub_cell_type","group1"),
                  variable.name='Pathways',
                  value.name='Score')


colnames(data_melt)[1] <- "celltype"
data_melt <- as.data.frame(data_melt)

pathways <- feature_list

for (i in pathways) {
  data_melt1 <- subset(data_melt,subset = Pathways %in% i)
  
  data_melt2 <- subset(data_melt1,subset = (group1 %in% c("HF","HM")))
  data_melt2$Score <- as.numeric(data_melt2$Score)
  Data_summary <- summarySE(data_melt2, 
                            measurevar="Score", 
                            groupvars=c("celltype","group1"))
  p<-ggplot(data_melt2,aes(x=celltype,y=Score,fill=group1))+
    geom_split_violin()+
    theme_bw()+
    geom_point(data = Data_summary,aes(x= celltype, y= Score),pch=19,
               position=position_dodge(0.25),size= 1)+ #绘制均值为点图
    geom_errorbar(data = Data_summary,aes(ymin = Score-ci, ymax= Score+ci),
                  width= 0.05,
                  position= position_dodge(0.25),
                  color="black",
                  alpha = 0.8,
                  size= 0.5) + theme(axis.text.x = element_text(size = 8)) +
    scale_fill_manual(values = c(
      "#f38989","#815e99"
    ))+
    
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
          panel.grid = element_blank(),
          legend.position = "none")+
    stat_compare_means(aes(group = group1),
                       label = "p.signif",
                       method = "t.test",
                       label.y = max(data_melt2$Score),
                       hide.ns = T)+
    labs(title = i,
         x = NULL,
         y = NULL)
  p
  ggsave(p,filename= paste0("~/figure/Mac/05_AUCell/3_",i,"_H.pdf"), width = 3, height = 5)
}



#### Fig.5F ----
library(ComplexHeatmap)
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Macrophage")
seurat_obj <- NormalizeData(seurat_obj)
DefaultAssay(seurat_obj) <- "RNA"
genes <- c(
  "CCL3",
  "CCL4",
  "CXCL8",
  "OSM",
  "CCL3L1",
  "CSF1R",
  "CCL4L2",
  "NLRP3",
  "KDM6B",
  "HLA-DRB1",
  "CD74",
  "IFI6",
  "CXCL13",
  "CCL8",
  "CD86",
  "ADGRE2",
  "CSF2RA",
  "HIF1A",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DQA1",
  "HLA-DRA",
  "HLA-DRB1",
  "ITGAX",
  "PLSCR1",
  "POU2F2",
  "CCL19",
  "BTK",
  "CLNK",
  "CX3CR1",
  "GAB2",
  "JAK3",
  "NR4A3",
  "MALT1",
  "RASGRP1",
  "TNFSF18",
  "TNFSF4")

genes <- unique(genes)
seurat_obj <- subset(seurat_obj, features = genes)
expr <- AverageExpression(seurat_obj, group.by = "group1", assays = "RNA")[["RNA"]]
expr <- as.data.frame(expr)
dim(expr)
expr_z <- t(scale(t(expr)))
genes <- rownames(expr_z)

expr_z <- expr_z[, c("NF","NM","HF","HM")]

ann <- HeatmapAnnotation(
  Group = c(
    NF = "Normal",
    NM = "Normal",
    HF = "HCM",
    HM = "HCM"
  ),
  Sex = c(
    NF = "Female",
    NM = "Male",
    HF = "Female",
    HM = "Male"
  ),
  col = list(
    Group = c("Normal" = "#f9a341", "HCM" = "#af93c4"),
    Sex   = c("Female" = "#4dae47", "Male" = "#277fb8")
  )
)

pdf("~/03_figure/HCM/Mac/expr.pdf",width = 5,height = 7)
ComplexHeatmap::Heatmap(expr_z,
                        cluster_rows = F,
                        cluster_columns = F,
                        row_order = genes,
                        column_order = c("NF","NM","HF","HM"),
                        row_names_side = "left",
                        top_annotation = ann,
                        col = colorRampPalette(c("#67a4cc", "white", "#ec5051"))(100),)

dev.off()

#### Fig.5G ----
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type_leiden0.6 %in% "Macrophage")
seurat_obj <- NormalizeData(seurat_obj)
DefaultAssay(seurat_obj) <- "RNA"

gene <- "NLRP3"
expr_positive <- FetchData(seurat_obj, vars = gene)
cells_pos <- rownames(expr_positive)
meta <- seurat_obj@meta.data[cells_pos, , drop = FALSE]


expr_positive <- expr_positive[expr_positive[[gene]] > 0, , drop = FALSE]
cells_pos <- rownames(expr_positive)
meta_positive <- seurat_obj@meta.data[cells_pos, , drop = FALSE]

plot_df <- data.frame(expression = expr_positive[[gene]], 
                      group = meta_positive$group1)
plot_df$group <- factor(
  plot_df$group,
  levels = c("NF", "NM", "HF", "HM")
)

library(ggplot2)
library(ggpubr)
ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
  geom_violin(trim = FALSE,alpha=0.7) +
  geom_boxplot(width=0.1, alpha=1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(
      c("HF","HM"),
      c("NF","NM"),
      c("NF","HF"),
      c("NM","HM")
    )
  )+
  scale_fill_manual(
    values = c(
      NF = "#f9a341",
      NM = "#af93c4",
      HF = "#4dae47",
      HM = "#277fb8"
    )
  )+
  theme_classic() +
  labs(y = paste(gene, "expression"))
ggsave(paste0("~/03_figure/HCM/Mac/expr/",gene,"_expr.pdf"),width = 4.5,height = 4)

gene <- "HLA-DRB1"
expr_positive <- FetchData(seurat_obj, vars = gene)
cells_pos <- rownames(expr_positive)
meta <- seurat_obj@meta.data[cells_pos, , drop = FALSE]


expr_positive <- expr_positive[expr_positive[[gene]] > 0, , drop = FALSE]
cells_pos <- rownames(expr_positive)
meta_positive <- seurat_obj@meta.data[cells_pos, , drop = FALSE]

plot_df <- data.frame(expression = expr_positive[[gene]], 
                      group = meta_positive$group1)
plot_df$group <- factor(
  plot_df$group,
  levels = c("NF", "NM", "HF", "HM")
)

library(ggplot2)
library(ggpubr)
ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
  geom_violin(trim = FALSE,alpha=0.7) +
  geom_boxplot(width=0.1, alpha=1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(
      c("HF","HM"),
      c("NF","NM"),
      c("NF","HF"),
      c("NM","HM")
    )
  )+
  scale_fill_manual(
    values = c(
      NF = "#f9a341",
      NM = "#af93c4",
      HF = "#4dae47",
      HM = "#277fb8"
    )
  )+
  theme_classic() +
  labs(y = paste(gene, "expression"))
ggsave(paste0("~/03_figure/HCM/Mac/expr/",gene,"_expr.pdf"),width = 4.5,height = 4)


gene <- "HLA-DPB1"
expr_positive <- FetchData(seurat_obj, vars = gene)
cells_pos <- rownames(expr_positive)
meta <- seurat_obj@meta.data[cells_pos, , drop = FALSE]


expr_positive <- expr_positive[expr_positive[[gene]] > 0, , drop = FALSE]
cells_pos <- rownames(expr_positive)
meta_positive <- seurat_obj@meta.data[cells_pos, , drop = FALSE]

plot_df <- data.frame(expression = expr_positive[[gene]], 
                      group = meta_positive$group1)
plot_df$group <- factor(
  plot_df$group,
  levels = c("NF", "NM", "HF", "HM")
)

library(ggplot2)
library(ggpubr)
ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
  geom_violin(trim = FALSE,alpha=0.7) +
  geom_boxplot(width=0.1, alpha=1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(
      c("HF","HM"),
      c("NF","NM"),
      c("NF","HF"),
      c("NM","HM")
    )
  )+
  scale_fill_manual(
    values = c(
      NF = "#f9a341",
      NM = "#af93c4",
      HF = "#4dae47",
      HM = "#277fb8"
    )
  )+
  theme_classic() +
  labs(y = paste(gene, "expression"))
ggsave(paste0("~/03_figure/HCM/Mac/expr/",gene,"_expr.pdf"),width = 4.5,height = 4)

