#### Fig.3A ----
library(ggalluvial)
library(tidyverse)
library(cols4all)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
subcell_color <- c(
  'Cardiomyocyte_I' = "#f9b769",
  'Cardiomyocyte_II' = "#af93c4",
  'Cardiomyocyte_III' = "#67a4cc"
)


seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Cardiomyocyte")
pdf("~/03_figure/HCM/Car/01_umap1.pdf",width = 6, height = 4)
DimPlot(seurat_obj, 
        reduction = "umap", 
        group.by='cell_type_leiden0.6',
        label = T,
        cols = subcell_color,
        repel = TRUE,raster = FALSE)
dev.off()


#### Fig.3B ----
library(Seurat)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("~/03_figure/HCM/Car/09_亚群鉴定/")


seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
Idents(seurat_obj) <- "cell_type"
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Cardiomyocyte")
seurat_obj <- NormalizeData(seurat_obj)

Idents(seurat_obj) <- "cell_type_leiden0.6"
all_mk <-FindAllMarkers(seurat_obj)
write.csv(all_mk,"~/02_work_data/HCM_NF/ScRNA.marker/Car_sub.csv")


top_mk <- c("PDK4","MYH6","LPL","RYR2", "USP13",
            # CarI
            "EBF1","MYL2","CRYAB","DCN","GSN",
            #  CarII 
            "NPPA","RTN4","MYL4","XIRP2","SLIT2"
            # CarIII 
)

desired_order <- c("Cardiomyocyte_I", "Cardiomyocyte_II", "Cardiomyocyte_III")

p <- DotPlot(seurat_obj,
             features = top_mk,
             scale.by ='size')+
  scale_size(range = c(1, 6))+
  theme( axis.text.x.bottom =element_text(angle =45,hjust =1,vjust =1))+
  scale_y_discrete(limits = desired_order)
p
ggsave("marker.pdf",plot = p, width = 9,height = 3.5)


#### Fig.3C ----
library(ggalluvial)
library(tidyverse)
library(cols4all)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
subcell_color <- c(
  'Cardiomyocyte_I' = "#f9b769",
  'Cardiomyocyte_II' = "#af93c4",
  'Cardiomyocyte_III' = "#67a4cc"
)
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
ggsave(file.path("~/03_figure/HCM/Car/01_umap_ratio/01_ratio2.pdf"),
       ggplot2::last_plot(),
       height= 2,
       width= 8)



#### Fig.3D ----
library(GSVA)
library(msigdbr)
library(dplyr)
library(stringr)
library(Seurat)

category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
Dataset <- Dataset %>% filter(str_detect(Dataset$gs_name,"GOBP_"))
Dataset$gs_name <- Dataset$gs_name|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
a <- as.data.frame(unique(Dataset$gs_name))

seurat.obj <- readRDS("~/data/seurat_obj_group1.rds")
seurat.obj <- subset(seurat.obj,subset = cell_type %in% "Cardiomyocyte")
seurat.obj <- NormalizeData(seurat.obj)
expr <- seurat.obj@assays$RNA@data

feature_list <- c("Tricarboxylic acid metabolic process",
                  "Pyruvate metabolic process",
                  "Electron transport chain",
                  "Cellular respiration",
                  "Lipid catabolic process",
                  "Regulation of fatty acid oxidation",
                  "Fatty acid beta oxidation",
                  "Lipid oxidation",
                  "Energy derivation by oxidation of organic compounds",
                  "Fatty acid catabolic process",
                  "Atp metabolic process",
                  "Glucose catabolic process",
                  "Oxidative phosphorylation",
                  "Fatty acyl coa catabolic process",
                  "Lipid homeostasis")

feature_list[!(feature_list %in% Dataset$gs_name)]

dataset <- Dataset %>% filter(Dataset$gs_name %in% feature_list)
length(unique(dataset$gs_name))

geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name)

param <- gsvaParam(exprData = expr,
                   geneSets = geneSets,
                   kcdf = "Gaussian")

gsva_mat <- gsva(param, verbose = TRUE)

write.csv(gsva_mat,"~/data/GSVA/gsva_go_matrix_car.csv")

library(limma)

gsva_mat <- read.csv("~/data/GSVA/gsva_go_matrix_car.csv")
rownames(gsva_mat) <- gsva_mat$X
gsva_mat <- gsva_mat[,-1]
colnames(gsva_mat) <- gsub("\\.", "-", colnames(gsva_mat))

seurat.obj <- readRDS("~/data/seurat_obj_group1.rds")
seurat.obj <- subset(seurat.obj,subset = cell_type %in% "Cardiomyocyte")
seurat.obj <- NormalizeData(seurat.obj)

metadata <- seurat.obj@meta.data
metadata$cell_id <- rownames(metadata)

## HCM:Normal
group_list <- metadata[,c("cell_id","disease")]
group_list <- group_list[match(colnames(gsva_mat), group_list$cell_id),]
colnames(group_list) <- c("cell_id","group")
design <- model.matrix(~0+factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_mat)

contrast.matrix <- makeContrasts(HCM-NF, levels = design)

fit <- lmFit(gsva_mat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(diff)

write.csv(diff,"~/data/GSVA/gsva_diff_car_HCM_NF.csv")

# vis
diff <- read.csv("~/data/GSVA/gsva_diff_car_HCM_NF.csv")
diff$group <- ifelse( diff$logFC > 0 & diff$adj.P.Val < 0.05 ,"up" ,
                      ifelse(diff$logFC < 0 & diff$P.Value < 0.05 ,"down","noSig"))
diff2 <- diff %>% arrange(t)

colnames(diff2)[1] <- "ID"
diff2$ID <- factor(diff2$ID, levels = diff2$ID)
limt = max(abs(diff2$t))

ggbarplot(diff2,x="ID",y="t",fill="group",
          palette = c("down" = "#008020",
                      "noSig" = "#cccccc",
                      "up" = "#08519C"),
          sort.val = "asc",
          rotate = T,
          ylab = "t value of GSVA score",
          xlab = "GO pathways")

ggsave("~/figure/Car/11_GSVA_bar_HCM.pdf",height = 6,width = 9)


#### Fig.3E ----
library(limma)

gsva_mat <- read.csv("~/data/GSVA/gsva_go_matrix_car.csv")
rownames(gsva_mat) <- gsva_mat$X
gsva_mat <- gsva_mat[,-1]
colnames(gsva_mat) <- gsub("\\.", "-", colnames(gsva_mat))

seurat.obj <- readRDS("~/data/seurat_obj_group1.rds")
seurat.obj <- subset(seurat.obj,subset = cell_type %in% "Cardiomyocyte")
seurat.obj <- NormalizeData(seurat.obj)

metadata <- seurat.obj@meta.data
metadata$cell_id <- rownames(metadata)

## HF:NF
group_list <- metadata[,c("cell_id","group1")]
group_list <- group_list[match(colnames(gsva_mat), group_list$cell_id),]
colnames(group_list) <- c("cell_id","group")
design <- model.matrix(~0+factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_mat)

contrast.matrix <- makeContrasts(HF-NF, levels = design)

fit <- lmFit(gsva_mat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(diff)

write.csv(diff,"~/data/GSVA/gsva_diff_car_HF.csv")

## HM:NM
contrast.matrix <- makeContrasts(HM-NM, levels = design)

fit <- lmFit(gsva_mat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(diff)

write.csv(diff,"~/data/GSVA/gsva_diff_car_HM.csv")


## vis
diff_HF <- read.csv("~/data/GSVA/gsva_diff_car_HF.csv")
diff_HM <- read.csv("~/data/GSVA/gsva_diff_car_HM.csv")

diff_HF <- diff_HF[,c("X","t","P.Value","adj.P.Val")]
diff_HM <- diff_HM[,c("X","t","P.Value","adj.P.Val")]
colnames(diff_HF) <- c("Pathway","HF_t","HF_P.Value","HF_adj.P.Val")
colnames(diff_HM) <- c("Pathway","HM_t","HM_P.Value","HM_adj.P.Val")
diff <- merge(diff_HF,diff_HM,by="Pathway")
diff <- diff[diff$Pathway!="Mitochondrial electron transport cytochrome c to oxygen",]


diff$Pathway <- factor(diff$Pathway,levels = c("Tricarboxylic acid metabolic process",
                                               "Pyruvate metabolic process",
                                               "Electron transport chain",
                                               "Cellular respiration",
                                               "Lipid catabolic process",
                                               "Regulation of fatty acid oxidation",
                                               "Fatty acid beta oxidation",
                                               "Lipid oxidation",
                                               "Energy derivation by oxidation of organic compounds",
                                               "Fatty acid catabolic process",
                                               "Atp metabolic process",
                                               "Glucose catabolic process",
                                               "Oxidative phosphorylation",
                                               "Fatty acyl coa catabolic process",
                                               "Lipid homeostasis"))
diff <- diff[order(diff$Pathway), ]
diff <- diff[order(-rank(diff$Pathway)), ]

diff2 <- diff[,c("Pathway","HF_t","HM_t")]
ann <- diff[,c("Pathway","HF_adj.P.Val","HM_adj.P.Val")]

rownames(diff2) <- diff2$Pathway
diff2 <- diff2[,-1]
rownames(ann) <- ann$Pathway
ann <- ann[,-1]
ann2 <- -log10(ann)
ann2 <- round(ann2, 2)
ann2 <- as.matrix(ann2)

col_fun <- colorRampPalette(c("#2c7bb6", "white"))

pdf("~/figure/Car/11_GSVA_heatmap_FM.pdf",height = 4.5,width = 5)
pheatmap(diff2,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = ann2,
         col = col_fun(100)
)

dev.off()


#### Fig.3F ----
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
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Cardiomyocyte")
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

feature_list <- c(
  "Pyruvate metabolic process",
  "Tricarboxylic acid metabolic process")

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

pathways <- c("Pyruvate metabolic process",
              "Tricarboxylic acid metabolic process")

for (i in pathways) {
  data_melt1 <- subset(data_melt,subset = Pathways %in% i)
  data_melt2 <- data_melt1
  data_melt2$Score <- as.numeric(data_melt2$Score)
  data_melt2$celltype <- "Cardiomyocyte"
  Data_summary <- summarySE(data_melt2, 
                            measurevar="Score", 
                            groupvars=c("celltype","group1"))
  data_melt2$group1 <- factor(
    data_melt2$group1,
    levels = c("NF", "NM", "HF", "HM")
  )
  Data_summary$group1 <- factor(
    Data_summary$group1,
    levels = c("NF", "NM", "HF", "HM")
  )
  p <- ggplot(data_melt2,aes(x=group1,y=Score,fill=group1))+
    geom_violin()+
    theme_bw()+
    geom_point(data = Data_summary,aes(x= group1, y= Score),pch=19,
               position=position_dodge(0.25),size= 1)+ #绘制均值为点图
    geom_errorbar(data = Data_summary,aes(ymin = Score-ci, ymax= Score+ci),
                  width= 0.05,
                  position= position_dodge(0.25),
                  color="black",
                  alpha = 0.8,
                  size= 0.5) + theme(axis.text.x = element_text(size = 8)) +
    scale_fill_manual(values = c(
      "#f9b769","#af93c4",
      "#4dae47","#277fb8"
    ))+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
          panel.grid = element_blank())+
    geom_signif(comparisons = list(c("HF","HM"),
                                   c("NF","NM"),
                                   c("HF","NF"),
                                   c("HM","NM")
    ),
    map_signif_level = TRUE,
    step_increase = 0.1,
    test = "t.test",
    textsize = 4)+
    labs(title = i)
  
  ggsave(p,filename= paste0("~/figure/Car/03_AUCell/",i,".pdf"), width = 5.5, height = 4.5)
}



#### Fig.3G ----
library(ComplexHeatmap)

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Cardiomyocyte")
seurat_obj <- NormalizeData(seurat_obj)
DefaultAssay(seurat_obj) <- "RNA"

genes <- c("CD36",
           "CPT1A",
           "ACADVL",
           "ACSL1",
           "HADHB",
           "HADHA",
           "ECH1",
           "SLC5A1",
           "HK1",
           "ENO1",
           "CS",
           "IDH2",
           "OGDH",
           "SLC25A4",
           "CKMT2",
           "ATP1B1",
           "ATP2A2",
           "CHCHD3",
           "VDAC1",
           "SPG7",
           "NDUFS3",
           "VDAC2",
           "UQCRC1",
           "LETM1",
           "DNAJC11",
           "MRPS36",
           "TOMM40",
           "BCL2L1",
           "UQCRH",
           "NDUFA6",
           "NDUFS8",
           "NDUFS6",
           "FOXRED1",
           "COX5A",
           "NDUFB9",
           "NDUFA12",
           "HADHB",
           "ACAA2",
           "NDUFA5",
           "ETFDH",
           "NDUFB1")

genes <- unique(genes)
expr_z <- expr_z[rownames(expr_z)%in% genes,]

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


pdf("~/03_figure/HCM/Car/expr.pdf",width = 5,height = 8)
ComplexHeatmap::Heatmap(expr_z,
                        cluster_rows = F,
                        cluster_columns = F,
                        row_order = genes,
                        column_order = c("NF","NM","HF","HM"),
                        row_names_side = "left",
                        top_annotation = ann,
                        col = colorRampPalette(c("#67a4cc", "white", "#ec5051"))(100),)

dev.off()



#### Fig.3H ----
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Cardiomyocyte")
seurat_obj <- NormalizeData(seurat_obj)
DefaultAssay(seurat_obj) <- "RNA"

gene <- "CPT1A"
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

ggsave(paste0("~/03_figure/HCM/Car/expr/",gene,"_expr.pdf"),width = 4.5,height = 4)


gene <- "ACADVL"
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

ggsave(paste0("~/03_figure/HCM/Car/expr/",gene,"_expr.pdf"),width = 4.5,height = 4)

gene <- "ACSL1"
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

ggsave(paste0("~/03_figure/HCM/Car/expr/",gene,"_expr.pdf"),width = 4.5,height = 4)


gene <- "SLC25A4"
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

ggsave(paste0("~/03_figure/HCM/Car/expr/",gene,"_expr.pdf"),width = 4.5,height = 4)


gene <- "HK1"
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

ggsave(paste0("~/03_figure/HCM/Car/expr/",gene,"_expr.pdf"),width = 4.5,height = 4)


gene <- "ENO1"
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

ggsave(paste0("~/03_figure/HCM/Car/expr/",gene,"_expr.pdf"),width = 4.5,height = 4)








