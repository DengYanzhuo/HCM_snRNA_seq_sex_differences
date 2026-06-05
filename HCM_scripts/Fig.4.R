#### Fig.4A ----
library(ggalluvial)
library(tidyverse)
library(cols4all)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
subcell_color <- c(
  'Fibroblast_I' = "#277fb8",
  'Fibroblast_II' = "#b05a28",
  'Activated_fibroblast' = "#f48521",
  'Fibroblast' = "#277fb8"
)

seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Fibroblast")
seurat_obj$cell_type_leiden0.6 <- as.character(seurat_obj$cell_type_leiden0.6)
seurat_obj$cell_type_leiden0.6[seurat_obj$cell_type_leiden0.6 %in% c("Fibroblast_I","Fibroblast_II")] <- "Fibroblast"
pdf("~/03_figure/HCM/Fib/01_umap.pdf",width = 6, height = 4)
DimPlot(seurat_obj, 
        reduction = "umap", 
        group.by='cell_type_leiden0.6',
        label = T,
        cols = subcell_color,
        repel = TRUE,raster = FALSE)
dev.off()


#### Fig.4B ----
library(Seurat)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(dplyr)
setwd("~/03_figure/HCM/Fib/02_亚群鉴定/")


seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
Idents(seurat_obj) <- "cell_type"
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Fibroblast")
seurat_obj <- NormalizeData(seurat_obj)

seurat_obj$cell_type_leiden0.6 <- as.character(seurat_obj$cell_type_leiden0.6)
seurat_obj$cell_type_leiden0.6[seurat_obj$cell_type_leiden0.6 %in% c("Fibroblast_I","Fibroblast_II")] <- "Fibroblast"

Idents(seurat_obj) <- "cell_type_leiden0.6"
all_mk <-FindAllMarkers(seurat_obj)
write.csv(all_mk,"~/02_work_data/HCM_NF/ScRNA.marker/Fib_sub-改.csv")
all_mk <- read.csv("~/02_work_data/HCM_NF/ScRNA.marker/Fib_sub-改.csv")
top_mk <- all_mk %>%
  group_by(cluster) %>%
  top_n(wt = avg_log2FC, n =5)%>%
  pull(gene)
top_mk <- c("POSTN","FAP","COL1A1","MEOX1","TGFB1",
            # Activated_fibroblast
            "DCN","GSN","APOD","ABCA8","COL15A1",
            # Fibroblast
            
)

desired_order <- c("Activated_fibroblast", "Fibroblast")

p <- DotPlot(seurat_obj,
             features = top_mk,
             scale.by ='size')+
  scale_size(range = c(1, 6))+
  theme( axis.text.x.bottom =element_text(angle =45,hjust =1,vjust =1))+
  scale_y_discrete(limits = desired_order)
p
ggsave("~/03_figure/HCM/Fib/02_亚群鉴定/marker-改.pdf",plot = p, width = 8,height = 4)




#### Fig.4C ----
library(ggalluvial)
library(tidyverse)
library(cols4all)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
subcell_color <- c(
  'Fibroblast_I' = "#277fb8",
  'Fibroblast_II' = "#b05a28",
  'Activated_fibroblast' = "#f48521",
  'Fibroblast' = "#277fb8"
)

seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Fibroblast")
seurat_obj$cell_type_leiden0.6 <- as.character(seurat_obj$cell_type_leiden0.6)
seurat_obj$cell_type_leiden0.6[seurat_obj$cell_type_leiden0.6 %in% c("Fibroblast_I","Fibroblast_II")] <- "Fibroblast"


cellratio<-prop.table(table(seurat_obj$cell_type_leiden0.6,seurat_obj$group1),margin=2) 
cellratio<-as.data.frame(cellratio)
colnames(cellratio) <- c("cell_type_leiden0.6","sex","ratio")
cellratio <- subset(cellratio,ratio != 0)

p <- ggplot(cellratio, aes(x = sex, y = ratio, fill = cell_type_leiden0.6)) +
  geom_bar(position = "fill", stat="identity", color = 'white', alpha = 5, width = 0.95) +
  scale_fill_manual(values = subcell_color) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  coord_flip() 
p

pp <- ggplot(cellratio, aes(x = sex, y = ratio, fill = cell_type_leiden0.6,
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
ggsave(file.path("~/03_figure/HCM/Fib/01_umap_ratio/01_ratio.pdf"),
       ggplot2::last_plot(),
       height= 2,
       width= 6)



#### Fig.4D ----
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
seurat.obj <- subset(seurat.obj,subset = cell_type %in% "Fibroblast")
seurat.obj <- subset(seurat.obj,subset = cell_type_leiden0.6 %in% "Activated_fibroblast")
seurat.obj <- NormalizeData(seurat.obj)
expr <- AverageExpression(seurat.obj, group.by = "group1", assays = "RNA")[["RNA"]]

feature_list <- c("Canonical wnt signaling pathway",
                  "Collagen biosynthetic process",
                  "Cgmp mediated signaling",
                  "Collagen metabolic process",
                  "Fibroblast activation",
                  "Fibroblast growth factor receptor signaling pathway",
                  "Cellular response to insulin like growth factor stimulus",
                  "Response to transforming growth factor beta",
                  "Extracellular matrix assembly",
                  "Cell matrix adhesion",
                  "Adherens junction assembly",
                  "Adherens junction organization",
                  "Cell substrate adhesion",
                  "Cell matrix adhesion",
                  "Cell substrate junction organization",
                  "Transforming growth factor beta receptor signaling pathway",
                  "Small gtpase mediated signal transduction",
                  "Notch signaling pathway",
                  "Collagen fibril organization",
                  "Platelet derived growth factor receptor signaling pathway",
                  "Response to platelet derived growth factor",
                  "Positive regulation of extracellular matrix assembly",
                  "Positive regulation of cell substrate junction organization",
                  "Regulation of cellular response to transforming growth factor beta stimulus"
)

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

write.csv(gsva_mat,"~/data/GSVA/gsva_go_matrix_fiba.csv")

library(limma)

gsva_mat <- read.csv("~/data/GSVA/gsva_go_matrix_fiba.csv")
rownames(gsva_mat) <- gsva_mat$X
gsva_mat <- gsva_mat[,-1]
colnames(gsva_mat) <- gsub("\\.", "-", colnames(gsva_mat))

seurat.obj <- readRDS("~/data/seurat_obj_group1.rds")
seurat.obj <- subset(seurat.obj,subset = cell_type %in% "Fibroblast")
seurat.obj <- subset(seurat.obj,subset = cell_type_leiden0.6 %in% "Activated_fibroblast")
seurat.obj <- NormalizeData(seurat.obj)

metadata <- seurat.obj@meta.data
metadata$cell_id <- rownames(metadata)
group_list <- metadata[,c("cell_id","group1")]
group_list <- group_list[match(colnames(gsva_mat), group_list$cell_id),]
colnames(group_list) <- c("cell_id","group")

design <- model.matrix(~0+factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_mat)

contrast.matrix <- makeContrasts(HF-HM, levels = design)

fit <- lmFit(gsva_mat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(diff)

write.csv(diff,"~/data/GSVA/gsva_diff_fiba.csv")

diff <- read.csv("~/data/GSVA/gsva_diff_fiba.csv")
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

ggsave("~/figure/Fib/11_GSVA_bar.pdf",height = 6,width = 9)


#### Fig.4E ----
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
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Fibroblast")
seurat_obj <- subset(seurat_obj,subset = cell_type_leiden0.6 %in% "Activated_fibroblast")
seurat_obj <- NormalizeData(seurat_obj)

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
  "Calcium ion transport",
  "Calcium mediated signaling",
  "Cellular response to calcium ion",
  "Canonical wnt signaling pathway",
  "Collagen biosynthetic process",
  "Cgmp mediated signaling",
  "Collagen metabolic process",
  "Fibroblast activation",
  "Fibroblast growth factor receptor signaling pathway",
  "Cellular response to insulin like growth factor stimulus",
  "Wound healing",
  "Response to transforming growth factor beta",
  "Extracellular matrix assembly",
  "Cellular respiration",
  "Cell chemotaxis",
  "Cell matrix adhesion",
  "Cellular response to vascular endothelial growth factor stimulus",
  "Wound healing",
  "Collagen metabolic process",
  "Cell junction assembly",
  "Cell substrate adhesion",
  "Platelet derived growth factor receptor signaling pathway",
  "Calcium mediated signaling",
  "Transforming growth factor beta receptor signaling pathway",
  "Collagen fibril organization",
  "Notch signaling pathway")
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
  p <- ggplot(data_melt2, aes(x = group1, y = Score, fill = group1)) +
    geom_violin() + 
    theme_bw()+
    geom_point(data = Data_summary,aes(x= group1, y= Score),pch=19,
               position=position_dodge(0.25),size= 1)+
    geom_errorbar(data = Data_summary,aes(ymin = Score-ci, ymax= Score+ci),
                  width= 0.05,
                  position= position_dodge(0.25),
                  color="black",
                  alpha = 0.8,
                  size= 0.5)+
    scale_fill_manual(values = c(
      "#f48521","#277fb8"
    ))+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
          panel.grid = element_blank())+
    geom_signif(comparisons = list(c("HF","HM")),
                map_signif_level = TRUE,
                test = "t.test",
                tip_length = 0.01,
                textsize = 4)+
    labs(title = i,
         y = NULL,
         x = "Group")
  
  ggsave(p,filename= paste0("~/figure/Fib/07_AUCell/",i,"_H.pdf"), width = 4, height = 5)
}



#### Fig.4F ----
library(ComplexHeatmap)
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Fibroblast")
seurat_obj <- NormalizeData(seurat_obj)
DefaultAssay(seurat_obj) <- "RNA"

genes <- c("EGR1",
           "PTN",
           "ITGBL1",
           "FGF14",
           "AEBP1",
           "FAM155A",
           "FAP",
           "SLC44A5",
           "TSHZ2",
           "COL1A2",
           "COL1A1",
           "JAZF1",
           "ACTA2",
           "RUNX1",
           "LEF1",
           "RUNX2",
           "GLIS1",
           "TEAD4",
           "GLI2",
           "RELB",
           "RELA",
           "TEAD1",
           "ANKRD1",
           "COL3A1",
           "FERMT2",
           "ITGA8",
           "ITGB1",
           "PDGFD",
           "TGFB1I1",
           "TGFBRAP1")

genes <- unique(genes)

seurat_obj <- subset(seurat_obj, features = genes)
expr <- AverageExpression(seurat_obj, group.by = "group1", assays = "RNA")[["RNA"]]
expr <- as.data.frame(expr)
dim(expr)
expr_z <- t(scale(t(expr)))

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

genes <- rownames(expr_z)
pdf("~/03_figure/HCM/Fib/expr_fib.pdf",width = 5,height = 6)
ComplexHeatmap::Heatmap(expr_z,
                        cluster_rows = F,
                        cluster_columns = F,
                        row_order = genes,
                        column_order = c("NF","NM","HF","HM"),
                        row_names_side = "left",
                        top_annotation = ann,
                        col = colorRampPalette(c("#67a4cc", "white", "#ec5051"))(100),)

dev.off()
