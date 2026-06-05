#### Fig.1A ----
library(Seurat)
library(ggplot2)
library(patchwork)

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- 
  NormalizeData(seurat_obj,
                normalization.method = "LogNormalize")

seurat_obj$group1 <- ""
seurat_obj$group1[seurat_obj$sex == 'female' & seurat_obj$disease == 'HCM'] <- 'HF'  # HCM female
seurat_obj$group1[seurat_obj$sex == 'male' & seurat_obj$disease == 'HCM'] <- 'HM'  # HCM male
seurat_obj$group1[seurat_obj$sex == 'female' & seurat_obj$disease == 'NF'] <- 'NF'  # normal female
seurat_obj$group1[seurat_obj$sex == 'male' & seurat_obj$disease == 'NF'] <- 'NM'  # normal male

seurat_obj$cell_type <- ""
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Adipocyte')] <- 'Adipocyte'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Cardiomyocyte_I','Cardiomyocyte_II','Cardiomyocyte_III')] <- 'Cardiomyocyte'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Endocardial')] <- 'Endocardial'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Endothelial_I','Endothelial_II','Endothelial_III')] <- 'Endothelial'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Epicardial')] <- 'Epicardial'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Fibroblast_I','Fibroblast_II','Activated_fibroblast')] <- 'Fibroblast'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Lymphatic_endothelial')] <- 'Lymphatic_endothelial'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Lymphocyte')] <- 'Lymphocyte'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Macrophage','Proliferating_macrophage')] <- 'Macrophage'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Mast_cell')] <- 'Mast_cell'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Neuronal')] <- 'Neuronal'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('Pericyte_I','Pericyte_II')] <- 'Pericyte'
seurat_obj$cell_type[seurat_obj$cell_type_leiden0.6 %in% c('VSMC')] <- 'VSMC'

cell_color <- c(
  'Adipocyte' = "#a4cde1"
  ,'Cardiomyocyte' = "#277fb8"
  ,'Endocardial' = "#549da3"
  ,'Endothelial' = "#96cb8f"
  ,'Epicardial' = "#4dae47"
  ,'Fibroblast' = "#b79973"
  ,'Lymphatic_endothelial' = "#f38989"
  ,'Lymphocyte' = "#ec5051"
  ,'Macrophage' = "#f9b769"
  ,'Mast_cell' = "#d4a6a8"
  ,'Neuronal' = "#8660a8"
  ,'Pericyte' = "#af93c4"
  ,'VSMC' = "#f6f28f"
)

DimPlot(seurat_obj, 
        reduction = "umap", 
        group.by='cell_type',
        label = T,
        cols = cell_color,
        repel = TRUE,raster = FALSE)

ggsave(file.path("~/03_figure/HCM/大群/02_umap_celltype_data.pdf"),
       ggplot2::last_plot(),
       height= 6,
       width= 8)

#### Fig.1B ----
seurat_objHF <- subset(seurat_obj, group1 == "HF")
seurat_objHM <- subset(seurat_obj, group1 == "HM")
seurat_objNF <- subset(seurat_obj, group1 == "NF")
seurat_objNM <- subset(seurat_obj, group1 == "NM")

p1 <- DimPlot(seurat_objHF, 
              reduction = "umap", 
              group.by='cell_type',
              label = F,
              cols = cell_color,
              repel = TRUE,raster = FALSE)+
  labs(title = "HF")+
  theme(legend.position = "none")
p1

p2 <- DimPlot(seurat_objHM, 
              reduction = "umap", 
              group.by='cell_type',
              label = F,
              cols = cell_color,
              repel = TRUE,raster = FALSE)+
  labs(title = "HM")+
  theme(legend.position = "none")
p2

p3 <- DimPlot(seurat_objNF, 
              reduction = "umap", 
              group.by='cell_type',
              label = F,
              cols = cell_color,
              repel = TRUE,raster = FALSE)+
  labs(title = "NF")+
  theme(legend.position = "none")
p3

p4 <- DimPlot(seurat_objNM, 
              reduction = "umap", 
              group.by='cell_type',
              label = F,
              cols = cell_color,
              repel = TRUE,raster = FALSE)+
  labs(title = "NM")
p4

p1 | p2 | p3 | p4

ggsave(file.path("~/03_figure/HCM/大群/02_umap_group_data.pdf"),plot = last_plot(),height= 4,width= 18)


#### Fig.1C ----
library(Seurat)
library(ggplot2)
library(ggalluvial)

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
Idents(seurat_obj) <- "cell_type"

cellratio<-prop.table(table(seurat_obj$cell_type,seurat_obj$group1),margin=2)
cellratio<-as.data.frame(cellratio)
colnames(cellratio) <- c("cell_type","sex","ratio")

cellratio$cell_type <- factor(cellratio$cell_type,levels = unique(cellratio$cell_type))

cell_color <- c(
  'Adipocyte' = "#a4cde1"
  ,'Cardiomyocyte' = "#277fb8"
  ,'Endocardial' = "#549da3"
  ,'Endothelial' = "#96cb8f"
  ,'Epicardial' = "#4dae47"
  ,'Fibroblast' = "#b79973"
  ,'Lymphatic_endothelial' = "#f38989"
  ,'Lymphocyte' = "#ec5051"
  ,'Macrophage' = "#f9b769"
  ,'Mast_cell' = "#d4a6a8"
  ,'Neuronal' = "#8660a8"
  ,'Pericyte' = "#af93c4"
  ,'VSMC' = "#f6f28f"
)

p <- ggplot(cellratio, aes(x = sex, y = ratio, fill = cell_type)) +
  geom_bar(position = "fill", stat="identity", color = 'white', alpha = 5, width = 0.95) +
  scale_fill_manual(values = cell_color) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  coord_flip() 
p

pp <- ggplot(cellratio, aes(x = sex, y = ratio, fill = cell_type,
                            stratum = cell_type, alluvium = cell_type)) +
  scale_fill_manual(values = cell_color) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  coord_flip() 

p2 <- pp +
  geom_col(width = 0.6,
           color = 'white', size = 0.5) + #同时赋予白色描边，适当加粗
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0,
            color = 'white', size = 0.5)+ #同时赋予白色描边，适当加粗
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
p2

ggsave(file.path("~/03_figure/HCM/大群/03_柱状图.pdf"),
       plot =p2,
       height= 3,
       width= 15)



#### Fig.1D ----

library(Seurat)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)


## DEGs analysis
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")

## HCM (HF vs HM)

seurat_objHCM <- subset(seurat_obj, group == "HCM")

Idents(seurat_objHCM) <- "cell_type"

cell.types <- c("Cardiomyocyte","Adipocyte","Macrophage","Lymphatic_endothelial",
                "Fibroblast","VSMC","Endothelial","Endocardial","Pericyte","Neuronal",
                "Mast_cell","Lymphocyte","Epicardial")

for (cell_type in cell.types){
  seurat <- subset(seurat_objHCM,idents = cell_type)
  
  if (!("HF" %in% seurat$group1) | !("HM" %in% seurat$group1)) {
    message("Skipping ", cell_type, ": Missing 'HF' or 'HM' group")
    next
  }
  
  Idents(seurat) <- "group1"
  
  seurat <- NormalizeData(seurat)
  scRNA.markers <- FindMarkers(seurat,
                               ident.1 = "HF",
                               ident.2 = "HM")
  scRNA.markers$cell_type <- cell_type
  scRNA.markers$markers <- rownames(scRNA.markers)
  
  write.table(
    scRNA.markers,
    file = paste0("~/02_work_data/HCM_NF/ScRNA.marker/ScRNA.marker_H/",gsub(" ", "_",cell_type), ".txt"),
    quote = F,
    sep = ",",
    row.names = F
  )
  
  deg <- scRNA.markers
  head(deg)
  logFC.cutoff <- 0.25
  p_val_adj.cutoff <- 0.05
  
  deg$change <-
    as.factor(ifelse(
      deg$p_val_adj < p_val_adj.cutoff &
        abs(deg$avg_log2FC) > logFC.cutoff,
      ifelse(deg$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
      'NO'
    ))
  table(deg$change)
  
  write.table(
    deg,
    file = paste0("~/02_work_data/HCM_NF/findCYgene/deg_H/",gsub(" ", "_", cell_type), ".txt"),
    quote = F,
    sep = ",",
    row.names = F
  )
  
}


## NF (NF vs NM)
seurat_objNF <- subset(seurat_obj, group == "NF")

Idents(seurat_objNF) <- "cell_type"

cell.types <- c("Cardiomyocyte","Adipocyte","Macrophage","Lymphatic_endothelial",
                "Fibroblast","VSMC","Endothelial","Endocardial","Pericyte","Neuronal",
                "Mast_cell","Lymphocyte","Epicardial")

for (cell_type in cell.types){
  seurat <- subset(seurat_objNF,idents = cell_type)
  if (!("NF" %in% seurat$group1) | !("NM" %in% seurat$group1)) {
    message("Skipping ", cell_type, ": Missing 'NF' or 'NM' group")
    next
  }
  
  Idents(seurat) <- "group1"
  
  seurat <- NormalizeData(seurat)
  scRNA.markers <- FindMarkers(seurat,
                               ident.1 = "NF",
                               ident.2 = "NM")
  scRNA.markers$cell_type <- cell_type
  scRNA.markers$markers <- rownames(scRNA.markers)
  
  write.table(
    scRNA.markers,
    file = paste0("~/02_work_data/HCM_NF/ScRNA.marker/ScRNA.marker_N/",gsub(" ", "_",cell_type), ".txt"),
    quote = F,
    sep = ",",
    row.names = F
  )
  
  deg <- scRNA.markers
  head(deg)
  logFC.cutoff <- 0.25
  p_val_adj.cutoff <- 0.05
  
  deg$change <-
    as.factor(ifelse(
      deg$p_val_adj < p_val_adj.cutoff &
        abs(deg$avg_log2FC) > logFC.cutoff,
      ifelse(deg$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
      'NO'
    ))
  table(deg$change)
  
  write.table(
    deg,
    file = paste0("~/02_work_data/HCM_NF/findCYgene/deg_N/",gsub(" ", "_", cell_type), ".txt"),
    quote = F,
    sep = ",",
    row.names = F
  )
  
}


## female (HF vs NF)
seurat_objfemale <- subset(seurat_obj, sex == "female")

Idents(seurat_objfemale) <- "cell_type"

cell.types <- c("Cardiomyocyte","Adipocyte","Macrophage","Lymphatic_endothelial",
                "Fibroblast","VSMC","Endothelial","Endocardial","Pericyte","Neuronal",
                "Mast_cell","Lymphocyte","Epicardial")

for (cell_type in cell.types){
  seurat <- subset(seurat_objfemale,idents = cell_type)
  
  if (!("HF" %in% seurat$group1) | !("NF" %in% seurat$group1)) {
    message("Skipping ", cell_type, ": Missing 'HF' or 'NF' group")
    next
  }
  
  Idents(seurat) <- "group1"
  
  seurat <- NormalizeData(seurat)
  scRNA.markers <- FindMarkers(seurat,
                               ident.1 = "HF",
                               ident.2 = "NF")
  scRNA.markers$cell_type <- cell_type
  scRNA.markers$markers <- rownames(scRNA.markers)
  
  write.table(
    scRNA.markers,
    file = paste0("~/02_work_data/HCM_NF/ScRNA.marker/ScRNA.marker_F/",gsub(" ", "_",cell_type), ".txt"),
    quote = F,
    sep = ",",
    row.names = F
  )
  
  deg <- scRNA.markers
  head(deg)
  logFC.cutoff <- 0.25
  p_val_adj.cutoff <- 0.05
  
  deg$change <-
    as.factor(ifelse(
      deg$p_val_adj < p_val_adj.cutoff &
        abs(deg$avg_log2FC) > logFC.cutoff,
      ifelse(deg$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
      'NO'
    ))
  table(deg$change)
  
  write.table(
    deg,
    file = paste0("~/02_work_data/HCM_NF/findCYgene/deg_F/",gsub(" ", "_", cell_type), ".txt"),
    quote = F,
    sep = ",",
    row.names = F
  )
  
}

## male (HM vs NM)

seurat_objmale <- subset(seurat_obj, sex == "male")

Idents(seurat_objmale) <- "cell_type"

cell.types <- c("Cardiomyocyte","Adipocyte","Macrophage","Lymphatic_endothelial",
                "Fibroblast","VSMC","Endothelial","Endocardial","Pericyte","Neuronal",
                "Mast_cell","Lymphocyte","Epicardial")

for (cell_type in cell.types){
  seurat <- subset(seurat_objmale,idents = cell_type)
  
  if (!("HM" %in% seurat$group1) | !("NM" %in% seurat$group1)) {
    message("Skipping ", cell_type, ": Missing 'HM' or 'NM' group")
    next
  }
  
  Idents(seurat) <- "group1"
  
  seurat <- NormalizeData(seurat)
  scRNA.markers <- FindMarkers(seurat,
                               ident.1 = "HM",
                               ident.2 = "NM")
  scRNA.markers$cell_type <- cell_type
  scRNA.markers$markers <- rownames(scRNA.markers)
  
  write.table(
    scRNA.markers,
    file = paste0("~/02_work_data/HCM_NF/ScRNA.marker/ScRNA.marker_M/",gsub(" ", "_",cell_type), ".txt"),
    quote = F,
    sep = ",",
    row.names = F
  )
  
  deg <- scRNA.markers
  head(deg)
  logFC.cutoff <- 0.25
  p_val_adj.cutoff <- 0.05
  
  deg$change <-
    as.factor(ifelse(
      deg$p_val_adj < p_val_adj.cutoff &
        abs(deg$avg_log2FC) > logFC.cutoff,
      ifelse(deg$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
      'NO'
    ))
  table(deg$change)
  
  write.table(
    deg,
    file = paste0("~/02_work_data/HCM_NF/findCYgene/deg_M/",gsub(" ", "_", cell_type), ".txt"),
    quote = F,
    sep = ",",
    row.names = F
  )
  
}

## rbind
outdir <- "~/02_work_data/HCM_NF/findCYgene/deg_H/"
markers1 <- read_delim(paste0(outdir, "Adipocyte.txt"),show_col_types=FALSE)
markers1$cell_type <- "Adipocyte"
markers2 <- read_delim(paste0(outdir, "Cardiomyocyte.txt"),show_col_types=FALSE)
markers2$cell_type <- "Cardiomyocyte"
markers3 <- read_delim(paste0(outdir, "Endocardial.txt"),show_col_types=FALSE)
markers3$cell_type <- "Endocardial"
markers4 <- read_delim(paste0(outdir, "Endothelial.txt"),show_col_types=FALSE)
markers4$cell_type <- "Endothelial"
markers5 <- read_delim(paste0(outdir, "Fibroblast.txt"),show_col_types=FALSE)
markers5$cell_type <- "Fibroblast"
markers6 <- read_delim(paste0(outdir, "Lymphatic_endothelial.txt"),show_col_types=FALSE)
markers6$cell_type <- "Lymphatic_endothelial"
markers7 <- read_delim(paste0(outdir, "Lymphocyte.txt"),show_col_types=FALSE)
markers7$cell_type <- "Lymphocyte"
markers8 <- read_delim(paste0(outdir, "Macrophage.txt"),show_col_types=FALSE)
markers8$cell_type <- "Macrophage"
markers9 <- read_delim(paste0(outdir, "Mast_cell.txt"),show_col_types=FALSE)
markers9$cell_type <- "Mast_cell"
markers10 <- read_delim(paste0(outdir, "Neuronal.txt"),show_col_types=FALSE)
markers10$cell_type <- "Neuronal"
markers11 <- read_delim(paste0(outdir, "Pericyte.txt"),show_col_types=FALSE)
markers11$cell_type <- "Pericyte"
markers12 <- read_delim(paste0(outdir, "VSMC.txt"),show_col_types=FALSE)
markers12$cell_type <- "VSMC"
markers13 <- read_delim(paste0(outdir, "Epicardial.txt"),show_col_types=FALSE)
markers13$cell_type <- "Epicardial"
deg <- rbind(markers1,markers2,markers3,markers4,markers5,markers6,markers7,markers8,markers9,markers10,markers11,markers12,markers13)

write.csv(deg,"~/02_work_data/HCM_NF/findCYgene/deg_H/rbind_H.csv")


outdir <- "~/02_work_data/HCM_NF/findCYgene/deg_N/"

markers1 <- read_delim(paste0(outdir, "Adipocyte.txt"),show_col_types=FALSE)
markers1$cell_type <- "Adipocyte"
markers2 <- read_delim(paste0(outdir, "Cardiomyocyte.txt"),show_col_types=FALSE)
markers2$cell_type <- "Cardiomyocyte"
markers3 <- read_delim(paste0(outdir, "Endocardial.txt"),show_col_types=FALSE)
markers3$cell_type <- "Endocardial"
markers4 <- read_delim(paste0(outdir, "Endothelial.txt"),show_col_types=FALSE)
markers4$cell_type <- "Endothelial"
markers5 <- read_delim(paste0(outdir, "Fibroblast.txt"),show_col_types=FALSE)
markers5$cell_type <- "Fibroblast"
markers6 <- read_delim(paste0(outdir, "Lymphatic_endothelial.txt"),show_col_types=FALSE)
markers6$cell_type <- "Lymphatic_endothelial"
markers7 <- read_delim(paste0(outdir, "Lymphocyte.txt"),show_col_types=FALSE)
markers7$cell_type <- "Lymphocyte"
markers8 <- read_delim(paste0(outdir, "Macrophage.txt"),show_col_types=FALSE)
markers8$cell_type <- "Macrophage"
markers9 <- read_delim(paste0(outdir, "Mast_cell.txt"),show_col_types=FALSE)
markers9$cell_type <- "Mast_cell"
markers10 <- read_delim(paste0(outdir, "Neuronal.txt"),show_col_types=FALSE)
markers10$cell_type <- "Neuronal"
markers11 <- read_delim(paste0(outdir, "Pericyte.txt"),show_col_types=FALSE)
markers11$cell_type <- "Pericyte"
markers12 <- read_delim(paste0(outdir, "VSMC.txt"),show_col_types=FALSE)
markers12$cell_type <- "VSMC"

deg <- rbind(markers1,markers2,markers3,markers4,markers5,markers6,markers7,markers8,markers9,markers10,markers11,markers12)

write.csv(deg,"~/02_work_data/HCM_NF/findCYgene/deg_N/rbind_N.csv")


outdir <- "~/02_work_data/HCM_NF/findCYgene/deg_F/"

markers1 <- read_delim(paste0(outdir, "Adipocyte.txt"),show_col_types=FALSE)
markers1$cell_type <- "Adipocyte"
markers2 <- read_delim(paste0(outdir, "Cardiomyocyte.txt"),show_col_types=FALSE)
markers2$cell_type <- "Cardiomyocyte"
markers3 <- read_delim(paste0(outdir, "Endocardial.txt"),show_col_types=FALSE)
markers3$cell_type <- "Endocardial"
markers4 <- read_delim(paste0(outdir, "Endothelial.txt"),show_col_types=FALSE)
markers4$cell_type <- "Endothelial"
markers5 <- read_delim(paste0(outdir, "Fibroblast.txt"),show_col_types=FALSE)
markers5$cell_type <- "Fibroblast"
markers6 <- read_delim(paste0(outdir, "Lymphatic_endothelial.txt"),show_col_types=FALSE)
markers6$cell_type <- "Lymphatic_endothelial"
markers7 <- read_delim(paste0(outdir, "Lymphocyte.txt"),show_col_types=FALSE)
markers7$cell_type <- "Lymphocyte"
markers8 <- read_delim(paste0(outdir, "Macrophage.txt"),show_col_types=FALSE)
markers8$cell_type <- "Macrophage"
markers9 <- read_delim(paste0(outdir, "Mast_cell.txt"),show_col_types=FALSE)
markers9$cell_type <- "Mast_cell"
markers10 <- read_delim(paste0(outdir, "Neuronal.txt"),show_col_types=FALSE)
markers10$cell_type <- "Neuronal"
markers11 <- read_delim(paste0(outdir, "Pericyte.txt"),show_col_types=FALSE)
markers11$cell_type <- "Pericyte"
markers12 <- read_delim(paste0(outdir, "VSMC.txt"),show_col_types=FALSE)
markers12$cell_type <- "VSMC"

deg <- rbind(markers1,markers2,markers3,markers4,markers5,markers6,markers7,markers8,markers9,markers10,markers11,markers12)

write.csv(deg,"~/02_work_data/HCM_NF/findCYgene/deg_F/rbind_F.csv")


outdir <- "~/02_work_data/HCM_NF/findCYgene/deg_M/"
markers1 <- read_delim(paste0(outdir, "Adipocyte.txt"),show_col_types=FALSE)
markers1$cell_type <- "Adipocyte"
markers2 <- read_delim(paste0(outdir, "Cardiomyocyte.txt"),show_col_types=FALSE)
markers2$cell_type <- "Cardiomyocyte"
markers3 <- read_delim(paste0(outdir, "Endocardial.txt"),show_col_types=FALSE)
markers3$cell_type <- "Endocardial"
markers4 <- read_delim(paste0(outdir, "Endothelial.txt"),show_col_types=FALSE)
markers4$cell_type <- "Endothelial"
markers5 <- read_delim(paste0(outdir, "Fibroblast.txt"),show_col_types=FALSE)
markers5$cell_type <- "Fibroblast"
markers6 <- read_delim(paste0(outdir, "Lymphatic_endothelial.txt"),show_col_types=FALSE)
markers6$cell_type <- "Lymphatic_endothelial"
markers7 <- read_delim(paste0(outdir, "Lymphocyte.txt"),show_col_types=FALSE)
markers7$cell_type <- "Lymphocyte"
markers8 <- read_delim(paste0(outdir, "Macrophage.txt"),show_col_types=FALSE)
markers8$cell_type <- "Macrophage"
markers9 <- read_delim(paste0(outdir, "Mast_cell.txt"),show_col_types=FALSE)
markers9$cell_type <- "Mast_cell"
markers10 <- read_delim(paste0(outdir, "Neuronal.txt"),show_col_types=FALSE)
markers10$cell_type <- "Neuronal"
markers11 <- read_delim(paste0(outdir, "Pericyte.txt"),show_col_types=FALSE)
markers11$cell_type <- "Pericyte"
markers12 <- read_delim(paste0(outdir, "VSMC.txt"),show_col_types=FALSE)
markers12$cell_type <- "VSMC"
deg <- rbind(markers1,markers2,markers3,markers4,markers5,markers6,markers7,markers8,markers9,markers10,markers11,markers12)

write.csv(deg,"~/02_work_data/HCM_NF/findCYgene/deg_M/rbind_M.csv")


## visualisation

deg1 <- read.csv("~/02_work_data/HCM_NF/findCYgene/deg_H/rbind_H.csv")
deg1 <- deg1[,c(7,8,9)] 

deg12 <- subset(deg1,subset = deg1$change %in% "UP")
deg12$group <- "HF"
deg13 <- subset(deg1,subset = deg1$change %in% "DOWN")
deg13$group <- "HM"

deg2 <- read.csv("~/02_work_data/HCM_NF/findCYgene/deg_N/rbind_N.csv")
deg2 <- deg2[,c(7,8,9)]

deg22 <- subset(deg2,subset = deg2$change %in% "UP")
deg22$group <- "NF"
deg23 <- subset(deg2,subset = deg2$change %in% "DOWN")
deg23$group <- "NM"

deg <- rbind(deg12,deg13,deg22,deg23)
a <- deg %>% dplyr::count(change, cell_type,group,sort = TRUE)
deg.num <- a

deg.num$disease <- ""
deg.num$disease[which((deg.num$group == c("HF")) | (deg.num$group == c("HM")))]="HCM"
deg.num$disease[which((deg.num$group == c("NF")) | (deg.num$group == c("NM")))]="Normal"

deg.num$d <- ""
deg.num$d[which((deg.num$cell_type == c("Adipocyte")))]=1
deg.num$d[which((deg.num$cell_type == c("Cardiomyocyte")))]=2
deg.num$d[which((deg.num$cell_type == c("Endocardial")))]=3
deg.num$d[which((deg.num$cell_type == c("Endothelial")))]=4
deg.num$d[which((deg.num$cell_type == c("Fibroblast")))]=5
deg.num$d[which((deg.num$cell_type == c("Lymphocyte")))]=6
deg.num$d[which((deg.num$cell_type == c("Macrophage")))]=7
deg.num$d[which((deg.num$cell_type == c("Mast_cell")))]=8
deg.num$d[which((deg.num$cell_type == c("Neuronal")))]=9
deg.num$d[which((deg.num$cell_type == c("Pericyte")))]=10
deg.num$d[which((deg.num$cell_type == c("VSMC")))]=11
deg.num$d[which((deg.num$cell_type == c("Lymphatic_endothelial")))]=12

deg.num$d[which((deg.num$cell_type == c("Epicardial")))]=13

dd <- deg.num$d
dd <- as.numeric(dd)
deg.num <- cbind(deg.num, dd)
deg.num$n <- as.numeric(deg.num$n)

deg.num <- deg.num[deg.num$cell_type != 'Epicardial',]

disease_H <- dplyr::filter(deg.num, disease == "HCM") %>% 
  group_by(cell_type) %>% arrange(group) %>%
  mutate(pos = cumsum(n) - n / 2)

disease_N <- dplyr::filter(deg.num, disease == "Normal") %>% 
  group_by(cell_type) %>% arrange(group) %>%
  mutate(pos = cumsum(n) - n / 2)

Tmale <- paste0("Up-regulated DEGs in sex")
ct <- c("Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial",
        "Fibroblast", "Lymphocyte", "Macrophage", "Mast_cell",
        "Neuronal", "Pericyte", "VSMC","Lymphatic_endothelial")  # ,"Epicardial"

barwidth = 0.4
ggplot() +
  geom_bar(data = disease_H,
           mapping = aes(y=n,
                         x=dd,
                         color = group,
                         fill=group),
           stat="identity",
           position='stack',
           width = barwidth) +
  geom_bar(data = disease_N,
           mapping = aes(x = dd + barwidth + 0.05,
                         y = n,
                         color = group,
                         fill=group),
           stat="identity",
           position='stack',
           width = barwidth) +
  theme_bw()+
  labs(x="",y="",title=Tmale)+
  scale_fill_manual(values = c("#ca9e95",  "#92a7c4", "#ebd4cc", "#d4e4f2"))+
  scale_color_manual(values =c("#883745", "#395694", "#883745", "#395694"))+
  theme(axis.ticks.length = unit(0.5,'cm'),
        axis.text.x = element_text(angle=45,hjust = 1),
        plot.title = element_text(size=12, hjust = 0.5))+
  scale_x_continuous(breaks = c(1.25,2.25,3.25,4.25,5.25,6.25,7.25,8.25,9.25,10.25,11.25,12.25), 
                     labels = ct)

ggsave(file.path("~/03_figure/HCM/大群/05_差异基因堆叠柱状图_sex.pdf"),
       ggplot2::last_plot(),
       height= 4,
       width= 6)


#### Fig.1E ----

deg3 <- read.csv("~/02_work_data/HCM_NF/findCYgene/deg_F/rbind_F.csv")
deg3 <- deg3[,c(7,8,9)]

deg32 <- subset(deg3,subset = deg3$change %in% "UP")
deg32$group <- "HF"
deg33 <- subset(deg3,subset = deg3$change %in% "DOWN")
deg33$group <- "NF"

deg4 <- read.csv("~/02_work_data/HCM_NF/findCYgene/deg_M/rbind_M.csv")
deg4 <- deg4[,c(7,8,9)]

deg42 <- subset(deg4,subset = deg4$change %in% "UP")
deg42$group <- "HM"
deg43 <- subset(deg4,subset = deg4$change %in% "DOWN")
deg43$group <- "NM"


deg <- rbind(deg32,deg33,deg42,deg43)
a <- deg %>% dplyr::count(change, cell_type,group,sort = TRUE)
deg.num <- a

deg.num$sex <- ""
deg.num$sex[which((deg.num$group == c("HF")) | (deg.num$group == c("NF")))]="female"
deg.num$sex[which((deg.num$group == c("HM")) | (deg.num$group == c("NM")))]="male"

deg.num$d <- ""
deg.num$d[which((deg.num$cell_type == c("Adipocyte")))]=1
deg.num$d[which((deg.num$cell_type == c("Cardiomyocyte")))]=2
deg.num$d[which((deg.num$cell_type == c("Endocardial")))]=3
deg.num$d[which((deg.num$cell_type == c("Endothelial")))]=4
deg.num$d[which((deg.num$cell_type == c("Fibroblast")))]=5
deg.num$d[which((deg.num$cell_type == c("Lymphocyte")))]=6
deg.num$d[which((deg.num$cell_type == c("Macrophage")))]=7
deg.num$d[which((deg.num$cell_type == c("Mast_cell")))]=8
deg.num$d[which((deg.num$cell_type == c("Neuronal")))]=9
deg.num$d[which((deg.num$cell_type == c("Pericyte")))]=10
deg.num$d[which((deg.num$cell_type == c("VSMC")))]=11
deg.num$d[which((deg.num$cell_type == c("Lymphatic_endothelial")))]=12

dd <- deg.num$d
dd <- as.numeric(dd)
deg.num <- cbind(deg.num, dd)
deg.num$n <- as.numeric(deg.num$n)

disease_F <- dplyr::filter(deg.num, sex == "female") %>% 
  group_by(cell_type) %>% arrange(group) %>%
  mutate(pos = cumsum(n) - n / 2)

disease_M <- dplyr::filter(deg.num, sex == "male") %>% 
  group_by(cell_type) %>% arrange(group) %>%
  mutate(pos = cumsum(n) - n / 2)

Tmale <- paste0("Up-regulated DEGs in disease")
ct <- c("Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial",
        "Fibroblast", "Lymphocyte", "Macrophage", "Mast_cell",
        "Neuronal", "Pericyte", "VSMC","Lymphatic_endothelial")  # ,"Epicardial"

barwidth = 0.4
ggplot() +
  geom_bar(data = disease_F,
           mapping = aes(y=n,
                         x=dd,
                         color = group,
                         fill=group),
           stat="identity",
           position='stack',
           width = barwidth) +
  geom_bar(data = disease_M,
           mapping = aes(x = dd + barwidth + 0.05,
                         y = n,
                         color = group,
                         fill=group),
           stat="identity",
           position='stack',
           width = barwidth) +
  theme_bw()+
  labs(x="",y="",title=Tmale)+
  scale_fill_manual(values = c("#ca9e95",  "#92a7c4", "#ebd4cc", "#d4e4f2"))+
  scale_color_manual(values =c("#883745", "#395694", "#883745", "#395694"))+
  theme(axis.ticks.length = unit(0.5,'cm'),
        axis.text.x = element_text(angle=45,hjust = 1),
        plot.title = element_text(size=12, hjust = 0.5))+
  scale_x_continuous(breaks = c(1.25,2.25,3.25,4.25,5.25,6.25,7.25,8.25,9.25,10.25,11.25,12.25), 
                     labels = ct)

ggsave(file.path("~/03_figure/HCM/大群/05_差异基因堆叠柱状图_disease.pdf"),
       ggplot2::last_plot(),
       height= 4,
       width= 6)


#### Fig.1F ----
deg_H <- read.csv("~/02_work_data/HCM_NF/findCYgene/deg_H/rbind_H.csv")

gene_cell <- deg_H[deg_H$change == "UP",c("cell_type","markers") ]

Adipocyte <- gene_cell[gene_cell$cell_type == "Adipocyte",]
Cardiomyocyte <- gene_cell[gene_cell$cell_type == "Cardiomyocyte",]
Endothelial <- gene_cell[gene_cell$cell_type =="Endothelial",]
Endocardial <- gene_cell[gene_cell$cell_type =="Endocardial",]
Fibroblast <- gene_cell[gene_cell$cell_type =="Fibroblast",]
Lymphatic_endothelial <- gene_cell[gene_cell$cell_type =="Lymphatic_endothelial",]
Lymphocyte <- gene_cell[gene_cell$cell_type =="Lymphocyte",]
Macrophage <- gene_cell[gene_cell$cell_type =="Macrophage",]
Mast_cell <- gene_cell[gene_cell$cell_type =="Mast_cell",]
Neuronal <- gene_cell[gene_cell$cell_type =="Neuronal",]
Pericyte <- gene_cell[gene_cell$cell_type =="Pericyte",]
VSMC <- gene_cell[gene_cell$cell_type =="VSMC",]

H_UP <- cbind(Adipocyte$markers, Cardiomyocyte$markers, Endocardial$markers,
              Endothelial$markers, Fibroblast$markers, Lymphatic_endothelial$markers,
              Lymphocyte$markers, Macrophage$markers, Mast_cell$markers,
              Neuronal$markers, Pericyte$markers, VSMC$markers)

H_UP <- as.data.frame(H_UP)

colnames(H_UP) <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast", 
                    "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal",
                    "Pericyte","VSMC")

write.csv(H_UP, file = "~/02_work_data/HCM_NF/FC/H_UP.csv")

library(UpSetR)
datatry <-read.csv('~/02_work_data/HCM_NF/FC/H_UP.csv')
datatry <- datatry[,-1]

a <- fromList(datatry)
upset(fromList(datatry), order.by = "freq")
sets <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast", 
          "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal",
          "Pericyte","VSMC")

dev.off()

pdf("~/03_figure/HCM/大群/07_upset_HF.pdf", width=8.5, height = 6.5)
upset(fromList(datatry),
      order.by = "freq",
      nsets = 100,
      nintersects = 20,
      mb.ratio = c(0.7, 0.3),
      number.angles = 0,
      point.size = 1.8,
      line.size = 0.8,
      sets.x.label = "Set Size",
      main.bar.color = 'black',  #"#00B0D7"
      sets.bar.color = c("#a4cde1","#f6f28f","#96cb8f","#b79973",
                         "#f9b769","#549da3","#8660a8","#f38989",
                         "#d4a6a8","#af93c4","#277fb8","#ec5051"),
      matrix.color = "black",
      mainbar.y.label = "HF Intersection Sizes",
      text.scale = 1.25,
      shade.color = "black",
      shade.alpha = 0.1,
      
      queries = list(
        list(query=intersects, params=list("Adipocyte"), color="#a4cde1", active=T),
        list(query=intersects, params=list('Cardiomyocyte'), color="#277fb8", active=T),
        list(query=intersects, params=list('Endocardial'), color="#549da3", active=T),
        list(query=intersects, params=list('Endothelial'), color="#96cb8f", active=T),
        list(query=intersects, params=list('Fibroblast'), color="#b79973", active=T),
        list(query=intersects, params=list('Lymphatic_endothelial'), color="#f38989", active=T),
        list(query=intersects, params=list('Lymphocyte'), color="#ec5051", active=T),
        list(query=intersects, params=list('Macrophage'), color="#f9b769", active=T),
        list(query=intersects, params=list('Mast_cell'), color="#d4a6a8", active=T),
        list(query=intersects, params=list('Neuronal'), color="#8660a8", active=T),
        list(query=intersects, params=list('Pericyte'), color="#af93c4", active=T),
        list(query=intersects, params=list('VSMC'), color="#f6f28f", active=T)
      )
)

dev.off()



#### Fig.1G ----
deg_H <- read.csv("~/02_work_data/HCM_NF/findCYgene/deg_H/rbind_H.csv")

gene_cell <- deg_H[deg_H$change == "DOWN",c("cell_type","markers") ]

Adipocyte <- gene_cell[gene_cell$cell_type == "Adipocyte",]
Cardiomyocyte <- gene_cell[gene_cell$cell_type == "Cardiomyocyte",]
Endothelial <- gene_cell[gene_cell$cell_type =="Endothelial",]
Endocardial <- gene_cell[gene_cell$cell_type =="Endocardial",]
Fibroblast <- gene_cell[gene_cell$cell_type =="Fibroblast",]
Lymphatic_endothelial <- gene_cell[gene_cell$cell_type =="Lymphatic_endothelial",]
Lymphocyte <- gene_cell[gene_cell$cell_type =="Lymphocyte",]
Macrophage <- gene_cell[gene_cell$cell_type =="Macrophage",]
Mast_cell <- gene_cell[gene_cell$cell_type =="Mast_cell",]
Neuronal <- gene_cell[gene_cell$cell_type =="Neuronal",]
Pericyte <- gene_cell[gene_cell$cell_type =="Pericyte",]
VSMC <- gene_cell[gene_cell$cell_type =="VSMC",]

H_DOWN <- cbind(Adipocyte$markers, Cardiomyocyte$markers, Endocardial$markers,
                Endothelial$markers, Fibroblast$markers, Lymphatic_endothelial$markers,
                Lymphocyte$markers, Macrophage$markers, Mast_cell$markers,
                Neuronal$markers, Pericyte$markers, VSMC$markers)

H_DOWN <- as.data.frame(H_DOWN)

colnames(H_DOWN) <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast", 
                      "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal",
                      "Pericyte","VSMC")

write.csv(H_DOWN, file = "~/02_work_data/HCM_NF/FC/H_DOWN.csv")


datatry <-read.csv('~/02_work_data/HCM_NF/FC/H_DOWN.csv')
datatry <- datatry[,-1]

library(UpSetR)
a <- fromList(datatry)
upset(fromList(datatry), order.by = "freq")
sets <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast", 
          "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal",
          "Pericyte","VSMC")

dev.off()
dev.new()
pdf("~/03_figure/HCM/大群/07_upset_HM.pdf", width=8.5, height = 6.5)
upset(fromList(datatry),
      order.by = "freq",
      nsets = 100,
      nintersects = 20,
      mb.ratio = c(0.7, 0.3),
      number.angles = 0,
      point.size = 1.8,
      line.size = 0.8,
      sets.x.label = "Set Size",
      main.bar.color = 'black',  #"#00B0D7"
      sets.bar.color = c("#f9b769","#a4cde1","#277fb8","#f6f28f",
                         "#b79973","#96cb8f","#af93c4","#549da3",
                         "#f38989","#ec5051","#d4a6a8","#8660a8"),
      matrix.color = "black",
      mainbar.y.label = "HM Intersection Sizes",
      text.scale = 1.25,
      shade.color = "black",
      shade.alpha = 0.1,
      
      queries = list(
        list(query=intersects, params=list("Adipocyte"), color="#a4cde1", active=T),
        list(query=intersects, params=list('Cardiomyocyte'), color="#277fb8", active=T),
        list(query=intersects, params=list('Endocardial'), color="#549da3", active=T),
        list(query=intersects, params=list('Endothelial'), color="#96cb8f", active=T),
        list(query=intersects, params=list('Fibroblast'), color="#b79973", active=T),
        list(query=intersects, params=list('Lymphatic_endothelial'), color="#f38989", active=T),
        list(query=intersects, params=list('Lymphocyte'), color="#ec5051", active=T),
        list(query=intersects, params=list('Macrophage'), color="#f9b769", active=T),
        list(query=intersects, params=list('Mast_cell'), color="#d4a6a8", active=T),
        list(query=intersects, params=list('Pericyte'), color="#af93c4", active=T),
        list(query=intersects, params=list('VSMC'), color="#f6f28f", active=T)
      )
)

dev.off()


#### Fig.1H ----
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")

# Adipocyte
library(clusterProfiler)
library(org.Hs.eg.db)
datatry <-read.csv('~/02_work_data/HCM_NF/FC/H_UP.csv')
datatry <- datatry[,-1]

adipocyte <- datatry$Adipocyte
other <- unique(unlist(datatry[, -1]))
unique_adipocyte <- setdiff(adipocyte, other)
bp <- enrichGO(unique_adipocyte,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = "BP", 
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, "~/02_work_data/HCM_NF/intersection/GO/unique_Adipocyte_F.rds")

### HM
datatry <-read.csv('~/02_work_data/HCM_NF/FC/H_DOWN.csv')
datatry <- datatry[,-1]

adipocyte <- datatry$Adipocyte
other <- unique(unlist(datatry[, -1]))
unique_adipocyte <- setdiff(adipocyte, other)
bp <- enrichGO(unique_adipocyte,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = "BP", 
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, "~/02_work_data/HCM_NF/intersection/GO/unique_Adipocyte_M.rds")

pathway_F <- readRDS("~/02_work_data/HCM_NF/intersection/GO/unique_Adipocyte_F.rds")
pathway_M <- readRDS("~/02_work_data/HCM_NF/intersection/GO/unique_Adipocyte_M.rds")
pathway_F_unique <- pathway_F[pathway_F$Description %in% setdiff(pathway_F$Description, pathway_M$Description), ]
pathway_M_unique <- pathway_M[pathway_M$Description %in% setdiff(pathway_M$Description,pathway_F$Description), ]
pathway_F_unique <- pathway_F_unique[pathway_F_unique$pvalue<0.05,]
pathway_M_unique <- pathway_M_unique[pathway_M_unique$pvalue<0.05,]


sel_path_F <- c("triglyceride biosynthetic process",
                "regulation of cholesterol storage",
                "fatty-acyl-CoA metabolic process",
                "regulation of fatty acid beta-oxidation",
                "fatty acid homeostasis",
                "response to glucagon")

sel_F <- pathway_F %>%
  filter(Description %in% sel_path_F)
sel_F$group <- "HF"

sel_path_M <- c("activation of Janus kinase activity",
                "nitric-oxide synthase biosynthetic process",
                "positive regulation of mitochondrial membrane permeability",
                "negative regulation of fat cell differentiation",
                "glucose import",
                "lipid translocation")
sel_M <- pathway_M %>%
  filter(Description %in% sel_path_M)
sel_M$group <- "HM"

haha2 <- rbind(sel_F,sel_M)
haha2 <- haha2 %>% arrange(group)
pathway <- haha2$Description
haha2$Description <- factor(haha2$Description, levels = rev(pathway))
haha2$group <- factor(haha2$group, levels = c("HF", "HM"))


pdf("~/03_figure/HCM/大群/10_veen_GO/unique_Adipocyte_GOBP.pdf", width = 6, height = 4.5) 
catdotplot(haha2,
           x = group,
           y = Description,
           color = -log(pvalue),
           size = Count,
           dot_scale = 8.5
) +
  scale_size(range = c(3, 10)) +   
  
  scale_color_gradient2(low = "white", mid = "#FFED99", high ="#C00224") +
  guides(
    color = guide_colorbar(order = 1),  
    size = guide_legend(order = 2)      
  )+
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5,
    colour = "black",
    size = 12,
    
  ),
  axis.text.y = element_text(
    colour = "black",
    size = 12
  ),
  plot.title = element_text(hjust = 5, vjust = 0, size = 12),
  legend.position = "right",
  panel.spacing.x = unit(0, "pt")
  )
dev.off()


# Cardiomyocyte
datatry <-read.csv('~/02_work_data/HCM_NF/FC/H_UP.csv')
datatry <- datatry[,-1]

adipocyte <- datatry$Cardiomyocyte
other <- unique(unlist(datatry[, -2]))
unique_adipocyte <- setdiff(adipocyte, other)
bp <- enrichGO(unique_adipocyte,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = "BP", 
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, "~/02_work_data/HCM_NF/intersection/GO/unique_Cardiomyocyte_F.rds")

datatry <-read.csv('~/02_work_data/HCM_NF/FC/H_DOWN.csv')
datatry <- datatry[,-1]

adipocyte <- datatry$Cardiomyocyte
other <- unique(unlist(datatry[, -2]))
unique_adipocyte <- setdiff(adipocyte, other)
bp <- enrichGO(unique_adipocyte,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = "BP", 
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, "~/02_work_data/HCM_NF/intersection/GO/unique_Cardiomyocyte_M.rds")



pathway_F <- readRDS("~/02_work_data/HCM_NF/intersection/GO/unique_Cardiomyocyte_F.rds")
pathway_M <- readRDS("~/02_work_data/HCM_NF/intersection/GO/unique_Cardiomyocyte_M.rds")
pathway_F_unique <- pathway_F[pathway_F$Description %in% setdiff(pathway_F$Description, pathway_M$Description), ]
pathway_M_unique <- pathway_M[pathway_M$Description %in% setdiff(pathway_M$Description,pathway_F$Description), ]
pathway_F_unique <- pathway_F_unique[pathway_F_unique$pvalue<0.05,]
pathway_M_unique <- pathway_M_unique[pathway_M_unique$pvalue<0.05,]

sel_path_F <- c("cardiac muscle cell development",
                "cardiac muscle cell differentiation",
                "cardiac myofibril assembly",
                "positive regulation of fatty acid beta-oxidation",
                "negative regulation of mitochondrial membrane potential",
                "sarcomere organization")

sel_F <- pathway_F %>%
  filter(Description %in% sel_path_F)
sel_F$group <- "HF"

sel_path_M <- c("cell-matrix adhesion",
                "cell-cell junction assembly",
                "canonical Wnt signaling pathway",
                "cell-cell adhesion mediated by cadherin",
                "glycerolipid biosynthetic process",
                "cardiac left ventricle morphogenesis"
)
sel_M <- pathway_M %>%
  filter(Description %in% sel_path_M)
sel_M$group <- "HM"

haha2 <- rbind(sel_F,sel_M)
haha2 <- haha2 %>% arrange(group)
pathway <- haha2$Description
haha2$Description <- factor(haha2$Description, levels = rev(pathway))
haha2$group <- factor(haha2$group, levels = c("HF", "HM"))

pdf("~/03_figure/HCM/大群/10_veen_GO/unique_Cardiomyocyte_GOBP.pdf", width = 6, height = 4.5) 
catdotplot(haha2,
           x = group,
           y = Description,
           color = -log(pvalue),
           size = Count,
           dot_scale = 8.5
) +
  scale_size(range = c(3, 10)) +   
  
  scale_color_gradient2(low = "white", mid = "#FFED99", high ="#C00224") +
  guides(
    color = guide_colorbar(order = 1),  
    size = guide_legend(order = 2)      
  )+
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5,
    colour = "black",
    size = 12,
    
  ),
  axis.text.y = element_text(
    colour = "black",
    size = 12
  ),
  plot.title = element_text(hjust = 5, vjust = 0, size = 12),
  legend.position = "right",
  panel.spacing.x = unit(0, "pt")
  )
dev.off()


# Macrophage
datatry <-read.csv('~/02_work_data/HCM_NF/FC/H_UP.csv')
datatry <- datatry[,-1]

adipocyte <- datatry$Macrophage
other <- unique(unlist(datatry[, -8]))
unique_adipocyte <- setdiff(adipocyte, other)
bp <- enrichGO(unique_adipocyte,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = "BP", 
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, "~/02_work_data/HCM_NF/intersection/GO/unique_Macrophage_F.rds")

datatry <-read.csv('~/02_work_data/HCM_NF/FC/H_DOWN.csv')
datatry <- datatry[,-1]

adipocyte <- datatry$Macrophage
other <- unique(unlist(datatry[, -8]))
unique_adipocyte <- setdiff(adipocyte, other)
bp <- enrichGO(unique_adipocyte,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = "BP", 
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, "~/02_work_data/HCM_NF/intersection/GO/unique_Macrophage_M.rds")

pathway_F <- readRDS("~/02_work_data/HCM_NF/intersection/GO/unique_Macrophage_F.rds")
pathway_M <- readRDS("~/02_work_data/HCM_NF/intersection/GO/unique_Macrophage_M.rds")
pathway_F_unique <- pathway_F[pathway_F$Description %in% setdiff(pathway_F$Description, pathway_M$Description), ]
pathway_M_unique <- pathway_M[pathway_M$Description %in% setdiff(pathway_M$Description,pathway_F$Description), ]
pathway_F_unique <- pathway_F_unique[pathway_F_unique$pvalue<0.05,]
pathway_M_unique <- pathway_M_unique[pathway_M_unique$pvalue<0.05,]

sel_path_F <- c("immunoglobulin mediated immune response",
                "regulation of leukocyte mediated immunity",
                "immunoglobulin production",
                "MHC class II protein complex assembly",
                "lymphocyte activation involved in immune response",
                "complement receptor mediated signaling pathway")

sel_F <- pathway_F %>%
  filter(Description %in% sel_path_F)
sel_F$group <- "HF"

sel_path_M <- c("regulation of cAMP-mediated signaling",
                "cell-substrate junction assembly",
                "calcium ion transmembrane transport",
                "Notch signaling pathway",
                "ATP synthesis coupled electron transport",
                "response to angiotensin")
sel_M <- pathway_M %>%
  filter(Description %in% sel_path_M)
sel_M$group <- "HM"

haha2 <- rbind(sel_F,sel_M)
haha2 <- haha2 %>% arrange(group)
pathway <- haha2$Description
haha2$Description <- factor(haha2$Description, levels = rev(pathway))
haha2$group <- factor(haha2$group, levels = c("HF", "HM"))

pdf("~/03_figure/HCM/大群/10_veen_GO/unique_Macrophage_GOBP.pdf", width = 6, height = 4.5) 
catdotplot(haha2,
           x = group,
           y = Description,
           color = -log(pvalue),
           size = Count,
           dot_scale = 8.5
) +
  scale_size(range = c(3, 10)) +   
  
  scale_color_gradient2(low = "white", mid = "#FFED99", high ="#C00224") +
  guides(
    color = guide_colorbar(order = 1),  
    size = guide_legend(order = 2)      
  )+
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5,
    colour = "black",
    size = 12,
    
  ),
  axis.text.y = element_text(
    colour = "black",
    size = 12
  ),
  plot.title = element_text(hjust = 5, vjust = 0, size = 12),
  legend.position = "right",
  panel.spacing.x = unit(0, "pt")
  )
dev.off()


# Fibroblast
datatry <-read.csv('~/02_work_data/HCM_NF/FC/H_UP.csv')
datatry <- datatry[,-1]

adipocyte <- datatry$Fibroblast
other <- unique(unlist(datatry[, -5]))
unique_adipocyte <- setdiff(adipocyte, other)
bp <- enrichGO(unique_adipocyte,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = "BP", 
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, "~/02_work_data/HCM_NF/intersection/GO/unique_Fibroblast_F.rds")

datatry <-read.csv('~/02_work_data/HCM_NF/FC/H_DOWN.csv')
datatry <- datatry[,-1]

adipocyte <- datatry$Fibroblast
other <- unique(unlist(datatry[, -5]))
unique_adipocyte <- setdiff(adipocyte, other)
bp <- enrichGO(unique_adipocyte,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = "BP", 
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, "~/02_work_data/HCM_NF/intersection/GO/unique_Fibroblast_M.rds")

pathway_F <- readRDS("~/02_work_data/HCM_NF/intersection/GO/unique_Fibroblast_F.rds")
pathway_M <- readRDS("~/02_work_data/HCM_NF/intersection/GO/unique_Fibroblast_M.rds")
pathway_F_unique <- pathway_F[pathway_F$Description %in% setdiff(pathway_F$Description, pathway_M$Description), ]
pathway_M_unique <- pathway_M[pathway_M$Description %in% setdiff(pathway_M$Description,pathway_F$Description), ]
pathway_F_unique <- pathway_F_unique[pathway_F_unique$pvalue<0.05,]
pathway_M_unique <- pathway_M_unique[pathway_M_unique$pvalue<0.05,]
sel_path_F <- c("insulin-like growth factor receptor signaling pathway",
                "positive regulation of response to wounding",
                "second-messenger-mediated signaling",
                "cGMP-mediated signaling",
                "calcineurin-NFAT signaling cascade",
                "response to calcium ion")

sel_F <- pathway_F %>%
  filter(Description %in% sel_path_F)
sel_F$group <- "HF"

sel_path_M <- c("collagen fibril organization",
                "mitochondrial calcium ion homeostasis",
                "negative regulation of insulin receptor signaling pathway",
                "Rab protein signal transduction",
                "regulation of ATP biosynthetic process",
                "negative regulation of nitric oxide metabolic process")
sel_M <- pathway_M %>%
  filter(Description %in% sel_path_M)
sel_M$group <- "HM"

haha2 <- rbind(sel_F,sel_M)
haha2 <- haha2 %>% arrange(group)
pathway <- haha2$Description
haha2$Description <- factor(haha2$Description, levels = rev(pathway))
haha2$group <- factor(haha2$group, levels = c("HF", "HM"))


pdf("~/03_figure/HCM/大群/10_veen_GO/unique_Fibroblast_GOBP.pdf", width = 6, height = 4.5) 
catdotplot(haha2,
           x = group,
           y = Description,
           color = -log(pvalue),
           size = Count,
           dot_scale = 8.5
) +
  scale_size(range = c(3, 10)) +   
  
  scale_color_gradient2(low = "white", mid = "#FFED99", high ="#C00224") +
  guides(
    color = guide_colorbar(order = 1),  
    size = guide_legend(order = 2)      
  )+
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5,
    colour = "black",
    size = 12,
  ),
  axis.text.y = element_text(
    colour = "black",
    size = 12
  ),
  plot.title = element_text(hjust = 5, vjust = 0, size = 12),
  # axis.text.x = element_text(face = "italic"),
  legend.position = "right",
  panel.spacing.x = unit(0, "pt")
  )
dev.off()



#### Fig.1I ----
library(msigdbr)
library(Seurat)
library(clusterProfiler)
library(ggplot2)

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")

# gsea_H
seurat_obj <- subset(seurat_obj, group == "HCM")
Idents(seurat_obj) <- "cell_type"

cell.types <- c("Adipocyte","Cardiomyocyte","Macrophage","Lymphatic_endothelial",
                "Fibroblast","VSMC","Endothelial","Endocardial",
                "Pericyte","Neuronal","Mast_cell","Lymphocyte")

for (cell.type in cell.types) {
  seurat <- subset(seurat_obj, cell_type == cell.type)
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
    file = paste0("~/02_work_data/HCM_NF/gsea.input/gsea_H/",cell.type,"_input.csv"),
    quote = F,
    sep = ",",
    row.names = F
  )
}
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
# gsea_N
seurat_obj <- subset(seurat_obj, group == "NF")
Idents(seurat_obj) <- "cell_type"

cell.types <- c("Adipocyte","Cardiomyocyte","Macrophage","Lymphatic_endothelial",
                "Fibroblast","VSMC","Endothelial","Endocardial",
                "Pericyte","Neuronal","Mast_cell","Lymphocyte")

for (cell.type in cell.types) {
  seurat <- subset(seurat_obj, cell_type == cell.type)
  Idents(seurat) <- "group1"
  seurat <- NormalizeData(seurat)
  
  gsea.input <- FindMarkers(
    seurat,
    ident.1 = "NF",
    ident.2 = "NM",
    min.pct = 0,
    logfc.threshold = 0
  )
  gsea.input$celltype <- cell.type
  gsea.input$gene <- rownames(gsea.input)
  
  write.table(
    gsea.input,
    file = paste0("~/02_work_data/HCM_NF/gsea.input/gsea_N/",cell.type,"_input.csv"),
    quote = F,
    sep = ",",
    row.names = F
  )
}
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")

# gsea_F
seurat_objfemale <- subset(seurat_obj, sex == "female")
Idents(seurat_objfemale) <- "cell_type"

cell.types <- c("Cardiomyocyte","Adipocyte","Macrophage","Lymphatic_endothelial",
                "Fibroblast","VSMC","Endothelial","Endocardial",
                "Pericyte","Neuronal","Mast_cell","Lymphocyte")

# cell.types <- "Adipocyte"
for (cell.type in cell.types) {
  seurat <- subset(seurat_objfemale, cell_type == cell.type)
  Idents(seurat) <- "group1"
  seurat <- NormalizeData(seurat)
  
  gsea.input <- FindMarkers(
    seurat,
    ident.1 = "HF",
    ident.2 = "NF",
    min.pct = 0,
    logfc.threshold = 0
  )
  gsea.input$celltype <- cell.type
  gsea.input$gene <- rownames(gsea.input)
  
  write.table(
    gsea.input,
    file = paste0("~/02_work_data/HCM_NF/gsea.input/gsea_F/",cell.type,"_input.csv"),
    quote = F,
    sep = ",",
    row.names = F
  )
}
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")

# gsea_M
seurat_obj <- subset(seurat_obj, sex == "male")
Idents(seurat_obj) <- "cell_type"

cell.types <- c("Adipocyte","Cardiomyocyte","Macrophage","Lymphatic_endothelial",
                "Fibroblast","VSMC","Endothelial","Endocardial",
                "Pericyte","Neuronal","Mast_cell","Lymphocyte")

for (cell.type in cell.types) {
  seurat <- subset(seurat_obj, cell_type == cell.type)
  Idents(seurat) <- "group1"
  seurat <- NormalizeData(seurat)
  
  gsea.input <- FindMarkers(
    seurat,
    ident.1 = "HM",
    ident.2 = "NM",
    min.pct = 0,
    logfc.threshold = 0
  )
  gsea.input$celltype <- cell.type
  gsea.input$gene <- rownames(gsea.input)
  
  write.table(
    gsea.input,
    file = paste0("~/02_work_data/HCM_NF/gsea.input/gsea_M/",cell.type,"_input.csv"),
    quote = F,
    sep = ",",
    row.names = F
  )
}


## GSEA analysis
category <- c("C5")
genesets <- msigdbr(species = "Homo sapiens", category = category)
genesets <- subset(genesets, select = c("gs_name", "gene_symbol"))

cell.types <- c("Cardiomyocyte","Macrophage","Lymphatic_endothelial",
                "Fibroblast","Endothelial","Endocardial",
                "Pericyte","Neuronal","Mast_cell","Lymphocyte")
for (i in cell.types){
  gsea.input <-read.csv(paste0("~/02_work_data/HCM_NF/gsea.input/gsea_F/",i,"_input.csv"), sep = ",", header = T)
  rownames(gsea.input) <- gsea.input$gene
  gsea.input <- gsea.input[order(gsea.input$avg_log2FC, decreasing = T),]
  genelist <- structure(gsea.input$avg_log2FC, names = rownames(gsea.input))
  
  write.table(genelist, quote = F, sep = "\t", col.names = F, file = paste0("~/02_work_data/HCM_NF/gsea_data/gsea_F",i,"_F.rnk"))
  genelist <- genelist[which(genelist!=0)]
  res <- GSEA(genelist, TERM2GENE = genesets, eps = 0 ,pvalue= 1 ,minGSSize = 1)
  
  saveRDS(res, paste0("~/02_work_data/HCM_NF/gsea_data/gsea_F/",i,"_F.rds"))
  write.csv(res, paste0("~/02_work_data/HCM_NF/gsea_data/gsea_F/",i,"_F.csv"), row.names = F)
}

for (i in cell.types){
  gsea.input <-read.csv(paste0("~/02_work_data/HCM_NF/gsea.input/gsea_M/",i,"_input.csv"), sep = ",", header = T)
  rownames(gsea.input) <- gsea.input$gene
  gsea.input <- gsea.input[order(gsea.input$avg_log2FC, decreasing = T),]
  genelist <- structure(gsea.input$avg_log2FC, names = rownames(gsea.input))
  
  write.table(genelist, quote = F, sep = "\t", col.names = F, file = paste0("~/02_work_data/HCM_NF/gsea_data/gsea_M/",i,"_M.rnk"))
  
  genelist <- genelist[which(genelist!=0)] 
  res <- GSEA(genelist, TERM2GENE = genesets, eps = 0 ,pvalue= 1 ,minGSSize = 1)
  
  saveRDS(res, paste0("~/02_work_data/HCM_NF/gsea_data/gsea_M/",i,"_M.rds"))
  
  
  write.csv(res, paste0("~/02_work_data/HCM_NF/gsea_data/gsea_M/",i,"_M.csv"), row.names = F)
}


for (i in cell.types){
  gsea.input <-read.csv(paste0("~/02_work_data/HCM_NF/gsea.input/gsea_H/",i,"_input.csv"), sep = ",", header = T)
  rownames(gsea.input) <- gsea.input$gene
  gsea.input <- gsea.input[order(gsea.input$avg_log2FC, decreasing = T),]
  genelist <- structure(gsea.input$avg_log2FC, names = rownames(gsea.input))
  
  write.table(genelist, quote = F, sep = "\t", col.names = F, file = paste0("~/02_work_data/HCM_NF/gsea_data/gsea_H/",i,"_H.rnk"))
  
  genelist <- genelist[which(genelist!=0)] 
  res <- GSEA(genelist, TERM2GENE = genesets, eps = 0 ,pvalue= 1 ,minGSSize = 1)
  
  saveRDS(res, paste0("~/02_work_data/HCM_NF/gsea_data/gsea_H/",i,"_H.rds"))
  
  
  write.csv(res, paste0("~/02_work_data/HCM_NF/gsea_data/gsea_H/",i,"_H.csv"), row.names = F)
}

for (i in cell.types){
  gsea.input <-read.csv(paste0("~/02_work_data/HCM_NF/gsea.input/gsea_N/",i,"_input.csv"), sep = ",", header = T)
  rownames(gsea.input) <- gsea.input$gene
  gsea.input <- gsea.input[order(gsea.input$avg_log2FC, decreasing = T),]
  genelist <- structure(gsea.input$avg_log2FC, names = rownames(gsea.input))
  
  write.table(genelist, quote = F, sep = "\t", col.names = F, file = paste0("~/02_work_data/HCM_NF/gsea_data/gsea_N/",i,"_N.rnk"))
  
  genelist <- genelist[which(genelist!=0)] 
  res <- GSEA(genelist, TERM2GENE = genesets, eps = 0 ,pvalue= 1 ,minGSSize = 1)
  
  saveRDS(res, paste0("~/02_work_data/HCM_NF/gsea_data/gsea_N/",i,"_N.rds"))
  
  
  write.csv(res, paste0("~/02_work_data/HCM_NF/gsea_data/gsea_N/",i,"_N.csv"), row.names = F)
}

# HF:HM
res1 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Adipocyte_H.csv")
res1$cell_type <- "Adipocyte"
res1$group <- "HF:HM"
res2 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Cardiomyocyte_H.csv")
res2$cell_type <- "Cardiomyocyte"
res2$group <- "HF:HM"
res3 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Endothelial_H.csv")
res3$cell_type <- "Endothelial"
res3$group <- "HF:HM"
res4 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Fibroblast_H.csv")
res4$cell_type <- "Fibroblast"
res4$group <- "HF:HM"
res5 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Lymphatic_endothelial_H.csv")
res5$cell_type <- "Lymphatic_endothelial"
res5$group <- "HF:HM"
res6 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Lymphocyte_H.csv")
res6$cell_type <- "Lymphocyte"
res6$group <- "HF:HM"
res7 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Macrophage_H.csv")
res7$cell_type <- "Macrophage"
res7$group <- "HF:HM"
res8 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Mast_cell_H.csv")
res8$cell_type <- "Mast_cell"
res8$group <- "HF:HM"
res9 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Endocardial_H.csv")
res9$cell_type <- "Endocardial"
res9$group <- "HF:HM"
res10 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Neuronal_H.csv")
res10$cell_type <- "Neuronal"
res10$group <- "HF:HM"
res11 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/Pericyte_H.csv")
res11$cell_type <- "Pericyte"
res11$group <- "HF:HM"
res12 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_H/VSMC_H.csv")
res12$cell_type <- "VSMC"
res12$group <- "HF:HM"
resH <- rbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12)
write.csv(resH,"~/02_work_data/HCM_NF/gsea_data/gsea_H/resH.csv",row.names = F)

# NF:NM
res1 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Adipocyte_N.csv")
res1$cell_type <- "Adipocyte"
res1$group <- "NF:NM"
res2 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Cardiomyocyte_N.csv")
res2$cell_type <- "Cardiomyocyte"
res2$group <- "NF:NM"
res3 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Endothelial_N.csv")
res3$cell_type <- "Endothelial"
res3$group <- "NF:NM"
res4 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Fibroblast_N.csv")
res4$cell_type <- "Fibroblast"
res4$group <- "NF:NM"
res5 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Lymphatic_endothelial_N.csv")
res5$cell_type <- "Lymphatic_endothelial"
res5$group <- "NF:NM"
res6 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Lymphocyte_N.csv")
res6$cell_type <- "Lymphocyte"
res6$group <- "NF:NM"
res7 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Macrophage_N.csv")
res7$cell_type <- "Macrophage"
res7$group <- "NF:NM"
res8 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Mast_cell_N.csv")
res8$cell_type <- "Mast_cell"
res8$group <- "NF:NM"
res9 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Endocardial_N.csv")
res9$cell_type <- "Endocardial"
res9$group <- "NF:NM"
res10 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Neuronal_N.csv")
res10$cell_type <- "Neuronal"
res10$group <- "NF:NM"
res11 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/Pericyte_N.csv")
res11$cell_type <- "Pericyte"
res11$group <- "NF:NM"
res12 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_N/VSMC_N.csv")
res12$cell_type <- "VSMC"
res12$group <- "NF:NM"
resN <- rbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12)
write.csv(resN,"~/02_work_data/HCM_NF/gsea_data/gsea_N/resN.csv",row.names = F)

# HF:NF
res1 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Adipocyte_F.csv")
res1$cell_type <- "Adipocyte"
res1$group <- "HF:NF"
res2 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Cardiomyocyte_F.csv")
res2$cell_type <- "Cardiomyocyte"
res2$group <- "HF:NF"
res3 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Endothelial_F.csv")
res3$cell_type <- "Endothelial"
res3$group <- "HF:NF"
res4 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Fibroblast_F.csv")
res4$cell_type <- "Fibroblast"
res4$group <- "HF:NF"
res5 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Lymphatic_endothelial_F.csv")
res5$cell_type <- "Lymphatic_endothelial"
res5$group <- "HF:NF"
res6 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Lymphocyte_F.csv")
res6$cell_type <- "Lymphocyte"
res6$group <- "HF:NF"
res7 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Macrophage_F.csv")
res7$cell_type <- "Macrophage"
res7$group <- "HF:NF"
res8 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Mast_cell_F.csv")
res8$cell_type <- "Mast_cell"
res8$group <- "HF:NF"
res9 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Endocardial_F.csv")
res9$cell_type <- "Endocardial"
res9$group <- "HF:NF"
res10 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Neuronal_F.csv")
res10$cell_type <- "Neuronal"
res10$group <- "HF:NF"
res11 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/Pericyte_F.csv")
res11$cell_type <- "Pericyte"
res11$group <- "HF:NF"
res12 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_F/VSMC_F.csv")
res12$cell_type <- "VSMC"
res12$group <- "HF:NF"
resF <- rbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12)
write.csv(resF,"~/02_work_data/HCM_NF/gsea_data/gsea_F/resF.csv",row.names = F)

# HM:NM
res1 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Adipocyte_M.csv")
res1$cell_type <- "Adipocyte"
res1$group <- "HM:NM"
res2 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Cardiomyocyte_M.csv")
res2$cell_type <- "Cardiomyocyte"
res2$group <- "HM:NM"
res3 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Endothelial_M.csv")
res3$cell_type <- "Endothelial"
res3$group <- "HM:NM"
res4 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Fibroblast_M.csv")
res4$cell_type <- "Fibroblast"
res4$group <- "HM:NM"
res5 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Lymphatic_endothelial_M.csv")
res5$cell_type <- "Lymphatic_endothelial"
res5$group <- "HM:NM"
res6 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Lymphocyte_M.csv")
res6$cell_type <- "Lymphocyte"
res6$group <- "HM:NM"
res7 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Macrophage_M.csv")
res7$cell_type <- "Macrophage"
res7$group <- "HM:NM"
res8 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Mast_cell_M.csv")
res8$cell_type <- "Mast_cell"
res8$group <- "HM:NM"
res9 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Endocardial_M.csv")
res9$cell_type <- "Endocardial"
res9$group <- "HM:NM"
res10 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Neuronal_M.csv")
res10$cell_type <- "Neuronal"
res10$group <- "HM:NM"
res11 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/Pericyte_M.csv")
res11$cell_type <- "Pericyte"
res11$group <- "HM:NM"
res12 <- read.csv("~/02_work_data/HCM_NF/gsea_data/gsea_M/VSMC_M.csv")
res12$cell_type <- "VSMC"
res12$group <- "HM:NM"
resM <- rbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12)
write.csv(resM,"~/02_work_data/HCM_NF/gsea_data/gsea_M/resM.csv",row.names = F)

resALL <- rbind(resH,resN,resF,resM)
write.csv(resALL,"~/02_work_data/HCM_NF/gsea_data/resALL.csv",row.names = F)

## visualisation
adata <- read_csv("~/02_work_data/HCM_NF/gsea_data/resALL.csv")
a <- as.data.frame(str_split_fixed(adata$ID  ,"_",2))
b <- cbind(a,adata)
a <- subset(b,b$V1 %in% 'GOBP')
a <- a[,c(-1,-3)]

a$group <- factor(a$group, levels = c("HF:HM", "NF:NM","HF:NF","HM:NM"))

a$V2 <- tolower(a$V2)

# metabolism
meta_cell <- c("Adipocyte","Cardiomyocyte","Fibroblast","VSMC","Endothelial","Endocardial")

meta_description <- c("atp_metabolic_process",
                      "atp_synthesis_coupled_electron_transport",
                      "regulation_of_triglyceride_biosynthetic_process",
                      "ribose_phosphate_metabolic_process",
                      "acetyl_coa_metabolic_process",
                      "glucose_catabolic_process",
                      "fatty_acyl_coa_metabolic_process",
                      "citrate_metabolic_process",
                      "cellular_respiration",
                      "respiratory_electron_transport_chain",
                      "oxidative_phosphorylation",
                      "respiratory_system_development"
)
adataQQ <- filter(a, V2 %in% meta_description & cell_type %in% meta_cell)
adataQQ <- adataQQ[adataQQ$pvalue<0.05,]

p1 <- adataQQ |>
  ggplot(aes(
    x = group,
    y = V2,
    color = NES,
    size = -log10(pvalue)
  )) +
  geom_point(shape = 18) +
  facet_grid(. ~ cell_type) + 
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing.x = unit(0, "pt")
  ) +
  coord_fixed() +
  guides(x = guide_axis(angle = 45)) + 
  scale_color_gradient2(low = "#0571b0", 
                        mid = "white", high = "#ca0020")
p1
ggsave( "~/03_figure/HCM/大群/11_GSEA/11_GSEA_metabolism.pdf", plot = p1, height = 6, width = 9)

# remodeling
remodel_cell <- c("Cardiomyocyte","Endothelial","Endocardial",
                  "Fibroblast","VSMC","Adipocyte",
                  "Pericyte")
remodel_description <- c(
  "calcium_ion_transport",
  "regulation_of_anatomical_structure_size",
  "cell_fate_commitment",
  "cardiac_muscle_cell_contraction",
  "cell_growth",
  "tissue_migration",
  "stem_cell_differentiation",
  "regulation_of_calcium_ion_transport",
  "regulation_of_actin_filament_based_process",
  "positive_regulation_of_calcium_ion_transport",
  "cell_substrate_junction_organization",
  "cell_junction_assembly",
  "fibroblast_activation",
  "fibroblast_growth_factor_receptor_signaling_pathway",
  "response_to_fibroblast_growth_factor"
)


adataQQ <- filter(a, V2 %in% remodel_description & cell_type %in% remodel_cell)
adataQQ <- adataQQ[adataQQ$pvalue<0.05,]

p1 <- adataQQ |>
  ggplot(aes(
    x = group,
    y = V2,
    color = NES,
    size = -log10(pvalue)
  )) +
  geom_point(shape = 18) +
  facet_grid(. ~ cell_type) + 
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing.x = unit(0, "pt")
  ) +
  coord_fixed() +
  guides(x = guide_axis(angle = 45)) + 
  scale_color_gradient2(low = "#0571b0", 
                        mid = "white", high = "#ca0020")
p1
ggsave( "~/03_figure/HCM/大群/11_GSEA/11_GSEA_remodel.pdf", plot = p1, height = 5, width = 10)


# inflammation
meta_cell <- c("Adipocyte","Fibroblast","Lymphocyte","Endocardial",
               "Lymphatic_endothelial","Macrophage","Mast_cell")
meta_description <- c(
  "immune_response_regulating_signaling_pathway",
  "lymphocyte_activation_involved_in_immune_response",
  "cell_activation_involved_in_immune_response",
  "leukocyte_mediated_immunity",
  "adaptive_immune_response",
  "t_cell_activation",
  "cytokine_mediated_signaling_pathway",
  "response_to_tumor_necrosis_factor",
  "activation_of_immune_response",
  "positive_regulation_of_interleukin_2_production",
  "response_to_interleukin_1"
)


adataQQ <- filter(a, V2 %in% meta_description & cell_type %in% meta_cell)
adataQQ <- adataQQ[adataQQ$pvalue<0.05,]

p1 <- adataQQ |>
  ggplot(aes(
    x = group,
    y = V2,
    color = NES,
    size = -log10(pvalue)
  )) +
  geom_point(shape = 18) +
  facet_grid(. ~ cell_type) + 
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing.x = unit(0, "pt")
  ) +
  coord_fixed() +
  guides(x = guide_axis(angle = 45)) + 
  scale_color_gradient2(low = "#0571b0", 
                        mid = "white", high = "#ca0020")
p1
ggsave( "~/03_figure/HCM/大群/11_GSEA/11_GSEA_immune.pdf", plot = p1, height = 6, width = 11)


