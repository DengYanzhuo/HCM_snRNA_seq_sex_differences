#### Fig.6A ----
library(CellChat)
library(patchwork)
library(circlize)
library(Seurat)
library(tidyverse)
options(future.globals.maxSize = 4000000000)
source("/home/gongfengcz/Retina/src/computing/custom_function.R")

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj, subset = cell_type %in% c("Cardiomyocyte","Macrophage","Adipocyte",
                                                           "Endocardial","Neuronal","Fibroblast",
                                                           "Endothelial","Lymphatic_endothelial",
                                                           "VSMC","Pericyte","Lymphocyte","Mast_cell"
))
Idents(seurat_obj) <- "cell_type"
seurat_obj <- NormalizeData(seurat_obj)

group1 <- "HF"
output.dir <-
  paste0("/home/dengyz/02_work_data/HCM_NF/cellchat/",group1, "/")
dir.create(output.dir, recursive = T)
seurat_obj <- seurat_obj[, seurat_obj[["group1"]] == group1]
names(table(Idents(seurat_obj)))
print(table(Idents(seurat_obj)))

expr <- seurat_obj@assays$RNA@data

data.input <- expr
dim(data.input)
data.input[1:4, 1:4]
meta <- as.data.frame(Idents(seurat_obj))
colnames(meta) <- "labels"

unique(meta$labels)
cellchat <-
  createCellChat(object = data.input,
                 meta = meta,
                 group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <-
  setIdent(cellchat, ident.use = "labels")
levels(cellchat@idents)
groupSize <-
  as.numeric(table(cellchat@idents))

CellChatDB <-
  CellChatDB.human
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <-
  subsetData(cellchat)
cellchat <-
  identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

saveRDS(cellchat, file = "~/02_work_data/HCM_NF/cellchat/HF/cellchat.rds")



seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj, subset = cell_type %in% c("Cardiomyocyte","Macrophage","Adipocyte",
                                                           "Endocardial","Neuronal","Fibroblast",
                                                           "Endothelial","Lymphatic_endothelial",
                                                           "VSMC","Pericyte","Lymphocyte","Mast_cell"
))
Idents(seurat_obj) <- "cell_type"
seurat_obj <- NormalizeData(seurat_obj)

group1 <- "HM"
output.dir <-
  paste0("/home/dengyz/02_work_data/HCM_NF/cellchat/",group1, "/")
dir.create(output.dir, recursive = T)
seurat_obj <- seurat_obj[, seurat_obj[["group1"]] == group1]
names(table(Idents(seurat_obj)))
print(table(Idents(seurat_obj)))

expr <- seurat_obj@assays$RNA@data

data.input <- expr
dim(data.input)
data.input[1:4, 1:4]
meta <- as.data.frame(Idents(seurat_obj))
colnames(meta) <- "labels"

unique(meta$labels)
cellchat <-
  createCellChat(object = data.input,
                 meta = meta,
                 group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <-
  setIdent(cellchat, ident.use = "labels")
levels(cellchat@idents)
groupSize <-
  as.numeric(table(cellchat@idents))

CellChatDB <-
  CellChatDB.human
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <-
  subsetData(cellchat)
cellchat <-
  identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

saveRDS(cellchat, file = "~/02_work_data/HCM_NF/cellchat/HM/cellchat.rds")

library(igraph, lib.loc = "/home/dengyz/R/x86_64-pc-linux-gnu-library/4.3")
library(CellChat)
library(patchwork)
library(Seurat)
library(cowplot)
library(tidyverse)
source("/home/dengyz/01_code/HCM/custom_function.R")

setwd("~/03_figure/HCM/cellchat/HF_HM/")

outdir <- "~/02_work_data/HCM_NF/cellchat/"

library(CellChat)
library(patchwork)
library(cowplot)
cellchat.control  <- readRDS("~/02_work_data/HCM_NF/cellchat/HM/cellchat.rds")
cellchat.case  <- readRDS("~/02_work_data/HCM_NF/cellchat/HF/cellchat.rds")

object.list <- list(HM = cellchat.control, 
                    HF = cellchat.case)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list),cell.prefix=TRUE)
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use =c("#67a4cc",'#f38989') )
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",color.use =c("#67a4cc",'#f38989'))
gg1 + gg2  
ggsave("1-配受体比较.pdf",width = 5,height = 3)



#### Fig.6B ----
cellchat.control  <- readRDS("~/02_work_data/HCM_NF/cellchat/HM/cellchat.rds")
cellchat.case  <- readRDS("~/02_work_data/HCM_NF/cellchat/HF/cellchat.rds")

object.list <- list(HM = cellchat.control, 
                    HF = cellchat.case)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list),cell.prefix=TRUE)
cellchat

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("4-网络图.pdf",height = 5,width = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count,
                   weight.scale = T,
                   label.edge= F,
                   edge.weight.max = weight.max[2],
                   edge.width.max = 5, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()


#### Fig.6C ----
cellchat.control  <- readRDS("~/02_work_data/HCM_NF/cellchat/HM/cellchat.rds")
cellchat.case  <- readRDS("~/02_work_data/HCM_NF/cellchat/HF/cellchat.rds")

object.list <- list(HM = cellchat.control, 
                    HF = cellchat.case)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list),cell.prefix=TRUE)
cellchat

pdf("2-相互作用强度.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()



#### Fig.6D ----
cellchat.control  <- readRDS("~/02_work_data/HCM_NF/cellchat/HM/cellchat.rds")
cellchat.case  <- readRDS("~/02_work_data/HCM_NF/cellchat/HF/cellchat.rds")

object.list <- list(HM = cellchat.control, 
                    HF = cellchat.case)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list),cell.prefix=TRUE)
cellchat

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
p <- gg1 + gg2
pdf("3-热图.pdf",p,width = 8,height = 4)
print(p)
dev.off()


#### Fig.6E ----
cellchat.control  <- readRDS("~/02_work_data/HCM_NF/cellchat/HM/cellchat.rds")
cellchat.case  <- readRDS("~/02_work_data/HCM_NF/cellchat/HF/cellchat.rds")

object.list <- list(HM = cellchat.control, 
                    HF = cellchat.case)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list),cell.prefix=TRUE)
cellchat

HF <- cellchat@netP$HF
HM <- cellchat@netP$HM

HF_signaling <- unique(HF$pathways)
HM_signaling <- unique(HM$pathways)

signalings <- union(HF_signaling,HM_signaling)


gg1 <- rankNet(cellchat,
               mode = "comparison", stacked = T,
               color.use = c("#67a4cc",'#f38989'),
               signaling = signalings,
               do.stat = F) +coord_flip()
gg2 <- rankNet(cellchat,
               mode = "comparison", stacked = F,
               color.use = c("#67a4cc",'#f38989'),
               signaling = signalings,
               do.stat = TRUE) +coord_flip()

gg1+gg2
ggsave("5-堆叠柱形图.pdf", plot = gg1+gg2, width = 10,height = 10)


#### Fig.6F ----
cellchat.control  <- readRDS("~/02_work_data/HCM_NF/cellchat/HM/cellchat.rds")
cellchat.case  <- readRDS("~/02_work_data/HCM_NF/cellchat/HF/cellchat.rds")

object.list <- list(HM = cellchat.control, 
                    HF = cellchat.case)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list),cell.prefix=TRUE)
cellchat

df.net <- subsetCommunication(cellchat,
                              signaling = c("LEP","CD99","PTPRM","CADM","NT","IGF","BMP","CD45",
                                            "NECTIN","APP","SEMA3","SEMA6","ADIPONECTIN","LAMININ",
                                            "COLLAGEN","FGF","TENASCIN","THBS","FN1","NOTCH",
                                            "EPHB","PROS","GAS","PTN","KIT"))
HF <- df.net$HF

table(HF$source)
table(HF$target)

HM <- df.net$HM

table(HM$source)
table(HM$target)

write.csv(HF,"~/02_work_data/HCM_NF/cellchat/HF_net.csv")
write.csv(HM,"~/02_work_data/HCM_NF/cellchat/HM_net.csv")

HF <- read.csv("~/02_work_data/HCM_NF/cellchat/HF_net.csv")
HM <- read.csv("~/02_work_data/HCM_NF/cellchat/HM_net.csv")
HF <- HF[,-1]
HM <- HM[,-1]
HF <- HF[HF$source %in% c("Cardiomyocyte","Fibroblast","Adipocyte","Endothelial","Macrophage"),]
HM <- HM[HM$source %in% c("Cardiomyocyte","Fibroblast","Adipocyte","Endothelial","Macrophage"),]

colnames(HF)[5] <- "HF_prob"
colnames(HF)[6] <- "HF_pval"
colnames(HM)[5] <- "HM_prob"
colnames(HM)[6] <- "HM_pval"

commu_diff_HF <- HF %>% 
  dplyr::full_join(HM) %>% 
  dplyr::select(source, target, ligand, receptor,
                interaction_name, interaction_name_2, 
                pathway_name, annotation, evidence, 
                HF_prob, HF_pval, 
                HM_prob, HM_pval) %>%
  dplyr::mutate(
    HF_prob = tidyr::replace_na(HF_prob, 0),
    HM_prob = tidyr::replace_na(HM_prob, 0)) %>%
  dplyr::filter(HF_prob > HM_prob)


commu_diff_HM <- HM %>% 
  dplyr::full_join(HF) %>% 
  dplyr::select(source, target, ligand, receptor,
                interaction_name, interaction_name_2, 
                pathway_name, annotation, evidence, 
                HF_prob, HF_pval, 
                HM_prob, HM_pval) %>%
  dplyr::mutate(
    HF_prob = tidyr::replace_na(HF_prob, 0),
    HM_prob = tidyr::replace_na(HM_prob, 0)) %>%
  dplyr::filter(HF_prob < HM_prob)

commu_diff_HF$diff <- commu_diff_HF$HF_prob - commu_diff_HF$HM_prob
commu_diff_HM$diff <- commu_diff_HM$HM_prob - commu_diff_HM$HF_prob


## Cardiomyocyte
pairLR.use <- as.data.frame(c(
  "LAMA2_ITGA1_ITGB1",
  "LAMA4_ITGA1_ITGB1",
  "COL4A2_ITGA1_ITGB1",
  "COL4A5_ITGA1_ITGB1",
  "FGF1_FGFR1",
  "SEMA3C_PLXND1",
  "GAS6_MERTK"
))
colnames(pairLR.use) <- 'interaction_name'
a <- unique(commu_diff_HF$pathway_name[commu_diff_HF$interaction_name %in% pairLR.use$interaction_name])
a


p <- netVisual_bubble(cellchat, 
                      sources.use = c(1), 
                      targets.use = c(1:10), 
                      comparison = c(1, 2),
                      angle.x = 45,
                      pairLR.use = pairLR.use
)
p <- p + labs(title = "Cardiomyocyte")
p

ggsave("6-气泡图-Cardiomyocyte.pdf",height = 4.2,width = 7)


## Fibroblast
pairLR.use <- as.data.frame(c(
  "LAMA2_ITGA1_ITGB1",
  "COL6A3_ITGA1_ITGB1",
  "THBS2_CD36",
  "LAMA2_ITGA6_ITGB1",
  "LAMA2_ITGA9_ITGB1",
  "LAMC1_ITGA1_ITGB1",
  "COL4A2_ITGA1_ITGB1",
  "APP_CD74"
))
colnames(pairLR.use) <- 'interaction_name'

c <- unique(commu_diff_HF$pathway_name[commu_diff_HF$interaction_name %in% pairLR.use$interaction_name])
c

p <- netVisual_bubble(cellchat, 
                      sources.use = 6 ,
                      targets.use = c(1:3,5:10),
                      comparison = c(1, 2),
                      angle.x = 45,
                      pairLR.use = pairLR.use
)
p <- p + labs(title = "Fibroblast")
p

ggsave("6-气泡图-Fibroblast.pdf",height = 4,width = 6.5)



#### Fig.6G ----
cellchat.control  <- readRDS("~/02_work_data/HCM_NF/cellchat/HM/cellchat.rds")
cellchat.case  <- readRDS("~/02_work_data/HCM_NF/cellchat/HF/cellchat.rds")

object.list <- list(HM = cellchat.control, 
                    HF = cellchat.case)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list),cell.prefix=TRUE)
cellchat

pathways.show <- c("FGF")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf(paste0("./8-弦图/",pathways.show,".pdf"),height = 10,width = 20)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "chord",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()

pathways.show <- c("LAMININ")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf(paste0("./8-弦图/",pathways.show,".pdf"),height = 10,width = 20)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "chord",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()

