#### Fig.S5A ----
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
i = 1
cellchat1<-netAnalysis_computeCentrality(object.list[[i]],slot.name="netP") # HM
cellchat2<-netAnalysis_computeCentrality(object.list[[i+1]],slot.name="netP") # HF

pdf("10-netAnalysis_signalingRole_scatter.pdf", height = 5, width = 10)
gg1<-netAnalysis_signalingRole_scatter(cellchat1)
gg2<-netAnalysis_signalingRole_scatter(cellchat2)
gg1+gg2
dev.off()


#### Fig.S5B ----
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



## Adipocyte
pairLR.use <- as.data.frame(c(
  "COL4A1_ITGA1_ITGB1",
  "COL4A2_ITGA1_ITGB1",
  "LAMA4_ITGA1_ITGB1",
  "IGF1_IGF1R",
  "SEMA3C_PLXND1",
  "FGF1_FGFR1",
  "FGF2_FGFR1",
  "LAMA4_DAG1",
  "APP_CD74"
))
colnames(pairLR.use) <- 'interaction_name'
b <- unique(commu_diff_HF$pathway_name[commu_diff_HF$interaction_name %in% pairLR.use$interaction_name])
b

p <- netVisual_bubble(cellchat, 
                      sources.use = 3,
                      targets.use = c(1:10),
                      comparison = c(1, 2),
                      angle.x = 45,
                      pairLR.use = pairLR.use,
)
p <- p + labs(title = "Adipocyte")
p

ggsave("6-气泡图-Adipocyte.pdf",height = 4,width = 7)


## Endothelial
pairLR.use <- as.data.frame(c(
  "COL4A2_ITGA1_ITGB1",
  "COL4A1_ITGA1_ITGB1",
  "COL4A2_ITGA9_ITGB1",
  "LAMA4_ITGA9_ITGB1",
  "LAMA4_ITGA6_ITGB1",
  "LAMA4_DAG1",
  "APP_CD74"
))
colnames(pairLR.use) <- 'interaction_name'

d <- unique(commu_diff_HF$pathway_name[commu_diff_HF$interaction_name %in% pairLR.use$interaction_name])
d


p <- netVisual_bubble(cellchat, 
                      sources.use = 7 ,
                      targets.use = c(1:3,5:10),
                      comparison = c(1, 2),
                      angle.x = 45,
                      pairLR.use = pairLR.use
)
p <- p + labs(title = "Endothelial")
p

ggsave("6-气泡图-Endothelial.pdf",height = 4,width = 6.5)



#### Fig.S5C ----
cellchat.control  <- readRDS("~/02_work_data/HCM_NF/cellchat/HM/cellchat.rds")
cellchat.case  <- readRDS("~/02_work_data/HCM_NF/cellchat/HF/cellchat.rds")

object.list <- list(HM = cellchat.control, 
                    HF = cellchat.case)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list),cell.prefix=TRUE)
cellchat

pathways.show <- c("COLLAGEN")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf(paste0("./8-弦图/",pathways.show,".pdf"),height = 10,width = 20)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "chord",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()

pathways.show <- c("IGF")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf(paste0("./8-弦图/",pathways.show,".pdf"),height = 10,width = 20)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "chord",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()



#### Fig.S5D ----
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
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
par(mfrow = c(1,2), xpd=TRUE)  #输出为几行几列图形
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show,
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",
                                                  names(object.list)[i]))
}
p <- ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
pdf(paste0("./8-热图/",pathways.show,".pdf"),height = 4,width = 8)
print(p)
dev.off()

pathways.show <- c("COLLAGEN") 
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
par(mfrow = c(1,2), xpd=TRUE)  #输出为几行几列图形
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show,
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",
                                                  names(object.list)[i]))
}
p <- ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
pdf(paste0("./8-热图/",pathways.show,".pdf"),height = 4,width = 8)
print(p)
dev.off()

pathways.show <- c("LAMININ") 
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
par(mfrow = c(1,2), xpd=TRUE)  #输出为几行几列图形
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show,
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",
                                                  names(object.list)[i]))
}
p <- ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
pdf(paste0("./8-热图/",pathways.show,".pdf"),height = 4,width = 8)
print(p)
dev.off()

pathways.show <- c("IGF") 
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
par(mfrow = c(1,2), xpd=TRUE)  #输出为几行几列图形
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show,
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",
                                                  names(object.list)[i]))
}
p <- ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
pdf(paste0("./8-热图/",pathways.show,".pdf"),height = 4,width = 8)
print(p)
dev.off()