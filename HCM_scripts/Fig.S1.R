#### Fig.S1A ----
library("xlsx")
library(VennDiagram)
library(venn)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

deg_F <- read.csv("~/02_work_data/HCM_NF/findCYgene/deg_F/rbind_F.csv")
gene_cell <- deg_F[deg_F$change == "UP",c("cell_type","markers") ]

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

F_UP <- cbind(Adipocyte$markers, Cardiomyocyte$markers, Endocardial$markers,
              Endothelial$markers, Fibroblast$markers, Lymphatic_endothelial$markers,
              Lymphocyte$markers, Macrophage$markers, Mast_cell$markers,
              Neuronal$markers, Pericyte$markers, VSMC$markers)
colnames(F_UP) <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast", 
                    "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal",
                    "Pericyte","VSMC")

write.xlsx(F_UP, file = "~/02_work_data/HCM_NF/FC/F_UP.xlsx")

deg_M <- read.csv("~/02_work_data/HCM_NF/findCYgene/deg_M/rbind_M.csv")
gene_cell <- deg_M[deg_M$change == "UP",c("cell_type","markers") ]

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

M_UP <- cbind(Adipocyte$markers, Cardiomyocyte$markers, Endocardial$markers,
              Endothelial$markers, Fibroblast$markers, Lymphatic_endothelial$markers,
              Lymphocyte$markers, Macrophage$markers, Mast_cell$markers,
              Neuronal$markers, Pericyte$markers, VSMC$markers)

colnames(M_UP) <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast", 
                    "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal",
                    "Pericyte","VSMC")

write.xlsx(M_UP, file = "~/02_work_data/HCM_NF/FC/M_UP.xlsx")

gene_cell_F <- deg_F[deg_F$change == "UP",c("cell_type","markers") ]
gene_cell_M <- deg_M[deg_M$change == "UP",c("cell_type","markers") ]

cell_type <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast", 
               "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal",
               "Pericyte","VSMC")

for (i in cell_type) {
  cell_F <- gene_cell_F[gene_cell_F$cell_type == i,"markers"]
  cell_M <- gene_cell_M[gene_cell_M$cell_type == i,"markers"]
  cell_F_unique <- setdiff(cell_F, cell_M) # 找female组特有上调差异基因
  cell_M_unique <- setdiff(cell_M, cell_F) # 找male组特有上调差异基因
  write.xlsx(cell_F_unique,paste0("~/02_work_data/HCM_NF/intersection/",i,"_inter_F.xlsx"))
  write.xlsx(cell_M_unique,paste0("~/02_work_data/HCM_NF/intersection/",i,"_inter_M.xlsx"))
  
}

datatry1 <-read.xlsx('~/02_work_data/HCM_NF/FC/F_UP.xlsx',sheetIndex = 1)
datatry2 <-read.xlsx('~/02_work_data/HCM_NF/FC/M_UP.xlsx',sheetIndex = 1)

# Adipocyte
dataForVennDiagram=list(HF_UP=datatry1$Adipocyte ,
                        HM_UP=datatry2$Adipocyte)
pdf("~/03_figure/HCM/大群/10_veen_GO/veen_Adipocyte.pdf",width = 4, height = 4)
venn(dataForVennDiagram, 
     zcolor = c("#CB9494", "#8690A7"), 
     opacity = 0.25, 
     cexil = 0.5, 
     cexsn = 0.5, 
     ellipse = FALSE,
     borders = FALSE, 
     box = FALSE,
     sncs = 1.5,
     ilcs = 1.8,
)

title("Adipocyte", adj = 0,line = -3, cex.main = 1.5,font.main = 1)
dev.off()

# Cardiomyocyte
pdf("~/03_figure/HCM/大群/10_veen_GO/veen_Cardiomyocyte.pdf",width = 4, height = 4)
dataForVennDiagram=list(HF_UP=datatry1$Cardiomyocyte ,
                        HM_UP=datatry2$Cardiomyocyte)
venn(dataForVennDiagram, 
     zcolor = c("#CB9494", "#8690A7"), 
     opacity = 0.25, 
     cexil = 0.5, 
     cexsn = 0.5, 
     ellipse = FALSE,
     borders = FALSE, 
     box = FALSE,
     sncs = 1.5,
     ilcs = 1.8,
)

title("Cardiomyocyte", adj = 0,line = -3, cex.main = 1.5,font.main = 1)
dev.off()

# Fibroblast
pdf("~/03_figure/HCM/大群/10_veen_GO/veen_Fibroblast.pdf",width = 4, height = 4)
dataForVennDiagram=list(HF_UP=datatry1$Fibroblast ,
                        HM_UP=datatry2$Fibroblast)
venn(dataForVennDiagram, 
     zcolor = c("#CB9494", "#8690A7"), 
     opacity = 0.25, 
     cexil = 0.5, 
     cexsn = 0.5, 
     ellipse = FALSE,
     borders = FALSE, 
     box = FALSE,
     sncs = 1.5,
     ilcs = 1.8,
)

title("Fibroblast", adj = 0,line = -3, cex.main = 1.5,font.main = 1)
dev.off()

# Macrophage
pdf("~/03_figure/HCM/大群/10_veen_GO/veen_Macrophage.pdf",width = 4, height = 4)
dataForVennDiagram=list(HF_UP=datatry1$Macrophage ,
                        HM_UP=datatry2$Macrophage)
venn(dataForVennDiagram, 
     zcolor = c("#CB9494", "#8690A7"), 
     opacity = 0.25, 
     cexil = 0.5, 
     cexsn = 0.5, 
     ellipse = FALSE,
     borders = FALSE, 
     box = FALSE,
     sncs = 1.5,
     ilcs = 1.8,
)

title("Macrophage", adj = 0,line = -3, cex.main = 1.5,font.main = 1)
dev.off()


## GO
cell_type <- c("Adipocyte","Cardiomyocyte","Fibroblast","Macrophage")
library(readxl)

# cell_type <- "Adipocyte"
for (i in cell_type){
  mydata <- read_excel(paste0("~/02_work_data/HCM_NF/intersection/",i,"_inter_F.xlsx"))
  mydata <- mydata[,-1]
  colnames(mydata) <- "INTER"
  bp1 <- enrichGO(mydata$INTER,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  term <- bp1@result
  saveRDS(term, paste0("~/02_work_data/HCM_NF/intersection/GO/",i,"_inter_F.rds"))
}


for (i in cell_type){
  mydata <- read_excel(paste0("~/02_work_data/HCM_NF/intersection/",i,"_inter_M.xlsx"))
  mydata <- mydata[,-1]
  colnames(mydata) <- "INTER"
  bp1 <- enrichGO(mydata$INTER,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  term <- bp1@result
  saveRDS(term, paste0("~/02_work_data/HCM_NF/intersection/GO/",i,"_inter_M.rds"))
}


# Adipocyte
term <- readRDS("~/02_work_data/HCM_NF/intersection/GO/Adipocyte_inter_F.rds")
go <- c(
  "regulation of small GTPase mediated signal transduction",
  "cell-matrix adhesion",
  "extracellular matrix organization",
  "cGMP-mediated signaling",
  "response to transforming growth factor beta",
  "energy homeostasis"
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)

ggplot(data = term,
       aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=0.5, fill= "#CB9494",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=6,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave("~/03_figure/HCM/大群/10_veen_GO/GO_F_Adipocyte.pdf",width = 8, height = 3)

# M
term <- readRDS("~/02_work_data/HCM_NF/intersection/GO/Adipocyte_inter_M.rds")
go <- c(
  "lipid localization",
  "transforming growth factor beta1 production",
  "collagen fibril organization",
  "cell-substrate junction organization",
  "regulation of insulin receptor signaling pathway",
  "insulin-like growth factor receptor signaling pathway"
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
ggplot(data = term,
       aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=0.5, fill= "#8690A7",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=6,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave("~/03_figure/HCM/大群/10_veen_GO/GO_M_Adipocyte.pdf",width = 8, height = 3)



# Cardiomyocyte
term <- readRDS("~/02_work_data/HCM_NF/intersection/GO/Cardiomyocyte_inter_F.rds")
go <- c(
  "muscle system process",
  "cGMP metabolic process",
  "cardiac muscle contraction",
  "actomyosin structure organization",
  "sarcomere organization",
  "cell growth"
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
ggplot(data = term,
       aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=0.5, fill= "#CB9494",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=6,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave("~/03_figure/HCM/大群/10_veen_GO/GO_F_Cardiomyocyte.pdf",width = 8, height = 3)

term <- readRDS("~/02_work_data/HCM_NF/intersection/GO/Cardiomyocyte_inter_M.rds")
go <- c(
  "cardiac muscle tissue development",
  "cardiac ventricle morphogenesis",
  "cardiocyte differentiation",
  "ephrin receptor signaling pathway",
  "cell-substrate junction organization",
  "collagen-activated signaling pathway"
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
ggplot(data = term,
       aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=0.5, fill= "#8690A7",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=6,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave("~/03_figure/HCM/大群/10_veen_GO/GO_M_Cardiomyocyte.pdf",width = 8, height = 3)


# Fibroblast
term <- readRDS("~/02_work_data/HCM_NF/intersection/GO/Fibroblast_inter_F.rds")
go <- c(
  "axonogenesis",
  "small GTPase mediated signal transduction",
  "extracellular matrix organization",
  "regulation of cellular response to growth factor stimulus",
  "Wnt signaling pathway",
  "transmembrane receptor protein serine/threonine kinase signaling pathway"
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
ggplot(data = term,
       aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=0.5, fill= "#CB9494",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=6,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave("~/03_figure/HCM/大群/10_veen_GO/GO_F_Fibroblast.pdf",width = 8, height = 3)


term <- readRDS("~/02_work_data/HCM_NF/intersection/GO/Fibroblast_inter_M.rds")
go <- c(
  "positive regulation of supramolecular fiber organization",
  "regulation of actin filament organization",
  "positive regulation of integrin-mediated signaling pathway",
  "cell-substrate junction organization",
  "ERK1 and ERK2 cascade",
  "positive regulation of MAP kinase activity"
  
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
ggplot(data = term,
       aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=0.5, fill= "#8690A7",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=6,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave("~/03_figure/HCM/大群/10_veen_GO/GO_M_Fibroblast.pdf",width = 8, height = 3)

# Macrophage
term <- readRDS("~/02_work_data/HCM_NF/intersection/GO/Macrophage_inter_F.rds")
go <- c(
  "peptide antigen assembly with MHC class II protein complex",
  "regulation of T cell activation",
  "leukocyte cell-cell adhesion",
  "leukocyte mediated immunity",
  "cytokine-mediated signaling pathway",
  "activation of immune response"
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
ggplot(data = term,
       aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=0.5, fill= "#CB9494",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=6,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave("~/03_figure/HCM/大群/10_veen_GO/GO_F_Macrophage.pdf",width = 8, height = 3)

term <- readRDS("~/02_work_data/HCM_NF/intersection/GO/Macrophage_inter_M.rds")
go <- c(
  "positive regulation of lymphocyte differentiation",
  "immune response-regulating cell surface receptor signaling pathway",
  "cellular response to insulin stimulus",
  "immune response-regulating signaling pathway",
  "regulation of Fc receptor mediated stimulatory signaling pathway",
  "cellular response to interleukin-6"
)
term <- subset(term,term$Description %in% go)
term$labelx=rep(0,nrow(term))
term$labely=seq(nrow(term),1)
ggplot(data = term,
       aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=0.5, fill= "#8690A7",width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = term$Description),size=6,hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave("~/03_figure/HCM/大群/10_veen_GO/GO_M_Macrophage.pdf",width = 8, height = 3)


#### Fig.S1B ----
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = group %in% "HCM")
seurat_obj <- NormalizeData(seurat_obj)

cell.types <- c("Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial", 
                "Fibroblast", "Lymphatic_endothelial", "Lymphocyte", 
                "Macrophage", "Mast_cell", "Neuronal", "Pericyte", "VSMC")
for (i in cell.types) {
  sub_seurat_obj <- subset(seurat_obj, cell_type == i)
  
  countexp_Seurat=sc.metabolism.Seurat(obj = sub_seurat_obj,
                                       method = "AUCell", 
                                       imputation =F, 
                                       #ncores = 2, 
                                       metabolism.type = "KEGG")
  
  saveRDS(countexp_Seurat, paste0("~/02_work_data/HCM_NF/scMetabolism/HCM/",i,"_Seurat_metabolism.rds"))
}

cell.types <- c("Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial", 
                "Fibroblast", "Lymphatic_endothelial", "Lymphocyte", 
                "Macrophage", "Mast_cell", "Neuronal", "Pericyte", "VSMC")

all_data <- list()

for (i in cell.types) {
  file_path <- paste0("~/02_work_data/HCM_NF/scMetabolism/HCM/", i, "_Seurat_metabolism.rds")
  obj <- readRDS(file_path)
  
  metabolism_data <- obj@assays$METABOLISM$score
  metabolism_data <- as.matrix(metabolism_data)
  
  meta_cells <- colnames(metabolism_data)
  new_meta_cells <- gsub("\\.1\\.(\\d+)$", "-1-\\1", meta_cells)
  colnames(metabolism_data) <- new_meta_cells
  
  obj@assays[["METABOLISM"]] <- NULL
  obj[["METABOLISM"]] <- CreateAssayObject(counts = metabolism_data)
  
  all_data[[i]] <- obj
  
}

merged_obj <- merge(all_data[[1]],all_data[[3]]) 
merged_obj1 <- merge(merged_obj,all_data[[4]])
merged_obj2 <- merge(merged_obj1,all_data[[5]])
merged_obj3 <- merge(merged_obj2,all_data[[6]])
merged_obj4 <- merge(merged_obj3,all_data[[7]])
merged_obj5 <- merge(merged_obj4,all_data[[8]])
merged_obj6 <- merge(merged_obj5,all_data[[9]])
merged_obj7 <- merge(merged_obj6,all_data[[10]])
merged_obj8 <- merge(merged_obj7,all_data[[11]])
merged_obj9 <- merge(merged_obj8,all_data[[12]])
merged_obj <- merge(merged_obj9,all_data[[2]])

metabolism_data_df <- as.data.frame(merged_obj@assays$METABOLISM@data)
merged_obj@assays[["METABOLISM"]] <- NULL
merged_obj@assays[["METABOLISM"]][["score"]] <- metabolism_data_df

saveRDS(merged_obj,"~/02_work_data/HCM_NF/scMetabolism/HCM/countexp.Seurat.rds")

countexp.Seurat <- readRDS("~/02_work_data/HCM_NF/scMetabolism/HCM/countexp.Seurat.rds")

meta <- c(2,10,15,20,54)
for (a in meta) {
  input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[a]
  metadata <- countexp.Seurat@meta.data
  metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
  
  input.parameter <- "group1"
  
  metadata[, input.parameter] <- as.character(metadata[, input.parameter])
  metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, ])
  gg_table <- c()
  for (i in 1:length(input.pathway)) {
    gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], 
                                      input.pathway[i], metabolism.matrix_sub[, i]))
  }
  
  gg_table_df <- as.data.frame(gg_table)
  gg_table_df$V3 <- as.numeric(gg_table_df$V3)
  gg_table_df$scale <- scale(gg_table_df$V3)
  
  gg_table_df$scaled_value <- as.numeric(gg_table_df$scale[, 1])
  mean_HM <- mean(gg_table_df[gg_table_df[,1] == "HM", 5], na.rm = TRUE)
  mean_HF <- mean(gg_table_df[gg_table_df[,1] == "HF", 5], na.rm = TRUE)
  
  comparison_text <- if(mean_HM > mean_HF) {
    "HM > HF"
  } else if(mean_HF > mean_HM) {
    "HF > HM"
  } else {
    "HM = HF"
  }
  ggplot(data = gg_table_df, aes(x = gg_table_df[, 1], y = gg_table_df[,5],
                                 fill = gg_table_df[, 1]))+ 
    geom_violin() + 
    scale_fill_manual(values = c('HM'='lightblue','HF'='pink'))+
    geom_boxplot(width=0.18,fill="white",alpha=0.8,outlier.shape = NA)+
    geom_signif(comparisons = list(c("HM", "HF")),
                map_signif_level = TRUE,
                test = "t.test")+
    ylab(input.pathway) + xlab(input.parameter) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    labs(fill = input.parameter) + 
    annotate("text", 
             x = 1.5,  # x轴中间位置（HM=1, HF=2）
             y = max(gg_table_df[,5], na.rm = TRUE) - 0.1,  # y轴最小值下方
             label = comparison_text,
             size = 4,
             color = "black") +NULL
  ggsave(paste0("~/03_figure/HCM/大群/13_metabolism_sign/all_vb_",input.pathway,".pdf"),plot = last_plot(),width = 4,height = 5)
}

#### Fig.S1C ----
countexp.Seurat <- readRDS("~/02_work_data/HCM_NF/scMetabolism/HCM/countexp.Seurat.rds")
countexp.Seurat$group2 <- paste0(countexp_Seurat$cell_type,"_",countexp_Seurat$sex)

group_color <- c(
  'Adipocyte_male' = "#a4cde1","Adipocyte_female" = "#7AB6D4"
  ,'Cardiomyocyte_male' = "#277fb8","Cardiomyocyte_female" = "#1F6595"
  ,'Endocardial_male' = "#549da3","Endocardial_female" = "#447F84"
  ,'Endothelial_male' = "#96cb8f","Endothelial_female" = "#5FB054"
  ,'Fibroblast_male' = "#b79973", "Fibroblast_female" = "#9D7B51"
  ,'Lymphatic_endothelial_male' = "#f38989", "Lymphatic_endothelial_female" = "#EF5F5F"
  ,'Lymphocyte_male' = "#ec5051", "Lymphocyte_female" = "#C81616"
  ,'Macrophage_male' = "#f9b769", "Macrophage_female" = "#F69826"
  ,'Mast_cell_male' = "#d4a6a8", "Mast_cell_female" = "#C48689"
  ,'Neuronal_male' = "#8660a8", "Neuronal_female" = "#6F4D8D"
  ,'Pericyte_male' = "#af93c4", "Pericyte_female" = "#956FB1"
  ,'VSMC_male' = "#f6f28f", "VSMC_female" = "#EEE426"
)

#input.pathway <- c("Glycolysis / Gluconeogenesis",1
#                   "Citrate cycle (TCA cycle)",2
#                   "Oxidative phosphorylation",15
#                   "Pyruvate metabolism",10
#                   "Fatty acid degradation",20
#                   "Glutathione metabolism",54
#                   "Steroid hormone biosynthesis",22
#                   "Arachidonic acid metabolism",29
#                   "Arginine and proline metabolism",42
#                   "Synthesis and degradation of ketone bodies"21
#)

a <- 10

input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[a]

phenotype="group2"
ncol <- 1
input.parameter <- phenotype
metadata <- countexp.Seurat@meta.data
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
metadata[, input.parameter] <- as.character(metadata[, input.parameter])
metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, 
])
gg_table <- c()
for (i in 1:length(input.pathway)) {
  gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], 
                                    input.pathway[i], metabolism.matrix_sub[, i]))
}
gg_table <- data.frame(gg_table)
gg_table[, 3] <- as.numeric(as.character(gg_table[, 3]))

cell.types <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast",
                "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal","Pericyte","VSMC")

plot_list <- list()

for (celltype in cell.types) {
  select_group <- c(paste0(celltype, "_male"), paste0(celltype, "_female"))
  sub_gg_table <- gg_table[gg_table[, 1] %in% select_group, ]
  
  plot_data <- data.frame(
    Group = sub_gg_table[, 1],
    Value = sub_gg_table[, 3]
  )
  
  p <- ggplot(data = plot_data, aes(x = Group, y = Value, fill = Group)) + 
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = group_color) +
    geom_signif(comparisons = list(select_group),
                map_signif_level = TRUE,
                test = "t.test") +
    ylab(NULL) + xlab(NULL) + 
    ggtitle(celltype) +
    theme_bw() + 
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  plot_list[[celltype]] <- p
}
combined_plot <- wrap_plots(plot_list, nrow = 1)
combined_plot
ggsave(paste0("~/03_figure/HCM/大群/13_metabolism_改/Pyruvate_metabolism.pdf"),width = 24,height = 6,plot = combined_plot)



a <- 1

input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[a]

phenotype="group2"
ncol <- 1
input.parameter <- phenotype
metadata <- countexp.Seurat@meta.data
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
metadata[, input.parameter] <- as.character(metadata[, input.parameter])
metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, 
])
gg_table <- c()
for (i in 1:length(input.pathway)) {
  gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], 
                                    input.pathway[i], metabolism.matrix_sub[, i]))
}
gg_table <- data.frame(gg_table)
gg_table[, 3] <- as.numeric(as.character(gg_table[, 3]))

cell.types <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast",
                "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal","Pericyte","VSMC")

plot_list <- list()

for (celltype in cell.types) {
  select_group <- c(paste0(celltype, "_male"), paste0(celltype, "_female"))
  sub_gg_table <- gg_table[gg_table[, 1] %in% select_group, ]
  
  plot_data <- data.frame(
    Group = sub_gg_table[, 1],
    Value = sub_gg_table[, 3]
  )
  
  p <- ggplot(data = plot_data, aes(x = Group, y = Value, fill = Group)) + 
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = group_color) +
    geom_signif(comparisons = list(select_group),
                map_signif_level = TRUE,
                test = "t.test") +
    ylab(NULL) + xlab(NULL) + 
    ggtitle(celltype) +
    theme_bw() + 
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  plot_list[[celltype]] <- p
}
combined_plot <- wrap_plots(plot_list, nrow = 1)
combined_plot
ggsave(paste0("~/03_figure/HCM/大群/13_metabolism_改/Glycolysis.pdf"),width = 24,height = 6,plot = combined_plot)


a <- 15

input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[a]

phenotype="group2"
ncol <- 1
input.parameter <- phenotype
metadata <- countexp.Seurat@meta.data
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
metadata[, input.parameter] <- as.character(metadata[, input.parameter])
metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, 
])
gg_table <- c()
for (i in 1:length(input.pathway)) {
  gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], 
                                    input.pathway[i], metabolism.matrix_sub[, i]))
}
gg_table <- data.frame(gg_table)
gg_table[, 3] <- as.numeric(as.character(gg_table[, 3]))

cell.types <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast",
                "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal","Pericyte","VSMC")

plot_list <- list()

for (celltype in cell.types) {
  select_group <- c(paste0(celltype, "_male"), paste0(celltype, "_female"))
  sub_gg_table <- gg_table[gg_table[, 1] %in% select_group, ]
  
  plot_data <- data.frame(
    Group = sub_gg_table[, 1],
    Value = sub_gg_table[, 3]
  )
  
  p <- ggplot(data = plot_data, aes(x = Group, y = Value, fill = Group)) + 
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = group_color) +
    geom_signif(comparisons = list(select_group),
                map_signif_level = TRUE,
                test = "t.test") +
    ylab(NULL) + xlab(NULL) + 
    ggtitle(celltype) +
    theme_bw() + 
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  plot_list[[celltype]] <- p
}
combined_plot <- wrap_plots(plot_list, nrow = 1)
combined_plot
ggsave(paste0("~/03_figure/HCM/大群/13_metabolism_改/Oxidative_phosphorylation.pdf"),width = 24,height = 6,plot = combined_plot)


a <- 2

input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[a]

phenotype="group2"
ncol <- 1
input.parameter <- phenotype
metadata <- countexp.Seurat@meta.data
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
metadata[, input.parameter] <- as.character(metadata[, input.parameter])
metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, 
])
gg_table <- c()
for (i in 1:length(input.pathway)) {
  gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], 
                                    input.pathway[i], metabolism.matrix_sub[, i]))
}
gg_table <- data.frame(gg_table)
gg_table[, 3] <- as.numeric(as.character(gg_table[, 3]))

cell.types <- c("Adipocyte","Cardiomyocyte","Endocardial","Endothelial","Fibroblast",
                "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell","Neuronal","Pericyte","VSMC")

plot_list <- list()

for (celltype in cell.types) {
  select_group <- c(paste0(celltype, "_male"), paste0(celltype, "_female"))
  sub_gg_table <- gg_table[gg_table[, 1] %in% select_group, ]
  
  plot_data <- data.frame(
    Group = sub_gg_table[, 1],
    Value = sub_gg_table[, 3]
  )
  
  p <- ggplot(data = plot_data, aes(x = Group, y = Value, fill = Group)) + 
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = group_color) +
    geom_signif(comparisons = list(select_group),
                map_signif_level = TRUE,
                test = "t.test") +
    ylab(NULL) + xlab(NULL) + 
    ggtitle(celltype) +
    theme_bw() + 
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  plot_list[[celltype]] <- p
}
combined_plot <- wrap_plots(plot_list, nrow = 1)
combined_plot
ggsave(paste0("~/03_figure/HCM/大群/13_metabolism_改/TCA_cycle.pdf"),width = 24,height = 6,plot = combined_plot)


#### Fig.S1D ----
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
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Adipocyte")
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

feature_list <- c("Atp metabolic process",
                  "Glycolytic process through glucose 6 phosphate",
                  "Energy derivation by oxidation of organic compounds",
                  "Cellular respiration", 
                  "Acyl coa metabolic process",
                  "Glucose catabolic process",
                  "Triglyceride catabolic process",
                  "Triglyceride metabolic process",
                  "Fatty acid beta oxidation",
                  "Fatty acid catabolic process",
                  "Fatty acyl coa catabolic process",
                  "Regulation of fatty acid oxidation",
                  "Lipid catabolic process",
                  "Lipid homeostasis",
                  "Lipid oxidation",
                  "Positive regulation of lipid metabolic process",
                  "Positive regulation of lipid catabolic process")

feature_list[!(feature_list %in% Dataset$gs_name)]

dataset <- Dataset %>% filter(Dataset$gs_name %in% feature_list) # 835 genes
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

meta <- subset(meta,subset = sub_cell_type %in% 
                 c("Adipocyte"))


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

pathways <- c("Atp metabolic process",
              "Glycolytic process through glucose 6 phosphate",
              "Energy derivation by oxidation of organic compounds",
              "Cellular respiration", 
              "Acyl coa metabolic process",
              "Glucose catabolic process",
              "Triglyceride catabolic process",
              "Triglyceride metabolic process",
              "Fatty acid beta oxidation",
              "Fatty acid catabolic process",
              "Fatty acyl coa catabolic process",
              "Regulation of fatty acid oxidation",
              "Lipid catabolic process",
              "Lipid homeostasis",
              "Lipid oxidation",
              "Positive regulation of lipid metabolic process",
              "Positive regulation of lipid catabolic process")

for (i in pathways) {
  data_melt1 <- subset(data_melt,subset = Pathways %in% i)
  data_melt2 <- subset(data_melt1,subset = (group1 %in% c("HF","HM")))
  data_melt2$Score <- as.numeric(data_melt2$Score)
  Data_summary <- summarySE(data_melt2, 
                            measurevar="Score", 
                            groupvars=c("celltype","group1"))
  
  mean_HM <- mean(data_melt2[data_melt2[,2] == "HM", 4], na.rm = TRUE)
  mean_HF <- mean(data_melt2[data_melt2[,2] == "HF", 4], na.rm = TRUE)
  
  comparison_text <- if(mean_HM > mean_HF) {
    "HM > HF"
  } else if(mean_HF > mean_HM) {
    "HF > HM"
  } else {
    "HM = HF"
  }
  
  p <- ggplot(data_melt2, aes(x = group1, y = Score, fill = group1)) +
    geom_boxplot(outlier.shape = NA) + 
    theme_bw()+
    scale_fill_manual(values = c(
      "#ebd4cc","#d4e4f2"
    ))+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
          panel.grid = element_blank(),
          legend.position = "none")+
    geom_signif(comparisons = list(c("HF","HM")),
                map_signif_level = TRUE,
                test = "t.test",
                tip_length = 0.01,
                textsize = 4)+
    labs(title = i,
         y = NULL,
         x = "Group")+
    annotate("text", 
             x = 1.5,
             y = max(data_melt2[,4], na.rm = TRUE),
             label = comparison_text,
             size = 4,
             color = "black")
  p
  ggsave(p,filename= paste0("~/figure/Adi/AUCell/",i,".pdf"), width = 3, height = 5)
}


#### Fig.S1E ----
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

feature_list <- c("Atp metabolic process",
                  "Glycolytic process through glucose 6 phosphate",
                  "Energy derivation by oxidation of organic compounds",
                  "Cellular respiration", 
                  "Acyl coa metabolic process",
                  "Glucose catabolic process",
                  "Triglyceride catabolic process",
                  "Triglyceride metabolic process",
                  "Fatty acid beta oxidation",
                  "Fatty acid catabolic process",
                  "Fatty acyl coa catabolic process",
                  "Regulation of fatty acid oxidation",
                  "Lipid catabolic process",
                  "Lipid homeostasis",
                  "Lipid oxidation",
                  "Positive regulation of lipid metabolic process",
                  "Positive regulation of lipid catabolic process")

feature_list[!(feature_list %in% Dataset$gs_name)]

dataset <- Dataset %>% filter(Dataset$gs_name %in% feature_list) # 835 genes
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

pathways <- c("Atp metabolic process",
              "Glycolytic process through glucose 6 phosphate",
              "Energy derivation by oxidation of organic compounds",
              "Cellular respiration", 
              "Acyl coa metabolic process",
              "Glucose catabolic process",
              "Triglyceride catabolic process",
              "Triglyceride metabolic process",
              "Fatty acid beta oxidation",
              "Fatty acid catabolic process",
              "Fatty acyl coa catabolic process",
              "Regulation of fatty acid oxidation",
              "Lipid catabolic process",
              "Lipid homeostasis",
              "Lipid oxidation",
              "Positive regulation of lipid metabolic process",
              "Positive regulation of lipid catabolic process")

for (i in pathways) {
  data_melt1 <- subset(data_melt,subset = Pathways %in% i)
  data_melt2 <- subset(data_melt1,subset = (group1 %in% c("HF","HM")))
  data_melt2$Score <- as.numeric(data_melt2$Score)
  Data_summary <- summarySE(data_melt2, 
                            measurevar="Score", 
                            groupvars=c("celltype","group1"))
  
  mean_HM <- mean(data_melt2[data_melt2[,2] == "HM", 4], na.rm = TRUE)
  mean_HF <- mean(data_melt2[data_melt2[,2] == "HF", 4], na.rm = TRUE)
  
  comparison_text <- if(mean_HM > mean_HF) {
    "HM > HF"
  } else if(mean_HF > mean_HM) {
    "HF > HM"
  } else {
    "HM = HF"
  }
  
  p <- ggplot(data_melt2, aes(x = group1, y = Score, fill = group1)) +
    geom_boxplot(outlier.shape = NA) + 
    theme_bw()+
    scale_fill_manual(values = c(
      "#ca9e95","#92a7c4"
    ))+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
          panel.grid = element_blank(),
          legend.position = "none")+
    geom_signif(comparisons = list(c("HF","HM")),
                map_signif_level = TRUE,
                test = "t.test",
                tip_length = 0.01,
                textsize = 4)+
    labs(title = i,
         y = NULL,
         x = "Group")+
    annotate("text", 
             x = 1.5,
             y = max(data_melt2[,4], na.rm = TRUE),
             label = comparison_text,
             size = 4,
             color = "black")
  p
  ggsave(p,filename= paste0("~/figure/Car/AUCell/",i,".pdf"), width = 3, height = 5)
}



