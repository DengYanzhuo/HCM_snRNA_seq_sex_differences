#### Fig.2A ----

options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) 
library(data.table) 
library(Seurat)
setwd("~/figure/大群/16_WGCNA/")
rm(list = ls())
load(file = "~/data/WGCNA/step1_input2.Rdata")
load(file = "~/data/WGCNA/step2_power_value2.Rdata")
load(file = "~/data/WGCNA/step3_genes_modules2.Rdata")

datTraits$group <- as.factor(datTraits$group)
design <- model.matrix(~0+datTraits$group)
colnames(design) <- levels(datTraits$group)

MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs <- orderMEs(MES0)
moduleTraitCor <- cor(MEs,design,use = "p")

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                     signif(moduleTraitPvalue,1),")")
dim(textMatrix) <- dim(moduleTraitCor)
save(design, file = "~/data/WGCNA/step4_design.Rdata")

# top panel
design_t <- t(design)
moduleTraitCor_t <- t(moduleTraitCor)
MEs_t <- as.data.frame(t(MEs))
textMatrix_t <- t(textMatrix)

pdf("step4_Module-trait-relationship_heatmap.pdf",width = 20,height = 6)
labeledHeatmap(Matrix = moduleTraitCor_t,
               yLabels = colnames(design),
               xLabels = names(MEs),
               xSymbols = names(MEs),
               colorLabels = F,
               colors = colorRampPalette(c("#105eb7", "white","#d7131a"))(50),
               textMatrix = textMatrix_t,
               setStdMargins = F,
               cex.text = 1.2,
               xLabelsPosition = "top",
               xLabelsAdj = 0.5,
               cex.lab.x = 0.7,
               cex.lab.y = 2,
               xLabelsAngle = 0,
               zlim = c(-1,1),
               main = "Module-trait relationships")

dev.off()

# middle panel
gene_module <- data.frame(gene=colnames(datExpr),
                          module=moduleColors)

gene <- read.csv("~/data/deg/rbind_H.csv")
gene <- gene[gene$p_val_adj<0.05,]
gene <- gene[abs(gene$avg_log2FC)>0.25,]

cell_types <- unique(gene$cell_type)
modules <- unique(gene_module$module)

result_df <- data.frame(module = character(),
                        cell_type = character(),
                        pvalue = numeric(),
                        stringsAsFactors = FALSE)

for (i in cell_types) {
  sub_gene <- gene[gene$cell_type == i,]
  cell_type_genes <- sub_gene$markers
  
  all_genes <- colnames(datExpr)
  cell_type_genes <- intersect(cell_type_genes,all_genes)
  for (j in modules){
    module_genes <- gene_module[gene_module$module == j,]
    module_genes <- module_genes$gene
    
    a <- length(intersect(cell_type_genes,module_genes))  # in celltype & in module
    b <- length(module_genes)-a  # only in module
    c <- length(cell_type_genes)-a  # only in celltype
    d <- length(all_genes) - (a + b + c)  # neither
    
    contingency_table <- matrix(c(a, b, c, d), nrow = 2, ncol = 2, byrow = TRUE)
    rownames(contingency_table) <- c("In_module", "Not_in_module")
    colnames(contingency_table) <- c("In_cell_type", "Not_in_cell_type")
    
    fisher_result <- fisher.test(contingency_table, alternative = "greater")
    or_value <- as.numeric(fisher_result$estimate)
    result_df <- rbind(result_df, 
                       data.frame(module = j,
                                  cell_type = i,
                                  pvalue = fisher_result$p.value,
                                  stringsAsFactors = FALSE)) 
  }
  
}

result_df <- result_df[result_df$cell_type != "Epicardial",]
library(dplyr)
result_df <- result_df %>%
  group_by(cell_type) %>%
  mutate(
    FDR = p.adjust(pvalue, method = "BH"),
  ) %>%
  ungroup()
result_df$log10FDR <- -log10(result_df$FDR)

library(ComplexHeatmap)
library(circlize)
library(dplyr)

heatmap_data <- result_df %>%
  select(cell_type, module, log10FDR) %>%
  pivot_wider(names_from = module, values_from = log10FDR) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

star_matrix <- result_df %>%
  select(cell_type, module, FDR) %>%
  mutate(star = ifelse(FDR < 0.05, "*", "")) %>%
  select(cell_type, module, star) %>%
  pivot_wider(
    names_from = module,
    values_from = star,
    values_fill = ""
  ) %>%
  as.data.frame() %>% 
  column_to_rownames("cell_type") %>% 
  as.matrix()

row_order <- c("Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial", "Fibroblast",
               "Lymphatic_endothelial","Lymphocyte","Macrophage","Mast_cell",
               "Neuronal","Pericyte","VSMC")
row_order <- row_order[row_order %in% rownames(heatmap_data)]

col_order <- c("red", "cyan", "tan", "turquoise", "blue", 
               "salmon", "green", "black", "lightcyan", "grey60", 
               "purple", "pink", "brown", "greenyellow", "yellow", 
               "magenta", "midnightblue","grey")
col_order <- col_order[col_order %in% colnames(heatmap_data)]

heatmap_data_ordered <- heatmap_data[row_order, col_order]
star_matrix_ordered <- star_matrix[row_order, col_order]

library(ComplexHeatmap)


pdf("step5_Module-celltype-relationship_heatmap2.pdf",width = 9,height = 4.5)

Heatmap(
  heatmap_data_ordered,
  name = "-log10(FDR)",
  
  col = colorRampPalette(c("white", "#8bc96d", "#5c9e43"))(100),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(star_matrix_ordered[i, j] == "*") {
      grid.text("*", x, y, gp = gpar(fontsize = 12))
    }
  },
  
  row_order = 1:nrow(heatmap_data_ordered), 
  column_order = 1:ncol(heatmap_data_ordered),
  show_row_names = TRUE,
  row_names_side = "left",
  
  show_column_names = T,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  border = "grey60",
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5),
    labels_gp = gpar(fontsize = 5),
    legend_height = unit(10.5, "cm") 
  )
  
)
dev.off()


# bottom panel
library(data.table)
library(dplyr)

gwas_data <- fread('~/data/GWAS/gwas-catalog-download-associations-alt-full.tsv', quote = "", fill = TRUE, check.names = TRUE) %>% 
  dplyr::filter(P.VALUE < 9e-6)

uniprot <- fread('~/data/GWAS/uniprotkb_taxonomy_id_9606_2025_11_26.tsv')
mapping <- uniprot %>%
  dplyr::filter(Reviewed == 'reviewed') %>% 
  tidyr::separate_rows(`Gene Names`) %>% 
  dplyr::rename(GeneName = `Gene Names`)

g_data <- as.data.frame(gwas_data$DISEASE.TRAIT)
g_data <- g_data %>%
  add_count(`gwas_data$DISEASE.TRAIT`)
g2_data <- as.data.frame(unique(g_data$`gwas_data$DISEASE.TRAIT`))


disease_list <- c("Left_ventricular_mass",
                  "Cardiac_hypertrophy",
                  "Left_ventricle_systolic_dysfunction",
                  "Left_ventricle_wall_thickness",
                  "Hypertrophic_cardiomyopathy")

Left_ventricular_mass <- c("Left ventricular mass",
                           "Left ventricular mass to end-diastolic volume ratio (MTAG)",
                           "Left ventricular mass (MTAG)",
                           "Left ventricular mass to end-diastolic volume ratio",
                           "Left ventricular mass x antihypertensive treatment interaction",
                           "Left ventricular mass indexed by body surface area",
                           "Left ventricular concentricity (left ventricular mass / left ventricular end diastolic volume)",
                           "Left ventricular mass (indexed to body surface area)")

Left_ventricle_systolic_dysfunction <- c("Left ventricular ejection fraction",
                                         "Left ventricular ejection fraction (MTAG)",
                                         "Left  ventricle systolic dysfunction",
                                         "Left ventricular global longitudinal strain",
                                         "Left ventricular global longitudinal strain (MTAG)")

Left_ventricle_wall_thickness <- c("Left ventricle wall thickness",
                                   "Left ventricular end systole anterolateral wall thickness",
                                   "Left ventricular end systole inferoseptal wall thickness",
                                   "Left ventricular end systole anterior wall thickness",
                                   "Left ventricular end systole inferior wall thickness",
                                   "Left ventricular end systole anterospetal wall thickness",
                                   "Left ventricular end systole inferolateral wall thickness",
                                   "Left ventricular end diastole inferoseptal wall thickness",
                                   "Left ventricular end diastole inferior wall thickness",
                                   "Left ventricular end diastole anterolateral wall thickness",
                                   "Left ventricular end diastole anterior wall thickness",
                                   "Left ventricular end diastole anterospetal wall thickness",
                                   "Left ventricular wall thickness (maximal)",
                                   "Left ventricular wall thickness (mean)"
)
Hypertrophic_cardiomyopathy <- c("Hypertrophic cardiomyopathy (MTAG)",
                                 "Hypertrophic cardiomyopathy",
                                 "Hypertrophic cardiomyopathy (sarcomere negative)",
                                 "Hypertrophic cardiomyopathy (sarcomere positive)",
                                 "Hypertrophic cardiomyopathy (ICD-10 coded)",
                                 "Hypertrophic cardiomyopathy (sarcomere-positive)",
                                 "Hypertrophic cardiomyopathy (sarcomere-negative)"
)

gwas_genelist <- list()
gene_df <- data.frame(group = character(),
                      dis = character(),
                      gene = character())

for (disease in disease_list){
  phenotypes <- get(disease)
  for(dis_cate in phenotypes){
    sub_df <- gwas_data %>%
      dplyr::filter(DISEASE.TRAIT == dis_cate)
    gwas_genelist[[disease]][[dis_cate]] <- sub_df
    
    human_genes <- data.frame(
      GeneName = unique(c(unlist(strsplit(sub_df$REPORTED.GENE.S., ", ")),
                          unlist(strsplit(sub_df$MAPPED_GENE, ",\\s*"))))
    ) %>% 
      dplyr::inner_join(mapping) %>% 
      dplyr::select(`Gene Names (primary)`) %>% 
      pull() %>% 
      unique()
    gene_df <- rbind(gene_df,
                     data.frame(group = disease,
                                dis = dis_cate,
                                gene = human_genes))
  }}

table(gene_df$group)

table(gene_df$dis)

count <- as.data.frame(table(gene_df$dis))

dcount <- as.data.frame(table(gene_df$group))

disease_phenotypes <- rbind(
  data.frame(
    Disease_Group = "Left_ventricular_mass",
    Phenotype = Left_ventricular_mass,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Disease_Group = "Left_ventricle_systolic_dysfunction",
    Phenotype = Left_ventricle_systolic_dysfunction,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Disease_Group = "Left_ventricle_wall_thickness",
    Phenotype = Left_ventricle_wall_thickness,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Disease_Group = "Hypertrophic_cardiomyopathy",
    Phenotype = Hypertrophic_cardiomyopathy,
    stringsAsFactors = FALSE
  )
)

colnames(count) <- c("Phenotype","gene_count")
count <- merge(count,disease_phenotypes,by="Phenotype")

write.csv(gene_df,"~/data/GWAS/HCM_GWAS.csv")

gene_module <- data.frame(gene=colnames(datExpr),
                          module=moduleColors)

gene <- read.csv("~/data/GWAS/HCM_GWAS.csv")
disease <- unique(gene$group)
modules <- unique(gene_module$module)

result_df <- data.frame(module = character(),
                        disease = character(),
                        pvalue = numeric(),
                        stringsAsFactors = FALSE)

for (i in disease) {
  sub_gene <- gene[gene$group == i,]
  disease_genes <- sub_gene$gene
  
  all_genes <- colnames(datExpr)
  disease_genes <- intersect(disease_genes,all_genes) 
  for (j in modules){
    module_genes <- gene_module[gene_module$module == j,]
    module_genes <- module_genes$gene
    
    a <- length(intersect(disease_genes,module_genes))  # in celltype & in module
    b <- length(module_genes)-a  # only in module
    c <- length(disease_genes)-a  # only in celltype
    d <- length(all_genes) - (a + b + c)  # neither
    
    contingency_table <- matrix(c(a, b, c, d), nrow = 2, ncol = 2, byrow = TRUE)
    rownames(contingency_table) <- c("In_module", "Not_in_module")
    colnames(contingency_table) <- c("In_disease", "Not_in_disease")
    
    fisher_result <- fisher.test(contingency_table, alternative = "greater")
    or_value <- as.numeric(fisher_result$estimate)
    result_df <- rbind(result_df, 
                       data.frame(module = j,
                                  disease = i,
                                  pvalue = fisher_result$p.value,
                                  stringsAsFactors = FALSE)) 
  }
  
}

library(dplyr)
library(tidyr)
result_df <- result_df %>%
  group_by(disease) %>%
  mutate(
    FDR = p.adjust(pvalue, method = "BH"),
  ) %>%
  ungroup()
result_df$log10FDR <- -log10(result_df$FDR)

library(ComplexHeatmap)
library(circlize)
library(dplyr)

result_df <- as.data.frame(result_df)
heatmap_data <- result_df %>%
  dplyr::select(disease, module, log10FDR) %>%
  pivot_wider(names_from = module, values_from = log10FDR) %>%
  column_to_rownames("disease") %>%
  as.matrix()

star_matrix <- result_df %>%
  dplyr::select(disease, module, FDR) %>%
  dplyr::mutate(star = ifelse(FDR < 0.05, "*", "")) %>%
  dplyr::select(disease, module, star) %>%
  tidyr::pivot_wider(
    names_from = module,
    values_from = star,
    values_fill = ""
  ) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("disease") %>%
  as.matrix()
col_order <- c("red", "cyan", "tan", "turquoise", "blue", 
               "salmon", "green", "black", "lightcyan", "grey60", 
               "purple", "pink", "brown", "greenyellow", "yellow", 
               "magenta", "midnightblue","grey")
col_order <- col_order[col_order %in% colnames(heatmap_data)]

heatmap_data_ordered <- heatmap_data[, col_order]
star_matrix_ordered <- star_matrix[, col_order]


pdf("step18_Module-GWAS-relationship_heatmap2.pdf",width = 9,height = 2.5)

Heatmap(
  heatmap_data_ordered,
  name = "-log10(FDR)",
  
  col = colorRampPalette(c("white", "#a4cde1", "#277fb8"))(100),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(star_matrix_ordered[i, j] == "*") {
      grid.text("*", x, y, gp = gpar(fontsize = 12))
    }
  },
  
  row_order = 1:nrow(heatmap_data_ordered), 
  column_order = 1:ncol(heatmap_data_ordered),
  
  show_row_names = TRUE,
  row_names_side = "left",
  
  show_column_names = T,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  border = "grey60",
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5),
    labels_gp = gpar(fontsize = 5),
    legend_height = unit(4, "cm")
  )
  
)
dev.off()


#### Fig.2B ----
load(file = '~/data/WGCNA/step1_input.Rdata')
load(file = "~/data/WGCNA/step2_power_value.Rdata")
load(file = "~/data/WGCNA/step3_genes_modules.Rdata")
load(file = "~/data/WGCNA/step4_design.Rdata")

OrgDb = "org.Hs.eg.db"
genetype = "SYMBOL"

choose_module <- c("brown","salmon","midnightblue","purple","tan") 

library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

gene_module <- data.frame(gene=colnames(datExpr),
                          module=moduleColors)
tmp <- bitr(gene_module$gene,
            fromType = genetype,
            toType = "ENTREZID",
            OrgDb = OrgDb )
gene_module_entrz <- merge(tmp,
                           gene_module, 
                           by.x=genetype, 
                           by.y="gene")

choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]

formula_res <- compareCluster(
  ENTREZID~module,
  data = choose_gene_module_entrz,
  fun = "enrichGO",
  OrgDb = OrgDb,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

save(gene_module, formula_res, file="~/data/WGCNA/step6_module_GO_term.Rdata")
write.csv(formula_res@compareClusterResult,
          file="~/data/WGCNA/step6_module_GO_term.csv")

load("~/data/WGCNA/step6_module_GO_term.Rdata")
module_term <- formula_res@compareClusterResult
description <- c(
  #brown
  "muscle cell development",
  "muscle system process",
  "mitochondrial transmembrane transport",
  "cardiac muscle cell action potential",
  "energy derivation by oxidation of organic compounds",
  
  #midnightblue
  "extracellular matrix organization",
  "regulation of integrin-mediated signaling pathway",
  "regulation of plasma membrane organization",
  "external encapsulating structure organization",
  "extracellular structure organization",
  
  #salmon
  "extracellular matrix organization",
  "morphogenesis of a branching epithelium",
  "smoothened signaling pathway",
  "positive regulation of cytosolic calcium ion concentration",
  "positive regulation of chondrocyte differentiation",
  
  #tan
  "small GTPase-mediated signal transduction",
  "cell-matrix adhesion",
  "cellular response to vascular endothelial growth factor stimulus",
  "calcium ion homeostasis",
  "positive regulation of Ras protein signal transduction"
  
)
term_selected <- module_term[module_term$Description %in% description,]

library(dplyr)
term_selected <- term_selected %>% 
  arrange(Cluster, p.adjust)

formula_res@compareClusterResult <- term_selected


dotp <- dotplot(formula_res,
                showCategory=10,
                includeAll = TRUE,
                label_format=90)
dotp
ggsave(dotp,filename= "step6_module_GO_term.pdf",
       width = 9.5,
       height = 6.5)


#### Fig.2C ----
load(file = '~/data/WGCNA/step1_input2.Rdata')
load(file = "~/data/WGCNA/step2_power_value2.Rdata")
load(file = "~/data/WGCNA/step3_genes_modules2.Rdata")
load(file = "~/data/WGCNA/step4_design2.Rdata")

gene_module <- data.frame(gene=colnames(datExpr),
                          module=moduleColors)

gene <- read.csv("~/data/GWAS/HCM_GWAS.csv")
disease <- unique(gene$group)
modules <- unique(gene_module$module)

j <- "brown"
i <- "Hypertrophic_cardiomyopathy"
sub_gene <- gene[gene$group == i,]
module_gene <- gene_module[gene_module$module == j,]
tgene <- intersect(sub_gene$gene,module_gene$gene)
print(tgene)
# "HSPB7"   "PDE3A"   "SRL"     "GMPR"    "MLIP"    "TMEM182" "ADPRHL1" "RNF207"  "TTN"  "RBFOX1" 

i <- "Left_ventricle_systolic_dysfunction"
sub_gene <- gene[gene$group == i,]
module_gene <- gene_module[gene_module$module == j,]
tgene <- intersect(sub_gene$gene,module_gene$gene)
print(tgene)
# "TTN"    "ALPK3"  "HSPB7"  "TRIM63" "GATA4"

seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Cardiomyocyte")
seurat_obj <- NormalizeData(seurat_obj)
DefaultAssay(seurat_obj) <- "RNA"
genes <- c("HSPB7",
           "PDE3A",
           "SRL",
           "GMPR",
           "MLIP",
           "ADPRHL1",
           "RBFOX1",
           "ALPK3",
           "HSPB7",
           "TRIM63",
           "GATA4" )
genes <- unique(genes)
seurat_obj <- subset(seurat_obj, features = genes)

Idents(seurat_obj) <- "group"
data <- FindMarkers(
  seurat_obj,
  ident.1 = "HCM",
  ident.2 = "NF",
  min.pct = 0,
  logfc.threshold = 0
)
data$gene <- rownames(data)
write.csv(data,"~/02_work_data/HCM_NF/gsea.input/gsea_H/GWAS_Cardiomyocyte_HCM_vs_NF.csv")

H_data <- read.csv("~/02_work_data/HCM_NF/gsea.input/gsea_H/Cardiomyocyte_input.csv")
N_data <- read.csv("~/02_work_data/HCM_NF/gsea.input/gsea_N/Cardiomyocyte_input.csv")
H_N_data <- read.csv("~/02_work_data/HCM_NF/gsea.input/gsea_H/GWAS_Cardiomyocyte_HCM_vs_NF.csv")

H_data <- H_data[H_data$gene %in% genes,]
H_data$group <- "HF: HM"
H_data <- H_data[,c("p_val","avg_log2FC","p_val_adj","gene","group")]
N_data <- N_data[N_data$gene %in% genes,]
N_data$group <- "NF: NM"
N_data <- N_data[,c("p_val","avg_log2FC","p_val_adj","gene","group")]
H_N_data <- H_N_data[H_N_data$gene %in% genes,]
H_N_data$group <- "HCM: Normal"
H_N_data <- H_N_data[,c("p_val","avg_log2FC","p_val_adj","gene","group")]
data2 <- rbind(H_data,N_data,H_N_data)
data2$`-log10(p_val)` <- -log10(data2$p_val)
data2 <- data2[data2$p_val_adj<0.05,]
data2$group <- factor(data2$group, 
                      levels = c("HF: HM", "NF: NM","HCM: Normal"))
max_finite <- max(data2$`-log10(p_val)`[is.finite(data2$`-log10(p_val)`)])
data2$`-log10(p_val)`[is.infinite(data2$`-log10(p_val)`)] <- max_finite + 1
source("/home/gongfengcz/Retina/src/visualization/custom_plot_function.R")
ggplot(data2, aes(x = group, y = gene)) +
  geom_point(aes(size = `-log10(p_val)`, 
                 color = avg_log2FC),
             shape = 19) +
  scale_color_gradient2(
    low = "#277fb8",
    mid = "white",
    high = "#ec5051",
    midpoint = 0,
    name = "log2FC"
  ) +
  scale_size_continuous(
    range = c(5, 10),
    name = "-log10(pvalue)"
  )+
  theme_cat()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "plain"),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    legend.box = "vertical"
  )
ggsave("~/03_figure/HCM/大群/16_WGCNA/step18_gene_log2FC_all.pdf",height = 5,width = 4)
