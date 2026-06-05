#### Fig.S4A ----
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Cardiomyocyte")
subcell_color <- c(
  'Cardiomyocyte_I' = "#f9b769",
  'Cardiomyocyte_II' = "#af93c4",
  'Cardiomyocyte_III' = "#67a4cc"
)


ratio<-prop.table(table(seurat_obj$cell_type_leiden0.6,seurat_obj$biosample_id),margin=2)
ratio<-as.data.frame(ratio)
colnames(ratio) <- c("celltype","biosample_id","ratio")
metadata <- seurat_obj@meta.data[, c("biosample_id", "group1")]
ratio$group1 <- metadata$group1[match(ratio$biosample_id, metadata$biosample_id)]

CarI <- ratio[ratio$celltype == "Cardiomyocyte_I", ]
CarI$group1 <-factor(CarI$group1,levels = c("NF","NM","HF","HM"))
ggplot(CarI, aes(x = CarI$group1, y = CarI$ratio))+ 
  labs(y="cell_ratio",x= "group",title = "Cardiomyocyte_I")+ 
  geom_boxplot(aes(fill = group1),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("HF" = "#ec5051",
                               "HM" = "#67a4cc",
                               "NF" = "#f38989",
                               "NM" = "#a4cde1"))+theme_classic() + 
  stat_compare_means(
    comparisons = list(
      c("NM", "HM"),
      c("NF", "HF"),
      c("NM", "NF"),
      c("HM", "HF")
      
    ),
    method = "wilcox.test",
    label = "p.signif"
  )

ggsave("./CarI_boxplot.pdf",width = 4,height = 3.5)


CarII <- ratio[ratio$celltype == "Cardiomyocyte_II", ]
CarII$group1 <-factor(CarII$group1,levels = c("NF","NM","HF","HM"))
ggplot(CarII, aes(x = CarII$group1, y = CarII$ratio))+ 
  labs(y="cell_ratio",x= "group",title = "Cardiomyocyte_II")+ 
  geom_boxplot(aes(fill = group1),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("HF" = "#ec5051",
                               "HM" = "#67a4cc",
                               "NF" = "#f38989",
                               "NM" = "#a4cde1"))+theme_classic() + 
  stat_compare_means(
    comparisons = list(
      c("NM", "HM"),
      c("NF", "HF"),
      c("NM", "NF"),
      c("HM", "HF")
      
    ),
    method = "wilcox.test",
    label = "p.signif"
  )

ggsave("./CarII_boxplot.pdf",width = 4,height = 3.5)


CarIII <- ratio[ratio$celltype == "Cardiomyocyte_III", ]
CarIII$group1 <-factor(CarIII$group1,levels = c("NF","NM","HF","HM"))
ggplot(CarIII, aes(x = CarIII$group1, y = CarIII$ratio))+ 
  labs(y="cell_ratio",x= "group",title = "Cardiomyocyte_III")+ 
  geom_boxplot(aes(fill = group1),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("HF" = "#ec5051",
                               "HM" = "#67a4cc",
                               "NF" = "#f38989",
                               "NM" = "#a4cde1"))+theme_classic() + 
  stat_compare_means(
    comparisons = list(
      c("NM", "HM"),
      c("NF", "HF"),
      c("NM", "NF"),
      c("HM", "HF")
      
    ),
    method = "wilcox.test",
    label = "p.signif"
  )

ggsave("./CarIII_boxplot.pdf",width = 4,height = 3.5)



#### Fig.S4B ----
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Fibroblast")
seurat_obj$cell_type_leiden0.6[seurat_obj$cell_type_leiden0.6 %in% c("Fibroblast_I","Fibroblast_II")] <- "Fibroblast"

ratio<-prop.table(table(seurat_obj$cell_type_leiden0.6,seurat_obj$biosample_id),margin=2)
ratio<-as.data.frame(ratio)
colnames(ratio) <- c("celltype","biosample_id","ratio")
metadata <- seurat_obj@meta.data[, c("biosample_id", "group1")]
ratio$group1 <- metadata$group1[match(ratio$biosample_id, metadata$biosample_id)]

Fib <- ratio[ratio$celltype == "Fibroblast", ]
Fib$group1 <-factor(Fib$group1,levels = c("NF","NM","HF","HM"))
ggplot(Fib, aes(x = Fib$group1, y = Fib$ratio))+ 
  labs(y="cell_ratio",x= "group",title = "Fibroblast")+ 
  geom_boxplot(aes(fill = group1),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("HF" = "#ec5051",
                               "HM" = "#67a4cc",
                               "NF" = "#f38989",
                               "NM" = "#a4cde1"))+theme_classic() + 
  stat_compare_means(
    comparisons = list(
      c("NM", "HM"),
      c("NF", "HF"),
      c("NM", "NF"),
      c("HM", "HF")
      
    ),
    method = "wilcox.test",
    label = "p.signif"
  )

ggsave("./Fib_boxplot.pdf",width = 4,height = 3.5)


Fiba <- ratio[ratio$celltype == "Activated_fibroblast", ]
Fiba$group1 <-factor(Fiba$group1,levels = c("NF","NM","HF","HM"))
ggplot(Fiba, aes(x = Fiba$group1, y = Fiba$ratio))+ 
  labs(y="cell_ratio",x= "group",title = "Activated_fibroblast")+ 
  geom_boxplot(aes(fill = group1),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("HF" = "#ec5051",
                               "HM" = "#67a4cc",
                               "NF" = "#f38989",
                               "NM" = "#a4cde1"))+theme_classic() + 
  stat_compare_means(
    comparisons = list(
      c("NM", "HM"),
      c("NF", "HF"),
      c("NM", "NF"),
      c("HM", "HF")
      
    ),
    method = "wilcox.test",
    label = "p.signif"
  )

ggsave("./FibA_boxplot.pdf",width = 4,height = 3.5)





#### Fig.S4C ----
seurat_obj <- readRDS("~/02_work_data/HCM_NF/raw_data/seurat_obj_group1.rds")
seurat_obj <- subset(seurat_obj,subset = cell_type %in% "Macrophage")

ratio<-prop.table(table(seurat_obj$cell_type_leiden0.6,seurat_obj$biosample_id),margin=2)
ratio<-as.data.frame(ratio)
colnames(ratio) <- c("celltype","biosample_id","ratio")
metadata <- seurat_obj@meta.data[, c("biosample_id", "group1")]
ratio$group1 <- metadata$group1[match(ratio$biosample_id, metadata$biosample_id)]

Mac <- ratio[ratio$celltype == "Macrophage", ]
Mac$group1 <-factor(Mac$group1,levels = c("NF","NM","HF","HM"))
ggplot(Mac, aes(x = Mac$group1, y = Mac$ratio))+ 
  labs(y="cell_ratio",x= "group",title = "Macrophage")+ 
  geom_boxplot(aes(fill = group1),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("HF" = "#ec5051",
                               "HM" = "#67a4cc",
                               "NF" = "#f38989",
                               "NM" = "#a4cde1"))+theme_classic() + 
  stat_compare_means(
    comparisons = list(
      c("NM", "HM"),
      c("NF", "HF"),
      c("NM", "NF"),
      c("HM", "HF")
      
    ),
    method = "wilcox.test",
    label = "p.signif"
  )

ggsave("./Mac_boxplot.pdf",width = 4,height = 3.5)


Macp <- ratio[ratio$celltype == "Proliferating_macrophage", ]
Macp$group1 <-factor(Macp$group1,levels = c("NF","NM","HF","HM"))
ggplot(Macp, aes(x = Macp$group1, y = Macp$ratio))+ 
  labs(y="cell_ratio",x= "group",title = "Proliferating_macrophage")+ 
  geom_boxplot(aes(fill = group1),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("HF" = "#ec5051",
                               "HM" = "#67a4cc",
                               "NF" = "#f38989",
                               "NM" = "#a4cde1"))+theme_classic() + 
  stat_compare_means(
    comparisons = list(
      c("NM", "HM"),
      c("NF", "HF"),
      c("NM", "NF"),
      c("HM", "HF")
      
    ),
    method = "wilcox.test",
    label = "p.signif"
  )

ggsave("./Macp_boxplot.pdf",width = 4,height = 3.5)
