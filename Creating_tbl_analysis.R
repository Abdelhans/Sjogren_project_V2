setwd("D:/data_umr_1227/Par_pureté_cyto_digitale")
library(dplyr)
library(DESeq2)
process_data <- function(gene_file, meta_file, res_cybers, count_file, meta_output_file, count_output_file,ctrl_meta_output_file, ctrl_count_output_file, cybersort_output_file, nomalize_file,SFnomalize_file) {
  
  # Étape 1: Charger les données
  cat("Étape 1: Chargement des données\n")
  gene <- read.csv(gene_file, header = TRUE)
  meta_matrix <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
  #meta_matrix <- meta_matrix[!grepl("^32151043_|^32160439_|^32150952_|^32151647_|^32151039_|^32140205_|^32151707_", meta_matrix$ID), ]
  #meta_matrix <- meta_matrix[!grepl("^32152113_|^32150099_|^32161284_|^32160215_|^32151137_|^32150762_|^32160509|^32151668|^32152294|^32150940|^32152162", meta_matrix$ID), ]
  count_matrix <- read.table(count_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  cybersort <- read.csv(res_cybers, header = TRUE)

  
  #somme des differents composantes des Lymphocytes T et B pour avoir le Total 
  cybersort <- cybersort %>%mutate(LT_purity = (T.cells.CD8 + T.cells.CD4.naive + T.cells.CD4.memory.resting +   T.cells.CD4.memory.activated + T.cells.follicular.helper + T.cells.gamma.delta) * 100) %>%select(-T.cells.CD8, -T.cells.CD4.naive, -T.cells.CD4.memory.resting, -T.cells.CD4.memory.activated, -T.cells.follicular.helper, -T.cells.gamma.delta)
  cybersort <- cybersort %>%mutate(LB_purity = (B.cells.naive + B.cells.memory) * 100) %>%select(-B.cells.naive, -B.cells.memory)
  cybersort$Mixture <- gsub("^X", "",cybersort$Mixture)
  cybersort$M_purity <- sprintf("%.2f", cybersort$Monocytes * 100)
  cybersort$N_purity <- sprintf("%.2f", cybersort$Neutrophils * 100)
  cybersort <- subset(cybersort, select = c("Mixture", "LT_purity", "LB_purity", "M_purity","N_purity"))
  
  # Étape 2: Modifier meta_matrix
  cat("Étape 2: Modification de meta_matrix\n")
  meta_matrix$Lineage_B <- meta_matrix$`B lineage`
  meta_matrix <- merge(meta_matrix, cybersort, by.x = 'ID', by.y= 'Mixture', all.x = TRUE)
 
  meta_matrix <- subset(meta_matrix, CONDITION %in% c("SJS", "CTRL") & GENDER %in% c("Female") & PURITY_CYTOMETRY >= 90.0 , select = c("ID", "GENDER", "AGE", "CELL_TYPE", "CONDITION","PURITY_CYTOMETRY","LT_purity", "LB_purity", "M_purity","N_purity", "T_cell_tot","Lineage_B"))
  meta_matrix <- meta_matrix[order(meta_matrix$ID), ]
  meta_matrix <- meta_matrix[order(meta_matrix$CONDITION, decreasing = FALSE), ] 
  meta_matrix$GROUPE <- cut(meta_matrix$AGE, breaks = c(0, 45, Inf), labels = c("Adult", "Old"))
  

  ctrl_meta_matrix <- subset(meta_matrix, CONDITION %in% c("CTRL") & GENDER %in% c("Female") & PURITY_CYTOMETRY >= 90.0 , select = c("ID", "GENDER", "AGE", "CELL_TYPE", "CONDITION","PURITY_CYTOMETRY","LT_purity", "LB_purity", "M_purity","N_purity", "T_cell_tot","Lineage_B"))
  ctrl_meta_matrix <- ctrl_meta_matrix[order(ctrl_meta_matrix$ID), ]
  ctrl_meta_matrix <- ctrl_meta_matrix[order(ctrl_meta_matrix$CONDITION, decreasing = FALSE), ]
  ctrl_meta_matrix$GROUPE <- cut(ctrl_meta_matrix$AGE, breaks = c(0, 45, Inf), labels = c("Adult", "Old"))
  # Étape 3: Filtrer count_matrix
  cat("Étape 3: Filtrage de count_matrix\n")
  count_matrix <- count_matrix[, colnames(count_matrix) %in% meta_matrix$ID]
  ctrl_count_matrix <- count_matrix[, colnames(count_matrix) %in% ctrl_meta_matrix$ID]
  
  
  designe <- subset(meta_matrix, select = c("CONDITION", "GROUPE", "ID"))
  designe$ID <- rownames(designe)
  designe <- subset(designe, select = -ID)
  designe$GROUPE <- factor(designe$GROUPE)
  designe$CONDITION <- factor(designe$CONDITION)
  
  
  dss <- DESeqDataSetFromMatrix(countData = count_matrix, colData = designe, design = ~ GROUPE + CONDITION)
  normalized_vsd = vst(dss, blind = FALSE)
  normalized_vsd_df <- as.data.frame(assay(normalized_vsd))
  
  
  sss <- estimateSizeFactors(dss)
  norm_sizefactor <- counts(sss, normalized = TRUE)
  norm_sizefactor_df <- as.data.frame(norm_sizefactor)
  
  
  # Étape 4: Écrire les fichiers de sortie
  cat("Étape 4: Écriture des fichiers de sortie\n")
  write.csv(meta_matrix, meta_output_file, row.names = FALSE)
  write.csv(count_matrix, count_output_file, row.names = TRUE)
  write.csv(ctrl_count_matrix, ctrl_count_output_file, row.names = TRUE)
  write.csv(ctrl_meta_matrix, ctrl_meta_output_file, row.names = FALSE)
  write.csv(cybersort, cybersort_output_file, row.names = TRUE)
  write.csv(normalized_vsd_df, nomalize_file, row.names = TRUE)
  write.csv(norm_sizefactor_df, SFnomalize_file, row.names = TRUE)
  cat("Traitement terminé.\n")
}

#process_data("mart_export.csv", "METADATA_SORTED_CELL_MCPcounter.tsv","CIBERSORTx_Job6_Results.csv", "ALL_SORTED_CELLS_MATRIX_ENSG.tsv", "meta_matrix.csv", "count_matrix.csv", "ctrl_meta_matrix.csv", "ctrl_count_matrix.csv","cybersort_collapse.csv", "normalized_vst_2.csv", "normalized_sizefactor_1.csv")
