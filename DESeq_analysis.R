setwd("D:/data_umr_1227/Par_pureté_cyto_digitale")

DESeq_analysis <- function(meta_file, count_file, gene_file, rescyber, ctrl_meta, ctrl_count, digital_purity ,cell_type) {
  
  library(DESeq2)
  library(pheatmap)
  library(dplyr)
  library(RColorBrewer)
  library(ggplot2)
  library(ggrepel)
  
  # Définir le chemin du dossier DESeq 
  
  print("Step 0 : Le dossier n'existe pas ,creation")
  folder_name <- paste0(cell_type, "_DESeq")
  folder_path <- file.path(getwd(), folder_name)
  
  if (!file.exists(folder_path)) {
    dir.create(folder_path)
  }
  # Charger les données
    
  print("Step 1 : Chargement des données")
  meta_matrix <- read.csv("meta_matrix.csv", header = TRUE)
  count_matrix <- read.csv("count_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
  gene_matrix <- read.csv("mart_export.csv", header = TRUE)
  cybersort <- read.csv("cybersort_collapse.csv", header = TRUE)
  ctrl_meta_matrix <- read.csv("ctrl_meta_matrix.csv", header = TRUE)
  ctrl_count_matrix <- read.csv("ctrl_count_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)

  print("Step 2 : Filtrage des métadonnées pour le type cellulaire spécifié :")
  # Filtrer les métadonnées pour le type cellulaire spécifié
  cell_meta_matrix <- meta_matrix[meta_matrix$CELL_TYPE == cell_type & meta_matrix[[digital_purity]] >= 90.0, c('ID', 'CONDITION', 'GROUPE')]
  ctrl_meta_matrix <- ctrl_meta_matrix[ctrl_meta_matrix$CELL_TYPE == cell_type & ctrl_meta_matrix[[digital_purity]] >= 90.0, c('ID', 'GROUPE')]
  
  print("Step 3 : Filtrage de la matrice de comptage pour les ID correspondant au type cellulaire :")
  # Filtrer la matrice de comptage pour les ID correspondant au type cellulaire
  cell_count_matrix <- count_matrix[, colnames(count_matrix) %in% cell_meta_matrix$ID]
  ctrl_count_matrix <- ctrl_count_matrix[, colnames(ctrl_count_matrix) %in% ctrl_meta_matrix$ID]
  # Vérifier les dimensions des matrices
  if (ncol(cell_count_matrix) == 0) {
    stop("Aucun ID correspondant trouvé pour le type cellulaire spécifié.")}
  if (ncol(ctrl_count_matrix) == 0) {
    stop("Aucun ID correspondant trouvé pour le type cellulaire spécifié pour le controle.")}
  
  print("Step 4: Mettre les variables en facteur :")
  
  # Mettre les variables en facteur
  
  cell_meta_matrix <- cell_meta_matrix[order(cell_meta_matrix$ID), ]
  if (!all(colnames(cell_count_matrix) == cell_meta_matrix$ID)) {
    stop("Les noms de colonnes de cell_count_matrix ne correspondent pas aux ID de cell_meta_matrix")}
  
  ctrl_meta_matrix <- ctrl_meta_matrix[order(ctrl_meta_matrix$ID), ]
  if (!all(colnames(ctrl_count_matrix) == ctrl_meta_matrix$ID)) {
    stop("Les noms de colonnes de ctrl_count_matrix ne correspondent pas aux ID de ctrl_meta_matrix")}
  
  cell_meta_matrix$GROUPE <- factor(cell_meta_matrix$GROUPE)
  cell_meta_matrix$CONDITION <- factor(cell_meta_matrix$CONDITION)
  ctrl_meta_matrix$GROUPE <- factor(ctrl_meta_matrix$GROUPE)
  
  print("Step 5: Créer un objet DESeq :")

  # Créer un objet DESeq
  dss <- DESeqDataSetFromMatrix(countData = cell_count_matrix, colData = cell_meta_matrix, design = ~ GROUPE + CONDITION)
  dss_ctrl <- DESeqDataSetFromMatrix(countData = ctrl_count_matrix, colData = ctrl_meta_matrix, design = ~ GROUPE)
  print("Step 6:  Fixation de la condition de référence :")
  
  # Fixer la condition de référence
  
  dss$CONDITION <- relevel(dss$CONDITION, ref = "CTRL")
  dss_ctrl$GROUPE <- relevel(dss_ctrl$GROUPE, ref = "Adult")

  print(" Step 7: Run process DESeq :")
  
  # Lancer le process DESeq
  des <- DESeq(dss)
  des_ctrl <- DESeq(dss_ctrl)

  print(" Step 8: Résultats Deseq :")
 
  # Résultats DESeq
  results <- results(des)
  results_ctrl <- results(des_ctrl)

  print("Step 9: objet DESeq2 into  dataframe :")
 
  results <- as.data.frame(results)
  results_ctrl <- as.data.frame(results_ctrl)
  print("Step 10: Filter résultats by  p-value adj :")
  
  # Filtrer les résultats par p-value ajustée
  
  DEG <- results %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5)
  DEG_ctrl <- results_ctrl %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5)

  print(" Step 11 : Filter result by log fold change :")
  
  # Filtrer les résultats par log fold change
  UP <- DEG%>%filter(DEG$log2FoldChange > 0.5)
  DOWN <- DEG%>%filter(DEG$log2FoldChange < -0.5)
  
  print("Step 12 :  add gene names :")
  
  # Ajouter les noms de gènes
  DEG_ctrl$Gene.stable.ID <- rownames(DEG_ctrl)
  results$Gene.stable.ID <- rownames(results)
  DEG$Gene.stable.ID <- rownames(DEG)
  UP$Gene.stable.ID <- rownames(UP)
  DOWN$Gene.stable.ID <- rownames(DOWN)
  
  print(" Step 13: find gene names  :")
 
    process_dataframe <- function(data_norm, gene_name) {
  merged_df <- merge(data_norm, gene_name, by = "Gene.stable.ID", all.x = TRUE)
  merged_df$Gene.name[is.na(merged_df$Gene.name)] <- merged_df$Gene.stable.ID[is.na(merged_df$Gene.name)]
  return(merged_df)
}
 
  DEG_ctrl <- process_dataframe(DEG_ctrl, gene_matrix)
  results<- process_dataframe(results, gene_matrix)
  DEG <- process_dataframe(DEG, gene_matrix)
  UP <- process_dataframe(UP, gene_matrix)
  DOWN <- process_dataframe(DOWN, gene_matrix)
  
  print("Step 14: Vérification et comparaison des gènes différentiellement exprimés")

  if (nrow(DEG_ctrl) == 0) {
    print("Aucun gène différentiellement exprimé pour les groupes d'âge.")
    genes_supprimes <- character(0)  # Créer un vecteur vide pour éviter les erreurs
  } else {
    genes_supprimes <- DEG_ctrl$Gene.stable.ID[DEG_ctrl$Gene.stable.ID %in% DEG$Gene.stable.ID]
    if (length(genes_supprimes) == 0) {
      print("Aucun des gènes différentiellement exprimés pour l'âge n'a été supprimé.")
    } else {
      print("Les gènes différentiellement exprimés pour l'âge suivants ont été supprimés:")
      print(genes_supprimes)
    }
  }
  # Continuer avec les gènes différentiellement exprimés
  print(" Step 15 : saving results  :")
  # Supprimer les gènes de DEG_ctrl de DEG
  DEG <- DEG[!DEG$Gene.stable.ID %in% genes_supprimes, ]
  DEG %>% arrange(padj)
  UP %>% arrange(padj)
  DOWN %>% arrange(padj)
  
  normalized_vsd = vst(dss, blind = FALSE)
  normalized_vsd_df <- as.data.frame(assay(normalized_vsd))
  normalized_vsd_df$Gene.stable.ID <- rownames(normalized_vsd_df)
  normalized_vsd_df <- merge(normalized_vsd_df, gene_matrix, by = "Gene.stable.ID", all.x = TRUE)
  normalized_vsd_df$Gene.name[is.na(normalized_vsd_df$Gene.name)] <- normalized_vsd_df$Gene.stable.ID[is.na(normalized_vsd_df$Gene.name)]
  # Identification des lignes où Gene.name est vide ou NA
  empty_gene_names <- is.na(normalized_vsd_df$Gene.name) | normalized_vsd_df$Gene.name == ""
  # Compteur pour créer des noms uniques
  counter <- 1
  # Boucle pour remplacer chaque valeur vide par un nom unique
  normalized_vsd_df$Gene.name[empty_gene_names] <- sapply(seq(sum(empty_gene_names)), function(x) {
    paste0("no_gene_name_", counter + x - 1)
  })
 
  # Rendre les valeurs de la colonne Gene.name uniques
  normalized_vsd_df$Gene.name <- make_unique(normalized_vsd_df$Gene.name)
  rownames(normalized_vsd_df) <- normalized_vsd_df$Gene.name
  normalized_vsd_df <- subset(normalized_vsd_df, select = -Gene.name)
  normalized_vsd_df <- subset(normalized_vsd_df, select = -Gene.stable.ID)
  counts_sum <- rowSums(normalized_vsd_df[, -1])
  normalized_vsd_df <- normalized_vsd_df[counts_sum != 0, ]

  write.csv(normalized_vsd_df, file = file.path(folder_path, paste0("VST_normalized_", cell_type, "_DESeq.csv")), row.names = TRUE)
  write.csv(results, file = file.path(folder_path, paste0("resultat_", cell_type, "_DESeq.csv")), row.names = TRUE)
  write.csv(DEG, file = file.path(folder_path, paste0("DEG_", cell_type, "_DESeq.csv")), row.names = TRUE)
  write.csv(UP, file = file.path(folder_path, paste0("UP_", cell_type, "_DESeq.csv")), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.csv(DOWN,file = file.path(folder_path, paste0("DOWN_", cell_type, "_DESeq.csv")), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.csv(DEG_ctrl, file = file.path(folder_path, paste0("DEGctrl_", cell_type, "_DESeq.csv")), row.names = FALSE)
  write.csv(genes_supprimes, file = file.path(folder_path, paste0("DEGgroupe_supprimés_", cell_type, "_DESeq.csv")), row.names = FALSE)
  write.csv(colData(des), file = file.path(folder_path, paste0("colData_", cell_type, "_DESeq.csv")), row.names = TRUE)
  write.csv(counts(des), file = file.path(folder_path, paste0("countData_", cell_type, "_DESeq.csv")), row.names = TRUE)
  write.csv(colData(des_ctrl), file = file.path(folder_path, paste0("ctrlcolData_", cell_type, "_DESeq.csv")), row.names = TRUE)
  write.csv(counts(des_ctrl), file = file.path(folder_path, paste0("ctrlcountData_", cell_type, "_DESeq.csv")), row.names = TRUE)

  print("Step 16: Création du volcanoplot :")
  
  # Ajouter une colonne pour la couleur du volcanoplot
  results$Significance <- "Not Significant"
  results$Significance[results$padj < 0.05 & results$log2FoldChange > 0.5] <- "Upregulated"
  results$Significance[results$padj < 0.05 & results$log2FoldChange < -0.5] <- "Downregulated"
  
  # Créer le volcanoplot
  volcano_plot <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
    theme_minimal() +
    geom_text_repel(data = subset(results, padj < 0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)),
                    aes(label = Gene.name),
                    size = 3) +
    labs(title = paste("Volcano Plot of Differentially Expressed Genes in", cell_type),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted p-value") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") + theme(plot.background = element_rect(fill = "white"),  # Assurer un fond blanc
        panel.background = element_rect(fill = "white"))
  
  # Afficher le volcanoplot
  print(volcano_plot)
   
  # Sauvegarder le volcanoplot
  ggsave(file.path(folder_path, paste0(cell_type, "_volcanoplot.png")), plot = volcano_plot, width = 10, height = 8, bg = "white")

  print(" Step 17 :  Estimation des comptes avec facteurs de tailles")

  sss <- estimateSizeFactors(dss)
  norm_sizefactor <- counts(sss, normalized = TRUE)
  norm_sizefactor_df <- as.data.frame(norm_sizefactor)

  norm_sizefactor_df$Gene.stable.ID <- rownames(norm_sizefactor_df)
  norm_sizefactor_df <- merge(norm_sizefactor_df, gene_matrix, by = "Gene.stable.ID", all.x = TRUE)
  norm_sizefactor_df$Gene.name[is.na(norm_sizefactor_df$Gene.name)] <- norm_sizefactor_df$Gene.stable.ID[is.na(norm_sizefactor_df$Gene.name)]
  # Identification des lignes où Gene.name est vide ou NA
  empty_gene_names <- is.na(norm_sizefactor_df$Gene.name) | norm_sizefactor_df$Gene.name == ""
  # Compteur pour créer des noms uniques
  counter <- 1
  # Boucle pour remplacer chaque valeur vide par un nom unique
  norm_sizefactor_df$Gene.name[empty_gene_names] <- sapply(seq(sum(empty_gene_names)), function(x) {
    paste0("no_gene_name_", counter + x - 1)
  })
  norm_sizefactor_df <- subset(norm_sizefactor_df, select = -Gene.stable.ID)
  # Rendre les valeurs de la colonne Gene.name uniques
  norm_sizefactor_df$Gene.name <- make_unique(norm_sizefactor_df$Gene.name)
  rownames(norm_sizefactor_df) <- norm_sizefactor_df$Gene.name
  norm_sizefactor_df <- subset(norm_sizefactor_df, select = -Gene.name)
  counts_sum <- rowSums(norm_sizefactor_df[, -1])
  norm_sizefactor_df <- norm_sizefactor_df[counts_sum != 0, ]
  write.csv(norm_sizefactor_df, file = file.path(folder_path, paste0("Sizefactornorm_", cell_type, "_DESeq.csv")), row.names = TRUE)


  print("Step 18: Visualisation de l'Analyse en Composantes Principales (ACP) :")
  
  # Visualisation de l'ACP
  vsd <- vst(dss, blind = FALSE)
  pca_data <- plotPCA(vsd, intgroup = "CONDITION", returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = CONDITION)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle("PCA of DESeq2 Normalized Counts") +
    theme_minimal()
  
  ggsave(file.path(folder_path, paste0(cell_type, "_PCA_plot.png")), plot = pca_plot, width = 10, height = 8, bg = "white")
  
  print("Step 19: Création de la matrice des distances de similarité :")
  
  # Création de la matrice des distances de similarité
  dist_matrix <- dist(t(assay(vsd)))
  dist_matrix_df <- as.data.frame(as.matrix(dist_matrix))
  write.csv(dist_matrix_df, file = file.path(folder_path, paste0("DistanceMatrix_", cell_type, "_DESeq.csv")), row.names = TRUE)
  
  print("Step 20: Visualisation du heatmap :")
  
  # Visualisation du heatmap
  heatmap_plot <- pheatmap(dist_matrix, 
                           clustering_distance_rows = dist_matrix, 
                           clustering_distance_cols = dist_matrix, 
                           color = colorRampPalette(rev(brewer.pal(n = 9, name = "Blues")))(255),
                           main = paste("Heatmap of Sample Distances for", cell_type))
  
  # Sauvegarder le heatmap
  save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
    stopifnot(!missing(filename))
    pdf(filename, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  save_pheatmap_pdf(heatmap_plot, file.path(folder_path, paste0("Heatmap_", cell_type, "_DESeq.pdf")))
  
  print("Analyse DESeq2 terminée.")
  
  }
##DESeq_analysis ("meta_matrix", "count_matrix", "mart_export", "CIBERSORTx_Job6_Results.csv", "ctrl_meta_matrix.csv", "ctrl_count_matrix.csv","LB_purity" , "BLymphocytes")
#"T_cell_tot","Lineage_B","Monocytic lineage","Neutrophils", "cybersort_Monocytes