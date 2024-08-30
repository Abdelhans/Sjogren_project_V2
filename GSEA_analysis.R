setwd("D:/data_umr_1227/Par_pureté_cyto_digitale")

library(dplyr)
library(ggplot2)
library(fgsea)
library(readxl)
library(ggrepel)
library(stringr)
library(data.table)
library(pheatmap)

GSEA_analysis <- function(folder_path) {
  
  # Déclaration de la fonction: Adjacency matrix to list 
  matrix_to_list <- function(pws){
    pws.l <- list()
    for (pw in colnames(pws)) {
      pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
    }
    return(pws.l)
  }
  
  # Déclaration de la fonction: prepare_gmt 
  prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
    # Read in gmt file
    gmt <- gmtPathways(gmt_file)
    hidden <- unique(unlist(gmt))
    
    # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
    mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                  nrow = length(hidden), ncol = length(gmt))
    for (i in 1:dim(mat)[2]){
      mat[,i] <- as.numeric(hidden %in% gmt[[i]])
    }
    
    #Subset to the genes that are present in our data to avoid bias
    hidden1 <- intersect(genes_in_data, hidden)
    mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
    # And get the list again
    final_list <- matrix_to_list(mat) # for this we use the function we previously defined
    
    if(savefile){
      saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
    }
    
    #print('Wohoo! .gmt conversion successfull!:)')
    return(final_list)
  }
  
  cat("Etape 1: Extraction des noms de fichiers et création du nouveau dossier...\n")
  # Extraction des noms de fichiers et création du nouveau dossier
  dossier_nom <- basename(folder_path)
  nom_sans_DESeq <- gsub("_DESeq.*", "", dossier_nom)
  nouveau_dossier <- paste0(nom_sans_DESeq, "_GSEA")
  chemin_nouveau_dossier <- file.path(dirname(folder_path), nouveau_dossier)
  if (!file.exists(chemin_nouveau_dossier)) {
    dir.create(chemin_nouveau_dossier)
  }
  output_file <- file.path(chemin_nouveau_dossier, "output_log.txt")
  # Ouvrir le fichier de sortie pour rediriger les sorties
  sink(output_file, append = TRUE)
  
  cat("Etape 2: Chargement des résultats DESeq...\n")
  # Chargement des résultats DESeq
  deseq_output <- read.csv(list.files(folder_path, pattern="^resultat_", full.names=TRUE))
  DGE <- read.csv(list.files(folder_path, pattern="^DEG_", full.names=TRUE), row.names = 1)
  normalized_data <- read.csv(list.files(folder_path, pattern="^VST_", full.names=TRUE), check.names = FALSE,row.names = 1)
  normalized_data <- normalized_data[rownames(normalized_data)%in%DGE$Gene.stable.ID,]
  condition_data <- read.csv(list.files(folder_path, pattern="^colData_", full.names=TRUE), check.names = FALSE, row.names = 1)
  # Extraire les nouveaux noms de colonnes depuis condition_data
  new_colnames <- condition_data$CONDITION
  # Remplacement des noms de colonnes dans normalized_data
  colnames(normalized_data) <- new_colnames
  normalized_data$Gene.stable.ID <- rownames(normalized_data)
  # Extraire les colonnes 'Gene.stable.ID' et 'Gene.name' de DGE
  DGE_reduit <- DGE[, c("Gene.stable.ID", "Gene.name")]
  
  # Fusionner DGE_reduit avec normalized_data
  normalized_data <- merge(DGE_reduit, normalized_data, by = "Gene.stable.ID")
  
  # Identifier les noms de gènes en double
  dup_indices <- which(duplicated(normalized_data$Gene.name) | duplicated(normalized_data$Gene.name, fromLast = TRUE))
  
  # Remplacer les valeurs en double par les valeurs correspondantes de 'Gene.stable.ID'
  normalized_data$Gene.name[dup_indices] <- normalized_data$Gene.stable.ID[dup_indices]
  
  # Définir les noms de ligne de normalized_data avec la colonne 'Gene.name'
  rownames(normalized_data) <- normalized_data$Gene.name
  
  # Supprimer les colonnes 'Gene.stable.ID' et 'Gene.name'
  normalized_data <- subset(normalized_data, select = -c(Gene.stable.ID, Gene.name))
  
  # Calculer les z-scores par ligne dans normalized_data
  z_scores <- apply(normalized_data, 1, function(x) (x - mean(x)) / sd(x))
  z_scores <- t(z_scores)
  
  write.csv(z_scores, file = file.path(chemin_nouveau_dossier, "zscore.csv"), quote = FALSE)
  
  cat("Etape 3: Ajout de la colonne 'diffexpressed'...\n")
  # Ajout de la colonne diffexpressed
  deseq_output <- deseq_output %>% mutate(diffexpressed = case_when(
    log2FoldChange > 0.5 & padj < 0.05 ~ 'UP',
    log2FoldChange < -0.5 & padj < 0.05 ~ 'DOWN',
    padj > 0.05 ~ 'NO'
  ))
  
  cat("Etape 4: Récupération des fichiers .gmt...\n")
  # Récupération des fichiers .gmt
  gmt_files <- list.files(pattern = '.gmt', full.names = TRUE)
  
  cat("Etape 5: Appel de la fonction prepare_gmt...\n")
  # Appel de la fonction prepare_gmt avec le premier fichier .gmt trouvé
  genes <- prepare_gmt(gmt_files[1], deseq_output$Gene.name, savefile = FALSE)
  
  cat("Etape 6: Calcul des classements des rangs des gènes...\n")
  # Calcul des classements des rangs des gènes
  rankings <- sign(deseq_output$log2FoldChange) * (-log10(deseq_output$pvalue))
  names(rankings) <- deseq_output$Gene.name
  rankings <- sort(rankings, decreasing = TRUE)
  
  cat("Etape 7: Exécution de fgsea...\n")
  # Exécution de fgsea
  GSEAresult <- fgsea(pathways = genes,
                      stats = rankings,
                      scoreType = 'std',
                      minSize = 10,
                      maxSize = 500,
                      nproc = 1)
  GSEAresult_sorted <- GSEAresult %>% arrange(padj)  
  pval=0.05
  
  msg <- GSEAresult_sorted %>% dplyr::filter(padj < !!pval) %>% arrange(desc(NES))
  
  message(paste("Number of signficant gene sets =", nrow(msg)))
  
  cat("Etape 8: Écriture des résultats bruts dans un fichier...\n")
  
  # Ajout de la colonne des valeurs absolues de NES
  GSEAresult_sorted$absNES <- abs(GSEAresult_sorted$NES)
  
  # Écriture des résultats bruts dans un fichier
  fwrite(GSEAresult_sorted, file = file.path(chemin_nouveau_dossier, "result_all_GSEA.csv"),  quote = FALSE)
  
  cat("Etape 9: Écriture des résultats filtrés suivant le PADJ..\n")
  # Écriture des résultats filtrés suivant le PADJ
  fwrite(GSEAresult_sorted[GSEAresult_sorted$padj < 0.05], file = file.path(chemin_nouveau_dossier, "GSEA_DEG.csv"),  quote = FALSE)
  
  #Écriture des résultats filtrés suivant le NES
  
  NES_HIGH <- GSEAresult_sorted[GSEAresult_sorted$absNES > 1.5, ]
  fwrite(NES_HIGH, file = file.path(chemin_nouveau_dossier, "GSEA_HIGH_NES.csv"), quote = FALSE)
  
  
  cat("Etape 10: Parcours des voies biologiques significativement enrichies...\n")
  
  # Parcours des voies biologiques significativement enrichies (pour les 4 premières voies)
  for (i in 1:min(4, nrow(GSEAresult_sorted))) {
    titre <- GSEAresult_sorted[i, "pathway"]
    titre <- as.character(titre)  # Assurez-vous que 'titre' est une chaîne de caractères
    
    # Génération des graphiques d'enrichissement
    plot <- plotEnrichment(genes[[titre]], rankings) + 
      labs(title = titre)
    chemin_fichier <- file.path(chemin_nouveau_dossier, paste0("plot_", i, "_", str_replace_all(titre, "/", "_"), ".png"))
    ggsave(chemin_fichier, plot, width = 10, height = 5)  # Taille arbitraire
    
    # Heatmap
    genes_dans_la_voie <- unlist(GSEAresult_sorted[i, "leadingEdge"])
    genes_a_afficher <- intersect(rownames(z_scores), genes_dans_la_voie)
    
    if (length(genes_a_afficher) > 1) {
      sous_ensemble_expression <- z_scores[genes_a_afficher, ]
      
      annotation_col <- data.frame(
        Group = ifelse(grepl("^CTRL", colnames(sous_ensemble_expression)), "CTRL", 
                       ifelse(grepl("^SJS", colnames(sous_ensemble_expression)), "SJS", NA))
      )
      rownames(annotation_col) <- colnames(sous_ensemble_expression)
      
      # Réorganiser les colonnes par ordre alphabétique des groupes
      ordered_cols <- order(annotation_col$Group)
      sous_ensemble_expression <- sous_ensemble_expression[, ordered_cols]
      annotation_col <- annotation_col[ordered_cols, , drop = FALSE]
      
      chemin_fichier_heatmap <- file.path(chemin_nouveau_dossier, paste0("heatmap_", i, "_", gsub("/", "_", titre), ".png"))
      
      # Vérifier les données avant de créer la heatmap
      if (ncol(sous_ensemble_expression) > 1 && nrow(sous_ensemble_expression) > 1) {
        pheatmap(sous_ensemble_expression, scale = "row", clustering_distance_rows = "euclidean", 
                 clustering_distance_cols = "euclidean", show_rownames = TRUE, show_colnames = FALSE,
                 main = titre, annotation_col = annotation_col, 
                 legend = TRUE, legend_labels = c("Low", "High"), legend_title = "Z-scores")
        dev.copy(png, filename = chemin_fichier_heatmap, width = 1000, height = 500)  # Sauvegarder l'image
        dev.off()  # Fermer le périphérique graphique
      } else {
        cat("Pas assez de données pour générer la heatmap pour la voie:", titre, "\n")
      }
    } else {
      cat("Pas assez de gènes à afficher pour la heatmap pour la voie:", titre, "\n")
    }
  }
  
  
  cat("Etape 11: visualisation des différentes voies enrichies...\n")
  
  top_terms <- 10
  GSEAresult_sorted_top <- GSEAresult_sorted %>%
    arrange(desc(absNES)) %>%
    head(top_terms)
  
  # Troncature des noms de termes pour un affichage clair
  GSEAresult_sorted_top$pathway <- sapply(GSEAresult_sorted_top$pathway, function(x) {
    if(nchar(x) > 40) {
      return(paste0(substr(x, 1, 37), "..."))
    } else {
      return(x)
    }
  })
  
  # Créer le graphique
  plot <- ggplot(GSEAresult_sorted_top, aes(x = reorder(pathway, -padj), y = NES, fill = size)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Top Enriched Pathways",
         x = "Pathway",
         y = "Normalized Enrichment Score",
         fill = "Gene Count") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA), # Fond blanc pour le graphique
          panel.background = element_rect(fill = "white", color = NA)) # Fond blanc pour le panneau
  
  # Chemin pour sauvegarder le graphique
  chemin_fichier_barplot <- file.path(chemin_nouveau_dossier, "enriched_pathways_barplot.png")
  
  # Sauvegarde du graphique
  ggsave(chemin_fichier_barplot, plot, width = 10, height = 5)  # Taille arbitraire
  
  # Affichage du chemin où le fichier est sauvegardé
  cat("Le graphique a été sauvegardé dans le fichier :", chemin_fichier_barplot, "\n")
  
  cat("Etape 12: voies avec le plus grand NES ...\n")
  
  # Afficher les noms des 2 premières voies avec le plus grand NES
  top_voies <- GSEAresult_sorted[order(-GSEAresult_sorted$absNES), ]
  top_voies <- top_voies[1:4, ]  # Sélectionner les deux premières voies
  
  cat(top_voies$pathway, "\n\n")
  
  # Déterminer si elles sont régulées positivement ou négativement
  for (i in 1:nrow(top_voies)) {
    titre <- as.character(top_voies[i, "pathway"])
    NES <- top_voies[i, "NES"]
    
    if (NES > 0) {
      cat(paste("Voie:", titre, "- Régulée positivement (NES =", NES, ")\n"))
    } else {
      cat(paste("Voie:", titre, "- Régulée négativement (NES =", NES, ")\n"))
    }
  }
}
  
cat("Toutes les étapes ont été complétées avec succès!\n")
# Utilisation de la fonction avec le chemin du dossier
#GSEA_analysis("D:/data_umr_1227/Par_pureté_cyto_digitale/BLymphocytes_DESeq")
