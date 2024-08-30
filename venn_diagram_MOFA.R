# Charger les packages nécessaires
library(dplyr)
library(VennDiagram)
library(readr)
library(grid)

# Fonction principale
venn_diagram_MOFA <- function(path, factor_num, prefix) {
  # Construire les sous-dossiers basés sur les types cellulaires
  cell_types <- c("BLymphocytes", "TLymphocytes", "Monocytes", "Neutrophils")
  
  # Construire le nom des fichiers basés sur les paramètres
  file_names <- paste0(prefix, "_contributors_", factor_num, "_", cell_types, ".csv")
  
  # Lire les fichiers et extraire les colonnes "Contributeur"
  data_list <- lapply(seq_along(cell_types), function(i) {
    file_path <- file.path(path, cell_types[i], file_names[i])
    if (file.exists(file_path)) {
      # Lire le fichier et sélectionner la colonne "Contributeur"
      read_csv(file_path) %>% pull(Contributeur)
    } else {
      stop(paste("File not found:", file_path))
    }
  })
  
  # Nommer les éléments de la liste pour les utiliser dans le diagramme de Venn
  names(data_list) <- cell_types
  
  # Trouver les gènes en commun
  common_genes <- Reduce(intersect, data_list)
  
  # Sauvegarder les gènes en commun dans un fichier CSV
  common_genes_file <- file.path(path, paste0(prefix, "_", factor_num, "_common_genes.csv"))
  write_csv(data.frame(Contributeur = common_genes), common_genes_file)
  
  # Définir les couleurs pour chaque type cellulaire
  venn_colors <- c("BLymphocytes" = "red", "TLymphocytes" = "blue", "Monocytes" = "green", "Neutrophils" = "purple")
  
  # Créer le diagramme de Venn
  venn_plot <- venn.diagram(
    x = data_list,
    category.names = cell_types,
    filename = NULL,
    output = TRUE,
    fill = venn_colors[cell_types],  # Utiliser les couleurs définies
    cex = 2,  # Taille du texte
    cat.cex = 2,  # Taille du texte des catégories
    cat.pos = 0,  # Position des catégories
    cat.dist = 0.05,  # Distance des labels des catégories du diagramme
    resolution = 300,  # Résolution de l'image
    imagetype = "png",  # Type de l'image
    height = 3000,  # Hauteur de l'image
    width = 3000  # Largeur de l'image
  )
  
  # Sauvegarder le diagramme de Venn
  venn_file <- file.path(path, paste0(prefix, "_", factor_num, "_venn.png"))
  png(venn_file, width = 3000, height = 3000, res = 300)
  grid.draw(venn_plot)
  dev.off()
  
  # Créer un tableau récapitulatif
  all_genes <- unique(unlist(data_list))
  summary_table <- sapply(data_list, function(genes) {
    # Compter le nombre d'occurrences de chaque gène
    sapply(all_genes, function(gene) sum(gene == genes))
  })
  
  summary_table <- as.data.frame(summary_table)
  rownames(summary_table) <- all_genes
  summary_table <- cbind(Gene = rownames(summary_table), summary_table)
  
  # Sauvegarder le tableau récapitulatif dans un fichier CSV
  summary_table_file <- file.path(path, paste0(prefix, "_", factor_num, "_summary_table.csv"))
  write_csv(summary_table, summary_table_file)
}

# Exemple d'utilisation
# venn_diagram_MOFA("D:/data_umr_1227/Par_pureté_cyto_digitale/MOFA_analysis_26", 12, "neg")
