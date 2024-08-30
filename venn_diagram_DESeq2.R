library(VennDiagram)
library(dplyr)

# Définir la fonction principale
venn_diagram_DESeq2 <- function(base_path) {
  
  # Créer un dossier pour les fichiers de sortie
  output_dir <- file.path(base_path, "DEG_DESeq2")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Récupérer les sous-dossiers avec l'extension _DESeq
  dossiers <- list.dirs(base_path, recursive = FALSE, full.names = TRUE)
  dossiers_deseq <- dossiers[grep("_DESeq$", dossiers)]
  
  # Initialiser une liste pour stocker les gènes de chaque type cellulaire
  genes_list <- list()
  
  # Parcourir chaque dossier et récupérer les fichiers CSV commençant par DEG_
  for (dossier in dossiers_deseq) {
    csv_files <- list.files(dossier, pattern = "^DEG_.*\\.csv$", full.names = TRUE)
    
    for (csv_file in csv_files) {
      # Lire le fichier CSV
      data <- read.csv(csv_file)
      
      # Filtrer les DEGs en fonction des critères (padj < 0.05 et |log2FoldChange| > 0.5)
      DEGs <- data %>%
        filter(padj < 0.05 & abs(log2FoldChange) > 0.5) %>%
        pull(Gene.name)  # Extraire les noms des gènes
      
      # Ajouter les gènes à la liste
      type_cellulaire <- basename(dossier)
      if (is.null(genes_list[[type_cellulaire]])) {
        genes_list[[type_cellulaire]] <- DEGs
      } else {
        genes_list[[type_cellulaire]] <- unique(c(genes_list[[type_cellulaire]], DEGs))
      }
    }
  }
  
  # Créer le diagramme de Venn pour les 4 types cellulaires avec des couleurs différentes
  venn.plot <- venn.diagram(
    x = genes_list,
    category.names = names(genes_list),
    filename = NULL,
    output = TRUE,
    fill = c("red", "green", "blue", "yellow"), # Différentes couleurs pour chaque type cellulaire
    cex = 2, # Taille de police pour les étiquettes
    cat.cex = 1, # Taille de police pour les catégories
    cat.pos = 0, # Position des catégories
    margin = 0.1, # Marge autour du diagramme
    main = "Venn Diagram of DEGs Across Cell Types"
  )
  
  # Sauvegarder le diagramme de Venn dans un fichier PNG
  venn_file <- file.path(output_dir, "venn_diagram_DEGs.png")
  png(filename = venn_file, width = 2000, height = 2000, res = 300) # Taille augmentée
  grid.draw(venn.plot)
  dev.off()
  
  # Calculer et afficher les intersections
  intersections <- calculate.overlap(genes_list)
  
  # Sauvegarder les intersections dans un fichier texte
  intersections_file <- file.path(output_dir, "intersections_DEGs.txt")
  sink(intersections_file)
  print(intersections)
  sink()
  
  # Identifier et sauvegarder les gènes uniques pour chaque type cellulaire
  for (type_cellulaire in names(genes_list)) {
    unique_genes <- setdiff(genes_list[[type_cellulaire]], unlist(genes_list[names(genes_list) != type_cellulaire]))
    unique_genes_df <- data.frame(Gene.name = unique_genes)
    
    # Sauvegarder les gènes uniques dans un nouveau fichier CSV
    unique_genes_file <- file.path(output_dir, paste0(type_cellulaire, "_unique_DEGs.csv"))
    write.csv(unique_genes_df, file = unique_genes_file, row.names = FALSE)
  }
  
  #return(intersections)
}

# Appel de la fonction avec le chemin de base
#base_path <- "D:/data_umr_1227/Par_pureté_cyto_digitale"  # Remplacez par votre chemin de base
#intersections <- venn_diagram_DESeq2(base_path)
