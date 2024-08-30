setwd("D:/data_umr_1227/Par_pureté_cyto_digitale")

# Chargement des librairies nécessaires
library(enrichR)
library(ggplot2)
library(dplyr)

# Fonction d'analyse EnrichR
EnrichR_analysis <- function(folder_path) {
  print("Étape 1 : Chargement et configuration de enrichR...")
  # Charger enrichR
  # Listes des bases de données enrichR disponibles
  listEnrichrSites()
  setEnrichrSite("Enrichr")
  db_enrichR <- c("GO_Biological_Process_2023")
  
  # Extraire le nom du dossier d'entrée
  dossier_nom <- basename(folder_path)
  
  # Extraire la première partie du nom avant "_DESeq"
  nom_sans_DESeq <- gsub("_DESeq.*", "", dossier_nom)
  
  # Créer le nom du nouveau dossier
  nouveau_dossier <- paste0(nom_sans_DESeq, "_Enrichir")
  
  # Chemin complet du nouveau dossier
  chemin_nouveau_dossier <- file.path(dirname(folder_path), nouveau_dossier)
  
  # Créer le nouveau dossier s'il n'existe pas déjà
  if (!file.exists(chemin_nouveau_dossier)) {
    dir.create(chemin_nouveau_dossier)
  }
  
  print("Étape 2 : Lecture des fichiers du dossier...")
  # Création des objets enrichR à tester
  DEG <- read.csv(list.files(folder_path, pattern="^DEG_", full.names=TRUE))
  UP <- read.csv(list.files(folder_path, pattern="^UP_", full.names=TRUE))
  DOWN <- read.csv(list.files(folder_path, pattern="^DOWN_", full.names=TRUE))
  
  # Création des objets enrichR à tester
  print("Étape 3 : Création des objets enrichR à tester...")
  # Processus pour enrichr_DEG
  tryCatch({
    enrichr_DEG <- enrichr(DEG$Gene.name, db_enrichR)
    enrichr_DEG <- do.call(rbind, enrichr_DEG)
    
    # Sélection des top termes enrichis
    top_terms <- 10
    enrichr_DEG_sorted <- enrichr_DEG %>%
      arrange(Adjusted.P.value) %>%
      head(top_terms)
    
    # Troncature des noms de termes pour un affichage clair
    enrichr_DEG_sorted$Term <- sapply(enrichr_DEG_sorted$Term, function(x) substr(x, 1, 40))
    
    # Création du graphique pour DEG
    plot_1_DEG <- ggplot(enrichr_DEG_sorted, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value), fill = Combined.Score)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = "Top Enriched Terms for DEG",
           x = "Term",
           y = "-log10 Adjusted p-value",
           fill = "Combined Score") +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA))
  }, error = function(e) {
    print(paste("Erreur lors de l'exécution de enrichr_DEG:", e))
    plot_1_DEG <- NULL
  })
  
  # Processus pour enrichr_UP
  tryCatch({
    enrichr_UP <- enrichr(UP$Gene.name, db_enrichR)
    enrichr_UP <- do.call(rbind, enrichr_UP)
    
    # Sélection des top termes enrichis
    top_terms <- 10
    enrichr_UP_sorted <- enrichr_UP %>%
      arrange(Adjusted.P.value) %>%
      head(top_terms)
    
    # Troncature des noms de termes pour un affichage clair
    enrichr_UP_sorted$Term <- sapply(enrichr_UP_sorted$Term, function(x) substr(x, 1, 40))
    
    # Création du graphique pour UP
    plot_2_UP <- ggplot(enrichr_UP_sorted, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value), fill = Combined.Score)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = "Top Enriched Terms for UP",
           x = "Term",
           y = "-log10 Adjusted p-value",
           fill = "Combined Score") +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA))
  }, error = function(e) {
    print(paste("Erreur lors de l'exécution de enrichr_UP:", e))
    plot_2_UP <- NULL
  })
  
  # Processus pour enrichr_DOWN
  tryCatch({
    enrichr_DOWN <- enrichr(DOWN$Gene.name, db_enrichR)
    enrichr_DOWN <- do.call(rbind, enrichr_DOWN)
    
    # Sélection des top termes enrichis
    top_terms <- 10
    enrichr_DOWN_sorted <- enrichr_DOWN %>%
      arrange(Adjusted.P.value) %>%
      head(top_terms)
    
    # Troncature des noms de termes pour un affichage clair
    enrichr_DOWN_sorted$Term <- sapply(enrichr_DOWN_sorted$Term, function(x) substr(x, 1, 40))
    
    # Création du graphique pour DOWN
    plot_3_DOWN <- ggplot(enrichr_DOWN_sorted, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value), fill = Combined.Score)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = "Top Enriched Terms for DOWN",
           x = "Term",
           y = "-log10 Adjusted p-value",
           fill = "Combined Score") +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA))
  }, error = function(e) {
    print(paste("Erreur lors de l'exécution de enrichr_DOWN:", e))
    plot_3_DOWN <- NULL
  })
  
  print("Étape 4 : Enregistrement des résultats...")
  # Étape 5 : Enregistrement des graphiques
  tryCatch({
    if (!is.null(plot_1_DEG)) {
      ggsave(file.path(chemin_nouveau_dossier, "EnrichR_DEG_BP_2023.png"), plot = plot_1_DEG, width = 10, height = 5)
    }
    if (!is.null(plot_2_UP)) {
      ggsave(file.path(chemin_nouveau_dossier, "enrichr_UP_BP_2023.png"), plot = plot_2_UP, width = 10, height = 5)
    }
    if (!is.null(plot_3_DOWN)) {
      ggsave(file.path(chemin_nouveau_dossier, "enrichr_DOWN_BP_2023.png"), plot = plot_3_DOWN, width = 10, height = 5)
    }
  }, error = function(e) {
    print(paste("Erreur lors de l'enregistrement des graphiques:", e))
  })
  
  print("Étape 6 : Sauvegarde des objets enrichR en CSV...")
  tryCatch({
    write.csv(enrichr_DEG, file = file.path(chemin_nouveau_dossier, "EnrichR_DEG.csv"), row.names = FALSE, quote = FALSE)
    write.csv(enrichr_UP, file = file.path(chemin_nouveau_dossier, "Enrichr_UP.csv"), row.names = FALSE, quote = FALSE)
    write.csv(enrichr_DOWN, file = file.path(chemin_nouveau_dossier, "Enrichr_DOWN.csv"), row.names = FALSE, quote = FALSE)
  }, error = function(e) {
    print(paste("Erreur lors de la sauvegarde des objets enrichR en CSV:", e))
  })
  
  print("Toutes les sorties ont été enregistrées avec succès dans le nouveau dossier.")
}

# Exemple d'utilisation
# EnrichR_analysis("D:/data_umr_1227/Par_pureté_cyto_digitale/BLymphocytes_DESeq")
