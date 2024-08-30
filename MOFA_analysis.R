MOFA_analysis <- function(folder_path, number_factors_to_start) {
  cat("Step 1 : Chargement des packages necessaires ...\n")
  
  library(reshape2)
  library(DESeq2)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(MOFA2)
  library(reticulate)
  library(enrichR)
  library(vroom)
  
  cat("Step 2 : Chargements des différentes fichiers necessaire pour l'analyse  ...\n")
  
  count_data_BLymphocytes <- read.csv(file.path(folder_path, "BLymphocytes_DESeq", "VST_normalized_BLymphocytes_DESeq.csv"), check.names = FALSE, row.names = 1)

  count_data_TLymphocytes <- read.csv(file.path(folder_path, "TLymphocytes_DESeq", "VST_normalized_TLymphocytes_DESeq.csv"), check.names = FALSE, row.names = 1)

  count_data_Monocytes <- read.csv(file.path(folder_path, "Monocytes_DESeq", "VST_normalized_Monocytes_DESeq.csv"), check.names = FALSE, row.names = 1)

  count_data_Neutrophils <- read.csv(file.path(folder_path, "Neutrophils_DESeq","VST_normalized_Neutrophils_DESeq.csv" ), check.names = FALSE, row.names = 1)

  meta_data <- vroom(file.path(folder_path, "METADATA_SORTED_CELL_MCPcounter.tsv"), delim = '\t')

  data_clinical_num <- vroom(file.path(folder_path, "CLINICAL_NUMERIC_STRD_SJS.tsv"), delim = '\t')
  data_clinical_num_signif <- read.csv(file.path(folder_path, "data_clinical_num_significatif.csv"), check.names = FALSE, row.names = 1)

  data_clinical_bin <- vroom(file.path(folder_path, "CLINICAL_BINARY_STRD_SJS.tsv"), delim = '\t')
  
  cat("Step 3 : preparation des fichiers pour obtenir des tableau pour chaque type cellulaire  ...\n")
  
  # Obtenir les valeurs uniques de la colonne sample_type
  types_cellulaires <- unique(meta_data$SAMPLE_TYPE)
  # Créer une liste pour stocker les tableaux séparés
  liste_tables <- list()
  # Boucle à travers chaque type cellulaire
  for (type_cellulaire in types_cellulaires) {
    # Filtrer les données pour le type cellulaire actuel
    table_separee <- meta_data %>%
      filter(SAMPLE_TYPE == type_cellulaire & CONDITION %in% c("SJS", "CTRL")) %>%
      select(ID, CONDITION)
    # Définir les noms de lignes basés sur les valeurs de la colonne ID
    rownames(table_separee) <- table_separee$ID
    # mise en facteur des variables pour l'analyse DESeq2
    
    table_separee$CONDITION <- factor(table_separee$CONDITION)
    table_separee$ID <- factor(table_separee$ID)
    # Assigner la table séparée à un nom de variable correspondant au type cellulaire
    nom_variable <- paste0("meta_", gsub(" ", "", type_cellulaire)) # Supprimer les espaces pour créer un nom de variable valide
    assign(nom_variable, table_separee)
    # Ajouter la table séparée à la liste
    liste_tables[[nom_variable]] <- table_separee
  }

  cat("Step 4 : calcule de la variance de chaque genes et selection des 5000  gènes avec la plus grande  ...\n")
  
  types_cellulaires <- c("BLymphocytes", "TLymphocytes", "Monocytes", "Neutrophils")
  count_data_all <- list(BLymphocytes = count_data_BLymphocytes,
                         TLymphocytes = count_data_TLymphocytes,
                         Monocytes = count_data_Monocytes,
                         Neutrophils = count_data_Neutrophils)
  
  # Listes pour stocker les résultats
  liste_top_5000_genes <- list()
  liste_normalized_counts <- list()
  
  # Boucle à travers chaque type cellulaire
  for (type_cellulaire in types_cellulaires) {
    
    # Calculer la variance pour chaque gène
    variance_genes <- apply(count_data_all[[type_cellulaire]], 1, var)
    
    # Trier les gènes par ordre décroissant de leur variance
    genes_sorted <- names(sort(variance_genes, decreasing = TRUE))
    
    # Sélectionner les 5000 premiers gènes avec la variance la plus élevée
    top_5000_genes <- genes_sorted[1:5000]
    
    # Ajouter les 5000 premiers gènes à la liste
    liste_top_5000_genes[[paste0("top_5000_genes_", type_cellulaire)]] <- top_5000_genes
    
    # Créer un nouveau data frame avec uniquement les gènes sélectionnés
    normalized_count <- count_data_all[[type_cellulaire]][top_5000_genes, ]
    
    # Ajouter le nouveau data frame à la liste
    nom_variable_normalized_count <- paste0("count_reduce_", type_cellulaire)
    assign(nom_variable_normalized_count, normalized_count)
    liste_normalized_counts[[nom_variable_normalized_count]] <- normalized_count
  }
  cat("Step 5 : Modifications des noms des colonnes des tableaux des comptes normalisées  ...\n")
  
  for (type_cellulaire in types_cellulaires) {
    nom_variable_normalized_count <- paste0("count_reduce_", gsub(" ", "", type_cellulaire))
    if (exists(nom_variable_normalized_count)) {
      normalized_count <- get(nom_variable_normalized_count)
      colnames(normalized_count) <- gsub("^(.*)_.*", "\\1", colnames(normalized_count))
      assign(nom_variable_normalized_count, normalized_count)
    }
  }
  
  cat("Step 6 : Modifications des noms des lignes des tableaux de nos metadonnées   ...\n")
  
  for (type_cellulaire in types_cellulaires) {
    nom_variable_meta <- paste0("meta_", gsub(" ", "", type_cellulaire))
    if (exists(nom_variable_meta)) {
      meta <- get(nom_variable_meta)
      rownames(meta) <- gsub("^(.*)_.*", "\\1", rownames(meta))
      meta$ID <- gsub("^(.*)_.*", "\\1", meta$ID)
      assign(nom_variable_meta, meta)
    }
  }
  
  cat("Step 7 : chargement de la fonction permettant de transformer les donner en tableau long format   ...\n")
  
  transformer_etapes <- function(normalized_count, meta, type_cellulaire) {
    # Ajouter la colonne feature
    normalized_count$feature <- rownames(normalized_count)
    # Convertir en format long
    Long_format <- melt(normalized_count, id.vars = "feature", variable.name = "sample", value.name = "value")
    # Ajouter la vue
    Long_format$view <- type_cellulaire
    # Convertir la colonne value en numérique
    Long_format$value <- as.numeric(Long_format$value)
    # Fusionner avec les métadonnées
    Long_format <- merge(Long_format, meta, by.x = "sample", by.y = "ID")
    # Réinitialiser les noms de lignes
    rownames(Long_format) <- NULL
    
    return(Long_format)}
  
  
  cat("Step 8 : chargement de la fonction permettant de transformer les donner en tableau long format et fusions des 4 tableaux  ...\n")
  
  # Créer une liste pour stocker les résultats de chaque type cellulaire
  liste_long_format <- list()
  
  # Boucle à travers chaque type cellulaire
  for (type_cellulaire in types_cellulaires) {
    nom_variable_normalized_count <- paste0("count_reduce_", gsub(" ", "", type_cellulaire))
    nom_variable_meta <- paste0("meta_", gsub(" ", "", type_cellulaire))
    
    # Vérifier si les objets existent
    if (exists(nom_variable_normalized_count) && exists(nom_variable_meta)) {
      normalized_count <- get(nom_variable_normalized_count)
      meta <- get(nom_variable_meta)
      
      # Appliquer les étapes de transformation
      long_format <- transformer_etapes(normalized_count, meta, type_cellulaire)
      
      # Ajouter le résultat à la liste
      liste_long_format[[paste0("Long_format_", gsub(" ", "", type_cellulaire))]] <- long_format
    }
  }
  
  # Combiner tous les résultats en un seul dataframe
  Long_format_all <- do.call(rbind, liste_long_format)
  
  cat("Step 9.0 : Ajouts des données cliniques aux  tableaux de nos metadonnées ...\n")
  
  rename_columns <- function(data) {
    new_colnames <- sapply(colnames(data), function(col) {
      # Supprimer les caractères spéciaux des noms de colonnes
      col <- gsub("[^A-Za-z0-9_]", "", col)
      # Récupérer les 10 premières lettres de chaque nom de colonne
      new_name <- substr(col, 1, 10)
      return(new_name)
    })
    # Attribuer les nouveaux noms de colonnes
    colnames(data) <- new_colnames
    return(data)
  }
  
  data_clinical_bin  <- rename_columns(data_clinical_bin)
  # Transformation en donnée numeric 
  data_clinic_bin  <- as.data.frame(lapply(data_clinical_bin, function(x) {
    x <- as.character(x)
    x[x == "No"] <- 0
    x[x == "Yes"] <- 1
    x[x == "Unknown"] <- 2
    return(as.numeric(as.character(x)))
  }))
  
  
  meta_all_celltype <- bind_rows(meta_BLymphocytes, meta_TLymphocytes,meta_Monocytes, meta_Neutrophils)
  meta_all_celltype_unique <- meta_all_celltype %>% distinct() %>% arrange(desc(ID)) 
  data_clinic_bin <- data_clinic_bin[data_clinic_bin$ID %in% meta_all_celltype_unique$ID,]
  
  Long_format_all <- Long_format_all %>% mutate(sample = as.numeric(as.character(sample)))
  
  # noms des colonnes sont dupliqués, utilisez make.unique pour les rendre uniques
  colnames(data_clinical_num_signif) <- make.unique(colnames(data_clinical_num_signif))

  Long_format_all <- Long_format_all %>% left_join(data_clinical_num_signif, by = c("sample" = "ID"))
    
  cat("Step 9.1 : creation de l'objet MOFA et apercu des données  ...\n")
  
  # Créer un objet MOFA
  MOFAobject <- create_mofa(Long_format_all)
  
  # Générer le plot de l'aperçu des données MOFA
  apercu_des_données <- plot_data_overview(MOFAobject)
  
  cat("Step 10 : Initialisation des options pour le modèle MOFA...\n")
  
  data_opts <- get_default_data_options(MOFAobject)
  
  # Initialisation des options du modèle avec un seul facteur
  num_factors <- number_factors_to_start
  result <- FALSE
  
  # Boucle pour incrémenter le nombre de facteurs jusqu'à ce que la condition soit vraie
  while (!result) {
    cat(paste0("Running MOFA with ", num_factors, " factors...\n"))
    
    # Définir le nombre de facteurs
    model_opts <- get_default_model_options(MOFAobject)
    model_opts$num_factors <- num_factors
    
    train_opts <- get_default_training_options(MOFAobject)
    train_opts$convergence_mode <- "fast"
    train_opts$maxiter <- 5000
    train_opts$seed <- 42
    
    # Préparation du modèle
    MOFAobject <- prepare_mofa(object = MOFAobject, data_options = data_opts, model_options = model_opts, training_options = train_opts)
    
    cat("Step 11 : Création du dossier de sauvegarde des résultats pour le modèle...\n")
    
    nouveau_dossier <- file.path(folder_path, paste0("MOFA_analysis_", num_factors))
    if (!file.exists(nouveau_dossier)) {
      dir.create(nouveau_dossier)
    }
    
    cat("Step 12 : Run du modèle MOFA...\n")
    
    reticulate::use_condaenv("C:/Users/utilisateur/miniconda3")
    MOFAmodel <- run_mofa(MOFAobject, use_basilisk = FALSE, outfile = file.path(nouveau_dossier, "modelmofa.hdf5"))
    # Assigner le modèle MOFA dans l'environnement global
    assign("MOFAmodel", MOFAmodel, envir = .GlobalEnv)
    
    cat("Step 13 : Vérification de la variance expliquée par le modèle...\n")
    
    var_explained <- as.data.frame(get_variance_explained(MOFAmodel)$r2_per_factor[[1]])
    var_explained$Total <- rowSums(var_explained)
    
    # Vérification si la dernière ligne est inférieure à 1%
    last_row_total <- var_explained[nrow(var_explained), "Total"]
    result <- last_row_total < 1
    
    if (!result) {
      cat(paste0("Variance expliquée par ", num_factors, " facteurs est supérieure à 1%.\n"))
      num_factors <- num_factors + 1
    } else {
      cat(paste0("Variance expliquée par ", num_factors, " facteurs est inférieure à 1%. Arrêt du processus.\n"))
    }
    
    cat("Step 14 : Plot de la variance expliquée, de la matrice de corrélation et de l'estimation de la variance...\n")
    
    # Plot de la variance expliquée
    variance_plot <- plot_variance_explained(MOFAmodel, x = "view", y = "factor", plot_total = TRUE)[[2]]
    ggsave(filename = file.path(nouveau_dossier, "variance_explained_plot.png"), plot = variance_plot, device = "png", width = 12, height = 8)
    
    correlation_matrix <- file.path(nouveau_dossier, "correlation_matrix.png")
    png(filename = correlation_matrix, width = 4000, height = 3000, units = "px", res = 800)
    plot_factor_cor(MOFAmodel)
    dev.off()
    
    estimations_de_variance <- plot_variance_explained(MOFAmodel)
    ggsave(filename = file.path(nouveau_dossier, "estimations_de_variance.png"), plot = estimations_de_variance, device = "png", width = 12, height = 8)
    ggsave(filename = file.path(nouveau_dossier, "data_overview_plot.png"), plot = apercu_des_données, device = "png", width = 12, height = 8)
    
    cat("Step 15 : Obtention des poids des facteurs pour separation en fonctions de nos groupes...\n")
    
    distribution_facto <- as.data.frame(get_expectations(MOFAmodel, "Z")$single_group)
    distribution_facto$ID <- rownames(distribution_facto)  
    
    # Filtrer les données en fonction des conditions et conserver uniquement les colonnes 
    meta_data_subset <- meta_data %>% filter(CONDITION %in% c("SJS", "CTRL")) %>% select(ID, CONDITION)
    meta_data_subset$ID <- gsub("^(.*)_.*", "\\1", meta_data_subset$ID)
    meta_data_subset <- distinct(meta_data_subset, ID, .keep_all = TRUE)
    distribution_facto <- merge(distribution_facto, meta_data_subset, by = "ID")
    
    # Séparer les groupes
    CTRL <- distribution_facto[distribution_facto$CONDITION == 'CTRL', ]
    SJS  <- distribution_facto[distribution_facto$CONDITION == 'SJS', ]
    
    # Obtenir les noms des facteurs
    factors <- colnames(distribution_facto)[grep("Factor", colnames(distribution_facto))]
    
    # Initialiser un dataframe pour stocker les résultats
    results_df <- data.frame(Facteur = character(),
                             Pvalue = numeric(),
                             Padjusted = numeric(),
                             Significativite = character(),
                             Mean = numeric(),
                             stringsAsFactors = FALSE)
    
    cat("Step 16 : Test de Student entre les deux groupes...\n")
    
    # Effectuer le test t pour chaque facteur
    results_df <- data.frame()
    for (factor_name in factors) {
      ttest <- t.test(CTRL[[factor_name]], SJS[[factor_name]])
      significant <- ifelse(p.adjust(ttest$p.value, method = "BH") < 0.05, "Significatif", "Non significatif")
      mean_value <- mean(distribution_facto[[factor_name]])
      results_df <- rbind(results_df, data.frame(Facteur = factor_name,
                                                 Pvalue = ttest$p.value,
                                                 Padjusted = p.adjust(ttest$p.value, method = "BH"),
                                                 Significativite = significant,
                                                 Mean = mean_value))
    }
    
    # Sauvegarder les résultats dans un fichier CSV
    write.csv(results_df, file = file.path(nouveau_dossier, "Student_t_test.csv"), row.names = FALSE)
    write.csv(distribution_facto, file = file.path(nouveau_dossier, "distribution_tableau.csv"), row.names = FALSE)
    
    tryCatch({
      # Vérifier s'il y a des facteurs significatifs
      facteurs_significatifs <- as.numeric(gsub("[^0-9]", "", results_df$Facteur[results_df$Significativite == "Significatif"]))
      
      # Assigner les facteurs significatifs dans l'environnement global
      assign("facteurs_significatifs", facteurs_significatifs, envir = .GlobalEnv)
      
      if (length(facteurs_significatifs) == 0) {
        cat("Pas de facteur significatif trouvé.\n")
      } else {
        cat("Step 17 : Plots des facteurs significatifs ...\n")
        
        # Vérifier si MOFAmodel est valide avant de l'utiliser
        if (exists("MOFAmodel", envir = .GlobalEnv) && !is.null(get("MOFAmodel", envir = .GlobalEnv))) {
          cat("Lancement des tracés des facteurs significatifs...\n")
          
          # Récupérer MOFAmodel de l'environnement global
          MOFAmodel <- get("MOFAmodel", envir = .GlobalEnv)
          
          # Tracer les facteurs significatifs
          Scatterplots_significatif <- plot_factors(MOFAmodel, factor = facteurs_significatifs, color_by = "CONDITION")
          ggsave(filename = file.path(nouveau_dossier, "Scatterplotsignif.png"), plot = Scatterplots_significatif, device = "png", width = 12, height = 8)
          
          Scatterplots <- plot_factor(MOFAmodel, factor = facteurs_significatifs, color_by = "CONDITION", add_violin = TRUE, dodge = TRUE, violin_alpha = 0.25)
          ggsave(filename = file.path(nouveau_dossier, "Scatterplots.png"), plot = Scatterplots, device = "png", width = 12, height = 8)
          
          cat("Les tracés ont été sauvegardés avec succès.\n")
        } else {
          cat("L'objet 'MOFAmodel' est NULL ou introuvable dans l'environnement global. Impossible de tracer les facteurs.\n")
        }
      }
    }, error = function(e) {
      cat("Une erreur s'est produite :\n")
      cat(e$message, "\n")
    })
    
    
    cat("Step 18 : Chargement de la fonction pour traiter les contributeurs et traitements ...\n")
    
    traiter_contributeurs <- function(contributors, facteur, view, sous_dossier) {
      top_contributors_df <- lapply(seq_along(contributors), function(i) {
        df <- data.frame(Contributeur = rownames(contributors[[i]]),
                         Value = as.vector(contributors[[i]]),
                         abs_Value = as.vector(abs(contributors[[i]])))
        return(df)
      })
      
      top_contributors_df <- do.call(rbind, top_contributors_df)
      
      # Sélectionner les 5% des meilleurs contributeurs
      cutoff_index <- ceiling(nrow(top_contributors_df) * 0.05)
      top_5_percent_contributors <- top_contributors_df[order(top_contributors_df$abs_Value, decreasing = TRUE), ][1:cutoff_index, ]
      top_5_percent_contributors <- top_5_percent_contributors[abs(top_5_percent_contributors$abs_Value) > 0.5, ]
      
      # Séparer les valeurs initiales en valeurs positives et négatives
      positive_contributors <- top_5_percent_contributors[top_5_percent_contributors$Value >= 0, ]
      negative_contributors <- top_5_percent_contributors[top_5_percent_contributors$Value < 0, ]
      
      # Renommer les colonnes
      colnames(positive_contributors)[colnames(positive_contributors) == "Value"] <- paste0("Factor", facteur)
      colnames(negative_contributors)[colnames(negative_contributors) == "abs_Value"] <- paste0("Factor", facteur, "_abs")
      
      # Nom du fichier CSV pour les contributeurs
      csv_name <- paste("top_5_percent_contributors_", facteur, "_", view, ".csv", sep = "")
      csv_pos <- paste("pos_contributors_", facteur, "_", view, ".csv", sep = "")
      csv_neg <- paste("neg_contributors_", facteur, "_", view, ".csv", sep = "")
      
      # Enregistrer les contributeurs dans des fichiers CSV
      write.csv(top_5_percent_contributors, file = file.path(sous_dossier, csv_name), row.names = FALSE)
      write.csv(positive_contributors, file = file.path(sous_dossier, csv_pos), row.names = FALSE)
      write.csv(negative_contributors, file = file.path(sous_dossier, csv_neg), row.names = FALSE)
    }
    
    cat("Step 19 : Différents plots restants...\n")
    
    # Liste des types cellulaires
    views <- c("BLymphocytes", "TLymphocytes", "Monocytes", "Neutrophils")
    
    for (view in views) {
      # Créer le sous-dossier s'il n'existe pas
      sous_dossier <- file.path(nouveau_dossier, gsub(" ", "", view))
      if (!file.exists(sous_dossier)) {
        dir.create(sous_dossier)
      }
      
      # Récupérer les contributeurs pour chaque facteur significatif
      if (length(facteurs_significatifs) == 0) {
        cat("Pas de facteur significatif trouvé pour la vue", view, "\n")
      } else {
        for (facteur in facteurs_significatifs) {
          # Récupérer les contributeurs
          contributors <- get_weights(MOFAmodel, view = view, factors = facteur, scale = TRUE)
          
          # Modifier les noms des lignes des facteurs pour chaque contributeur dans la liste
          contributors <- lapply(contributors, function(x) {
            if (!is.null(dim(x))) {
              rownames(x) <- gsub("_.*", "", rownames(x))
            }
            return(x)
          })
          
          # Filtrer les colonnes correspondant aux facteurs
          facteur_columns <- grep("^Factor[0-9]+$", colnames(contributors[[1]]))
          contributors <- lapply(contributors, function(x) x[, facteur_columns, drop = FALSE])
          
          # Appeler la fonction pour traiter les contributeurs
          traiter_contributeurs(contributors, facteur, view, sous_dossier)
          
          # Tracer le heatmap pour le facteur significatif et le type cellulaire actuel
          plot_heatmap <- plot_data_heatmap(MOFAmodel,
                                            view = view,
                                            factor = facteur,
                                            features = 20,
                                            cluster_rows = TRUE,
                                            cluster_cols = TRUE,
                                            show_rownames = TRUE,
                                            show_colnames = TRUE,
                                            denoise = TRUE,
                                            annotation_samples = "CONDITION")
          
          # Nom du fichier pour le heatmap
          heatmap_filename <- paste("heatmap_", facteur, "_", view, ".png", sep = "")
          
          # Enregistrer le heatmap dans le sous-dossier
          ggsave(file.path(sous_dossier, heatmap_filename), plot_heatmap, width = 12, height = 7)
          
          plot5 <- plot_weights(MOFAmodel, 
                                view = view, 
                                factor = facteur, 
                                nfeatures = 10, 
                                scale = TRUE)
          
          # Enregistrer le graphique dans 
          plot_filename <- paste("Poids_des_caracteristiques_", facteur, "_", view, ".png", sep = "")
          ggsave(file.path(sous_dossier, plot_filename), plot5, width = 12, height = 7)
          
          plot6 <- plot_data_scatter(MOFAmodel, view =  view,
                                     factor = facteur,
                                     features = 10, 
                                     add_lm = TRUE,  
                                     color_by = "CONDITION")
          # Enregistrer le graphique dans 
          plot_name <- paste("estimationderEgressionlineaire_", facteur, "_", view, ".png", sep = "")
          ggsave(file.path(sous_dossier, plot_name), plot6, width = 12, height = 7)
          
          plot7 <- plot_top_weights(MOFAmodel,
                                    view = view, 
                                    factor = facteur, 
                                    nfeatures = 20 )
          # Enregistrer le graphique dans 
          plot7_name <- paste("TOP_poids_des_caractéristiques_", facteur, "_", view, ".png", sep = "")
          ggsave(file.path(sous_dossier, plot7_name), plot7, width = 12, height = 7)
        }
      }
    }
  }
  
  cat("Step 20 : Résumer des résultats...\n")
  tryCatch({
  # Liste des dossiers correspondant aux analyses MOFA
  dossiers_MOFAs <- list.files(pattern = "^MOFA_analysis_[0-9]+$")
  
  # Initialiser une liste pour les résumés
  summary_list <- list()
  
  for (dossier_MOFA in dossiers_MOFAs) {
    # Récupérer le chiffre dans le nom du dossier
    analyse <- as.numeric(gsub("MOFA_analysis_", "", dossier_MOFA))
    
    # Lire le fichier Student_t_test.csv
    student_test_file <- file.path(dossier_MOFA, "Student_t_test.csv")
    if (file.exists(student_test_file)) {
      student_test_data <- read.csv(student_test_file)
      
      # Compter le nombre de facteurs significatifs
      nb_facteurs_significatifs <- sum(student_test_data$Significativite == "Significatif")
      
      # Trouver la valeur la plus élevée du p-value ajusté et la plus faible
      pvalues_adjusted <- student_test_data$Padjusted[student_test_data$Significativite == "Significatif"]
      high <- ifelse(length(pvalues_adjusted) > 0, max(pvalues_adjusted), NA)
      low <- ifelse(length(pvalues_adjusted) > 0, min(pvalues_adjusted), NA)
      
      # Calculer la somme des p-values ajustées
      total <- ifelse(length(pvalues_adjusted) > 0, sum(pvalues_adjusted), NA)
      
      # Créer un dataframe pour le résumé
      summary_df <- data.frame(Analyse = analyse, 
                               Facteur_Significatif = nb_facteurs_significatifs, 
                               High = high, 
                               Low = low,
                               Total = total)
      
      # Ajouter le dataframe au résumé
      summary_list[[dossier_MOFA]] <- summary_df
    } else {
      cat("Le fichier Student_t_test.csv n'a pas été trouvé dans le dossier", dossier_MOFA, "\n")
    }
  }
  
  # Concaténer tous les dataframes de la liste summary_list en un seul dataframe
  summary_df <- do.call(rbind, summary_list)
  
  # Trier le dataframe en fonction du nombre de facteurs significatifs (ordre décroissant)
  summary_df <- summary_df[order(summary_df$Facteur_Significatif, decreasing = TRUE), ]
  
  # Trouver le plus grand nombre de facteurs significatifs
  max_facteurs_significatifs <- max(summary_df$Facteur_Significatif, na.rm = TRUE)
  
  # Trier uniquement les lignes ayant le plus grand nombre de facteurs significatifs selon la colonne Low
  choix_df <- summary_df[summary_df$Facteur_Significatif == max_facteurs_significatifs, ]
  choix_df <- choix_df[order(choix_df$Low, na.last = TRUE), ]
  
  # Enregistrer le dataframe dans un fichier CSV
  write.csv(summary_df, file = "summary_results.csv", row.names = FALSE)
  write.csv(choix_df, file = "model_choice.csv", row.names = FALSE)
  })
  
  cat("Step 21 : Analyse EnrichR...\n")
  
  tryCatch({
    
    # Initialiser Enrichr
    listEnrichrSites()
    setEnrichrSite("Enrichr")
    db_enrichR <- c("GO_Biological_Process_2023")
    
    
    # Fonction pour troncature des noms de termes pour un affichage clair
    truncate_names <- function(name, max_length = 40) {
      if (nchar(name) > max_length) {
        return(paste0(substr(name, 1, max_length - 3), "..."))
      } else {
        return(name)
      }
    }
    
    # Chemin vers le dossier du facteur sélectionné
    facteur_selectionne <- choix_df$Analyse[1]
    dossier_facteur <- paste0("MOFA_analysis_", facteur_selectionne)
    enrichr_folder <- file.path(dossier_facteur, "EnrichR_analyse")
    
    # Créer le dossier enrichr_analyse s'il n'existe pas
    if (!dir.exists(enrichr_folder)) {
      dir.create(enrichr_folder, recursive = TRUE)
    }
    
    # Liste des types cellulaires
    types_cellulaires <- c("BLymphocytes", "TLymphocytes", "Monocytes", "Neutrophils")
    
    # Parcourir les sous-dossiers pour chaque type cellulaire
    for (type_cellulaire in types_cellulaires) {
      sous_dossier <- file.path(dossier_facteur, type_cellulaire)
      
      if (dir.exists(sous_dossier)) {
        # Lire tous les fichiers CSV non vides commençant par "pos_" ou "neg_"
        csv_files <- list.files(sous_dossier, pattern = "^(pos_|neg_).*\\.csv$", full.names = TRUE)
        
        for (csv_file in csv_files) {
          if (nrow(read.csv(csv_file)) > 1) {
            # Lire le fichier CSV
            data <- read.csv(csv_file)
            
            # Extraire les contributeurs (nom des gènes)
            gene_list <- data$Contributeur
            # Créer l'objet Enrichr pour chaque fichier
            enrichr_result <- enrichr(gene_list, db_enrichR)
            enrichr_result <- do.call(rbind, enrichr_result)
            
            # Sélectionner les top termes enrichis
            top_terms <- 10
            enrichr_result_sorted <- enrichr_result %>%
              arrange(Adjusted.P.value) %>%
              head(top_terms)
            
            # Troncature des noms de termes pour un affichage clair
            enrichr_result_sorted$Term <- sapply(enrichr_result_sorted$Term, truncate_names)
            
            # Créer le graphique
            plot <- ggplot(enrichr_result_sorted, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value), fill = Combined.Score)) +
              geom_bar(stat = "identity") +
              coord_flip() +
              labs(title = "Top Enriched Terms",
                   x = "Term",
                   y = "-log10 Adjusted p-value",
                   fill = "Combined Score") +
              theme_minimal() +
              theme(plot.background = element_rect(fill = "white", color = NA),
                    panel.background = element_rect(fill = "white", color = NA))
            
            # Nom de la variable Enrichr
            var_name <- gsub("\\.csv$", "", basename(csv_file))
            
            # Chemin pour sauvegarder le graphique
            chemin_fichier_barplot <- file.path(enrichr_folder, paste0("enriched_terms_barplot_", type_cellulaire, "_", var_name, ".png"))
            
            # Sauvegarde du graphique
            ggsave(chemin_fichier_barplot, plot, width = 10, height = 5)
            # Sauvegarder l'objet Enrichr
            write.csv(enrichr_result, file = file.path(enrichr_folder, paste0("enrichr_", var_name, ".csv")), row.names = FALSE)
          }
        }
      }
    }
  
  })
  
  cat("Step 22 :  Nettoyage de l'environnement  ...\n")
  
  tryCatch({
  # Nettoyer l'environnement
  if (exists("MOFAmodel", envir = .GlobalEnv)) {
    rm(MOFAmodel, envir = .GlobalEnv)
    cat("MOFAmodel a été supprimé de l'environnement global.\n")
  }
  
  if (exists("facteurs_significatifs", envir = .GlobalEnv)) {
    rm(facteurs_significatifs, envir = .GlobalEnv)
    cat("facteurs_significatifs a été supprimé de l'environnement global.\n")
  }
  })
}

cat("Toutes les étapes ont été complétées avec succès!\n")
# Utilisation de la fonction avec le chemin du dossier
#MOFA_analysis("D:/data_umr_1227/Par_pureté_cyto_digitale", 29)