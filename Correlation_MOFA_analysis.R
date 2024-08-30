Correlaton_MOFA_analysis <- function(folder_path) {
  
  # Charger les bibliothèques nécessaires
  library(dplyr)
  library(readr)
  library(corrr)
  library(vroom)
  library(ggplot2)
  
  data_clinical_num <- vroom(file.path(folder_path, "CLINICAL_NUMERIC_STRD_SJS.tsv"), delim = '\t')
  data_clinical_bin <- vroom(file.path(folder_path, "CLINICAL_BINARY_STRD_SJS.tsv"), delim = '\t')

  # Étape 1: Ouvrir le fichier choice.csv et récupérer la valeur de la colonne L4ANALYSE
  choice_file <- file.path(folder_path, "model_choice.csv")
  if (!file.exists(choice_file)) {
    stop("Le fichier choice.csv n'existe pas dans le répertoire spécifié.")
  }
  
  choice_data <- read_csv(choice_file)
  analyse_number <- choice_data$Analyse
  
  # Vérifiez si plusieurs valeurs sont présentes, si oui, sélectionnez la première
  if (length(analyse_number) > 1) {
    analyse_number <- analyse_number[1]
  }
  
  # Étape 2: Ouvrir le dossier MOFA_analyse suivi du numéro récupéré
  dossier_analyse <- file.path(folder_path, paste0("MOFA_analysis_", analyse_number))
  
  # Vérifier si le dossier existe
  if (!dir.exists(dossier_analyse)) {
    stop("Le dossier ", dossier_analyse, " n'existe pas.")
  }
  
  # Étape 3: Charger le fichier distribution.csv
  distribution_file <- file.path(dossier_analyse, "distribution_tableau.csv")
  
  # Vérifier si le fichier existe
  if (!file.exists(distribution_file)) {
    stop("Le fichier ", distribution_file, " n'existe pas.")
  }
  
  distribution_data <- read_csv(distribution_file)
  
  # Étape 4: Charger le fichier Student_t_test.csv et récupérer les facteurs significatifs
  student_test_file <- file.path(dossier_analyse, "Student_t_test.csv")
  
  # Vérifier si le fichier existe
  if (!file.exists(student_test_file)) {
    stop("Le fichier ", student_test_file, " n'existe pas.")
  }
  
  student_test_data <- read_csv(student_test_file)
  
  # Récupérer les noms des colonnes des facteurs significatifs
  significant_factors <- student_test_data %>%
    filter(Significativite == "Significatif") %>%
    pull(Facteur)
  
  # Ajouter la colonne "ID" à la liste des colonnes à conserver
  columns_to_keep <- c("ID", significant_factors)
  
  # Filtrer les colonnes du fichier distribution.csv pour ne garder que les facteurs significatifs et ID
  filtered_distribution_data <- distribution_data %>% select(all_of(columns_to_keep))
  
  filtered_distribution_data <- semi_join(filtered_distribution_data, data_clinical_num, by = "ID")
  filtered_clinical_data <- semi_join(data_clinical_num, filtered_distribution_data, by = "ID")
  filtered_clinical_binary_data <- semi_join(data_clinical_bin, filtered_distribution_data, by = "ID")
  
  # Supprimer les colonnes qui ont uniquement des valeurs NA
  filtered_clinical_data <- filtered_clinical_data %>% select_if(~ !all(is.na(.)))
  filtered_clinical_binary_data <- filtered_clinical_binary_data %>% select_if(~ !all(is.na(.)))
  
  # Fusionner les deux tableaux en utilisant la colonne ID comme clé de fusion pour les données binaires
  merged_data <- merge(filtered_clinical_binary_data, distribution_data[, c("ID", "CONDITION")], by = "ID", all.x = TRUE)
  
  # Créer une liste pour stocker les résultats
  results_list <- list()
  
  # Boucle sur chaque facteur (en supposant que les colonnes sont nommées comme "Factor1", "Factor2", etc.)
  for (factor in colnames(filtered_distribution_data)[-1]) {  # Exclure la colonne ID
    cat("Processing factor:", factor, "\n")  # Message de débogage
    temp_results <- list()
    
    for (variable in colnames(filtered_clinical_data)[-1]) {  # Exclure la colonne ID
      temp_data <- data.frame(
        ID = filtered_distribution_data$ID,
        Factor = filtered_distribution_data[[factor]],
        Variable = filtered_clinical_data[[variable]]
      )
      
      temp_data <- temp_data[complete.cases(temp_data), ]
      
      if (nrow(temp_data) >= 4 && var(temp_data$Variable, na.rm = TRUE) > 0) {
        correlation_test <- try(cor.test(temp_data$Variable, temp_data$Factor), silent = TRUE)
        
        if (!inherits(correlation_test, "try-error")) {
          temp_results[[variable]] <- data.frame(
            Facteur = factor,
            Variable = variable,
            Pvalue = correlation_test$p.value,
            Cor = correlation_test$estimate
          )
        }
      }
    }
    
    if (length(temp_results) > 0) {
      temp_results_df <- do.call(rbind, temp_results)
      temp_results_df$Padjusted <- p.adjust(temp_results_df$Pvalue, method = "fdr")
      temp_results_df$Significativite <- ifelse(temp_results_df$Pvalue <= 0.05, "Significatif", "Non significatif")
      results_list[[factor]] <- temp_results_df
    }
  }
  
  # Concaténer tous les résultats ajustés dans un seul dataframe
  adjusted_results_df <- do.call(rbind, results_list)
  variable_significative <- adjusted_results_df[adjusted_results_df$Significativite == "Significatif", ]
  write.csv(variable_significative, file = file.path(dossier_analyse, "variable_numeric_signif.csv"), row.names = FALSE)
  
  # Fonction pour récupérer les données des variables significatives
  retrieve_data <- function(data, variable_significative) {
    data_list <- list()
    for (i in 1:nrow(variable_significative)) {
      factor_name <- variable_significative[i, "Facteur"]
      variable_name <- variable_significative[i, "Variable"]
      factor_label <- factor_name  # Utiliser directement le nom du facteur
      variable_label <- variable_name  # Utiliser directement le nom de la variable
      temp_data <- data.frame(
        ID = filtered_distribution_data$ID,
        Factor = filtered_distribution_data[[factor_name]],
        Variable = filtered_clinical_data[[variable_name]]
      )
      temp_data <- temp_data[complete.cases(temp_data), ]
      if (nrow(temp_data) > 0) {
        data_list[[i]] <- list(data = temp_data, factor_label = factor_label, variable_label = variable_label)
      } else {
        cat("No complete cases for factor:", factor_name, "and variable:", variable_name, "\n")
      }
    }
    return(data_list)
  }
  
  # Récupérer les données des variables significatives avec les noms réels des facteurs et des variables
  significant_data_list <- retrieve_data(filtered_distribution_data, variable_significative)
  
  # Fonction pour calculer la p-value d'une régression linéaire
  calculate_pvalue <- function(data) {
    model <- lm(Variable ~ Factor, data = data)
    summary_model <- summary(model)
    p_value <- summary_model$coefficients[2, 4]
    return(p_value)
  }
  
  # Générer les scatterplots pour les variables significatives et les sauvegarder dans le dossier spécifié
  save_scatterplots <- function(data_list, output_dir) {
    plots <- list()
    for (i in 1:length(data_list)) {
      if (!is.null(data_list[[i]])) {
        p_value <- calculate_pvalue(data_list[[i]]$data)
        correlation <- cor(data_list[[i]]$data$Factor, data_list[[i]]$data$Variable)
        plot <- ggplot(data_list[[i]]$data, aes(x = Factor, y = Variable)) +
          geom_point() +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          labs(title = paste("Scatterplot for", data_list[[i]]$variable_label, "vs.", data_list[[i]]$factor_label, "\nP-value:", format(p_value, digits = 3),"\nCor:", format(correlation, digits = 3)),
               x = data_list[[i]]$factor_label, y = data_list[[i]]$variable_label)
        
        # Save each plot with the naming format: scatterplot_Factor_Variable_i.png
        factor_name <- gsub("[^A-Za-z0-9]", "", data_list[[i]]$factor_label)
        variable_name <- gsub("[^A-Za-z0-9]", "", data_list[[i]]$variable_label)
        plot_filename <- file.path(output_dir, paste0("scatterplot_", factor_name, "_", variable_name, "_", i, ".png"))
        ggsave(filename = plot_filename, plot = plot, width = 7, height = 5)
        
        plots[[i]] <- plot
      }
    }
    return(plots)
  }
  
  # Sauvegarder les scatterplots
  scatterplots <- save_scatterplots(significant_data_list, dossier_analyse)
  
  cat("Scatterplots saved successfully.\n")
  
}

# Appeler la fonction avec le chemin du dossier
#folder_path <- "D:/data_umr_1227/Par_pureté_cyto_digitale/DATA_work"
#Correlaton_MOFA_analysis(folder_path)
