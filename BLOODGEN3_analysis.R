
BloodGen_analysis <- function(folder_path) {
  # Définir le répertoire de travail
  setwd("D:/data_umr_1227/Par_pureté_cyto_digitale")
  
  # Charger les librairies nécessaires
  library(BloodGen3Module)
  library(SummarizedExperiment)
  library(dplyr)
  library(data.table)
  
  # Définir le chemin du nouveau dossier et créer le dossier s'il n'existe pas
  cat("Etape 1: Extraction des noms de fichiers et création du nouveau dossier...\n")
  dossier_nom <- basename(folder_path)
  nom_sans_DESeq <- gsub("_DESeq.*", "", dossier_nom)
  nouveau_dossier <- paste0(nom_sans_DESeq, "_BloodGen3")
  chemin_nouveau_dossier <- file.path(dirname(folder_path), nouveau_dossier)
  if (!file.exists(chemin_nouveau_dossier)) {
    dir.create(chemin_nouveau_dossier)
  }
  
  # Chargement des fichiers dans DESeq
  cat("Etape 2: Chargement des fichiers dans DESeq...\n")
  count_normalized <- read.csv(list.files(folder_path, pattern="^VST_normalized_", full.names=TRUE),  header = TRUE, check.names = FALSE,row.names = 1)
  gene_list_ensbl <- read.csv("D:/data_umr_1227/Par_pureté_cyto_digitale/mart_export.csv", header = TRUE, check.names = FALSE)
  coldata <- read.csv(list.files(folder_path, pattern="^colData_",  full.names=TRUE), row.names = 1)
  Module_listGen3 <- read.csv("D:/data_umr_1227/Par_pureté_cyto_digitale/Module_listGen3.csv", header = TRUE, check.names = FALSE)
  load("vst_RNAseq.RData")
  
  # Mise en format du tableau des counts normalized
  cat("Etape 3: Mise en format du tableau des counts normalized...\n")
  clean_colnames <- function(col_names) {
    cleaned_names <- gsub("^(\\d+)_.*", "\\1", col_names)
    return(cleaned_names)
  }
  clean_rownames <- function(row_names) {
    cleaned_names <- gsub("^(\\d+)_.*", "\\1", row_names)
    return(cleaned_names)
  }
  count_normalized$Gene.stable.ID <- rownames(count_normalized)
  count_normalized <- merge(count_normalized, gene_list_ensbl, by = "Gene.stable.ID")
  count_normalized <- count_normalized[!duplicated(count_normalized$Gene.name), ]
  #count_normalized <- aggregate(. ~ Gene.name, data = count_normalized, FUN = sum)
  rownames(count_normalized) <- count_normalized$Gene.name
  count_normalized <- count_normalized[, !colnames(count_normalized) %in% c("Gene.name", "Gene.stable.ID")]
  count_normalized <- count_normalized[rownames(count_normalized) %in% Module_listGen3$Gene, ]
  colnames(count_normalized) <- clean_colnames(colnames(count_normalized))
  count_normalized <- count_normalized[order(rownames(count_normalized)),]
  #count_normalized <- count_normalized[, !colnames(count_normalized) %in% "32152162", drop = FALSE]
  # Liste des colonnes à supprimer
  #cols_to_remove <- c("32152162", "autre_colonne", "encore_une_autre_colonne")
  #count_normalized <- count_normalized[, !colnames(count_normalized) %in% cols_to_remove, drop = FALSE]
  
  # Filtrer les colonnes et les lignes communes entre count_normalized et PreciseSADS_RNASeq
  cat("Etape 4: Filtrage des colonnes et des lignes communes...\n")
  common_rows <- rownames(PreciseSADS_RNASeq) %in% rownames(count_normalized)
  common_cols <- colnames(PreciseSADS_RNASeq) %in% colnames(count_normalized)
  filtered_precise <- PreciseSADS_RNASeq[common_rows, common_cols]
  filtered_precise_df <- as.data.frame(filtered_precise)
  extra_cols <- setdiff(colnames(count_normalized), colnames(filtered_precise))
  
  # Si des colonnes supplémentaires sont présentes dans count_normalized, les ajouter à filtered_precise
  if (length(extra_cols) > 0) {
    extra_cols_data <- count_normalized[common_rows, extra_cols]
    filtered_precise <- cbind(filtered_precise, extra_cols_data)
  }
  
  common_rows <- intersect(rownames(filtered_precise_df), rownames(count_normalized))
  common_cols <- intersect(colnames(filtered_precise_df), colnames(count_normalized))
  
  # Mettre à jour les valeurs de count_normalized avec celles de filtered_precise_df
  for (row in common_rows) {
    for (col in common_cols) {
      count_normalized[row, col] <- filtered_precise_df[row, col]
    }
  }
  
  colnames(count_normalized) <- paste0("X", colnames(count_normalized))
  
  cat("Etape 5: Mise en format du fichier coldata...\n")
  rownames(coldata) <- clean_rownames(rownames(coldata))
  rownames(coldata) <- paste0("X", rownames(coldata))
  coldata <- coldata[rownames(coldata) %in% colnames(count_normalized), ]
  
  cat("Etape 6: Vérification des fichiers...\n")
  if(all(colnames(count_normalized) == rownames(coldata))) {
    cat("Les fichiers sont conformes.\n")
  } else {
    cat("Les fichiers ne sont pas conformes.\n")
    stop("Les fichiers ne sont pas conformes.")
  }
  
  cat("Etape 7: Transformation en matrice...\n")
  objet_summrize <- SummarizedExperiment(assays=list(count_normalized), colData=coldata)
  #print(class(objet_summrize))
  
  cat("Etape 8: Création de l'objet SummarizedExperiment et exécution de Group...\n")
  Group_df <- Groupcomparison(objet_summrize,
                              sample_info = NULL,
                              FC = 1.3,
                              pval = 0.05,
                              FDR = TRUE,
                              Group_column = "CONDITION",              
                              Test_group = "SJS",
                              Ref_group = "CTRL")
  
  gridplot(Group_df, 
           cutoff = 10, 
           Ref_group = "CTRL")
  
  cat("Etape 9: Création de l'objet SummarizedExperiment et exécution de Individual...\n")
  objet_summrize_id <- SummarizedExperiment(assays=list(count_normalized), colData=coldata)
  
  Individual_df = Individualcomparison(objet_summrize_id, sample_info = NULL, FC = 0.5, DIFF = 1, Group_column = "CONDITION", Ref_group = "CTRL")
  
  fingerprintplot(Individual_df,
                  sample_info = NULL,
                  cutoff = 15,
                  rowSplit= TRUE ,
                  Group_column= "Group_test",
                  show_ref_group = FALSE, 
                  Ref_group = "Control",
                  Aggregate = NULL,
                  filename = chemin_nouveau_dossier,
                  height = 38,
                  width = 17)
  
  
  cat("Etape 10: Sauvegarde des images\n")
  file.rename(".pdf", file.path(chemin_nouveau_dossier, "blood_Group.pdf"))
  pdf_files <- list.files(pattern = "_BloodGen3\\.pdf$", full.names = TRUE)
  if (length(pdf_files) > 0) {
    file.rename(pdf_files, file.path(chemin_nouveau_dossier, basename(pdf_files)))
  }
  #file.rename("BLymphocytes_BloodGen3.pdf", file.path(chemin_nouveau_dossier, "blood_individual.pdf"))
  
  cat("Toutes les étapes ont été complétées avec succès!\n")
}
#BloodGen_analysis("D:/data_umr_1227/Par_pureté_cyto_digitale/TLymphocytes_DESeq")

                              