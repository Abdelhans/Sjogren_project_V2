# Description du Projet

Ce projet contient contient les scripts utilisés pour mon projet de stage de Master2 mené au UMR1227- LBAI de Brest. Le projet se concentre sur les changements dans le transcriptome de plusieurs cellules imùmunitaires dans la maladie de Sjögren par le biais de l’analyse transcriptomique, multi-omique et en réseau.

# Prérequis

Avant d'utiliser ces scripts, assurez-vous d'avoir installé les logiciels et paquets suivants :

### Liste des Packages R
- R (version 4.0 ou supérieure)
- RStudio (facultatif mais recommandé)
- **DESeq2** : Outil pour l'analyse des gènes différentiellement exprimés.
- **pheatmap** : Génération de heatmaps pour la visualisation des données.
- **dplyr** : Manipulation de données en R.
- **RColorBrewer** : Palettes de couleurs pour les graphiques.
- **ggplot2** : Création de graphiques élégants et complexes.
- **ggrepel** : Amélioration du placement des étiquettes dans les graphiques ggplot2.
- **fgsea** : Analyse de l'enrichissement des ensembles de gènes (Gene Set Enrichment Analysis).
- **readxl** : Lecture de fichiers Excel.
- **stringr** : Manipulation de chaînes de caractères.
- **data.table** : Manipulation efficace de grandes ensembles de données.
- **VennDiagram** : Création de diagrammes de Venn.
- **readr** : Lecture de données rectangulaires telles que CSV.
- **grid** : Graphiques de bas niveau (utilisé par pheatmap).
- **reshape2** : Transformation des données pour ggplot2.
- **purrr** : Programmation fonctionnelle pour manipulation de listes et vecteurs.
- **tibble** : Cadres de données modernes.
- **MOFA2** : Analyse factorielle multi-omique.
- **reticulate** : Intégration de Python dans R.
- **enrichR** : Enrichissement fonctionnel pour les listes de gènes.
- **vroom** : Lecture ultra-rapide des données tabulaires.
- **wgcna** : Outil pour l'analyse de reseau.
  
# Installation des Packages R

Vous pouvez installer tous les packages nécessaires en exécutant le script ci-dessous dans R ou RStudio :

```r
install.packages(c("pheatmap", "dplyr", "RColorBrewer", "ggplot2", "ggrepel", "fgsea", "readxl", "stringr", "data.table", "VennDiagram", "readr", "grid", "reshape2", "purrr", "tibble", "reticulate", "enrichR", "vroom"))

# Certains packages doivent être installés depuis Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "MOFA2"))
```
# Installation

### Option 1 : Cloner le dépôt GitHub

Clonez ce dépôt sur votre machine locale en utilisant la commande suivante :

```bash
git clone https://github.com/Abdelhans/Sjogren_project.git
```
### Option 2 : Télécharger les scripts manuellement
- Téléchargez les scripts individuellement depuis le dépôt GitHub.
- Enregistrez-les dans un dossier de votre choix sur votre machine locale.
- Ouvrez RStudio et chargez chaque script en utilisant le menu "File > Open File".
  
Lancer les Scripts R
Vous pouvez maintenant exécuter les scripts dans RStudio pour obtenir les résultats des analyses.

# Utilisation

### Script 1 : `Creating_tbl_analysis.R`
Ce script prépare les données en effectuant diverses transformations et filtrations nécessaires pour les analyses ultérieures. Il s'occupe des données brutes, ajuste les métadonnées et filtre les matrices de comptage.
#3# Fonction principale : `process_data` accomplit les tâches suivantes :
- Chargement des fichiers de données
- Filtrage des métadonnées et des matrices de comptage
- Normalisation des données
- Écriture des résultats dans des fichiers CSV
### Arguments : 
- gene_file: Fichier CSV contenant les informations sur les gènes
- meta_file: Fichier TSV avec les métadonnées
- res_cybers: Résultats de cybersort au format CSV
- count_file: Matrice de comptage au format TSV
- meta_output_file: Fichier CSV pour les métadonnées filtrées
- count_output_file: Fichier CSV pour la matrice de comptage filtrée
- ctrl_meta_output_file: Fichier CSV pour les métadonnées du contrôle
- ctrl_count_output_file: Fichier CSV pour la matrice de comptage du contrôle
- cybersort_output_file: Fichier CSV pour les résultats de cybersort
- nomalize_file: Fichier CSV pour les données normalisées
- SFnomalize_file: Fichier CSV pour les facteurs de taille normalisés
### Exécution
Pour exécuter ce script, assurez-vous que vous avez les fichiers d'entrée nécessaires et ajustez le chemin de travail (`setwd`) en fonction de l'emplacement de vos données. Ensuite, exécutez la fonction `process_data` avec les paramètres appropriés :
```r
source("Creating_tbl_analysis.R")
process_data(gene_file = "path/to/mart_export.csv",meta_file = "path/to/METADATA_SORTED_CELL_MCPcounter.tsv",res_cybers = "path/to/CIBERSORTx_Job6_Results.csv",count_file = "path/to/ALL_SORTED_CELLS_MATRIX_ENSG.tsv",meta_output_file = "path/to/meta_matrix.csv",count_output_file = "path/to/count_matrix.csv",ctrl_meta_output_file = "path/to/ctrl_meta_matrix.csv",ctrl_count_output_file = "path/to/ctrl_count_matrix.csv",cybersort_output_file = "path/to/cybersort_collapse.csv",nomalize_file = "path/to/normalized_vst_2.csv",SFnomalize_file = "path/to/normalized_sizefactor_1.csv")
```

### Script 2 : `DESeq_analysis.R`
Ce script réalise une analyse DESeq2 sur les données filtrées :
- Filtrage des données selon le type cellulaire spécifié
- Création d'objets DESeq2
- Analyse DESeq2 et sauvegarde des résultats
- Création et sauvegarde de graphiques (volcanoplot, PCA, heatmap)
### Arguments
- meta_file: Fichier CSV avec les métadonnées
- count_file: Fichier CSV avec la matrice de comptage
- gene_file: Fichier CSV avec les informations sur les gènes
- rescyber: Résultats de cybersort au format CSV
- ctrl_meta: Fichier CSV pour les métadonnées du contrôle
- ctrl_count: Fichier CSV pour la matrice de comptage du contrôle
- digital_purity: Nom de la colonne indiquant la pureté digitale
- cell_type: Type de cellule à analyser
### Exécution
```r
source("DESeq_analysis.R")
DESeq_analysis("meta_matrix.csv", "count_matrix.csv", "mart_export.csv", "cybersort_collapse.csv", "ctrl_meta_matrix.csv", "ctrl_count_matrix.csv", "LB_purity", "BLymphocytes")
```

### Script 3 : `venn_diagram_DESeq2.R`
Le script venn_diagram_DESeq2.R permet de générer un diagramme de Venn pour comparer les gènes différentiellement exprimés (DEGs) obtenus à partir des analyses DESeq2 effectuées sur différents types cellulaires. Le script calcule également les intersections entre les ensembles de DEGs et identifie les gènes uniques à chaque type cellulaire.
### Arguments : 
base_path: Le chemin de base où se trouvent les sous-dossiers contenant les résultats DESeq2.
### Exécution
```r
source("venn_diagram_DESeq2.R")
venn_diagram_DESeq2("D:/data_umr_1227/Par_pureté_cyto_digitale")
```

### Script 4 : `GSEA_analysis.R`

Le script GSEA_analysis.R permet de réaliser une analyse GSEA (Gene Set Enrichment Analysis) à partir des résultats d'une analyse DESeq2. Il calcule les scores d'enrichissement normalisés pour des ensembles de gènes, produit des heatmaps pour les voies biologiques significativement enrichies, et génère des graphiques pour visualiser les voies enrichies.
### Arguments : 
folder_path: Le chemin vers le dossier contenant les résultats DESeq2 et les fichiers nécessaires.
### Exécution
```r
source("GSEA_analysis.R")
GSEA_analysis("D:/data_umr_1227/Par_pureté_cyto_digitale/Type_cell_DESeq")
```
### Script 5 : `MOFA_analysis.R`

Ce script analyse les données de différents types cellulaires en utilisant MOFA2 pour l'analyse d'intégration des doonées multi omiques(réduction de dimensionnalité). Le script principal, `MOFA_analysis.R`, effectue l'analyse de variance, prépare les données, entraîne un modèle MOFA, et génère des visualisations.7

### Arguments : 
folder_path: Le chemin vers le dossier contenant les résultats DESeq2 et les fichiers nécessaires.

### Exécution
```r
source("R/MOFA_analysis.R")
MOFA_analysis("path/to/data", number_factors_to_start = 3)
```

### Script 6 : `venn_diagram_MOFA.R`

Ce script analyse les contributeurs pour différents types cellulaires en utilisant des fichiers CSV et génère un diagramme de Venn ainsi qu'un tableau récapitulatif des gènes. Il est conçu pour être utilisé avec des fichiers de résultats spécifiques.

### Arguments : 
- path : Chemin vers le dossier contenant les sous-dossiers pour chaque type cellulaire.
- factor_num : Numéro du facteur pour les fichiers de contributeurs.
- prefix : Préfixe des noms de fichiers des contributeurs.

# Exécution
```r
source("venn_diagram_MOFA.R")
venn_diagram_MOFA("D:/data_umr_1227/Par_pureté_cyto_digitale/MOFA_analysis_26", 12, "neg")
```

### Script 7 : `Correlation_MOFA_analysis.R`

Ce script analyse la corrélation entre les facteurs et les variables cliniques à partir de fichiers CSV et TSV. Il génère des scatterplots pour les variables significatives et effectue des tests de corrélation pour les facteurs et variables cliniques.

### Arguments : 

- folder_path : Chemin vers le dossier contenant les fichiers nécessaires à l'analyse.
  
### Exécution

```r
source("Correlaton_MOFA_analysis.R")
Correlaton_MOFA_analysis("D:/data_umr_1227/Par_pureté_cyto_digitale/MOFA_analysis_26")
```

### Script 8 : `WGCNA_analysis.Rmd`

Ce script réaliser une analyse mWGCNA pour identifier les modules de gènes co-exprimés dans nos données pour plusieurs type cellulaire . Le script inclut la préparation des données, la construction du réseau de co-expression, l'identification des modules, et la génération de graphiques pour visualiser les résultats. Le script est conçu pour être flexible et efficace, en permettant de charger directement les objets R préexistants afin de réduire le temps d'exécution.

L'analyse WGCNA peut être gourmande en ressources et en temps. Le temps nécessaire dépend des facteurs suivants :
- Nombre de gènes : Plus le nombre de gènes est élevé, plus l'analyse peut être longue.
- Nombre d'échantillons : Une plus grande quantité d'échantillons peut rallonger le temps d'exécution.
- Complexité des Calculs : Les étapes de transformation, construction du réseau, et détection des modules impliquent des calculs intensifs.
  
Le script WGCNA_analysis.R est conçu pour optimiser l'exécution et offrir une flexibilité maximale :

Les objets R créés lors des analyses précédentes sont chargés directement depuis des fichiers RData. Cela permet d'éviter la reconstruction complète des objets à chaque exécution, réduisant ainsi le temps d'analyse. Des lignes de commande commentées montrent les étapes pour générer les objets R utilisés dans l'analyse. Ces lignes peuvent être décommentées pour reproduire entièrement le processus.

### Arguments : 
- folder_path : Chemin vers le dossier contenant les fichiers nécessaires à l'analyse.
# Exécution
Si vous souhaitez reproduire les objets R, décommentez les lignes de code correspondantes et exécutez-les.

```r
source("WGCNA_analysis.Rmd")
WGCNA_analysis.Rmd("D:/data_umr_1227/Par_pureté_cyto_digitale")
```
### Script 9 : `machin_learning_wgcna.py`

Description :

Ce script Python effectue une analyse approfondie des gènes identifiés avec wgcna en utilisant des techniques de classification pour prédire des étiquettes cliniques : malade ou temoin basées sur le niveau d'expression de ces génes. Il se compose de plusieurs étapes clés, allant de la préparation des données à l'évaluation des modèles de machine learning. Le script est conçu pour être modulable et permet une analyse efficace des données.

Fonctionnalités du Script :

- Chargement et Préparation des Données 
- Prétraitement des Données
- Normalisation des données
- Construction et Évaluation des Modèles de Classification
- Entraînement et évaluation d'un modèle de régression logistique et d'un modèle de forêt aléatoire (Random Forest) sur les données.
- Calcul et affichage des performances des modèles à l'aide des métriques d'accuracy et de classification.
- Génération des matrices de confusion pour visualiser la performance des modèles.
- Visualisation des Résultats
- Tracé des courbes ROC pour évaluer la performance des modèles en termes de taux de faux positifs et de taux de vrais positifs.
- Affichage des importances des caractéristiques pour le modèle Random Forest.
- la génération d'un rapport PDF contenant les résultats et les visualisations.
  
Exécution :

Pour exécuter ce script, assurez-vous d'avoir installé les bibliothèques Python nécessaires, telles que pandas, sklearn, matplotlib, et seaborn. Le script suppose également que les chemins vers les fichiers de données sont correctement définis et que les fichiers requis sont présents aux emplacements spécifiés.

Arguments :

--data_path : Chemin vers le fichier de données génétiques.

--meta_data_path : Chemin vers le fichier des métadonnées.

--sno_list_path : Chemin vers le fichier de la liste des gènes d'intérêt.

```r
python machin_learning_wgcna.py --data "D:/data_umr_1227/Par_pureté_cyto_digitale/BLymphocytes_DESeq/Sizefactornorm_all_sampleBLymphocytes_DESeq.csv" --meta 'D:/data_umr_1227/Par_pureté_cyto_digitale/METADATA_SORTED_CELL_MCPcounter.tsv' --sno "D:/data_umr_1227/Par_pureté_cyto_digitale/WGCNA_analysis/Sno_gene_list.csv"


```

# Contribution
Les contributions sont les bienvenues ! Si vous avez des suggestions ou des améliorations, n'hésitez pas à créer une pull request ou à ouvrir une issue.
# Auteurs
Abdel YAYA-OYE AKIBOU 

# Remerciements
Merci à tous ceux qui ont contribué à ce projet, en particulier à Cristian IPERI pour son soutien et sa patience.
