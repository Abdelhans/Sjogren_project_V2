import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import argparse
import os

# Charger les données
def load_data(data_path, meta_data_path, sno_list_path):
    data_LB = pd.read_csv(data_path)
    data_LB.rename(columns={'Unnamed: 0': 'Gene'}, inplace=True)
    meta_data = pd.read_csv(meta_data_path, sep='\t')
    data_sno_list = pd.read_csv(sno_list_path)
    data_sno_list.rename(columns={'Unnamed: 0': ''}, inplace=True)
    return data_LB, meta_data, data_sno_list

# Filtrer les données
def filter_data(data_LB, data_sno_list):
    unknowns_values = data_sno_list['x'].values
    return data_LB[data_LB['Gene'].isin(unknowns_values)]

# Transposer le DataFrame
def transpose_dataframe(df):
    df_transposed = df.transpose()
    new_header = df_transposed.iloc[0]  # La première ligne comme en-tête
    df_transposed = df_transposed[1:]  # Retirer la première ligne du DataFrame
    df_transposed.columns = new_header  # Définir les nouvelles colonnes
    df_transposed.reset_index(inplace=True)
    df_transposed.rename(columns={'index': 'Gene'}, inplace=True)
    return df_transposed

# Ajouter les labels
def add_label(df, meta_data):
    id_to_condition = meta_data.set_index('ID')['CONDITION'].to_dict()
    df['LABEL'] = df['Gene'].map(id_to_condition)
    return df

# Remplacer les labels
def replace_labels(df, column_name='LABEL'):
    df[column_name] = df[column_name].replace({'SJS': 1, 'CTRL': 0})
    return df

# Modifier les colonnes pour garder les valeurs numériques
def modify_column_to_numeric(df, column_name='Gene'):
    df[column_name] = df[column_name].apply(lambda x: ''.join(filter(str.isdigit, x)))
    return df

# Génération des graphiques et du PDF
def generate_pdf(report, cm_log_reg_fig, cm_rf_fig, roc_fig, feature_importances_fig, output_path):
    c = canvas.Canvas(output_path, pagesize=letter)
    width, height = letter

    # Ajouter le titre du rapport
    c.setFont("Helvetica", 12)
    c.drawString(30, height - 40, "Rapport de Classification")
    
    # Ajouter le texte du rapport
    text = c.beginText(40, height - 60)
    text.setFont("Helvetica", 10)
    text.textLines(report)
    c.drawText(text)

    # Ajouter chaque graphique à une nouvelle page
    c.showPage()
    c.drawImage(cm_log_reg_fig, 0, height - 500, width=width, height=500)
    c.showPage()
    c.drawImage(cm_rf_fig, 0, height - 500, width=width, height=500)
    c.showPage()
    c.drawImage(roc_fig, 0, height - 500, width=width, height=500)
    c.showPage()
    c.drawImage(feature_importances_fig, 0, height - 500, width=width, height=500)

    c.save()

def main(data_path, meta_data_path, sno_list_path):
    # Chargement et prétraitement des données
    data_LB, meta_data, data_sno_list = load_data(data_path, meta_data_path, sno_list_path)
    filtered_BLymphocytes = filter_data(data_LB, data_sno_list)
    filtered_BLymphocytes_transposed = transpose_dataframe(filtered_BLymphocytes)
    filtered_BLymphocytes_label = add_label(filtered_BLymphocytes_transposed, meta_data)
    filtered_BLymphocytes_label = replace_labels(filtered_BLymphocytes_label)
    final_data_BLymphocytes = modify_column_to_numeric(filtered_BLymphocytes_label)

    # Préparation des données pour les modèles
    X = final_data_BLymphocytes.drop(columns=['Gene', 'LABEL'])
    Y = final_data_BLymphocytes['LABEL']

    # Diviser les données en ensemble d'entraînement et de test
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3, random_state=42)

    # Standardiser les données
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Initialiser et entraîner le modèle de régression logistique
    log_reg = LogisticRegression(max_iter=1000, solver='saga')
    log_reg.fit(X_train_scaled, Y_train)

    # Prédictions et évaluation
    Y_pred_log_reg = log_reg.predict(X_test_scaled)
    accuracy_log_reg = accuracy_score(Y_test, Y_pred_log_reg)
    report_log_reg = classification_report(Y_test, Y_pred_log_reg)
    cm_log_reg = confusion_matrix(Y_test, Y_pred_log_reg)

    # Initialiser et entraîner le modèle Random Forest
    rf = RandomForestClassifier(random_state=42)
    rf.fit(X_train_scaled, Y_train)
    Y_pred_rf = rf.predict(X_test_scaled)
    accuracy_rf = accuracy_score(Y_test, Y_pred_rf)
    report_rf = classification_report(Y_test, Y_pred_rf)
    cm_rf = confusion_matrix(Y_test, Y_pred_rf)

    # Calcul des courbes ROC
    Y_prob_log_reg = log_reg.predict_proba(X_test_scaled)[:, 1]
    fpr_log_reg, tpr_log_reg, _ = roc_curve(Y_test, Y_prob_log_reg)
    roc_auc_log_reg = auc(fpr_log_reg, tpr_log_reg)

    Y_prob_rf = rf.predict_proba(X_test_scaled)[:, 1]
    fpr_rf, tpr_rf, _ = roc_curve(Y_test, Y_prob_rf)
    roc_auc_rf = auc(fpr_rf, tpr_rf)

    # Calcul de l'importance des caractéristiques
    importances = rf.feature_importances_
    feature_importances = pd.Series(importances, index=X.columns).sort_values(ascending=False)

    # Création des graphiques
    plt.figure(figsize=(10, 7))
    sns.heatmap(cm_log_reg, annot=True, fmt='d', cmap='Blues', xticklabels=['Témoin', 'Malade'], yticklabels=['Témoin', 'Malade'])
    plt.title('Matrice de Confusion - Régression Logistique')
    plt.xlabel('Prédictions')
    plt.ylabel('Véritables')
    cm_log_reg_fig = "cm_log_reg.png"
    plt.savefig(cm_log_reg_fig)
    plt.close()

    plt.figure(figsize=(10, 7))
    sns.heatmap(cm_rf, annot=True, fmt='d', cmap='Blues', xticklabels=['Témoin', 'Malade'], yticklabels=['Témoin', 'Malade'])
    plt.title('Matrice de Confusion - Forêt Aléatoire')
    plt.xlabel('Prédictions')
    plt.ylabel('Véritables')
    cm_rf_fig = "cm_rf.png"
    plt.savefig(cm_rf_fig)
    plt.close()

    plt.figure(figsize=(10, 7))
    plt.plot(fpr_log_reg, tpr_log_reg, lw=2, label=f'Régression Logistique (AUC = {roc_auc_log_reg:.2f})')
    plt.plot(fpr_rf, tpr_rf, lw=2, label=f'Forêt Aléatoire (AUC = {roc_auc_rf:.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Taux de Faux Positifs')
    plt.ylabel('Taux de Vrais Positifs')
    plt.title('Courbes ROC')
    plt.legend(loc='lower right')
    roc_fig = "roc_curve.png"
    plt.savefig(roc_fig)
    plt.close()

    plt.figure(figsize=(10, 7))
    feature_importances.plot(kind='bar')
    plt.title('Importance des Caractéristiques - Random Forest')
    plt.xlabel('Caractéristiques')
    plt.ylabel('Importance')
    feature_importances_fig = "feature_importances.png"
    plt.savefig(feature_importances_fig)
    plt.close()

    # Création du rapport PDF
    report = f"Régression Logistique - Accuracy: {accuracy_log_reg}\n{report_log_reg}\n\n" \
             f"Forêt Aléatoire - Accuracy: {accuracy_rf}\n{report_rf}\n"

    base_name = os.path.basename(data_path)
    file_name = os.path.splitext(base_name)[0]
    output_pdf_path = f"result_machine_learning_{file_name}.pdf"

    generate_pdf(report, cm_log_reg_fig, cm_rf_fig, roc_fig, feature_importances_fig, output_pdf_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Cellular data analysis.')
    parser.add_argument('--data_path', required=True, help='Path to the data file.')
    parser.add_argument('--meta_data_path', required=True, help='Path to the metadata file.')
    parser.add_argument('--sno_list_path', required=True, help='Path to the sno list file.')

    args = parser.parse_args()

    main(args.data_path, args.meta_data_path, args.sno_list_path)



#python .\machin_learning_wgcna.py --data "D:\data_umr_1227\Par_pureté_cyto_digitale\BLymphocytes_DESeq\Sizefactornorm_all_sampleBLymphocytes_DESeq.csv" 
#--meta 'D:/data_umr_1227/Par_pureté_cyto_digitale/METADATA_SORTED_CELL_MCPcounter.tsv' 
#--sno "D:/data_umr_1227/Par_pureté_cyto_digitale/WGCNA_analysis/Sno_gene_list.csv"