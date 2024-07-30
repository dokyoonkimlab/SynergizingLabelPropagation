import numpy as np
from sklearn import metrics
import pandas as pd

def evaluate(drug_score_on_ppi, disease_contribution, truth_RA):
    computed_drug_score = drug_score_on_ppi @ disease_contribution
    
    fpr, tpr, _ = metrics.roc_curve(truth_RA, computed_drug_score, pos_label=1)
    auc_score = metrics.auc(fpr, tpr)
    
    # Approved drugs evaluation
    idx_truth_approved = np.where(truth_RA == 1)[0]
    fpr_approved, tpr_approved, _ = metrics.roc_curve(truth_RA[idx_truth_approved], computed_drug_score[idx_truth_approved], pos_label=1)
    auc_approved = metrics.auc(fpr_approved, tpr_approved)
    
    # Precision at top 10%
    sorted_indices = np.argsort(computed_drug_score)[::-1]
    sorted_drug_scores = computed_drug_score[sorted_indices]
    sorted_truth_RA = truth_RA[sorted_indices]
    
    top_10_percent_index = int(len(sorted_drug_scores) * 0.1)
    top_10_percent_truth_RA = sorted_truth_RA[:top_10_percent_index]
    
    precision_at_top_10_percent = np.sum(top_10_percent_truth_RA) / len(top_10_percent_truth_RA)
    
    num_drugs_related_to_RA = np.sum(truth_RA)
    total_num_drugs = len(truth_RA)
    baseline_precision_at_top_10_percent = num_drugs_related_to_RA / total_num_drugs
    
    # Recall at top 10%
    num_drugs_related_to_RA_at_top_10_percent = np.sum(top_10_percent_truth_RA)
    recall_at_top_10_percent = num_drugs_related_to_RA_at_top_10_percent / num_drugs_related_to_RA
    baseline_recall_at_top_10_percent = num_drugs_related_to_RA / total_num_drugs
    
    return {
        "auc_all": auc_score,
        "auc_approved": auc_approved,
        "precision_at_top_10_percent": precision_at_top_10_percent,
        "baseline_precision_at_top_10_percent": baseline_precision_at_top_10_percent,
        "recall_at_top_10_percent": recall_at_top_10_percent,
        "baseline_recall_at_top_10_percent": baseline_recall_at_top_10_percent
    }

def save_results(drug_score_on_ppi, filename='drug_score_RA_heatmap.csv'):
    df_save = pd.DataFrame(drug_score_on_ppi)
    df_save.to_csv(filename)
