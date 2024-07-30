from scripts.load_data import load_data
from scripts.ssl_on_diseases import calculate_ssl_diseases, calculate_disease_scores
from scripts.ssl_on_drugs import calculate_ssl_drugs, calculate_drug_scores
from scripts.evaluation import evaluate, save_results
import numpy as np

def main():
    # Load data
    print('dataLoaded')
    data = load_data()
    
    # Disease SSL
    inv_lap_DDN = calculate_ssl_diseases(data['disease_snp_matrix'])
    result_score_f = calculate_disease_scores(data['disease_snp_gene']['final_disease'], inv_lap_DDN)
    
    # Find RA index
    disease_list = data['disease_snp_gene']['final_disease']
    disease_list = np.array(disease_list)
    RA_index = np.where(disease_list == 'Phe_714.1')[0][0]  # RA
    
    # Normalize RA scores
    f_RA = result_score_f[:, RA_index]
    f_RA[f_RA == np.max(f_RA)] = 10
    second_max_value = np.partition(f_RA, -2)[-2]
    f_RA[f_RA != 10] = (f_RA[f_RA != 10] - np.min(f_RA)) / (second_max_value - np.min(f_RA))
    f_RA[RA_index] = 1
    
    # Drug SSL
    inv_lap_ppi_dr = calculate_ssl_drugs(data['ppi_net'], data['drug_target_matrix'])
    drug_score_on_ppi = calculate_drug_scores(data['drug_info']['drug_name'].to_numpy(), disease_list, data['disease_snp_gene']['final_gene'], inv_lap_ppi_dr, data['gene_list'], RA_index)
    
    # Evaluate
    results = evaluate(drug_score_on_ppi, f_RA, data['drug_ground_truth']['Rheumatoid_Arthritis'].to_numpy())
    
    # Save results
    save_results(drug_score_on_ppi)
    
    # Print results
    for key, value in results.items():
        print(f"{key}: {value:.4f}")

if __name__ == "__main__":
    main()
