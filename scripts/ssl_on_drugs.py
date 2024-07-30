import numpy as np

def calculate_ssl_drugs(ppi_net, DrT_matrix, mu=0.5):
    diag_ppi = np.diag(np.sum(ppi_net, 1))
    laplacian_ppi = diag_ppi - ppi_net
    inv_lap_ppi = np.linalg.inv(np.eye(len(ppi_net)) + mu * laplacian_ppi)
    
    # Combine PPI and drug target matrices
    DrNet = metrics.pairwise.cosine_similarity(DrT_matrix)
    upper_ = np.hstack([ppi_net, DrT_matrix.transpose()])
    lower_ = np.hstack([DrT_matrix, DrNet])
    ppi_Dr_net = np.vstack([upper_, lower_])
    
    diag_ppi_dr = np.diag(np.sum(ppi_Dr_net, 1))
    lap_ppi_dr = diag_ppi_dr - ppi_Dr_net
    inv_lap_ppi_dr = np.linalg.inv(np.eye(len(lap_ppi_dr)) + mu * lap_ppi_dr)
    
    return inv_lap_ppi_dr

def calculate_drug_scores(drug_name, final_disease, final_gene, inv_lap_ppi_dr, gene_list, RA_index):
    drug_score_on_ppi = np.zeros([len(drug_name), len(final_disease)])
    
    for i in range(len(final_disease)):
        target_idx = RA_index
        alter_idx = i
        
        target_gene = final_gene[target_idx]
        alter_gene = final_gene[alter_idx]
        
        union_gene = np.union1d(target_gene, alter_gene)
        union_gene_idx = np.where(np.isin(gene_list, union_gene))[0]
        
        y_ppi = np.zeros([len(gene_list), 1])
        y_ppi[union_gene_idx] = 1
        y_drug = np.zeros([len(drug_name), 1])
        
        y_concat = np.vstack([y_ppi, y_drug])
        
        f_ppi_drug = inv_lap_ppi_dr @ y_concat
        f_drug = f_ppi_drug[len(y_ppi):len(y_concat) + 1]
        f_drug = (f_drug - np.min(f_drug)) / (np.max(f_drug) - np.min(f_drug))
        f_drug = f_drug.reshape((-1, 1))
        
        drug_score_on_ppi[:, i] = f_drug.transpose()
        
    return drug_score_on_ppi
