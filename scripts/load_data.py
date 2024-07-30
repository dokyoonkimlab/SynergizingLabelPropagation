import gzip
import pickle
import pandas as pd

def load_pickle(file_path):
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def load_compressed_pickle(file_path):
    with gzip.open(file_path, 'rb') as f:
        return pickle.load(f)

def load_csv(file_path):
    return pd.read_csv(file_path, index_col=0)

# Load data
def load_data():
    data = {}
    data['disease_snp_gene'] = load_pickle('data/final_disease_snp_gene.pkl')
    data['disease_snp_matrix'] = load_compressed_pickle('data/final_disease_pruned_snp_matrix.pkl.gz')
    data['ppi_net'] = load_compressed_pickle('data/combined_ppi_net.pkl.gz')
    data['gene_list'] = load_csv('data/STRING_DB/used_full_protein_symbol.csv').to_numpy().squeeze()
    data['drug_target_matrix'] = load_compressed_pickle('data/drug_target_matrix.pkl.gz')
    data['drug_info'] = load_csv('data/drugbank/drugbank_human.csv')
    data['drug_ground_truth'] = load_csv('data/drug_ground_truth.csv')
    return data
