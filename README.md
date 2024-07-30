# Context-specific Drug Scoring for Precision Re-purposing via Synergizing Label Propagation with Heterogeneous Networks

## Abstract
While biological and pharmaceutical knowledge networks have significantly propelled drug repurposing efforts, reliance solely on these networks is insufficient for accurately addressing genetic and phenotypic variances. This limitation highlights the need for an integrative approach that leverages context-specific data to enhance the precision of drug repurposing. We introduce a network-based integrative drug scoring approach that synergistically incorporates data-driven and knowledge-driven networks without the necessity of their direct integration. We developed synergistic label propagation algorithms that facilitate information transfer from data-driven to knowledge-driven networks. To enable context-specific drug repurposing, we constructed a data-driven disease-disease association network utilizing European-specific genetic information and a knowledge-driven drug-target protein association network. In a proof-of-concept study, drug scoring was applied to identify candidate drugs for rheumatoid arthritis, asthma, and multiple sclerosis.

## Directory Structure
- `data/`: Contains all the necessary data files.
  - `final_disease_snp_gene.pkl`
  - `final_disease_pruned_snp_matrix.pkl`
  - `combined_ppi_net.pkl.gz`
  - `STRING_DB/used_full_protein_symbol.csv`
  - `drug_target_matrix.pkl.gz`
  - `drugbank/drugbank_human.csv`
  - `drug_ground_truth.csv`
- `scripts/`: Contains Python scripts for loading data, performing SSL, and evaluating results.
  - `load_data.py`
  - `ssl_on_diseases.py`
  - `ssl_on_drugs.py`
  - `evaluation.py`
- `README.md`: Project description and instructions.
- `main.py`: Main script to run the project.

## Setup
1. Install necessary packages:
    ```sh
    pip install -r requirements.txt
    ```

2. Ensure that all data files are in the `data/` directory as shown in the directory structure.

## Usage
Run the main script to perform SSL and evaluation:
```sh
python main.py
