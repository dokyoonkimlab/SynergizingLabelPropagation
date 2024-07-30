import numpy as np
from sklearn import metrics

def calculate_ssl_diseases(DpS_matrix, mu=0.5):
    similarity_DDN = metrics.pairwise.cosine_similarity(DpS_matrix)
    diag_DDN = np.diag(np.sum(similarity_DDN, 1))
    laplacian_DDN = diag_DDN - similarity_DDN
    inv_lap_DDN = np.linalg.inv(np.eye(len(DpS_matrix)) + mu * laplacian_DDN)
    return inv_lap_DDN

def calculate_disease_scores(final_disease, inv_lap_DDN):
    result_score_f = np.zeros([len(final_disease), len(final_disease)])
    for i in range(len(final_disease)):
        y = np.zeros(len(final_disease))
        y[i] = 1
        f = inv_lap_DDN @ y
        result_score_f[:, i] = f
    return result_score_f
