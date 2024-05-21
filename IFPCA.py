### This file contain codes to implement IF-PCA 

from sklearn.cluster import SpectralClustering

import numpy as np
import scipy
from scipy.optimize import linear_sum_assignment
from scipy.linalg import svd
from scipy.special import erf
from scipy.stats import norm
from sklearn.metrics import accuracy_score
from sklearn.cluster import KMeans

def normal_cdf(x,mu,std):
    return 0.5 * (1 + erf((x - mu) / (std * np.sqrt(2))))

def estimate_Gaussian_parameter (g_samples):
    return np.mean(g_samples), np.std(g_samples)


def calculate_single_gene_gaussian_ks (sample):
    '''
    This function assumes data follows Gaussian distribution.
    It calcualte KS statistics by fitting a Gaussian distribution,
    then compare ECDF to the CDF for fitted distribution.
    '''

    n = len(sample)
    mu,std = estimate_Gaussian_parameter (sample)

    sorted_sample = sorted(sample)
    
    # Calculate the empirical distribution function (EDF)
    edf = [i / n for i in range(1, n + 1)]
    
    # Calculate the cumulative distribution function (CDF) using the provided function
    cdf = [normal_cdf(x,mu,std) for x in sorted_sample]
    
    # Compute the absolute differences between EDF and CDF
    differences = [abs(edf_i - cdf_i) for edf_i, cdf_i in zip(edf, cdf)]
    
    # Identify the maximum absolute difference as the KS statistic
    ks_statistic = np.sqrt(n)*max(differences)
    return ks_statistic

#Get HC point,
def HCT(x, alpha, n,param_size):
    x = (x-np.mean(x))/np.std(x)
    m = round(param_size*alpha)
    pr = np.zeros(param_size)
    pr = 1 - norm.cdf(x)
    pr = np.sort(pr)
    s = (1/param_size)*(1 + np.arange(param_size))
    HC = np.zeros(m)
    HC = np.divide(np.sqrt(param_size)*(s[0:m] - pr[0:m]),np.sqrt(s[0:m] + np.maximum(0, np.sqrt(n)*(s[0:m] - pr[0:m]))))
    ind = np.where(HC == np.max(HC))[0]
    L = min(ind)
    return L


def get_IFPCA_acc (X,y,K,effron=True, manifold = False):
    ''' 
    This function inputs n by p matrix X, number of clusters K, and true label.
    Output accuracy score after performing IF-PCA.
    The effron parameter determines whether one wants to perform the effron null step.
    The manifold parameter checks if the method is part of the IFPCA+ or not.
    If manifold = True, implement IF-PCA with no normalization on X and a small modification to PCA .
    '''
    if not manifold:
        # Normalize each column
        mean = np.mean(X, axis=0)
        std = np.std(X, axis=0)
        X_normalized = (X - mean) / std
        X=X_normalized

    print("Running IF steps")
    n, param_size = np.shape(X)  
    phi_nj = [calculate_single_gene_gaussian_ks(X[:,i]) for i in range(param_size)]
    phi_nj=np.array(phi_nj)
    norm_phi = phi_nj
    if effron:
        norm_phi = (phi_nj-np.mean(phi_nj))/np.std(phi_nj)
    k = HCT(norm_phi, 0.5, n, param_size)
    #get top k in k-s score
    sort_phi = np.argsort(norm_phi)
    ind_select_phi= sort_phi[-k:]
    X_selected = X[:,ind_select_phi]
    IF_accs = []
    if manifold:
        print("Running PCA (for manifold)")
        for _ in range(5):
            acc = get_acc_complex_PCA(X_selected,y, K)
            IF_accs.append(acc)
        print(IF_accs)
        avg_acc = sum(IF_accs) / len(IF_accs)
    else:
        print("Running PCA")
        for _ in range(5):
            acc = get_acc_PCA(X_selected,y, K)
            IF_accs.append(acc)
        print(IF_accs)
        avg_acc = sum(IF_accs) / len(IF_accs)
    return avg_acc



############ Helper functions

def best_clustering_error(true_labels, predicted_labels):
    """
    Calculates the best clustering error    across all permutations of cluster assignments.
    """

    # Find the best matching between true labels and predicted labels
    cost_matrix = np.zeros((len(np.unique(true_labels)), len(np.unique(predicted_labels))))
    for i, true_label in enumerate(np.unique(true_labels)):
        for j, predicted_label in enumerate(np.unique(predicted_labels)):
            common_elements = np.sum((true_labels == true_label) & (predicted_labels == predicted_label))
            cost_matrix[i, j] = -common_elements  # negative since linear_sum_assignment finds the minimum cost

    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Map predicted labels to true labels based on the best matching
    mapped_labels = np.zeros_like(predicted_labels)
    for true_label, col_index in zip(np.unique(true_labels), col_ind):
        mapped_labels[predicted_labels == np.unique(predicted_labels)[col_index]] = true_label

    # Calculate accuracy using mapped labels
    accuracy = accuracy_score(true_labels, mapped_labels)
    return accuracy

def get_acc_complex_PCA (data, y, K):
    U, s, VT = svd(data)
    K0 = K
    if K<4:
        K0=4
    PCA_input = U[:,0:K0]
    Cluster = KMeans(n_clusters=K)
    Cluster.fit(PCA_input)
    y_pred = Cluster.predict(PCA_input)
    acc = best_clustering_error(y_pred,y)
    return acc

def get_acc_PCA (data, y, K):
    U, s, VT = svd(data)        
    PCA_input = U[:,0:K-1]
    Cluster = KMeans(n_clusters=K)
    Cluster.fit(PCA_input)
    y_pred = Cluster.predict(PCA_input)
    acc = best_clustering_error(y_pred,y)
    return acc



