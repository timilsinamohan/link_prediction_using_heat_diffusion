
__author__ = 'mohan'
import pandas as pd
import networkx as nx
from networkx.algorithms.bipartite.matrix import biadjacency_matrix
import numpy as np
from sklearn.metrics import precision_recall_curve, auc
from sklearn.preprocessing import normalize
import random
from sklearn import metrics
import time


def compute_score(TG):

    iteration = 6
    lambda_diff = 1.0
    #recon_mat = np.zeros([m,n])
    print "for iteration:",iteration
    #preparing R matrix ####
    gamma = 0.85
    H = B.copy()
    #print H
    g = (1.0/n * np.ones(n))
    g = g.reshape(n,1)
    R = gamma * H + (1.00-gamma) * g
    I = np.eye(n,n,dtype=np.float64)
    V = I + (lambda_diff/iteration) * R
    state_matrix = TG.copy()


    for j in xrange(iteration):
        state_matrix_new = V.dot(state_matrix.T).T
        state_matrix = state_matrix_new.copy()

    return state_matrix


def innerfold(IDX,m,n):
    mask_idx = np.unravel_index(IDX, (m, n))
    tumor_gene_relation_copy = matrix.copy()
    target_idx = np.unravel_index(IDX, (m, n))
    #print "Inside INNERFOLD \n",tumor_gene_relation_copy

    for i in range(len(mask_idx[0])):
        tumor_gene_relation_copy[mask_idx[0][i], mask_idx[1][i]] = 0

    score = compute_score(tumor_gene_relation_copy)

    score = normalize(score, norm='l2')

    Ground_Truth = matrix.copy()
    Ground_Truth = np.array(Ground_Truth)



    #print Ground_Truth[target_idx],"\n",score[target_idx]
    prec, recall, _ = precision_recall_curve(Ground_Truth[target_idx],score[target_idx])
    print "AUC-PR", auc(recall, prec)

    fpr, tpr, threshold = metrics.roc_curve(Ground_Truth[target_idx],score[target_idx])
    #fpr, tpr, threshold = metrics.roc_curve(actual, predicted)
    roc_auc = metrics.auc(fpr, tpr)

    print "AUC-ROC score:",roc_auc

    return roc_auc


df1 = pd.read_csv("tumor_interaction_gene.txt",sep = ",")
tumor_nodes_col = df1["tumor"]
gene_nodes_col =df1["gene"]

#print tumor_nodes_col,gene_nodes_col
edgelist = zip(tumor_nodes_col,gene_nodes_col)
B = nx.Graph()
B.add_nodes_from(tumor_nodes_col, bipartite=0)
B.add_nodes_from(gene_nodes_col, bipartite=1)
B.add_edges_from(edgelist)
tum_nodes = set(n for n,d in B.nodes(data=True) if d['bipartite']==0)
gene_nodes =set(n for n,d in B.nodes(data=True) if d['bipartite']==1)

#print tum_nodes, gene_nodes
matrix = biadjacency_matrix(B, row_order= tum_nodes, column_order=gene_nodes)
matrix = matrix.A

print "Tumor Gene Shape:",matrix.shape

m = matrix.shape[0]
n = matrix.shape[1]

G = nx.read_edgelist("fusion_scores.txt", delimiter= " ",
                     nodetype=str,
                     data=(('weight',float),))

print "total number of nodes and edges:",G.number_of_nodes(),G.number_of_edges()


Gene_Gene_Adj_mat = nx.adjacency_matrix(G, nodelist=gene_nodes,weight='none')

A = np.array(Gene_Gene_Adj_mat.todense(), dtype=np.float64)

weight_matrix = nx.attr_matrix(G, edge_attr='weight', rc_order=gene_nodes)
weight_matrix = np.array(weight_matrix)

heat_matrix = np.zeros([n,n])
#print "Heat Matrix Creation started:"
gene_list = list(gene_nodes)
G= nx.from_numpy_matrix(A)

print "Heat Matrix filling started:"
for i in range(n):
    for j in range(n):
         if A[j,i] == 1.0:

            heat_matrix[i,j] = weight_matrix[j,i]/G.degree(j)

         if (i==j):

            if G.degree(i):
                heat_matrix[i,j] = (-1.0 / G.degree(i)) * sum(weight_matrix[i,:])


print "Heat Matrix Completed:"

B = heat_matrix.copy()

FOLDS = 10
sz = m * n
IDX = list(range(sz))
fsz = int(sz/FOLDS)
np.random.shuffle(IDX)
offset = 0
print "Fold size",fsz
AUC_test = np.zeros(FOLDS)

for f in xrange(FOLDS):
    print "Fold:",f
    start_time = time.time()
    IDX1 = random.sample(xrange(sz),fsz)
    AUC_test[f] = innerfold(IDX1,m,n)
    print("--- %s seconds ---" % (time.time() - start_time))
    offset += fsz
print "Mean AUC-PR",AUC_test.mean(), AUC_test





