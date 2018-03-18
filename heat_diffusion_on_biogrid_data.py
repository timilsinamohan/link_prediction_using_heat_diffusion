
__author__ = 'mohan'
import pandas as pd
import networkx as nx
from networkx.algorithms.bipartite.matrix import biadjacency_matrix
import numpy as np
from sklearn.metrics import precision_recall_curve, auc
from sklearn.preprocessing import normalize
from sklearn import metrics
import time


def compute_score(TG):

    ##############################

    iteration = 6
    lambda_diff = 1.0
    #preparing R matrix ####
    gamma = 0.85

    H = A.copy()
    g = (1.0/n * np.ones(n))
    g = g.reshape(n,1)
    R = gamma * H + (1.00-gamma) * g
    I = np.eye(n,n,dtype=np.float64)
    V = (I + (lambda_diff/iteration) * R )

    state_matrix = TG.copy()
    print "computing started"
    for j in xrange(iteration):
        state_matrix = V.dot(state_matrix.T).T

        state_matrix = state_matrix.copy()


    return state_matrix


def innerfold(IDX,m,n):
    mask_idx = np.unravel_index(IDX, (m, n))
    tumor_gene_relation_copy = matrix.copy()
    target_idx = np.unravel_index(IDX, (m, n))


    for i in range(len(mask_idx[0])):
        tumor_gene_relation_copy[mask_idx[0][i], mask_idx[1][i]] = 0

    print "masking done:"
    score = compute_score(tumor_gene_relation_copy)
    print "shape of the final score:",score.shape

    score = normalize(score, norm='l2')

    Ground_Truth = matrix.copy()
    Ground_Truth = np.array(Ground_Truth)

    prec, recall, _ = precision_recall_curve(Ground_Truth[target_idx],score[target_idx])
    print "AUC-PR", auc(recall, prec)


    fpr, tpr, threshold = metrics.roc_curve(Ground_Truth[target_idx],score[target_idx])
    #fpr, tpr, threshold = metrics.roc_curve(actual, predicted)
    roc_auc = metrics.auc(fpr, tpr)

    print "AUC-ROC score:",roc_auc

    return roc_auc




####This is tumor-Gene Part
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

#print matrix.shape

m = matrix.shape[0]
n = matrix.shape[1]



G = nx.read_edgelist("physical_interaction_gene.txt", nodetype=str)
#GG.add_edges_from(edgelist)
Gene_Gene_Adj_mat = nx.adjacency_matrix(G, nodelist=gene_nodes)
print "information of this physical network:",nx.info(G)

A = np.array(Gene_Gene_Adj_mat.todense(), dtype=np.float64)

A = normalize(A, axis= 1, norm= 'l1')

A = A.transpose()

deg = []
for i in G.nodes():
    deg.append(G.degree(i))

np.fill_diagonal(A,-1*deg)



FOLDS = 10
sz = m * n
IDX = list(range(sz))
fsz = int(sz/FOLDS)
np.random.shuffle(IDX)
offset = 0
print "Fold size",fsz
AUC_test = np.zeros(FOLDS)

for f in xrange(FOLDS):
    #idx_test = IDX[offset:offset + fsz]
    IDX1 = IDX[offset:offset + fsz]
    start_time = time.time()
    AUC_test[f] = innerfold(IDX1,m,n)
    print("--- %s seconds ---" % (time.time() - start_time))

    offset += fsz
print "Mean AUC-PR",AUC_test.mean(), AUC_test




