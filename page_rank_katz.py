
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

def graph_from_biadjacency_matrix(M):

    print "Bi adjacency matrix started:"
    #print "am i different matrix:\n",M
    # Give names to the nodes in the two node sets
    U = [ "u{}".format(i) for i in range(M.shape[0]) ]
    V = [ "v{}".format(i) for i in range(M.shape[1]) ]

    # Create the graph and add each set of nodes
    G = nx.Graph()
    G.add_nodes_from(U, bipartite=0)
    G.add_nodes_from(V, bipartite=1)

    G.add_edges_from([(U[i], V[j]) for i, j in zip(*M.nonzero())])

    print "Bi Partite graph created:"
    return G,U



def katz_centrality_score(TG):
    iteration = 4
    state_matrix = TG.copy()

    Adj_mat = Gene_Gene_Adj_mat.copy()

    alpha = 0.15
    beta = 0.0001
    for j in range(iteration):
        state_matrix_next = Adj_mat.dot(alpha * state_matrix.T).T + beta
        state_matrix = state_matrix_next.copy()



    return state_matrix


def personalized_pageRank(TG):

    iteration = 6

    #preparing R matrix ####
    gamma = 0.85
    ####get the transition matrix####
    H = normalize(A.T, norm='l1', axis=0)


    g = (1.0/n * np.ones(n))
    g = g.reshape(n,1)
    R = gamma * H + (1.00-gamma) * g


    state_matrix = TG.copy()
    ####vectorized the operation of matrix multiplication#######

    for j in range(iteration):
        state_matrix_next = R.dot(state_matrix.T).T
        state_matrix = state_matrix_next.copy()

    return state_matrix



def innerfold(IDX,m,n):
    mask_idx = np.unravel_index(IDX, (m, n))
    tumor_gene_relation_copy = matrix.copy()
    target_idx = np.unravel_index(IDX, (m, n))

    for i in range(len(mask_idx[0])):
        tumor_gene_relation_copy[mask_idx[0][i], mask_idx[1][i]] = 0

    #score = katz_centrality_score(tumor_gene_relation_copy)
    score = personalized_pageRank(tumor_gene_relation_copy)

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

# G = nx.read_edgelist("fusion_scores.txt", delimiter= " ",
#                      nodetype=str,
#                      data=(('weight',float),))

G = nx.read_edgelist("physical_interaction_gene.txt", nodetype=str)
print "total number of nodes and edges:",G.number_of_nodes(),G.number_of_edges()


Gene_Gene_Adj_mat = nx.adjacency_matrix(G, nodelist=gene_nodes,weight='none')

A = np.array(Gene_Gene_Adj_mat.todense(), dtype=np.float64)

# weight_matrix = nx.attr_matrix(G, edge_attr='weight', rc_order=gene_nodes)
# weight_matrix = np.array(weight_matrix)

heat_matrix = np.zeros([n,n])

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





