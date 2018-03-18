
import networkx as nx
import pandas as pd
from networkx.algorithms.bipartite.matrix import biadjacency_matrix
import random
import numpy as np
from sklearn.metrics import precision_recall_curve, auc
from sklearn import metrics
from networkx.algorithms import bipartite
from sklearn.preprocessing import normalize
from scipy import sparse


def graph_from_biadjacency_matrix(M):


    ############################

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

    return G,U,V


def resource_allocation_index(TG_graph,gen,tumor_gene_relation_bi_adj_mat):

    G = bipartite.projected_graph(TG_graph, gen)
    A = nx.to_scipy_sparse_matrix(G)
    degrees = A.sum(axis=0)

    with np.errstate(divide='ignore'):
        weights = sparse.csr_matrix(1./(degrees)) ###to avoid divide by 0
    AA =  A.multiply(weights) * A.T

    tumor_gene_relation_bi_adj_mat = sparse.csr_matrix(tumor_gene_relation_bi_adj_mat)

    scores = np.dot(tumor_gene_relation_bi_adj_mat,AA)
    return scores


def adamic_adar(TG_graph,gen,tumor_gene_relation_bi_adj_mat):

    G = bipartite.projected_graph(TG_graph, gen)
    A = nx.to_scipy_sparse_matrix(G)
    degrees = A.sum(axis=0)

    with np.errstate(divide='ignore'):
        weights = sparse.csr_matrix(1./np.log(degrees)) ###to avoid divide by 0 
    AA =  A.multiply(weights) * A.T

    tumor_gene_relation_bi_adj_mat = sparse.csr_matrix(tumor_gene_relation_bi_adj_mat)

    scores = np.dot(tumor_gene_relation_bi_adj_mat,AA)
    return scores


def innerfold(IDX1,m,n):
    mask_idx = np.unravel_index(IDX1, (m, n))
    tumor_gene_relation_copy = matrix.copy()
    target_idx = np.unravel_index(IDX1, (m, n))

    for i in range(len(mask_idx[0])):
        tumor_gene_relation_copy[mask_idx[0][i], mask_idx[1][i]] = 0

    TG,tum,gen = graph_from_biadjacency_matrix(tumor_gene_relation_copy)
    #scores = adamic_adar(TG,gen,tumor_gene_relation_copy)
    scores = resource_allocation_index(TG,gen,tumor_gene_relation_copy)
    scores = scores.A

    #scores = normalize(scores.astype(np.float64), norm='l2')

    Ground_Truth = matrix.copy()
    Ground_Truth = np.array(Ground_Truth)
    prec, recall, _ = precision_recall_curve(Ground_Truth[target_idx],scores[target_idx])
    print "AUC-PR score:",auc(recall, prec)
    fpr, tpr, threshold = metrics.roc_curve(Ground_Truth[target_idx],scores[target_idx])
    roc_auc = metrics.auc(fpr, tpr)

    print "AUC-ROC score:",roc_auc

    return roc_auc,auc(recall, prec)

df1 = pd.read_csv("tumor_interaction_gene.txt",sep = ",")
tumor_nodes_col = df1["tumor"]
gene_nodes_col =df1["gene"]

#print tumor_nodes_col,gene_nodes_col
edgelist = zip(tumor_nodes_col,gene_nodes_col)
B = nx.DiGraph()
B.add_nodes_from(tumor_nodes_col, bipartite=0)
B.add_nodes_from(gene_nodes_col, bipartite=1)
B.add_edges_from(edgelist)
tum_nodes = set(n for n,d in B.nodes(data=True) if d['bipartite']==0)
gene_nodes =set(n for n,d in B.nodes(data=True) if d['bipartite']==1)


matrix = biadjacency_matrix(B, row_order= tum_nodes, column_order=gene_nodes)
matrix = matrix.A

m = matrix.shape[0]
n = matrix.shape[1]
FOLDS = 10
sz = m * n
fsz = int(sz/FOLDS)
#np.random.shuffle(IDX)
offset = 0
print "Fold size",fsz
AUC_test = np.zeros(FOLDS)
AUC_roc_test = np.zeros(FOLDS)

for f in xrange(FOLDS):
    print "Fold:",f
    IDX1 = random.sample(xrange(sz),fsz)
    roc,pr = innerfold(IDX1,m,n)
    AUC_test[f] = pr
    AUC_roc_test[f] = roc
    offset += fsz

print "Mean AUC-PR and standard deviation",AUC_test.mean(), AUC_test.std()
print "Mean AUC-ROC",AUC_roc_test.mean(), AUC_roc_test.std()
