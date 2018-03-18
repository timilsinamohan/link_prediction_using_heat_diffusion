
import networkx as nx
import pandas as pd
from networkx.algorithms.bipartite.matrix import biadjacency_matrix
import random
import numpy as np
from sklearn.metrics import precision_recall_curve, auc
from sklearn import metrics
from networkx.algorithms import bipartite



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


def jaccard_similarity(TG_graph,gen,tumor_gene_relation_bi_adj_mat):

    #for jaccard similarity nearest neighbors#

    G = bipartite.overlap_weighted_projected_graph(TG_graph, gen)
    print "calculating the link prediction scores",G.number_of_edges()
    A = nx.to_scipy_sparse_matrix(G)
    scores = np.dot(tumor_gene_relation_bi_adj_mat,A)
    return  scores




def nearest_neighbors_weighted_projection(TG_graph,gen,tumor_gene_relation_bi_adj_mat):
    #for nearest neighbors#

    #G = bipartite.weighted_projected_graph(TG_graph, gen)
    #for jaccard similarity nearest neighbors#

    #G = bipartite.overlap_weighted_projected_graph(TG_graph, gen)

    #for other baseline link prediction alogrithm#
    G = bipartite.projected_graph(TG_graph, gen)
    print "calculating the link prediction scores",G.number_of_edges()
    preds = nx.preferential_attachment(G, ebunch=G.edges())

    pred_scores_graph = nx.Graph()
    pred_scores_graph.add_nodes_from(gen)
    cnt = 0
    for u,v,p in preds:
        pred_scores_graph.add_edge(u,v,weight = p)

    print cnt



    print "sorted !"
    #G = pred_scores_graph
    A = nx.to_scipy_sparse_matrix(G)
    scores = np.dot(tumor_gene_relation_bi_adj_mat,A)
    return scores


def innerfold(IDX1,m,n):
    mask_idx = np.unravel_index(IDX1, (m, n))
    tumor_gene_relation_copy = matrix.copy()
    target_idx = np.unravel_index(IDX1, (m, n))

    for i in range(len(mask_idx[0])):
        tumor_gene_relation_copy[mask_idx[0][i], mask_idx[1][i]] = 0

    TG,tum,gen = graph_from_biadjacency_matrix(tumor_gene_relation_copy)
    scores = nearest_neighbors_weighted_projection(TG,gen,tumor_gene_relation_copy)
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
