import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx

G = nx.Graph()
edges = [('COMMD1','COMMD2'),('COMMD1','COMMD3'),('COMMD1','COMMD9'),('COMMD1','CCDC22'),('COMMD1','CCDC93'),('COMMD2','COMMD3'),('COMMD2','COMMD5'),('COMMD2','COMMD9'),('COMMD2','CCDC22'),('COMMD2','CCDC93'),('COMMD2','C16ORF62'),('COMMD3','COMMD6'),('COMMD3','CCDC22'),('COMMD3','CCDC93'),('COMMD3','C16ORF62'),('COMMD5','COMMD6'),('COMMD5','CCDC22'),('COMMD5','CCDC93'),('COMMD6','COMMD9'),('COMMD6','CCDC22'),('COMMD6','C16ORF62'),('CCDC22','CCDC93'),('CCDC22','C16ORF62'),('CCDC93','C16ORF62')]
G.add_edges_from(edges)
pos = nx.circular_layout(G)
labels = {}
li = ['COMMD1','COMMD2','COMMD3','COMMD5','COMMD6','COMMD9','CCDC22','CCDC93','C16ORF62']
for s in li:
    labels[s] = s

nx.draw(G, pos, node_size=1500, node_color='lightgreen', font_size=8, font_weight='bold')
nx.draw_networkx_labels(G, pos, labels, font_size=8)

plt.savefig('commd9.png')
