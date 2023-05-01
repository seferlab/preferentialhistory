import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx

G = nx.Graph()
G.add_edges_from([('ATF2','ATF3'),('ATF2','ATF4'),('ATF2','BATF'),('ATF2','CREB'),('ATF2','JUN'),('ATF4','CEBP'),('ATF6','XBP1'),('BATF','CEBP'),('BATF','CREB'),('BATF','E4BP4'),('BATF','JUN'),('BATF','OASIS'),('BATF','PAR'),('CEBP','PAR'),('CREB','XBP1'),('FOS','JUN')])

pos = nx.circular_layout(G)
labels = {}
li = ['ATF2','ATF3','ATF4','BATF','CREB','JUN','CEBP','ATF6','XBP1','E4BP4','OASIS','PAR','FOS']
for s in li:
    labels[s] = s

nx.draw(G, pos, node_size=1500, node_color='lightblue', font_size=10, font_weight='bold')
nx.draw_networkx_labels(G, pos, labels, font_size=10)

plt.savefig('bzip.png')
