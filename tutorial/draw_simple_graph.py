import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import sys

def read_input(FNAME):
    input_edges = []
    with open(FNAME) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            assert(len(li) == 2)
            input_edges.append((li[0],li[1]))
    return input_edges

fname = sys.argv[1]
ofname = sys.argv[2]
nG = int(sys.argv[3])
col = sys.argv[4]

#assert(int(dup) >= 1 and int(dup) <= nG)
#assert(int(anchor) >= 1 and int(anchor) <= nG)

#ofname = fname + '.png'

input_edges = read_input(fname)

G = nx.Graph()
G.add_edges_from(input_edges)
print(G.edges)

pos = nx.circular_layout(G)
#labels = {}
#li = ['ATF2','ATF3','ATF4','BATF','CREB','JUN','CEBP','ATF6','XBP1','E4BP4','OASIS','PAR','FOS']
#for s in li:
#    labels[s] = s

nx.draw(G, pos, with_labels=True, node_size=2500, node_color=col, font_size=30, font_weight='bold')
#nx.draw_networkx_labels(G, pos, labels, font_size=10)

plt.savefig(ofname)
