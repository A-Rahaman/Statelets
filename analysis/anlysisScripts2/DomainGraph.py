import matplotlib.pyplot as plt
import networkx as nx
import scipy.io
import numpy as np

allHC1 = scipy.io.loadmat('HC_eWeights.mat')
allHC  = allHC1['e_weights']
X,Y = allHC.shape
edges1 = scipy.io.loadmat('DOmainedges_string.mat')
edges = edges1['Edges']
#print(Y)
trans = np.zeros(Y)
for hc in range(0, Y):
    HC = allHC[0,hc]
    s1,s2 = HC.shape
    if s1>1:
        G = nx.Graph()
        #G.add_nodes_from([i for i in range(1,48)])
        G.add_nodes_from(["SC","AUD","VIS","SM","CC","DM","CB"])
        sizeX,sizeY = edges.shape
        for m in range(0, sizeX):
           G.add_edge(edges[m][0], edges[m][1], weight=HC[m])

        alledges = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] != 0]
        pos = nx.circular_layout(G)  # positions for all nodes

        fig= plt.figure()
        # nodes
        nx.draw_networkx_nodes(G, pos, node_size=200)

        # edges
        nx.draw_networkx_edges(G, pos, edgelist=alledges,
                       width=2,edge_cmap = plt.cm.jet,edge_vmin =-1,edge_vmax = 1)
        #sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = vmin, vmax=vmax))
        #sm._A = []
        #plt.colorbar(sm)
        #nx.draw_networkx_edges(G, pos, edgelist=esmall,width=6, alpha=0.5, edge_color='b', style='dashed')

        # labels
        nx.draw_networkx_labels(G, pos, font_size=12, font_family='sans-serif')
        trans[hc] = nx.transitivity(G)
        plt.axis('off')
        #name = "/home/users/mrahaman1/Statelet_V2/scripts/shapeletAlgorithmScripts/SZgraphs/graph_" + str(hc) + ".gpickle"
        #nx.write_gpickle(G, "SZ.gpickle")
        #nx.write_gpickle(G, name)
        plt.savefig('/home/users/mrahaman1/Statelet_V2/scripts/shapeletAlgorithmScripts/HCdomgraphs/graph_'+str(hc)+'.png')
        plt.show()
scipy.io.savemat('Domtrans_sz'+'.mat',mdict={'tsz':trans})
