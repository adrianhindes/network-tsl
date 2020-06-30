# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:55:35 2019

@author: hindesa
"""
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from networkTSLfunc import colourNetwork, TSL


stockStart = 100
# Construct Motifs

mu = 3

nodes = [0,1,2,3]
groups = {0,1}
resources = {2,3}
#0,1 social
#2,3 ecological
edgeConfigs = { 'direct1': [(0,1), (0,3), (1,2)], 
                'direct2': [(0,3), (1,2), (2,3)], 
                'direct3': [(0,1), (1,2), (2,3), (0,3)],
                'shared1': [(0,3), (1,3), (2,3)],
                'shared2': [(0,1), (0,3), (1,3)],
                'shared3': [(0,3), (1,3), (2,3), (0,1)],
                'mShared1': [(0,3), (0,2), (1,2), (1,3)],
                'mShared2': [(0,3), (0,2), (1,2), (1,3), (0,1)],
                'mShared3': [(0,3), (0,2), (1,2), (1,3), (2,3)],
                'mShared4': [(0,3), (0,2), (1,2), (1,3), (2,3), (0,1)],
                'asymmShare1': [(0,2), (0,3), (1,2)],
                'asymmShare2': [(0,2), (0,3), (1,2), (0,1)],
                'asymmShare3': [(0,2), (0,3), (1,2), (2,3)],
                'asymmShare4': [(0,2), (0,3), (1,2), (2,3), (0,1)],
                }

#Choose which motifs to compare
motif = 'direct1'

# Population of social nodes
# and initial cooperation fractions
pop = 50
initFrac1 = 0.7
initFrac2 = 0.3
frac = {0: initFrac1, 1: initFrac2}

G = nx.Graph()

G.add_nodes_from(nodes)

G.add_edges_from(edgeConfigs[motif])

# Square layout for motifs
pos = {0: [-0.4, 0.4], 1:[0.4,0.4], 2:[0.4,-0.4], 3:[-0.4,-0.4]}


for j in resources:
    G.node[j]['type'] = 'ecological'
    G.node[j]['stock'] = stockStart

    
# Populate social nodes and their level of cooperation
# Each social node has some population of extractors
# Preassign proportion of extractors within community that are cooperating
for j in groups:
    G.node[j]['type'] = 'social'
    G.node[j]['pop'] = pop


colours = colourNetwork(G, 0.5)

plt.figure(1)
nx.draw_networkx(G, pos=pos, node_color = colours, with_labels=True)
plt.title('Initial Motif ' + motif1)
plt.show()


k = 100 # number of parameter choices
possibleF = list(map(lambda x : x/100, list(range(100))))
initF1 = random.choices(possibleF, k=k)
initF2 = random.choices(possibleF, k=k)

finalFval1 = []
finalFval2 = []
for j in range(k):
    G.node[0]['fc'] = initF1[j]
    G.node[1]['fc'] = initF2[j]
    results1 = TSL(G, mu)
    finalFc1 = ((results1[2])[0])[-1]
    finalFc2 = ((results1[2])[1])[-1]
    
    finalFval1.append(finalFc1)
    finalFval2.append(finalFc2)
    


