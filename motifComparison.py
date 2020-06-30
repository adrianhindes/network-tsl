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
motif1 = 'direct1'
motif2 = 'asymmShare2'

# Population of social nodes
# and initial cooperation fractions
pop = 50
initFrac1 = 0.7
initFrac2 = 0.3
frac1 = {0: initFrac1, 1: initFrac2}
frac2 = {0: initFrac1, 1: initFrac2}

G1 = nx.Graph()
G2 = nx.Graph()
G1.add_nodes_from(nodes)
G2.add_nodes_from(nodes)
G1.add_edges_from(edgeConfigs[motif1])
G2.add_edges_from(edgeConfigs[motif2])

# Square layout for motifs
pos = {0: [-0.4, 0.4], 1:[0.4,0.4], 2:[0.4,-0.4], 3:[-0.4,-0.4]}


for j in resources:
    G1.nodes[j]['type'] = 'ecological'
    G1.nodes[j]['stock'] = stockStart
    
    G2.nodes[j]['type'] = 'ecological'
    G2.nodes[j]['stock'] = stockStart
    
# Populate social nodes and their level of cooperation
# Each social node has some population of extractors
# Preassign proportion of extractors within community that are cooperating
for j in groups:
    G1.nodes[j]['type'] = 'social'
    G1.nodes[j]['fc'] = frac1[j]
    G1.nodes[j]['pop'] = pop
    G2.nodes[j]['type'] = 'social'
    G2.nodes[j]['pop'] = pop
    G2.nodes[j]['fc'] = frac2[j]

G1.nodes[0]['fc'] = 0.8
G1.nodes[1]['fc'] = 0.35
#Plot initial motifs

colours = colourNetwork(G1, 0.5)

plt.figure(1)
nx.draw_networkx(G1, pos=pos, node_color = colours, with_labels=True)
plt.title('Initial Motif ' + motif1)
plt.show()

colours = colourNetwork(G2, 0.5)

plt.figure(2)
nx.draw_networkx(G2, pos=pos, node_color = colours, with_labels=True)
plt.title('Initial Motif ' + motif2)
plt.show()


results1 = TSL(G1, mu)
time = results1[0]
m1fc1 = (results1[2])[0]
m1fc2 = (results1[2])[1]

results2 = TSL(G2, mu)
m2fc1 = (results2[2])[0]
m2fc2 = (results2[2])[1]

plt.xlabel('Time')
plt.ylabel('Fraction of Cooperators')
plt.plot(time,m1fc1, label = 'Node 0 ' +motif1)
plt.plot(time,m2fc1, label = 'Node 0 ' +motif2)
plt.legend(loc='best')
plt.show()

plt.xlabel('Time')
plt.ylabel('Fraction of Cooperators')
plt.plot(time,m1fc2, label = 'Node 1 ' + motif1)
plt.plot(time,m2fc2, label = 'Node 1 ' + motif2)
plt.legend(loc='best')
plt.show()

