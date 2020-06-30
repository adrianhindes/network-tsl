# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:55:35 2019

@author: hindesa
"""
import random
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Random resource network
# Using ecological subsystem only version of TSL

#Assume n,m same for both social networks
#total no. agents
n = 20
#number of links
#May change depending on how network is generated, links are added
#if an isolated node is made
m = 30 

mu = 1.2
ec = 0.483/50. #level of effort (cooperators) #level of effort (defectors)
ed = mu*ec
Rmax = 250
c = random.sample(range(20, 60), n)
d, q = 50, 1

# Network 

G = nx.gnm_random_graph(n, m)


# Need network to have no isolated nodes
# If there is an isolated node, find it and link it the next node
for k in range(n):
    if len([j for j in G.neighbors(k)]) == 0:
        G.add_edge(k,(k+1))
    else:
        pass
m = G.number_of_edges() #Reassign in case rewire


#Populate resources with random levels of stock
stocks = random.sample(range(int(np.floor(Rmax/4)),Rmax), n)

for j in range(n):
    G.node[j]['node_size'] = stocks[j]
    
# Randomize coupling strengths between 0.1 and 0.5
deltas = random.sample(set(np.linspace(0.1, 0.5, num=50)), m)

k = 0
for u,v,l in G.edges(data=True):
    l['weight'] = deltas[k]
    k += 1


nx.draw_networkx(G, pos=None, node_size = stocks)

plt.title('Initial Network')
plt.show()
    

# Extraction function
def ext(f):
    E = n*(f*ec+(1-f)*ed)
    return E


# initial condition
# static fraction of cooperators
fc = 0.7
R = 200

t = 0
tEnd = 30 #end point
dt = 0.1 #time step
# Lists to store values to plot later
time = []
rLists = []
k = 0
while k < n:
    rLists.append([])
    k += 1



while t<tEnd:
    R = np.zeros(n)
    for k in range(n):
        R[k] = G.node[k]['node_size']

        dRself = c[k] - d*(R[k]/Rmax)**2 - q*ext(fc)*R[k]
        differences = [G.edges[j,k]['weight']*(G.node[j]['node_size'] - R[k])  for j in G.neighbors(k)]
        inOutFlow = sum(differences)
        dR = (dRself + inOutFlow)*dt
        
        R[k] += dR
        G.node[k]['node_size'] += dR

    #update quantities
    
    t += dt
    
    #append lists
    time.append(t)
    for k in range(n):
        rLists[k].append(R[k])


# Plot
plt.xlabel('Time')
plt.ylabel('Resource Stock')
for i in range(n):
    plt.plot(time, rLists[i])
plt.show()

stocks = nx.get_node_attributes(G, 'node_size')
sizes = [stocks[k] for k in stocks]
nx.draw_networkx(G, pos=None, node_size = sizes)

plt.title('Final Network')
plt.show()


