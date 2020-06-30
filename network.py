import sys
import matplotlib.pyplot as plt
from networkx import nx
import numpy as np
import random

n = 50 # 10 nodes
m = 110 # 20 edges

G = nx.gnm_random_graph(n, m)

mu = 4 # degree of cheating
ec = 0.483/n #level of effort (cooperators) #level of effort (defectors)
ed = mu*ec

#state of node
# if 0 = cooperator (blue)
# if 1 = defector (red)

states = np.zeros(n)

# Randomly assign cooperators and defectors
for k in range(n):
    G.node[k]['state'] = random.randint(0,1)
    
colours = []
    
for k in range(n):
    if G.node[k]['state'] == 0:
        colours.append('b')
    else: 
        colours.append('r')
        
def coopNeighbors(G,n):
    # Count cooperating neighbours
    # Returns proportion of neighbours who cooperate
    neighbours = [k for k in G.neighbors(n)]
    coops = 0
    for k in neighbours:
        if G.node[k]['state'] == 0:
            coops += 1
        else:
            pass
    return coops/len(neighbours)

def gompertz(f):
    #parameters (from paper)
    h = 0.34
    t = -150
    g = -10
    gomp = h * np.exp(t * np.exp(g*f))
    return gomp

def utilCoop(fc,R):
    E = ext(fc)
    uc = ec*((cobbdoug(E,R)/E)-w)
    return uc

def utilDefec(fc,R):
    E = ext(fc)
    ud = ed*((cobbdoug(E,R)/E)-w)-gompertz(fc)*(ed-ec)/ed
    return ud


    

nx.draw_networkx(G, pos=None, node_color=colours, with_labels=True, node_size=100)
plt.show()