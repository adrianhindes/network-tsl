# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 13:40:24 2020

@author: adria
"""

import random
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from tqdm import tqdm
from networkTSLfunc import *
from motifCount import *


nSoc = 20
nRes = 5
m = 50
def generateNet(nSoc,nRes, m, rRange = (50, 100), popRange = (30, 50), fcRange = (40, 60)):
    #Generate social-ecological network for use in TSL model
    n = nSoc + nRes
    G = nx.gnm_random_graph(n, m)

    # Allocate social and ecological nodes
    pop = set([j for j in range(n)])
    groups = set(random.sample(pop, nSoc))
    resources = pop.difference(groups)

    for k in range(n):
        if k in groups:
            G.nodes[k]['type'] = 'social'
        else:
            G.nodes[k]['type'] = 'ecological'
    # Need network to have no isolated nodes
    # If there is an isolated node, find it and link it to another node
    for k in range(n):
        if len([j for j in G.neighbors(k)]) == 0:
            neb = random.choice(tuple(pop.difference([k])))
            G.add_edge(k,neb)
        else:
            pass
  
    #Make sure every group is connected to at least one resource
    for g in groups:
        nebs = [j for j in G.neighbors(g)]
        types = [G.nodes[k]['type'] for k in nebs]
    
        if not('ecological' in types):
            r = random.choice(tuple(resources))
            G.add_edge(g,r)
        else:
            pass
    
    #Populate resources with random levels of stock
    Rmin, Rmax = rRange
    for j in resources:
        G.nodes[j]['stock'] = random.sample(range(Rmin,Rmax),1)[0]
    
    # Populate social nodes and their level of cooperation
    # Each social node has some population of extractors
    # Preassign proportion of extractors within community that are cooperating
    popMin, popMax = popRange

    fcMin, fcMax = fcRange
    
    for k in groups:
        G.nodes[k]['pop'] = random.sample(range(popMin,popMax),1)[0]
        G.nodes[k]['fc'] = random.sample(range(fcMin,fcMax),1)[0]/100.
        
    return G

G = generateNet(nSoc,nRes,m)

fcMin = 40
fcMax = 60

iters = 10

n = 0
groups = [k for k in G.nodes() if G.nodes[k]['type'] == 'social']

for n in tqdm(range(iters)):
    # Reassign starting values
    for k in groups:
        G.nodes[k]['fc'] =  random.sample(range(fcMin,fcMax),1)[0]/100.
    t, res = TSL(G)
    plt.plot(t,res)

