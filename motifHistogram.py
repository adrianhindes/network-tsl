# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:14:15 2020

@author: Adrian Hindes
Null model for motif correlation
"""
import random
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import itertools
import networkx.algorithms.isomorphism as iso
from tqdm import tqdm

#Assume n,m same for both social networks
#total no. communities
nSocial = 5
#total no. resource nodes
nResource = 5
n = nSocial + nResource

m = 15

reps = 20

motifsLinks = {
        'IA': [(0,3), (1,2)],
        'IB': [(0,1), (0,3), (1,2)], 
        'IC': [(0,3), (1,2), (2,3)], 
        'ID': [(0,1), (1,2), (2,3), (0,3)],
        'IIA': [(0,3), (1,3)],
        'IIB': [(0,1), (0,3), (1,3)],
        'IIC': [(0,3), (1,3), (2,3)],
        'IID': [(0,3), (1,3), (2,3), (0,1)],
        'IIIA': [(0,3), (0,2), (1,2), (1,3)],
        'IIIB': [(0,3), (0,2), (1,2), (1,3), (0,1)],
        'IIIC': [(0,3), (0,2), (1,2), (1,3), (2,3)],
        'IIID': [(0,3), (0,2), (1,2), (1,3), (2,3), (0,1)],
        'IVA': [],
        'IVB': [(0,1)],
        'IVC': [(2,3)],
        'IVD': [(0,1), (2,3)],
        'VA': [(0,2), (0,3), (1,2)],
        'VB': [(0,2), (0,3), (1,2), (0,1)],
        'VC': [(0,2), (0,3), (1,2), (2,3)],
        'VD': [(0,2), (0,3), (1,2), (2,3), (0,1)],
        'VIA': [(0,1), (0,3)],
        'VIB': [(0,1), (0,3), (2,3)],
        'VIC': [(0,2), (0,3), (0,1)],
        'VID': [(0,1), (0,3), (0,2), (2,3)],
        'VIIA': [(0,3)],
        'VIIB': [(0,3), (2,3)],
        'VIIC': [(0,3), (0,2)],
        'VIID': [(0,2), (0,3), (2,3)]
        }

motifs = {}

for k in motifsLinks.keys():
    Gr = nx.Graph()
    Gr.add_nodes_from([0,1,2,3])
    Gr.nodes[0]['type'] = 'social'
    Gr.nodes[1]['type'] = 'social'
    Gr.nodes[2]['type'] = 'ecological'
    Gr.nodes[3]['type'] = 'ecological'
    
    Gr.add_edges_from(motifsLinks[k])
    motifs[k] = Gr
    
nm = iso.categorical_node_match('type', ['social','resource'])

def mcounter(G, motifs):
    
    #initialize count
    mcount = dict(zip(motifsLinks.keys(), list(map(int, np.zeros(len(motifs))))))
    
    nodes = G.nodes()
    actors = [k for k in nodes if G.nodes[k]['type'] == 'social']
    resources = [k for k in nodes if G.nodes[k]['type'] == 'ecological']
    
    quads = list(itertools.product(*[actors, actors, resources, resources]))
    # Filter duplicates
    quads = [list(mot) for mot in quads if len(set(mot)) == 4]
    #Filter repeated entries
    quads = map(list, map(np.sort,quads))
    fquads = []
    [fquads.append(quad) for quad in quads if not fquads.count(quad)]
    
    for q in fquads:
        subg = G.subgraph(q)
        matches = [net for net in motifs if nx.is_isomorphic(subg,motifs[net],node_match=nm)]
        if not len(matches) == 1:
            raise Exception("No isomorphic motifs found for" + str(q))
        
        mcount[matches[0]] += 1
    
    return mcount


data = {k:[] for k in motifsLinks.keys()}

for s in tqdm(range(reps)):
    G = nx.gnm_random_graph(n, m)
    
    for k in range(n):
        if len([j for j in G.neighbors(k)]) == 0:
            neb = random.sample(range(n),1)[0]
            G.add_edge(k,neb)
        else: pass

    pop = set([j for j in range(n)])
    groups = set(random.sample(pop, nSocial))
    resources = pop.difference(groups)

    for k in range(n):
        if k in groups:
            G.nodes[k]['type'] = 'social'
        else:
            G.nodes[k]['type'] = 'ecological'
        
    for g in groups:
        nebs = [j for j in G.neighbors(g)]
        types = [G.nodes[k]['type'] for k in nebs]
    
        if not('ecological' in types):
            r = random.choice(tuple(resources))
            G.add_edge(g,r)
        else: pass
    
    count = mcounter(G,motifs)
    for mot in motifsLinks.keys():
        data[mot].append(count[mot])
    
