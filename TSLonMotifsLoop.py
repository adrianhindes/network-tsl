# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:55:35 2019

@author: hindesa
"""
import random
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from tqdm import tqdm
from networkTSLfunc import *
import pandas as pd

#Initial network conditions
fc0 = 0.6
res = 50
pop = 50

#Motif list
motifs = {
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
#Construct graphs
motifGs = {}
social = [0,1]
resources = [2,3]

for k in motifs.keys():
    Gr = nx.Graph()
    Gr.add_nodes_from([0,1,2,3])
    for j in social: 
        Gr.nodes[j]['type'] = 'social'
        Gr.nodes[j]['fc'] = fc0
        Gr.nodes[j]['pop'] = pop
    for j in resources: 
        Gr.nodes[j]['type'] = 'ecological'
        Gr.nodes[j]['stock'] = res

    Gr.add_edges_from(motifs[k])
    motifGs[k] = Gr

res = {key:[] for key in motifs.keys()}

for G in tqdm(motifGs.keys()):
    for fc in np.linspace(0.1,0.9,9):
        Gr = motifGs[G]
        for n in social:
            Gr.nodes[n]['fc'] = fc
        
        result = TSL(Gr, lam=0.5, mu=1.5)
        _, fcList = result
        res[G].append(fcList[-1])
    
resDf = pd.DataFrame.from_dict(res)
resDf = resDf.transpose()
#resDf.to_csv('TSLonMotifs.csv')
    
    

