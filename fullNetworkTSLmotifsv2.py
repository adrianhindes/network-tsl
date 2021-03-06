# -*- coding: utf-8 -*-#
import random
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import itertools
import networkx.algorithms.isomorphism as iso
from tqdm import tqdm

# Social parameters
mu = 1.2
ec = 0.483/50. #level of effort (cooperators) #level of effort (defectors)
ed = mu*ec
w = 15
lam = 0.6 # Social ostracism coupling

# Resource stock parameters
c, d, q = 50, 50, 1
# Resource link strength
delta = 0.2
# Max resource capacity
Rmax = 100

#Assume n,m same for both social networks
#total no. communities
nSocial = 20
#total no. resource nodes
nResource = 10
n = nSocial + nResource


#number of links
#May change depending on how network is generated, links are added
#if an isolated node is made
m = 100

# Network 

G = nx.gnm_random_graph(n, m)



# Allocate social and ecological nodes
pop = set([j for j in range(n)])
groups = set(random.sample(pop, nSocial))
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

m = G.number_of_edges() #Reassign in case rewire

#Populate resources with random levels of stock
Rmin = 50
for j in resources:
    G.nodes[j]['stock'] = random.sample(range(Rmin,Rmax),1)[0]
    
# Populate social nodes and their level of cooperation
# Each social node has some population of extractors
# Preassign proportion of extractors within community that are cooperating
popMin = 30
popMax = 50

fcMin = 40
fcMax = 60
for k in groups:
    G.nodes[k]['pop'] = random.sample(range(popMin,popMax),1)[0]
    G.nodes[k]['fc'] = random.sample(range(fcMin,fcMax),1)[0]/100.




# ----------------------
# Functions for TSL Model
# ----------------------



# Extraction function
def ext(f,pop):
    E = pop*(f*ec+(1-f)*ed)
    return E

# Ostracism function (Gompertz growth)
def gompertz(f):
    #parameters (from paper)
    h = 0.34
    t = -150
    g = -10
    gomp = h * np.exp(t * np.exp(g*f))
    return gomp


#Cobb-Douglas production function
def cobbdoug(E,R):
    gamma = 10
    a = 0.6
    b = 0.2
    product = gamma*(E**a)*(R**b)
    return product

    
#
def countCoops(G,k):
    #Calculate fraction of cooperators in the neighbourhood of a node
    #Get neighbouring nodes
    neighbours = [j for j in G.neighbors(k)]
    communities = groups.intersection(neighbours)
    
    if G.nodes[k]['type'] == 'social':
        nebCoops = [G.nodes[j]['fc'] for j in communities]
        fc = (sum(nebCoops) + G.nodes[k]['fc'])/(len(nebCoops) + 1)
        totalPop = sum([G.nodes[j]['pop'] for j in communities]) + G.nodes[k]['pop']
    else:
        if len(communities) == 0:
            fc = 0
            totalPop = 0
        else:
            nebCoops = [G.nodes[j]['fc'] for j in communities]
            fc = sum(nebCoops)/len(communities)
            totalPop = sum([G.nodes[j]['pop'] for j in communities])

    return (fc, totalPop)

population = 0
for k in groups:
    pop = G.nodes[k]['pop']
    population += pop
    
def countFrac(G):
    # List of cooperator populations
    coopPops = []
    for k in groups:
        coopPops.append(G.nodes[k]['fc']*G.nodes[k]['pop'])
        
    totalCoops = sum(coopPops)
    totalFc = totalCoops/population
    return totalFc
        

def utilNode(G, k):
    # Calculate payoffs for cooperators and defectors
    # within a given node in network
    fc = G.nodes[k]['fc']
    pop = G.nodes[k]['pop']
    fcNebs = countCoops(G,k)[0]
    
    # Node effort
    E = ext(fc,pop)
    
    # How much is this community extracting in total?
    extracting = resources.intersection([j for j in G.neighbors(k)])
    R = sum([G.nodes[j]['stock'] for j in extracting])
    
    #Payoffs
    pic = ec*((cobbdoug(E,R)/E)-w)
    pid = ed*((cobbdoug(E,R)/E)-w)
    H = (pid-pic)/pid
    
    #Utilities
    uc = pic
    ud = pid - H*(gompertz(fc)+lam*gompertz(fcNebs))
        
    return (uc,ud)

# -------------------------
# -------------------------
    

colours = []
    
for k in range(n):
    if G.nodes[k]['type'] == 'social':
        if G.nodes[k]['fc'] > 0.5:
            colours.append('b')
        else:
            colours.append('r')
    else: 
        colours.append('g')

    
t = 0
tEnd = 400 #end point
dt = 0.1 #time step
# Lists to store values to plot later
time = [0]



def TSLnetwork(G):
    resources = []
    groups = []
    
    for k in G.nodes():
        if G.nodes[k]['type'] == 
    
    tempR = []
    for k in resources:
        node = k
        stock = G.nodes[k]['stock']
        tempR.append((node, stock))
    
    rHistory = {node:[s] for (node, s) in tempR}

    tempF = []
    for k in groups:
        node = k
        frac = G.nodes[k]['fc']
        tempF.append((node, frac))

    fHistory = {node:[s] for (node, s) in tempF}

    fTotal = [countFrac(G)]
    while t<tEnd:
    
        # Update communities
        # List to store present cooperator fractions
        F = np.zeros(len(groups))
    
        for k in groups:
            F = G.nodes[k]['fc']
            
            Uc, Ud = utilNode(G,k)
        
            dfc = F*(1 - F)*(Uc - Ud)*dt
        
            F += dfc
            G.nodes[k]['fc'] += dfc
            fHistory[k].append(F)
    
        # Update resources
        # List to store present resource stocks
        for k in resources:
            R = G.nodes[k]['stock']
        
            # Get fraction of cooperators and no. extractors in neighbourhood
            fc, nebs = countCoops(G, k)

            dRself = c - d*(R/Rmax)**2 - q*ext(fc, nebs)*R
            nebResources = resources.intersection([j for j in G.neighbors(k)])
        
            differences = []
            for j in nebResources:
                diff = delta*(G.nodes[j]['stock']-R)
                G.nodes[j]['stock'] += -diff
                differences.append(diff)
        
            inOutFlow = sum(differences)
            dR = (dRself + inOutFlow)*dt
        
            R += dR
            G.nodes[k]['stock'] += dR
            rHistory[k].append(R)


        #update quantities
        fTotal.append(countFrac(G))
        t += dt
        time.append(t)
    

        

colours = []
    
for k in range(n):
    if G.nodes[k]['type'] == 'social':
        if G.nodes[k]['fc'] > 0.5:
            colours.append('b')
        else:
            colours.append('r')
    else: 
        colours.append('g')


        

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
    
    for q in tqdm(fquads):
        subg = G.subgraph(q)
        matches = [net for net in motifs if nx.is_isomorphic(subg,motifs[net],node_match=nm)]
        if not len(matches) == 1:
            raise Exception("No isomorphic motifs found for" + str(q))
        
        mcount[matches[0]] += 1
    
    return mcount

#count = mcounter(G,motifs)
#print(count)