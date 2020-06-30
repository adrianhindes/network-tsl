# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:55:35 2019

@author: hindesa
"""
import random
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx


# Social parameters
mu = 1.5
ec = 0.483/50. #level of effort (cooperators) #level of effort (defectors)
ed = mu*ec
w = 15
lam = 0.5 # Social ostracism coupling

# Resource stock parameters
c, d, q = 50, 50, 1
# Resource link strength
delta = 0.3
# Max resource capacity
Rmax = 250

# Construct Motifs


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
                'mediated1': [(0,1), (0,3), (2,3)],
                'mediated2': [(0,1), (0,3), (0,2), (2,3)]
                }

#Choose which motif to run
motif = 'asymmShare4'
G = nx.Graph()
G.add_nodes_from(nodes)
G.add_edges_from(edgeConfigs[motif])

# Square layout for motifs
pos = {0: [-0.4, 0.4], 1:[0.4,0.4], 2:[0.4,-0.4], 3:[-0.4,-0.4]}


#Populate resources with random levels of stock
Rmin = 50

for j in resources:
    G.node[j]['type'] = 'ecological'
    G.node[j]['stock'] = random.sample(range(Rmin,Rmax),1)[0]
    
# Populate social nodes and their level of cooperation
# Each social node has some population of extractors
# Preassign proportion of extractors within community that are cooperating
pop1 = 50
pop2 = 50
fc1 = 0.1
fc2 = 0.8

G.node[0]['type'] = 'social'
G.node[1]['type'] = 'social'
G.node[0]['pop'] = pop1
G.node[1]['pop'] = pop2
G.node[0]['fc'] = fc1
G.node[1]['fc'] = fc2




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

    
    
def countCoops(G,k):
    
    neighbours = [j for j in G.neighbors(k)]
    communities = groups.intersection(neighbours)
    
    if G.node[k]['type'] == 'social':
        nebCoops = [G.node[j]['fc'] for j in communities]
        fc = (sum(nebCoops) + G.node[k]['fc'])/(len(nebCoops) + 1)
        totalPop = sum([G.node[j]['pop'] for j in communities]) + G.node[k]['pop']
    else:
        if len(communities) == 0:
            fc = 0
            totalPop = 0
        else:
            nebCoops = [G.node[j]['fc'] for j in communities]
            fc = sum(nebCoops)/len(communities)
            totalPop = sum([G.node[j]['pop'] for j in communities])

    return (fc, totalPop)


def utilNode(G, k):
    # Calculate payoffs for cooperators and defectors
    # within a given node in network
    fc = G.node[k]['fc']
    pop = G.node[k]['pop']
    fcNebs = countCoops(G,k)[0]
    
    # Node effort
    E = ext(fc,pop)
    
    # How much is this community extracting in total?
    extracting = resources.intersection([j for j in G.neighbors(k)])
    R = sum([G.node[j]['stock'] for j in extracting])
    
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
    
for k in range(4):
    if G.node[k]['type'] == 'social':
        if G.node[k]['fc'] > 0.5:
            colours.append('b')
        else:
            colours.append('r')
    else: 
        colours.append('g')



nx.draw_networkx(G, pos=pos, node_color = colours, with_labels=True)

plt.title('Initial Network')
plt.show()
    
t = 0
tEnd = 200 #end point
dt = 0.1 #time step

# Lists to store values to plot later
time = [0]

tempR = []
for k in resources:
    stock = G.node[k]['stock']
    tempR.append((k, stock))
    
rHistory = {node:[s] for (node, s) in tempR}

tempF = []
for k in groups:
    node = k
    frac = G.node[k]['fc']
    tempF.append((node, frac))

fHistory = {node:[s] for (node, s) in tempF}



while t<tEnd:
    
    # Update communities
    # List to store present cooperator fractions
    F = np.zeros(len(groups))
    
    for k in groups:
        F = G.node[k]['fc']
        
        Uc, Ud = utilNode(G,k)
        
        dfc = F*(1 - F)*(Uc - Ud)*dt
        
        F += dfc
        G.node[k]['fc'] += dfc
        fHistory[k].append(F)
    
    # Update resources
    # List to store present resource stocks
    for k in resources:
        R = G.node[k]['stock']
        
        # Get fraction of cooperators and no. extractors in neighbourhood
        fc, nebs = countCoops(G, k)

        dRself = c - d*(R/Rmax)**2 - q*ext(fc, nebs)*R
        nebResources = resources.intersection([j for j in G.neighbors(k)])
        differences = [delta*(G.node[j]['stock'] - R)  for j in nebResources]
        inOutFlow = sum(differences)
        dR = (dRself + inOutFlow)*dt
        
        R += dR
        G.node[k]['stock'] += dR
        rHistory[k].append(R)


    #update quantities
    
    t += dt
    time.append(t)

        

colours = []
    
for k in range(4):
    if G.node[k]['type'] == 'social':
        if G.node[k]['fc'] > 0.5:
            colours.append('b')
        else:
            colours.append('r')
    else: 
        colours.append('g')

nx.draw_networkx(G, pos=pos, node_color = colours, with_labels=True)

plt.title('Final Network')
plt.show()

plt.xlabel('Time')
plt.ylabel('Resource Stock')
for i in resources:
    plt.plot(time, rHistory[i], label = 'Node ' + str(i))
plt.legend(loc='best')
plt.show()

plt.xlabel('Time')
plt.ylabel('Proportion of Cooperators')
for i in groups:
    plt.plot(time, fHistory[i], label= 'Node ' + str(i))
plt.legend(loc='best')
plt.show()

