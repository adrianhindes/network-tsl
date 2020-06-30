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
mu = 1.2
ec = 0.483/50. #level of effort (cooperators) #level of effort (defectors)
ed = mu*ec
w = 15
lam = 0.5 # Social ostracism coupling

# Resource stock parameters
c, d, q = 50, 50, 1
# Resource link strength
delta = 0.2
# Max resource capacity
Rmax = 250

#Assume n,m same for both social networks
#total no. communities
nSocial = 20
#total no. resource nodes
nResource = 10
n = nSocial + nResource


#number of links
#May change depending on how network is generated, links are added
#if an isolated node is made
m = 40

# Network 

G = nx.gnm_random_graph(n, m)

# Need network to have no isolated nodes
# If there is an isolated node, find it and link it the next node
for k in range(n):
    if len([j for j in G.neighbors(k)]) == 0:
        G.add_edge(k,(k+1))
    else:
        pass


# Allocate social and ecological nodes
pop = set([j for j in range(n)])
groups = set(random.sample(pop, nSocial))
resources = pop.difference(groups)

for k in range(n):
    if k in groups:
        G.nodes[k]['type'] = 'social'
    else:
        G.nodes[k]['type'] = 'ecological'


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

fcMin = 0
fcMax = 70
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

    
    
def countCoops(G,k):
    
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
        
layout = nx.spring_layout(G)
nx.draw_networkx(G, pos=layout, node_color = colours, with_labels=True, node_size=150)

plt.title('Initial Network')
plt.show()
    
t = 0
tEnd = 200 #end point
dt = 0.1 #time step
# Lists to store values to plot later
time = [0]

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
        differences = [delta*(G.nodes[j]['stock'] - R)  for j in nebResources]
        inOutFlow = sum(differences)
        dR = (dRself + inOutFlow)*dt
        
        R += dR
        G.nodes[k]['stock'] += dR
        rHistory[k].append(R)


    #update quantities
    
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

nx.draw_networkx(G, pos=layout, node_color = colours, with_labels=True, node_size=150)

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

