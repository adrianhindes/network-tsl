# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:55:35 2019

@author: hindesa
"""
import numpy as np
import networkx as nx



w = 15
lam = 0.5 # Social ostracism coupling

# Resource stock parameters
c, d, q = 50, 50, 1
# Resource link strength
delta = 0.3
# Max resource capacity
Rmax = 250

#Population
pop = 50


# ----------------------
# Functions for TSL Model
# ----------------------

def colourNetwork(G, cutoff):
    colours = []
    for k in range(nx.number_of_nodes(G)):
        if G.nodes[k]['type'] == 'social':
            if G.nodes[k]['fc'] > cutoff:
                colours.append('b')
            else:
                colours.append('r')
        else: 
            colours.append('g')
    return colours




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

    
    

# -------------------------
# -------------------------

def TSL(G, cheatMod):
    # G = input network 
    # (must be bimodal with 'type' = 'social' and 'ecological')
    # Runs TSL model on given network
    # Returns (time, resource data, coop fraction data)

    groups = [n for n in G.nodes() if G.nodes[n]['type'] == 'social']
    resources = [n for n in G.nodes() if G.nodes[n]['type'] == 'ecological']
    
    groups = set(groups)
    resources = set(resources)
    
    # Social parameters
    mu = cheatMod
    ec = 0.483/50. #level of effort (cooperators) #level of effort (defectors)
    ed = mu*ec

    # Populate social nodes and their level of cooperation
    # Each social node has some population of extractors
    # Preassign proportion of extractors within community that are cooperating

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
    
    # Define functions to count cooperators in a neighbourhood of a given node
    # Function dependent on network, need to know which nodes are social and
    # which are ecological
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
    
    # Calculate utility of a node in given network
    def ext(f,pop):
        E = pop*(f*ec+(1-f)*ed)
        return E
    
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
        


        
    return (time, rHistory, fHistory)
