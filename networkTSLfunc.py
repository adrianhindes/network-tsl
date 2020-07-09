# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:55:35 2019

@author: hindesa
"""
import numpy as np
import networkx as nx

w = 15
lam = 0.05 # Social ostracism coupling

# Resource stock parameters
c, d, q = 50, 50, 1
# Resource link strength
delta = 0.3
# Max resource capacity
Rmax = 250

# Social parameters
mu = 1.2
ec = 0.483/50. #level of effort (cooperators) #level of effort (defectors)
#ed = mu*ec

def ext(f,pop,mu):
    ed = mu*ec
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
    #Calculate fraction of cooperators in the neighbourhood of a node
    #Get neighbouring nodes
    groups = [n for n in G.nodes() if G.nodes[n]['type'] == 'social']

    neighbours = [j for j in G.neighbors(k)]
    communities = set(groups).intersection(neighbours)
    
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


def countFrac(G):
    # List of cooperator populations
    groups = [n for n in G.nodes() if G.nodes[n]['type'] == 'social']
    population = 0
    for k in groups:
        pop = G.nodes[k]['pop']
        population += pop
    

    coopPops = []
    for k in groups:
        coopPops.append(G.nodes[k]['fc']*G.nodes[k]['pop'])
        
    totalCoops = sum(coopPops)
    totalFc = totalCoops/population
    return totalFc

def payoffNode(G, k, mu):
    ed = mu*ec
    #Calculate the payoffs for cooperators and defectors in chosen node
    resources = [n for n in G.nodes() if G.nodes[n]['type'] == 'ecological']

    fc = G.nodes[k]['fc']
    pop = G.nodes[k]['pop']
    #Extraction effort
    E = ext(fc, pop, mu)
    #Local resource
    extracting = set(resources).intersection([j for j in G.neighbors(k)])
    R = sum([G.nodes[j]['stock'] for j in extracting])
    
    pic = ec*((cobbdoug(E,R)/E)-w)
    pid = ed*((cobbdoug(E,R)/E)-w)
    
    return (pic, pid)
    

def utilNode(G, k, lam, mu):
    # Calculate payoffs for cooperators and defectors
    # within a given node in network
    groups = [n for n in G.nodes() if G.nodes[n]['type'] == 'social']

    fc = G.nodes[k]['fc']
    neighbours = [j for j in G.neighbors(k)]
    nebs = set(groups).intersection(neighbours)
    nebs = list(nebs)
        
    pic, pid = payoffNode(G, k, mu)

    H = (pid-pic)/pid
    #Need to calculate H for external communities
    Hs = []
    for neb in nebs:
        nebc, nebd = payoffNode(G, neb, mu)
        nebH = (nebd-nebc)/nebd
        Hs.append(nebH)
        
    nebFcs = [G.nodes[j]['fc'] for j in nebs]
    gomps = list(map(gompertz, nebFcs))
    prods = [j*k for j in Hs for k in gomps]
    if len(prods) == 0: 
        avg = 0
    else: avg = sum(prods)/len(prods)
    #Utilities
    uc = pic
    ud = pid - H*(1-lam/2)*gompertz(fc)-(lam/2)*avg
        
    return (uc,ud)
    
    

# -------------------------
# -------------------------

def TSL(G,lam=0.5,c=50,mu=1.2):
    # G = input network 
    # (must be bimodal with 'type' = 'social' and 'ecological')
    # Runs TSL model on given network
    # Returns (time, resource data, coop fraction data)

    groups = [n for n in G.nodes() if G.nodes[n]['type'] == 'social']
    resources = [n for n in G.nodes() if G.nodes[n]['type'] == 'ecological']
    
    groups = set(groups)
    resources = set(resources)
    

    t = 0
    tEnd = 250 #end point
    dt = 0.1 #time step
    # Lists to store values to plot later
    time = [0]

    #rHistory = {node:[s] for (node, s) in tempR}


    #fHistory = {node:[s] for (node, s) in tempF}
    
    fTotal = [countFrac(G)]

    while t<tEnd:
    
        # Update communities
        # List to store present cooperator fractions
        F = np.zeros(len(groups))
    
        for k in groups:
            F = G.nodes[k]['fc']
        
            Uc, Ud = utilNode(G,k,lam=lam,mu=mu)
        
            dfc = F*(1 - F)*(Uc - Ud)*dt
        
            F += dfc
            G.nodes[k]['fc'] += dfc
            #fHistory[k].append(F)
    
        # Update resources
        # List to store present resource stocks
        for k in resources:
            R = G.nodes[k]['stock']
        
            # Get fraction of cooperators and no. extractors in neighbourhood
            fc, nebs = countCoops(G, k)

            dRself = c - d*(R/Rmax)**2 - q*ext(fc, nebs, mu)*R
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
            #rHistory[k].append(R)

        #update quantities
        fTotal.append(countFrac(G))
        t += dt
        time.append(t)
        
        
    return (time, fTotal)

