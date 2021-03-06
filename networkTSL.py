import sys
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import random
import networkx as nx

#total no. agents
n = 50
m = 120 # 20 edges


# Maximum amount of resource
Rmax = 500

# Social parameters

mu = 1.5 # degree of cheating
ec = 0.483/50. #level of effort (cooperators) #level of effort (defectors)
ed = mu*ec
w = 15 #cost of harvesting

# Resource parameters
# c inflow, d discharge, q catchability
c, d, q = 50, 50, 1

# Network 

#state of node


#fraction of cooperators initial
# Use this to probabilistically assign 
Fprob = 0.7

states = ['cooperator', 'defector']
prob = [Fprob, 1-Fprob]

G = nx.gnm_random_graph(n, m)

# Need network to have no isolated nodes
# If there is an isolated node, find it and link it the next node
for k in range(n):
    if len([j for j in G.neighbors(k)]) == 0:
        G.add_edge(k,(k+1))
    else:
        pass

# Randomly assign cooperators and defectors
for k in range(n):
    G.node[k]['state'] = random.choices(states, prob)[0]
    
colours = []
    
for k in range(n):
    if G.node[k]['state'] == 'cooperator':
        colours.append('b')
    else: 
        colours.append('r')
        
def coopNeighbors(G,k):
    # Count cooperating neighbours
    # Returns proportion of neighbours who cooperate
    neighbours = [j for j in G.neighbors(k)]
    coops = 0
    for j in neighbours:
        if G.node[j]['state'] == 'cooperator':
            coops += 1
        else:
            pass
    return coops/len(neighbours)


# Extraction function
def ext(f):
    E = n*(f*ec+(1-f)*ed)
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

    
    
def countCoops(G):
    count = 0
    for k in G.nodes():
        if G.node[k]['state'] == 'cooperator':
            count += 1
        else:
            pass
    return count/n


def utilNode(fc, R, k):
    # Calculate utility of given node in network
    # Depends on state of node
    E = ext(fc)
    
    if G.node[k]['state'] == 'cooperator':
        ui = ec*((cobbdoug(E,R)/E)-w)
    else:
        nc = coopNeighbors(G,k)
        ui = ed*((cobbdoug(E,R)/E)-w)-gompertz(nc)*(ed-ec)/ed
        
    return ui
        


# initial condition

#amount of resource available initial
R = 300

t = 0
tEnd = 1000 #end point
dt = 0.1 #time step
# Lists to store values to plot later
time = []
rList = []
fList = []

# Draw initial network
colours = []
    
for k in range(n):
    if G.node[k]['state'] == 'cooperator':
        colours.append('b')
    else: 
        colours.append('r')
nx.draw_networkx(G, pos=None, node_color=colours, with_labels=True, node_size=150)
plt.title('Initial Network')

plt.show()



while t<tEnd:
    F = countCoops(G)
    dR = (c - d*(R/Rmax)**2 - q*ext(F)*R)*dt
    
    #Choose random node to update
    updateNode = random.choice(range(n))
    compareNode = random.choice([k for k in  G.neighbors(updateNode)])
    
    # Calculate node utilities
    utilUpdate = utilNode(F, R, updateNode)
    utilCompare = utilNode(F, R, compareNode)
    
    delta = utilUpdate - utilCompare
    
    if delta<0:
        prob = np.abs(delta)/(np.abs(utilUpdate)+np.abs(utilCompare))
        states = [G.node[updateNode]['state'], G.node[compareNode]['state']]
        probs = [1-prob, prob]
        newState = random.choices(states, probs)[0]
        G.node[updateNode]['state'] = newState
    
    #update quantities
    R += dR
    t += dt
    
    #append lists
    time.append(t)
    rList.append(R)
    fList.append(F)


# Reset colours
colours = []
    
for k in range(n):
    if G.node[k]['state'] == 'cooperator':
        colours.append('b')
    else: 
        colours.append('r')

# Plot
plt.subplot(211)
plt.plot(time, fList, color='black')

plt.xlabel('Time')
plt.ylabel('Fraction of Cooperators')

plt.subplot(212)
plt.plot(time, rList, color='green')
plt.xlabel('Time')
plt.ylabel('Resource Stock')
plt.show()

nx.draw_networkx(G, pos=None, node_color=colours, with_labels=True, node_size=150)
plt.title('Final Network')

plt.show()
