import numpy as np
import matplotlib.pyplot as plt
import random
import networkx as nx

# Motif with separated social networks extracting from their own resource piles
# But resources are coupled
# Resource (2) is not being extracted

#Assume n,m same for both social networks
#total no. agents
n = 50
#number of links
#May change depending on how network is generated, links are added
#if an isolated node is made
m = 120 


# Maximum amount of resource
R1max = 200
R2max = 200

#amount of resource available initial
R1 = 150
R2 = 150

# Social parameters

mu = 1.5 # degree of cheating
ec = 0.483/50. #level of effort (cooperators) #level of effort (defectors)
ed = mu*ec
w = 15 #cost of harvesting

# Resource parameters
# c inflow, d discharge, q catchability
c1, d1 = 50, 50
c2, d2 = 20, 50
q = 1

delta = 0.2 # Resource leakage


# Network 

#state of node


#fraction of cooperators initial
# Use this to probabilistically assign
#Also assume 
F1prob, F2prob = 0.7, 0.4

states = ['cooperator', 'defector']
prob1 = [F1prob, 1-F1prob]
prob2 = [F2prob, 1-F2prob]

G1 = nx.gnm_random_graph(n, m)
G2 = nx.gnm_random_graph(n, m)

# Need network to have no isolated nodes
# If there is an isolated node, find it and link it the next node
for k in range(n):
    if len([j for j in G1.neighbors(k)]) == 0:
        G1.add_edge(k,(k+1))
    else:
        pass
    
    if len([j for j in G2.neighbors(k)]) == 0:
        G2.add_edge(k,(k+1))
    else:
        pass

# Randomly assign cooperators and defectors
for k in range(n):
    G1.node[k]['state'] = random.choices(states, prob1)[0]
    G2.node[k]['state'] = random.choices(states, prob2)[0]
    

        
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


def utilNode(fc, R, k, G):
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



t = 0
tEnd = 600 #end point
dt = 0.1 #time step
# Lists to store values to plot later
time = []
r1List = []
r2List = []
f1List = []
f2List = []
# Draw initial network

# Sort out colours
colours1 = []
colours2 = []
    
for k in range(n):
    if G1.node[k]['state'] == 'cooperator':
        colours1.append('b')
    else: 
        colours1.append('r')
        
    if G2.node[k]['state'] == 'cooperator':
        colours2.append('b')
    else: 
        colours2.append('r')

nx.draw_networkx(G1, pos=None, node_color=colours1, with_labels=True, node_size=150)
plt.title('Initial Network 1')
plt.show()

nx.draw_networkx(G2, pos=None, node_color=colours2, with_labels=True, node_size=150)
plt.title('Initial Network 2')

plt.show()



while t<tEnd:
    F1 = countCoops(G1)
    F2 = countCoops(G2)
    
    dR1 = (c1 - d1*(R1/R1max)**2 - q*ext(F1)*R1 + delta*(R2-R1))*dt
    dR2 = (c2 - d2*(R2/R2max)**2 - q*ext(F2)*R2 + delta*(R1-R2))*dt
    
    #Choose random node in each group to update
    for (G,F,R) in [(G1,F1,R1), (G2,F2,R2)]:
        updateNode = random.choice(range(n))
        compareNode = random.choice([k for k in  G.neighbors(updateNode)])
    
        # Calculate node utilities
        utilUpdate = utilNode(F, R, updateNode, G)
        utilCompare = utilNode(F, R, compareNode, G)
    
        delta = utilUpdate - utilCompare
    
        if delta<0:
            prob = np.abs(delta)/(np.abs(utilUpdate)+np.abs(utilCompare))
            states = [G.node[updateNode]['state'], G.node[compareNode]['state']]
            probs = [1-prob, prob]
            newState = random.choices(states, probs)[0]
            G.node[updateNode]['state'] = newState
    
    #update quantities
    R1 += dR1
    R2 += dR2
    t += dt
    
    #append lists
    time.append(t)
    r1List.append(R1)
    r2List.append(R2)
    f1List.append(F1)
    f2List.append(F2)


# Reset colours
colours1 = []
colours2 = []
    
for k in range(n):
    if G1.node[k]['state'] == 'cooperator':
        colours1.append('b')
    else: 
        colours1.append('r')
        
    if G2.node[k]['state'] == 'cooperator':
        colours2.append('b')
    else: 
        colours2.append('r')



# Plot
plt.plot(time, f1List, color='blue', label='Population 1')
plt.plot(time, f2List, color='red', label='Population 2')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Fraction of Cooperators')

plt.show()

plt.plot(time, r1List, color='lightgreen', label='Resource 1')
plt.plot(time, r2List, color='darkgreen', label='Resource 2')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Resource Stock')


plt.show()



nx.draw_networkx(G1, pos=None, node_color=colours1, with_labels=True, node_size=150)
plt.title('Final Network 1')
plt.show()

nx.draw_networkx(G2, pos=None, node_color=colours2, with_labels=True, node_size=150)
plt.title('Final Network 2')

plt.show()
