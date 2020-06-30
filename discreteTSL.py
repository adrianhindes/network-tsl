import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

#total no. agents
n = 50

# Maximum amount of resource
Rmax = 200

# Social parameters

mu = 4 # degree of cheating
ec = 0.483/n #level of effort (cooperators) #level of effort (defectors)
ed = mu*ec
w = 15 #cost of harvesting

# Resource parameters
# c inflow, d discharge, q catchability
c, d, q = 50, 50, 1


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

#Utility functions
def utilCoop(fc,R):
    E = ext(fc)
    uc = ec*((cobbdoug(E,R)/E)-w)
    return uc

def utilDefec(fc,R):
    E = ext(fc)
    ud = ed*((cobbdoug(E,R)/E)-w)-gompertz(fc)*(ed-ec)/ed
    return ud

# harvest-cooperator equations
def deriv(t,y):
    F, R = y
    dRdt = c - d*(R/Rmax)**2 - q*ext(F)*R
    dFdt = F*(1-F)*(utilCoop(F,R)-utilDefec(F,R))
    return dFdt, dRdt



# initial condition
#fraction of cooperators initial
F = 0.55
#amount of resource available initial
R = 100

t = 0
tEnd = 1000 #end point
dt = 0.1 #time step
# Lists to store values to plot later
time = []
fList = []
rList = []

while t<tEnd:
    dF = F*(1-F)*(utilCoop(F,R)-utilDefec(F,R))*dt
    dR = (c - d*(R/Rmax)**2 - q*ext(F)*R)*dt
    
    #update quantities
    F += dF
    R += dR
    t += dt
    
    #append lists
    time.append(t)
    fList.append(F)
    rList.append(R)




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
