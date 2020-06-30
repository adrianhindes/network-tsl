import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

#total no. agents
n = 50
#fraction of cooperators initial
fc0 = 0.7
#amount of resource available initial
R0 = 100
# Maximum amount of resource
Rmax = 200

# Social parameters

mu = np.linspace(2,4,num=20) # degree of cheating
ec = 0.483/n #level of effort (cooperators) #level of effort (defectors)
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
y0 = [fc0, R0]

# time points
t = np.linspace(0, 1000,num=1000)
interval = (t[0],t[-1])

# solve ODE

F = np.zeros(len(mu))
R = np.zeros(len(mu))

for j in range(len(mu)):
    ed = mu[j]*ec
    # Solve TSL model for each value of mu
    ret = spi.solve_ivp(deriv, interval, y0, t_eval=t)
    # Store final values for F and R (assuming convergence)
    F[j] = ret.y[0,-1]
    R[j] = ret.y[1,-1]


# Plot
plt.subplot(211)
plt.plot(F, mu, 'o', color='black')

plt.xlabel('Fraction Cooperators')
plt.ylabel('Degree of cheating mu')


plt.show()
