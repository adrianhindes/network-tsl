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

ec = 0.483/n #level of effort (cooperators)
ed = 1.826/n #level of effort (defectors)
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
    ud = ed*((cobbdoug(E,R))/E-w)-gompertz(fc)*(ed-ec)/ed
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


# solve ODE
interval = (t[0],t[-1])
ret = spi.solve_ivp(deriv, interval, y0, t_eval=t, method = 'RK45')

F = ret.y[0,:]
R = ret.y[1,:]

# Plot
plt.subplot(211)
plt.plot(t, F, 'b', lw=2)

plt.xlabel('Time')
plt.ylabel('Fraction Cooperators')

plt.subplot(212)
plt.plot(t,R,'r', lw=2)

plt.xlabel('Time')
plt.ylabel('Resource Stock')


plt.show()
