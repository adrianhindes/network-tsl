import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

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


# harvest-cooperator equations
def deriv(t, y):
    fc, R = y
    E = n*(fc*ec+(1-fc)*ed) # Extraction
    
    gamma = 10
    a = 0.6
    b = 0.2
    cobbdoug = gamma*(E**a)*(R**b)  
    
    h = 0.34
    t = -150
    g = -10
    gomp = h * np.exp(t * np.exp(g*fc))
    
    uc = ec*(cobbdoug/E)-w #Cooperator utility
    ud = ed*(cobbdoug/E)-w-gomp*(ed-ec)/ed #Defector utility
    
    dRdt = c - d*(R/Rmax)**2 - q*E*R
    dfcdt = fc*(1-fc)*(uc-ud)
    
    return [dRdt,dfcdt]


# initial condition
y0 = fc0, R0

# time points
t = np.linspace(0, 200,num=200)



ode =  spi.ode(deriv)

# BDF method suited to stiff systems of ODEs

t_end = 15000.
t_start = 1.
t_step = 1.
t_interval = np.arange(t_start, t_end, t_step)

ode.set_integrator('vode',nsteps=500,method='bdf')
ode.set_initial_value(y0,t_start)

ts = []
ys = []

while ode.successful() and ode.t < t_end:
    ode.integrate(ode.t + t_step)
    ts.append(ode.t)
    ys.append(ode.y)

t = np.vstack(ts)
F,R = np.vstack(ys).T



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


