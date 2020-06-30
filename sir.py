# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 12:21:28 2019
SEIR Model
Reinterpreted from MATLAB Code, MATH3501
@author: Adrian

Deterministic SIR Epidemic model
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population
N = 1000
# Initial number of infected and recovered
I0, R0 = 1, 0
# By def, all others are susceptible to infection
S0 = N - I0 - R0
# Contact rate and mean recovery rate (in 1/days)
beta, gamma = 0.2, 1./10
#grid of time points (in days)
t = np.linspace(0, 160,num=160)

# SIR Model 
def deriv(y, t, N, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

#initial conditions
y0 = S0, I0, R0
#Integrate SIR equations over timegrid
ret = odeint(deriv, y0, t, args=(N,beta,gamma))
S, I, R = ret.T

#plot data on three curves
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, axisbelow=True)
ax.plot(t, S/1000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I/1000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/1000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
ax.set_ylim(0,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
plt.show()