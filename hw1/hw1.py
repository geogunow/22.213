import numpy as np
from math import *
from copy import deepcopy as copy
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.integrate import quadrature

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size' : 14}
plt.rc('font', **font)

'''
Class object used to caclulate reactivity. It is defined by a set of input
<x> and output <y> values. Calling the calc function returns y(x0) where x0
is the requested input value.
'''
class interp_function(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def calc(self, x0):
        return np.interp(x0, self.x, self.y)
'''
A function that solves the point kinetics neutron diffusion problem
'''
def pke_solve(P0, lambdas, betas, gamma, rho, time):

    # correct time so initial time is 0
    time = time - time[0]
    
    # setup initial matrix
    L = len(lambdas)
    A0 = np.matrix( np.zeros( (L+1, L+1) ) )
    beta = sum(betas)
    A0[0, 0] = -beta/gamma
    for i in range(L):
        A0[0, i+1] = lambdas[i]
        A0[i+1, i+1] = -lambdas[i]
        A0[i+1, 0] = betas[i]/gamma


    # initialize solution
    T = len(time)
    N = np.matrix( np.zeros( (L+1, T) ) )

    # set initial conditions
    N[0, 0] = P0
    for i in range(L):
        N[i+1, 0] = P0 * betas[i] / (gamma * lambdas[i])

    # march along time steps
    for i, t in enumerate(time[1:]):
        
        # create appropriate matrix for the timestep
        A = copy(A0)
        A[0, 0] += rho(t) * beta / gamma

        # solve for the next time step
        dt = t - time[i]
        N[:, i+1] = expm( dt*A ) * N[:, i]

    return N

'''
A function that solves the inverse kinetics neutron diffusion problem
'''
def ike_solve(lambdas, betas, gamma, P, time):

    # correct time so initial time is 0
    time = time - time[0]
    T = len(time)
    
    # setup initial precursor arrays
    L = len(lambdas)
    C = np.zeros( (L, T) )
    beta = sum(betas)

    # set initial conditions
    for i in range(L):
        C[i, 0] = P(0) * betas[i] / (gamma * lambdas[i])

    # initialize reactivity
    rho = np.zeros(T)

    # march along time steps
    for n, t in enumerate(time[1:]):

        # add power component to reactivity
        dt = t - time[n]
        rho[n+1] = gamma/P(t) * (P(t) - P(time[n]))/dt + beta
        
        # calculate precurssor concentrations
        for i in range(L):
            
            # define integrand of production term
            def integrand(tp):
                return P(tp) * np.exp(-lambdas[i]*(t-tp))

            # calculate new precurssor concentrations
            e = 10**-5
            production =  betas[i] / gamma \
                    * quadrature(integrand, time[n], t, tol = e, rtol = e)[0]
            C[i, n+1] = C[i, n] * exp(-lambdas[i]*dt) + production

            # subtract precurssor contributions
            rho[n+1] -= gamma/P(t)*lambdas[i]*C[i, n+1]

    return rho/beta


'''
script to solve parts A & B
'''

# set physical parameters
half_lives = [55.6, 24.5, 16.3, 5.21, 2.37, 1.04, 0.424, 0.195]
betas = [0.000218, 0.001023, 0.000605, 0.00131, 0.00220, 0.00060, 0.000540,
        0.000152]
nu = 2.45
sigf = 0.05
v = 2200*100*sqrt(0.1/0.0253)
powden = 100.0

# declare transient conditions
eps = 10**-6
rho_x = [0, 5, 7, 27, 29 , 100]
rho_y = [0, 0, 0.1, 0.1, 0, 0]
rho_fun = interp_function(rho_x, rho_y)
time = np.linspace(0,100,1000)

# calculate derived physical parameters
P0 = powden / sigf / (1.602*10**-13*200)
gamma = 1.0 / (sigf*nu*v)
lambdas = log(2.0) / np.array(half_lives)
print "Beta = ", sum(betas)
print "Prompt Neutron Lifetime = ", gamma

# solve transient
N = pke_solve(P0, lambdas, betas, gamma, rho_fun.calc, time)

# get power and reactivity
P = N[0,:].A1
reactivity = rho_fun.calc(time)

# plot power
fig, ax1 = plt.subplots()
ax1.plot(time, P/P0, 'k-', lw=3)
ax1.set_xlabel('Time (sec)')
ax1.set_ylabel('Relative Power', color = 'k')
for t1 in ax1.get_yticklabels():
    t1.set_color('k')

# plot reactivity
ax2 = ax1.twinx()
ax2.plot(time, reactivity, 'r-', lw=3)
ax2.set_ylabel('Reactivity ($)', color = 'r')
for t1 in ax2.get_yticklabels():
    t1.set_color('r')

plt.show()

# plot precursors
L = len(betas)
for i in range(L):
    C = N[i+1,:].A1
    plt.plot(time, C/C[0], lw=3)
plt.xlabel('Time (sec)')
plt.ylabel('Relative Precursor Population')
plt.show()


'''
Script to solve part C
'''
# define desired power
time = [0, 2-eps, 2+eps, 10, 12, 20]
P = [1, 1, 100, 100, 10, 10]
power = interp_function(time, P).calc
time = np.linspace(0,20,1000)

# solve for reactivity
reactivity = ike_solve(lambdas, betas, gamma, power, time)

# plot power
fig, ax1 = plt.subplots()
P = power(time) / P[0]
ax1.plot(time, P, 'k-', lw=3)
ax1.plot([time[0]], [max(P)*1.01], 'wx')
ax1.set_xlabel('Time (sec)')
ax1.set_ylabel('Relative Power', color = 'k')
for t1 in ax1.get_yticklabels():
    t1.set_color('k')

# plot reactivity
ax2 = ax1.twinx()
ax2.plot(time, reactivity, 'r-', lw=3)
ax2.set_ylabel('Reactivity ($)', color = 'r')
for t1 in ax2.get_yticklabels():
    t1.set_color('r')

plt.show()
