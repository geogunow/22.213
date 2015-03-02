import numpy as np
from math import *
import matplotlib.pyplot as plt

'''
function for plotting fluxes, suspending program execution to display the plot
Arguments:
    delta:  spatial mesh used to construct the problem
    phi:    flux solution vector
    N:      number of mesh points in the problem
'''
def plot_fluxes(x, phi):
  
    # plot fluxes
    G = len(phi)
    for g in xrange(G):
        plt.plot(x, phi[g], 'x-')

    # create legend
    base = 'Flux '
    names = []
    for g in xrange(G):
        names.append(base + str(g+1))
    plt.legend(names)

    # label axes and show plot
    plt.xlabel('x (cm)')
    plt.ylabel('Normalized Flux')
    plt.show()

'''
Computes the dominance ratio of the matrix inv(A)*F
Arguments:
    A:  loss matrix
    F:  fission matrix
Output:
    dominance ratio
'''
def compute_dominance_ratio(A, F):
  
    print "Computing Dominance Ratio..."

    # compute full matrix (expensive)
    print "Inverting Matrix"
    M = np.linalg.inv(A) * F
    print "Matrix inversion complete"

    # compute eigenvalues
    eigs = np.real(np.linalg.eigvals(M))

    # sort list
    ordered = sorted(eigs)

    # compute dominance ratio
    eig1 = ordered[-1]
    eig2 = ordered[-2]
    dominance_ratio = eig2/eig1

    return dominance_ratio

'''
Calculates the L-2 norm of the error for the given cases relative to the
reference
Arguments:
    cases:      list of solution vectors to analyze
    reference:  reference solution vector
    N:          number of nodes involved in the solution
Output:
    A list of the RMS errors
'''
def calculate_RMS_errors(cases, reference, N):

    rms_errors = np.zeros( len(cases) )
    for i, case in enumerate(cases):
        error = 0
        for j, value in enumerate(case):
            if reference[j] != 0:
                error += (case[j]/reference[j] - 1)**2
        error /= N
        error = sqrt(error)
        rms_errors[i] = error

    return rms_errors

'''
Calculates error in scalars relative to reference
Arguments:
    cases:      list of values to analyze
    refernce:   refernce value
Output:
    A list of the relative errors
'''
def calculate_scalar_errors(cases, reference):
    return abs(cases - reference)/reference

'''
Computes the average parameter over nodes given a mesh refinement
Arguments:
    values:     list of values to average
    refinement: the step size over which values are averaged
Output:
    A numpy array of the averaged values
'''
def calculate_nodal_avg(values, refinement):
    
    L = len(values)
    avg_vals = np.zeros(L)
    
    for i in xrange(L):
        avg_vals[i] = sum(values[i*refinement:(i+1)*refinement]) / refinement

    return avg_vals

'''
plots two variables against the same x axis
'''
def dual_plot(x, y1, y2, xlabel='x', ylabel1='y1', ylabel2='y2'):
    
    # plot y1
    fig, ax1 = plt.subplots()
    ax1.plot(x, y1, 'kx-', lw=1)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1, color = 'k')
    for t1 in ax1.get_yticklabels():
        t1.set_color('k')

    # plot y2
    ax2 = ax1.twinx()
    ax2.plot(x, y2, 'rx-', lw=1)
    ax2.set_ylabel(ylabel2, color = 'r')
    for t1 in ax2.get_yticklabels():
        t1.set_color('r')
    plt.show()

'''
plots two variables against the same x axis, using logarithmic y axes
'''
def dual_log_plot(x, y1, y2, xlabel='x', ylabel1='y1', ylabel2='y2'):
    
    # plot y1
    fig, ax1 = plt.subplots()
    ax1.semilogy(x, y1, 'kx-', lw=1)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1, color = 'k')
    for t1 in ax1.get_yticklabels():
        t1.set_color('k')

    # plot y2
    ax2 = ax1.twinx()
    ax2.semilogy(x, y2, 'rx-', lw=1)
    ax2.set_ylabel(ylabel2, color = 'r')
    for t1 in ax2.get_yticklabels():
        t1.set_color('r')
    plt.show()

'''
plots two variables against the same logarithmic x axis, 
using logarithmic y axes
'''
def dual_log_log_plot(x, y1, y2, xlabel='x', ylabel1='y1', ylabel2='y2'):
    
    # plot y1
    fig, ax1 = plt.subplots()
    ax1.loglog(x, y1, 'kx-', lw=1)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1, color = 'k')
    for t1 in ax1.get_yticklabels():
        t1.set_color('k')

    # plot y2
    ax2 = ax1.twinx()
    ax2.loglog(x, y2, 'rx-', lw=1)
    ax2.set_ylabel(ylabel2, color = 'r')
    for t1 in ax2.get_yticklabels():
        t1.set_color('r')
    plt.show() 


'''
plots two variables against the same logarithmic x axis, 
using standard y axes
'''
def dual_logx_plot(x, y1, y2, xlabel='x', ylabel1='y1', ylabel2='y2'):
    
    # plot y1
    fig, ax1 = plt.subplots()
    ax1.semilogx(x, y1, 'kx-', lw=1)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1, color = 'k')
    for t1 in ax1.get_yticklabels():
        t1.set_color('k')

    # plot y2
    ax2 = ax1.twinx()
    ax2.semilogx(x, y2, 'rx-', lw=1)
    ax2.set_ylabel(ylabel2, color = 'r')
    for t1 in ax2.get_yticklabels():
        t1.set_color('r')
    plt.show() 
