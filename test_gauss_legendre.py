#!/usr/bin/python
"""
Test Gaussian quadrature integration of window functions.
"""
import numpy as np
import pylab as P
from scipy.special import legendre
import scipy.integrate

def weights(n, legendre):
    """
    Calculate the weight coefficients, w_i, for a Gauss-Legendre series on 
    [-1, 1], using the results of:
    http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
    """
    # Eq. 13
    w_i = 2. * (1. - xi**2.) / (n + 1)**2. / (leg[n+1](xi))**2.

def precompute():
    """
    x
    """
    pass

xi, wi = np.polynomial.legendre.leggauss(100)
x = np.linspace(-1., 1., 1000)

fn = lambda x: np.sinc(x)

for l in np.arange(20):
    
    gaussl = np.sum( wi*legendre(l)(xi) * fn(xi) )
    simps = scipy.integrate.simps(legendre(l)(x) * fn(x), x)
    
    if np.abs(gaussl) < 1e-14: gaussl = np.inf
    if np.abs(simps) < 1e-14: simps = 0.
    
    print "%2d %+8.8f" % (l, simps / gaussl)

