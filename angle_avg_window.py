#!/usr/bin/python
"""
Calculate angle-averaged tophat window function.
"""
import numpy as np
import scipy.integrate
from scipy.special import jn
import pylab as P

w = 10.
dr = 10.

# Integration:
#int k^2 dk mu dmu dphi

def window(k, theta_k, phi_k, x0=(0., 0., 0.)):
    """
    Fourier-space tophat window function.
    """
    x0, y0, z0 = x0
    kx = k * np.sin(theta_k) * np.cos(phi_k)
    ky = k * np.sin(theta_k) * np.sin(phi_k)
    kr = k * np.cos(theta_k)
    return np.sinc(kx * w) * np.sinc(ky * w) * np.sinc(kr * dr) \
         * np.cos(kx*x0 + ky*y0 + kr*z0)

def avg_window(k, nsamp=(100, 90), x0=(0., 0., 0.)):
    """
    Return the angle-averaged window function for a given k.
    """
    # Evaluate window fn. on a grid in theta_k, phi_k
    theta_samp, phi_samp = nsamp
    theta_k = np.linspace(0., np.pi, theta_samp)
    phi_k = np.linspace(0., 2.*np.pi, phi_samp)
    TH, PH = np.meshgrid(theta_k, phi_k)
    
    # Calculate window function for this k on the (theta, phi) grid
    _w = window(k, TH, PH, x0) * np.sin(TH) # Weight by diff. volume element
    y = [scipy.integrate.simps(_w[i], theta_k) for i in range(_w.shape[0])]
    return scipy.integrate.simps(y, phi_k)

def avg_window3d(k, nsamp=(100, 90), x0=(0., 0., 0.)):
    """
    Return the angle-averaged window function for a given k.
    """
    #FIXME
    # Evaluate window fn. on a grid in theta_k, phi_k
    theta_samp, phi_samp = nsamp
    theta_k = np.linspace(0., np.pi, theta_samp)
    phi_k = np.linspace(0., 2.*np.pi, phi_samp)
    TH, PH = np.meshgrid(theta_k, phi_k)
    
    # Calculate window function for this k on the (theta, phi) grid
    _w = window(k, TH, PH, x0) * np.sin(TH) # Weight by diff. volume element
    y = [scipy.integrate.simps(_w[i], theta_k) for i in range(_w.shape[0])]
    return scipy.integrate.simps(y, phi_k)

# Plotting
P.subplot(111)
kvals = np.logspace(-4., 0., 300)

#dr = 1.
#w = 1.
#P.plot(kvals, [avg_window(_k) for _k in kvals], 'r-', lw=1.8)
#w = 10.
#P.plot(kvals, [avg_window(_k) for _k in kvals], 'g-', lw=1.8)
#w = 100.
#P.plot(kvals, [avg_window(_k) for _k in kvals], 'b-', lw=1.8)

dr = 100.
#w = 1.
#P.plot(kvals, [avg_window(_k) for _k in kvals], 'r-', lw=1.8, alpha=0.5)
#w = 10.
#P.plot(kvals, [avg_window(_k) for _k in kvals], 'g-', lw=1.8, alpha=0.5)
w = 100.
P.plot(kvals, [avg_window(_k) for _k in kvals], 'b-', lw=1.8, alpha=0.5)
P.plot(kvals, [avg_window(_k, x0=(50.,0.,0.)) for _k in kvals], 'r-', lw=1.8, alpha=0.5)
P.plot(kvals, [avg_window(_k, x0=(100.,0.,0.)) for _k in kvals], 'g-', lw=1.8, alpha=0.5)

P.plot(kvals, [avg_window(_k, x0=(200.,0.,0.)) for _k in kvals], 'y-', lw=1.8, alpha=0.5)

P.plot(kvals, [avg_window(_k, x0=(300.,0.,0.)) for _k in kvals], 'k-', lw=1.8, alpha=0.5)

#kk = np.logspace(-4., 0., 500)
#P.plot(kk, 10.*np.pi*np.sinc(kk * dr), 'b-', lw=1.8)

P.xscale('log')

P.tight_layout()
P.show()
