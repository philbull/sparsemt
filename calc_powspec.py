#!/usr/bin/python
"""
Calculate the (correlated) power spectrum within some tophat window function 
on the sky.
"""
import numpy as np
import pylab as P
import scipy.interpolate
import scipy.integrate
import time

#import sys
#sys.path.append("src/")
import sparse

cumtrapz = scipy.integrate.cumtrapz

NQ = 400 #400 #300
NTHETA = 200
NPHI = 250

# Window function parameters
dr = 200. # Mpc
w = 200.
x0a = (0., 0., 0.)
x0b = (0., 0., 0.)


# FIXME: Would it make more sense to just precompute the F(k, k'; q) part? It's cosmology-independent, and can maybe be rescaled?

# FIXME: Does it matter if the patches of survey A and survey B overlap? They're sampling the same *aliased* density field...


def load_pk(fname="camb_pk_z0.dat"):
    """
    Load power spectrum from file.
    """
    # Load power spectrum from file
    _k, _pk = np.genfromtxt(fname).T
    # FIXME: Check conversion from h^-1 units
    h = 0.67
    _k *= h
    _pk *= h**3.
    
    # Interpolate in log-space and return interpolation function
    interp_pk = scipy.interpolate.interp1d( np.log(_k), np.log(_pk),
                                            bounds_error=False,
                                            fill_value=0. )
    pk = lambda k: np.exp( interp_pk(np.log(k)) )
    return pk

def sph_to_cart(kvec):
    """
    Convert spherical coordinates to Cartesian coordinates.
    """
    k, theta_k, phi_k = kvec
    kx = k * np.sin(theta_k) * np.cos(phi_k)
    ky = k * np.sin(theta_k) * np.sin(phi_k)
    kr = k * np.cos(theta_k)
    return kx, ky, kr


def window_tophat(kx, ky, kr, x0=(0., 0., 0.), sgn=1.):
    """
    Fourier-space tophat window function.
    """
    x0, y0, z0 = x0
    
    # FIXME: Missing factor of (Delta r) . w^2
    return np.sinc(kx * w) * np.sinc(ky * w) * np.sinc(kr * dr) \
         * np.cos(sgn * (kx*x0 + ky*y0 + kr*z0))
    #     * np.exp(sgn*1.j * (kx*x0 + ky*y0 + kr*z0)) # FIXME


def mode_cpl_fn(k, kp):
    """
    Calculate the mode-coupling function in (|k|, |k|') by integrating P(k) 
    over the Fourier-space window functions.
    """
    # Difference between k and q vectors
    kx, ky, kr = sph_to_cart(k)
    kxp, kyp, krp = sph_to_cart(kp)
    
    # Grid in q
    q_theta = np.linspace(0., np.pi, NTHETA)
    q_phi = np.linspace(0., 2.*np.pi, NPHI)
    QTH, QPHI = np.meshgrid(q_theta, q_phi)
    
    dth = q_theta[1] - q_theta[0]
    dphi = q_phi[1] - q_phi[0]
    
    # Loop over values of q
    wq = []
    for q in qgrid:
        Q = q * np.ones(QTH.shape)
        QX, QY, QR = sph_to_cart((Q, QTH, QPHI))
        
        # Construct 2D integrand
        ww = window_tophat(kx - QX, ky - QY, kr - QR, x0a, sgn=1.) \
           * window_tophat(kxp - QX, kyp - QY, krp - QR, x0b, sgn=-1.) \
           * np.sin(QTH) # FIXME
        
        wq.append( dth * dphi * np.sum(ww).real )
        #P.matshow(ww.real, cmap='RdBu', origin='lower', 
        #          extent=[np.min(q_theta), np.max(q_theta),
        #                  np.min(q_phi), np.max(q_phi)])
        #P.colorbar()
        #P.tight_layout()
        #P.show()
    
    return np.array(wq)

# Load power spectrum interpolator
pk = load_pk("camb_pk_z0.dat")
vol = dr * w**2.

qgrid = np.logspace(-3., np.log10(0.7), NQ) #0.2
kgrid = np.logspace(-3., np.log10(0.25), 40)
#kgrid = np.logspace(-3., np.log10(0.25), 10) # Reasonable kmax is 0.1

K, KP = np.meshgrid(kgrid, kgrid)
ZERO = np.zeros(K.shape).flatten()
#zero = np.zeros(kgrid.shape)
#one = np.ones(kgrid.shape)

# Loop over Fourier k (set all angles to zero)
print "Calculating F(k, k'; q) over grid..."
t0 = time.time()
Fkk = sparse.mode_cpl(qgrid, [K.flatten(), ZERO, ZERO], 
                             [KP.flatten(), ZERO, ZERO], w, dr)
#Fkk = sparse.mode_cpl(qgrid, [kgrid, zero, zero], 
#                             [kgrid, zero, zero], w, dr)
print "\tFinished in %d sec." % (time.time() - t0)


# Integrate over q for each grid element
print "Integrating over q..."
t0 = time.time()
Wkk = [ scipy.integrate.trapz(Fkk[:,i] * pk(qgrid) * qgrid**2., qgrid)
        for i in range(Fkk.shape[1]) ]
Wkk = np.array(Wkk).reshape(K.shape)
Wkk *= vol
#Wkk = vol * np.array(Wkk)
print "\tFinished in %1.1f sec." % (time.time() - t0)

"""
ax1 = P.subplot(211)
ax2 = P.subplot(212)
cols = ['r', 'y', 'g', 'b', 'm']
for i in range(Fkk.shape[1]):
    #y = scipy.integrate.cumtrapz(Fkk[:,i] * pk(qgrid) * qgrid**2., qgrid, initial=0.)
    ax1.plot(qgrid, Fkk[:,i] * pk(qgrid) * qgrid**2., color=cols[i], lw=1.8)
    
    ax2.plot(qgrid, scipy.integrate.cumtrapz(
                        Fkk[:,i] * pk(qgrid) * qgrid**2., 
                        qgrid, initial=0.),
             color=cols[i], lw=1.8)
    
    #_Fkk = mode_cpl_fn([kgrid[i], 0., 0.], [kgrid[i], 0., 0.])
    #P.plot(qgrid, _Fkk * pk(qgrid) * qgrid**2., color=cols[i], marker='.')
    
    ax2.axhline(pk(kgrid[i]) / (w**2. * dr)**1., color=cols[i], lw=1.5, alpha=0.5)
    ax1.axvline(kgrid[i], lw=2., color=cols[i], alpha=0.5)
    ax2.axvline(kgrid[i], lw=2., color=cols[i], alpha=0.5)

#P.subplot(111)
#P.plot(kgrid, Wkk, lw=1.8)
for ax in [ax1, ax2]:
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim((1e-12, 1e0))

P.show()
"""
"""
# Compare diagonals of power spectra
P.subplot(111)

P.plot(kgrid, pk(kgrid), 'k', lw=1.8)
P.plot(kgrid, Wkk, 'r--', lw=1.8)

P.xscale('log')
P.yscale('log')
P.ylim((1e2, 1e5))
P.tight_layout()
P.show()
"""

# Output to file
np.savetxt("wkk.dat", Wkk)
exit()


P.matshow(Wkk, extent=[np.log10(np.min(kgrid)), np.log10(np.max(kgrid)),
                       np.log10(np.min(kgrid)), np.log10(np.max(kgrid))])
P.colorbar()
P.show()

exit()

#Wkk = 

P.legend(loc='upper left')
P.xlabel("q [Mpc^-1]")
P.ylabel("W(k, k'; q)")

P.xscale('log')
P.yscale('log')
P.ylim((1e-10, 1e0))
P.tight_layout()
P.show()
