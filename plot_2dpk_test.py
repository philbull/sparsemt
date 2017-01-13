#!/usr/bin/python
"""
Plot diagonal part of power spectrum.
"""
import numpy as np
import pylab as P
import scipy.interpolate

# Load power spectrum
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
pk = load_pk()

# Get 2D power spectrum
#kgrid = np.logspace(-3., np.log10(0.1), 40)
kgrid = np.logspace(-3., np.log10(0.25), 40)
#pk2d = np.genfromtxt("wkk_200x200.dat")
#pk2dsmall = np.genfromtxt("wkk_200x200_smallk.dat")
#pk2dsmall2 = np.genfromtxt("wkk_200x200_smallk2.dat")
pk2d = np.genfromtxt("wkk_200x200_acc.dat")
pk2d400 = np.genfromtxt("wkk_400x400_acc.dat")
V = 200.**3.


# Compare diagonals of power spectra
P.subplot(111)
#P.plot(kgrid, pk(kgrid), 'k', lw=1.8)
#P.plot(kgrid, np.diag(pk2d), 'r--', lw=1.8)

P.plot(kgrid, np.diag(pk2d)/pk(kgrid) - 1., 'k-', lw=1.8)
P.plot(kgrid, np.diag(pk2d400)/pk(kgrid) - 1., 'r--', lw=1.8)
P.axhline(0., color='k')

P.xscale('log')
#P.yscale('log')
P.tight_layout()

"""
# Plot correlation matrix
def corrmat(mat):
    newmat = np.zeros(mat.shape)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            newmat[i,j] = mat[i,j] / np.sqrt(mat[i,i] * mat[j,j])
    return newmat

print corrmat(pk2d)[-4:,-4:]

#P.subplot(111)
P.matshow(corrmat(pk2d), cmap='RdBu', vmin=-1., vmax=1.)
P.colorbar()
"""
P.show()
