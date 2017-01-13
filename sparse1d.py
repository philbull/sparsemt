#!/usr/bin/python
"""
Construct Fourier-space covariance matrix for sparsely-sampled 1D survey.
"""
import numpy as np
import pylab as P
import scipy.interpolate

# Redefine sinc to exclude pi factor
sinc = lambda x: np.sinc(x / np.pi)

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
    pk = lambda k: np.exp( interp_pk(np.log(np.abs(k))) )
    return pk

def grid_resolution(w, x0):
    """
    Set resolution of integrator grid according to typical width of sinc 
    function or separation of patches, whichever is larger (as this sets the 
    max. frequency of oscillations in Fourier space).
    
    Returns number of q samples to use.
    """
    if len(x0) == 1:
        dx = w
    else:
        dx = np.max([w, np.max(np.diff(x0))])
    #NQ = int(1000 * (dx / 100.))
    NQ = int(400 * (dx / 100.))
    if NQ < 200: NQ = 200
    return NQ

def window(x, w, x0):
    """
    Return window function in real space, as a function of x.
    """
    win = np.zeros(x.shape)
    for _x0 in x0:
        win[np.where(np.logical_and(x >= _x0 - w, x <= _x0 + w))] = 1.
    return win

def modecpl_mat(q, k, kp, x0, x0p, w):
    """
    Construct mode-coupling matrix in 1D for two tophat window functions.
    """
    # Loop over patch origins
    x0_term = np.sum( [np.exp( -1.j * _x0 * (k - q) ) 
                       for _x0 in x0], axis=0 )
    x0p_term = np.sum( [np.exp( +1.j * _x0p * (kp - q) ) 
                       for _x0p in x0p], axis=0 )
    return sinc(w*(k-q)) * sinc(w*(kp-q)) * x0_term * x0p_term

def corrmat(m):
    """
    Calculate correlation matrix from a covariance matrix.
    """
    corr = m.copy()
    for i in range(corr.shape[0]):
        for j in range(corr.shape[1]):
            corr[i,j] = m[i,j] / np.sqrt(m[i,i] * m[j,j])
    return corr

def covmat(kvals, kpvals, x0, x0p):
    """
    Calculate covariance matrix in (k, k') by convolving input power spectrum 
    with window functions.
    """
    # Pre-calculate power spectrum on q grid
    qvals = np.logspace(-4., 0.5, NQ)
    qvals = np.concatenate((-qvals[::-1], qvals))
    pq = pk(np.abs(qvals))
    
    # Loop over (k, k') pairs
    cov = []
    i = -1
    K, KP = np.meshgrid(kvals, kpvals)
    for k, kp in zip(K.flatten(), KP.flatten()):
        i += 1
        #print "\tk = %2.2e, k' = %2.2e" % (k, kp)
        if i % 500 == 0: print "%d / %d" % (i, K.size)
        
        # Get mode-coupling matrix as a function of q, then integrate over P(q)
        y = modecpl_mat(qvals, k, kp, x0, x0p, w)
        cov.append( scipy.integrate.trapz(y * pq, qvals) )
    
    # Calculate effective volume
    veff = float(len(x0)) * 2. * w
    vol =  8.*np.pi * veff #1. #(nvol * 2.*w)**2. / (2.*(2.*np.pi)**2. * w * nvol)
    # FIXME: Where are some of these factors from?
    
    # Repackage into 2D arrays and return
    return 2.*np.pi * (2.*w)**2. * np.array(cov).reshape(K.shape) / vol
    

#FIXME: Shifting the origin with only one tophat doesn't give the same result!

# Define window function properties
w = 200.
x0 = [0.,]

CASE = 'A'

if CASE == 'A':
    # CASE A: Multiple small patches with large separation
    #w = 50.
    #w = 150.
    w = 50.
    #x0 = np.arange(-5., 6.) * 10.*w
    #x0 = np.arange(-5., 6.) * 3.5*w
    
    #x0 = [-1000., 1000.]
    #x0 = [-2050., -1150., -650., -350., -150., 0., 150., 350., 650., 1150., 2050.]
    x0 = [-1150., -650., -350., -150., 0., 150., 350., 650., 1150.]
    #x0 = [0.,]
    ##w = 50.*11.
    ##print np.diff(x0)
    print "dx = ", (np.max(x0) + w) - (np.min(x0) - w)
    #exit()
elif CASE == 'B':
    # CASE B: Uniform window over the *total* size of the survey
    w = 1020.
    x0 = [0.,]
elif CASE == 'C':
    # CASE C: Single window of the same size as the Case A multi-patch separation
    w = 200.
    x0 = [0.,]
elif CASE == 'D':
    # CASE D: Single window of the same size as a single patch of Case A
    w = 20.
    x0 = [0.,]

x0p = x0


# Set q grid resolution according to expected oscillation freq. of window fns.
NQ = grid_resolution(w, x0)
print "No. samples in q:", NQ

# Load power spectrum
pk = load_pk(fname="camb_pk_z0.dat")

# Define k-space grid
#k = np.logspace(-3., 0., 150) # 30
k = np.logspace(-3., 0., 500) # 30
cov = covmat(k, k, x0, x0p)

# Save covmat to file
np.savetxt("test_w50.cov", cov.real, header=" ".join(["%e" % _k for _k in k]))

#-------------------------------------------------------------------------------

# Plot power spectrum (diagonal)
#P.subplot(111)
#P.plot(k, pk(k), 'k-', lw=1.8)
#P.plot(k, np.diag(np.real(cov)), 'r--', lw=1.8)
#P.axhline(1.)
#P.xscale('log')
#P.yscale('log')
#P.ylim((1e1, 1e5))
#P.ylim((0., 2.))
#P.show()


# Plot signal correlation matrix
P.matshow(corrmat(cov.real), cmap='RdBu', vmin=-1., vmax=1.,
          extent=[np.log10(np.min(k)), np.log10(np.max(k)), 
                  np.log10(np.max(k)), np.log10(np.min(k))], 
          )

# Plot signal covariance matrix
#P.matshow(np.log10(np.abs(cov.real)), cmap='Greys', 
#          extent=[np.log10(np.min(k)), np.log10(np.max(k)), 
#                  np.log10(np.max(k)), np.log10(np.min(k))], 
#          )

P.colorbar()
P.grid()

P.show()
