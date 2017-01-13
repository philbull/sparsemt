#!/usr/bin/python
"""
More convenient Python wrappers for sparse sampling C functions.
"""
import numpy as np
import sparse

def mode_cpl(q, k, kp, w, dr):
    """
    Calculate mode-coupling function F(k, k'; q), which is the integral over 
    the angular parts of the convolution with the window function.
    """
    kk, theta_k, phi_k = k
    kkp, theta_kp, phi_kp = kp
    
    return sparse.mode_cpl_fn(
             q,
             np.atleast_1d(kk), np.atleast_1d(theta_k), np.atleast_1d(phi_k), 
             np.atleast_1d(kkp), np.atleast_1d(theta_kp), np.atleast_1d(phi_kp),
             w, dr).T

def symm_grid():
    """
    Return list of (k, k') pairs needed to evaluate the covariance, assuming 
    it is symmetric.
    """
    
def repack_grid():
    """
    Re-pack grid elements output by symm_grid() into a full square matrix.
    """
    
