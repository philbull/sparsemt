#!/usr/bin/python
"""
Test that FT of window function gives sensible results.
"""
import numpy as np
import pylab as P

sinc = lambda x: np.sinc(x / np.pi)

def window(k, w, x0):
    """
    FT of tophat window function covering [x0 - w, x0 + w].
    """
    return 2.*w * sinc(k*w) * (np.cos(k*x0) + 1.j*np.sin(k*x0))


w = 5.
k = np.linspace(-1., 1., 1000)

P.subplot(111)

P.plot(k, window(k, w, x0=0.).real, 'r-', lw=1.8)
P.plot(k, window(k, w=2.5, x0=-2.5) + window(k, w=2.5, x0=+2.5), 'b--', lw=1.8)

#P.plot(k, window(k, w/2., x0=-w/2.).real, 'y-', lw=1.8, alpha=0.5)

#P.plot(k, 2.*w*np.sinc(k*w/2.)*np.cos(k*w/2.), 'y--', lw=1.8)

#P.plot(k, np.sinc(k*w), 'r-')
#P.plot(k, np.sin(k*w) / (k*w), 'y--')

#P.plot(k, window(k, w=2.5, x0=-2.5).imag, 'm-', lw=1.8, alpha=0.5)
#P.plot(k, window(k, w=2.5, x0=+2.5), 'b-', lw=1.8, alpha=0.5)

#P.plot(k, window(k, w=2.5, x0=-2.5) + window(k, w=2.5, x0=+2.5), 'y-', lw=1.8)


P.grid()
P.tight_layout()
P.show()
