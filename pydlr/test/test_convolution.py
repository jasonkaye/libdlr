
import numpy as np

from pydlr import dlr, kernel
    
# -- Setup test system

N = 40
xi = -1.
beta = 1.
lamb = 10.
Nw = int(lamb)

n = np.arange(-Nw, Nw+1)
zeta = (1 - xi)/2
iwn = 1.j * np.pi/beta * (2*n + zeta) 

e1, e2 = 3., 3.3

def fermi_function(E, beta):
    p = np.argwhere(E > 0)
    m = np.argwhere(E <= 0)
    f = np.zeros_like(E)
    f[p] = np.exp(-beta*E[p]) / (1. + np.exp(-beta*E[p]))
    f[m] = 1. / (np.exp(beta*E[m]) + 1.)
    return f

# -- Solve test case with Legendre methods

from pypsitop.legendre_spectral import LegendreSpectral
from pypsitop.legendre_spectral_numba import free_greens_function

g1_naa = free_greens_function(np.array([[e1]]), beta, N, xi)
g2_naa = free_greens_function(np.array([[e2]]), beta, N, xi)

ls = LegendreSpectral(N)
tau_i = ls.get_tau_i(beta)
gg_naa = ls.convolution(g1_naa, g2_naa, beta, int(xi))

g1_i = np.squeeze(ls.to_collocation(g1_naa)).real
g2_i = np.squeeze(ls.to_collocation(g2_naa)).real
gg_i = np.squeeze(ls.to_collocation(gg_naa)).real

# -- DLR Matsubara convolution

d = dlr(lamb=lamb)

tau_l = d.get_tau(beta)
g1_l = d.tau_from_legendre(g1_naa[:, 0, 0])
g2_l = d.tau_from_legendre(g2_naa[:, 0, 0])

g1_x = d.dlr_from_tau(g1_l)
g2_x = d.dlr_from_tau(g2_l)

w_q = d.get_matsubara_frequencies(beta)
g1_q = d.matsubara_from_dlr(g1_x, beta)
g2_q = d.matsubara_from_dlr(g2_x, beta)

gg_q = g1_q * g2_q

gg_x = d.dlr_from_matsubara(gg_q, beta)
gg_l = d.tau_from_dlr(gg_x)
gg_i_ref = d.eval_dlr_tau(gg_x, tau_i, beta)

np.testing.assert_array_almost_equal(gg_i, gg_i_ref, decimal=9)

# -- DLR imaginary time convolution

gg_x_ref = d.convolution(g1_x, g2_x)
gg_l_ref = d.tau_from_dlr(gg_x_ref)

np.testing.assert_array_almost_equal(gg_l, gg_l_ref, decimal=9)

# --

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 9))

subp = [1, 1, 1]

plt.subplot(*subp); subp[-1] += 1
plt.plot(tau_i, g1_i, label='g1')
plt.plot(tau_i, g2_i, label='g2')
plt.plot(tau_i, gg_i, label='gg')

plt.plot(tau_l, g1_l, '.', label='g1 (dlr)')
plt.plot(tau_l, g2_l, '.', label='g2 (dlr)')
plt.plot(tau_l, gg_l, '.', label='gg (dlr)')

plt.plot(tau_l, gg_l_ref, 'x', label='gg (dlr conv)')

plt.legend(loc='best')

plt.show()
