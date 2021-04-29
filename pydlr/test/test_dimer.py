
import numpy as np

from pydlr import dlr
    
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
V = 4

E_aa = np.array([
    [e1, V],
    [V, e2],
    ])

def fermi_function(E, beta):
    p = np.argwhere(E > 0)
    m = np.argwhere(E <= 0)
    f = np.zeros_like(E)
    f[p] = np.exp(-beta*E[p]) / (1. + np.exp(-beta*E[p]))
    f[m] = 1. / (np.exp(beta*E[m]) + 1.)
    return f

E, U = np.linalg.eigh(E_aa)
rho = U.conj().T @ np.diag(fermi_function(E, beta)) @ U
n00 = rho[0, 0]

print(f'rho =\n {rho}')
print(f'n00 = {n00} (anal)')

# -- Solve test case with Legendre methods

from pypsitop.legendre_spectral import LegendreSpectral
from pypsitop.legendre_spectral_numba import free_greens_function

G_naa = free_greens_function(E_aa, beta, N, xi).real
G00_n = G_naa[:, 0, 0]

ls = LegendreSpectral(N)

G_iaa = ls.to_collocation(G_naa)
G00_i = G_iaa[:, 0, 0]
tau_i = ls.get_tau_i(beta)

n00_leg = -G00_i[-1].real
print(f'n00 = {n00_leg} (leg)')
np.testing.assert_array_almost_equal(n00, n00_leg)

T_nl = ls.to_matsubara(iwn, beta)
G_waa = np.tensordot(T_nl, G_naa, axes=(1,0))

G00_w = G_waa[:, 0, 0]
G00_w_anal = 1./ ( iwn - e1 - V**2/( iwn - e2 ) )
np.testing.assert_array_almost_equal(G00_w, G00_w_anal)

# -- Use DLR basis

d = dlr(lamb=lamb)

tau_l = d.get_tau(beta)
G_laa = d.tau_from_legendre(G_naa)
G00_l = G_laa[:, 0, 0] # Imaginary time Green's function on the DLR time points

G00_x = d.dlr_from_tau(G00_l)

np.testing.assert_array_almost_equal(G00_l, d.tau_from_dlr(G00_x))
np.testing.assert_array_almost_equal(G00_l, d.eval_dlr_tau(G00_x, tau_l, beta))

G00_i_ref = d.eval_dlr_tau(G00_x, tau_i, beta)

w_q = d.get_matsubara_frequencies(beta)
G00_q = d.matsubara_from_dlr(G00_x, beta)

np.testing.assert_array_almost_equal(G00_q, d.eval_dlr_freq(G00_x, w_q, beta))

G00_q_anal = 1./ ( w_q - e1 - V**2/( w_q - e2 ) )
np.testing.assert_array_almost_equal(G00_q, G00_q_anal)

G00_x_ref = d.dlr_from_matsubara(G00_q_anal, beta)

G00_l_ref = d.tau_from_dlr(G00_x_ref)
np.testing.assert_array_almost_equal(G00_l, G00_l_ref)

G00_q_ref = d.matsubara_from_dlr(G00_x_ref, beta)
np.testing.assert_array_almost_equal(G00_q, G00_q_ref)

# -- Viz

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 9))

subp = [4, 1, 1]

plt.subplot(*subp); subp[-1] += 1
plt.plot(tau_i, np.abs(G00_i - G00_i_ref), '.-')
plt.semilogy([], [])

plt.subplot(*subp); subp[-1] += 1
plt.plot(tau_i, G00_i, label=r'$G(\tau)$')
plt.plot(tau_l, G00_l, '.', label=r'$G(\tau_l)$ DLR points')
plt.xlabel(r'$\tau$')
plt.legend()

plt.subplot(*subp); subp[-1] += 1
l1 = plt.plot(iwn.imag, G00_w.real, '-', label='Re')
l2 = plt.plot(iwn.imag, G00_w.imag, '-', label='Im')
cre = l1[0].get_color()
cim = l2[0].get_color()
plt.plot(w_q.imag, G00_q.real, '.', color=cre)
plt.plot(w_q.imag, G00_q.imag, '.', color=cim)
plt.plot([], [], '.', color='gray', label='DLR')
plt.xlabel(r'$i\omega_n$')
plt.ylabel(r'$G(i\omega_n)$')
plt.legend()

plt.subplot(*subp); subp[-1] += 1
plt.plot(np.abs(G00_n), '.-', label='Leg')
plt.plot(np.abs(G00_x.real), '.-', label='Re DLR (tau)')
plt.plot(np.abs(G00_x.imag), '.-', label='Im DLR (tau)')

plt.plot(np.abs(G00_x_ref.real), 'x-', label='Re DLR (tau)')
plt.plot(np.abs(G00_x_ref.imag), 'x-', label='Im DLR (tau)')

plt.semilogy([], [])
plt.legend()
plt.xlabel(r'$n$')
plt.ylabel(r'$|G_n|$ expansion coeff')

plt.tight_layout()
plt.show()
