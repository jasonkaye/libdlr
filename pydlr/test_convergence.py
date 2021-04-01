"""
Solve for the density in a dimer problem
using analytical form in Matsuabara frequency space

Author: Hugo U.R. Strand (2021) 
"""


import numpy as np

from pydlr import dlr
    
# -- Setup test system

xi = -1.
beta = 100.

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

print(f'n00 = {n00:16.16E} (anal)')

# -- Use DLR basis
n_vec = []
dlr_vec = []
lamb_vec = [
    10., 100., 200., 400., 800., 1600., 3200.]

for lamb in lamb_vec:
    d = dlr(lamb=lamb)

    w_q = d.get_matsubara_frequencies(beta)
    delta_q = V**2/(w_q - e2)
    g_q = 1./(w_q - e1 - delta_q)
    g_x = d.dlr_from_matsubara(g_q, beta)
    g_beta = d.eval_dlr_tau(g_x, np.array([beta]), beta)
    n00_dlr = -g_beta[0].real
    
    print(f'n00 = {n00_dlr:16.16E} (dlr) lambda {lamb:2.2E}' + \
          f' ncoeff {d.rank:4d} error {np.abs(n00_dlr - n00):2.2E}')

    n_vec.append(n00_dlr)
    dlr_vec.append(d.rank)

n_vec = np.array(n_vec)

# -- Viz

import matplotlib.pyplot as plt

plt.figure(figsize=(6, 6))

subp = [1, 1, 1]
plt.subplot(*subp); subp[-1] += 1
plt.plot(lamb_vec, np.abs(n00 - n_vec), '.-')
#plt.semilogy([], [])
plt.loglog([], [])
plt.ylabel('Density error')
plt.xlabel(r'$\Lambda$')

plt.tight_layout()
plt.savefig('figure_density_convergence.pdf')
plt.show()
