""" Solving the SYK model using the DLR expansion

The non-linear problem is solved using both forward iteration
and a hybrid-Newton method.

Author: Hugo U.R. Strand (2021) """


import numpy as np
from scipy.optimize import root

from pydlr import dlr


def sigma_q_syk(g_q, U, d, beta):

    tau_l = d.get_tau(beta)
    tau_l_rev = beta - tau_l
    
    g_x = d.dlr_from_matsubara(g_q, beta)

    g_l = d.tau_from_dlr(g_x)
    g_l_rev = d.eval_dlr_tau(g_x, tau_l_rev, beta)
    
    sigma_l = U**2 * g_l**2 * g_l_rev

    sigma_q = d.matsubara_from_dlr(d.dlr_from_tau(sigma_l))

    return sigma_q


def solve_syk_fwd_iter(d, mu, beta=1., U=0.1, g0_l=None, max_iter=1000, tol=1e-12):

    print('='*72)
    print('SYK Forward iteration DLR solver')
    print('-'*72)
    print(f'mu = {mu}, U = {U}, beta = {beta}')
    print(f'lamb = {lamb}, n_dlr = {len(d.get_tau())}')
    print(f'max_iter = {max_iter}, tol = {tol}')
    print('='*72)
        
    if g0_l is not None: g_q = d.matsubara_from_dlr(d.dlr_from_tau(g0_l))
    else: g_q = d.free_greens_function_matsubara(np.array([[mu]]), beta)

    for iter in range(max_iter):
        sigma_q = sigma_q_syk(g_q, U, d, beta)
        g_q_old = g_q
        g_q = d.dyson_matsubara(np.array([[mu]]), sigma_q, beta)

        diff = np.max(np.abs(g_q - g_q_old))
        print(f'iter = {iter}, diff = {diff}')
        if diff < tol: break

    g_l = d.tau_from_dlr(d.dlr_from_matsubara(g_q, beta))
        
    return g_l


def solve_syk_root(d, mu, beta=1., U=0.1, g0_l=None, tol=1e-16, verbose=False):

    if verbose:
        print('='*72)
        print('SYK root DLR solver')
        print('-'*72)
        print(f'mu = {mu}, U = {U}, beta = {beta}')
        print(f'lamb = {lamb}, n_dlr = {len(d.get_tau())}')
        print(f'tol = {tol}')
        print('='*72)
        
    if g0_l is not None:
        g_l = g0_l[:, 0, 0]
    else:
        g_q = d.free_greens_function_matsubara(np.array([[mu]]), beta)
        g_l = d.tau_from_dlr(d.dlr_from_matsubara(g_q, beta)).real[:, 0, 0]
    
    def target_function(g_l):

        g_l = g_l.reshape((len(g_l), 1, 1))
        g_q = d.matsubara_from_dlr(d.dlr_from_tau(g_l), beta)        
        sigma_q = sigma_q_syk(g_q, U, d, beta)
        g_q_new = d.dyson_matsubara(np.array([[mu]]), sigma_q, beta)
        g_l_new = d.tau_from_dlr(d.dlr_from_matsubara(g_q_new, beta))

        return np.squeeze((g_l - g_l_new).real)

    sol = root(target_function, g_l, method='hybr', tol=tol)

    diff = np.max(np.abs(target_function(sol.x)))
    if verbose: print(f'nfev = {sol.nfev}, diff = {diff}')
    
    g_l = sol.x.reshape((len(g_l), 1, 1))

    return g_l

    
if __name__ == '__main__':

    U = 2.5
    mu0 = 0.
    beta = 10.0
    lamb = beta * 10

    d = dlr(lamb=lamb)
    tau_l = d.get_tau(beta)
    
    g_l_fwd = solve_syk_fwd_iter(d, mu0, beta=beta, U=U)
    g_l_root = solve_syk_root(d, mu0, beta=beta, U=U, verbose=True)

    mu_vec = np.linspace(-1, 1, num=100)
    density_vec = np.zeros_like(mu_vec)

    for idx, mu in enumerate(mu_vec):
        g_l = solve_syk_root(d, mu, beta=beta, U=U)
        g_x = d.dlr_from_tau(g_l)
        n = 1 + np.squeeze(d.eval_dlr_tau(g_x, np.array([beta]), beta).real)
        density_vec[idx] = n
        print(f'mu = {mu:+2.2E}, n = {n:+2.2E}')
    
    
    import matplotlib.pyplot as plt

    plt.figure(figsize=(6, 8))
    subp = [3, 1, 1]

    plt.subplot(*subp); subp[-1] += 1
    plt.plot(mu_vec, density_vec, '-')
    plt.xlabel(r'Chemical potential $\mu$')
    plt.ylabel(r'Density $n$')

    plt.subplot(*subp); subp[-1] += 1
    plt.title(r'SYK $\mu = ' + f'{mu0}' + r'$, $U = ' + f'{U}' + r'$, $\beta = ' + f'{beta}' + '$') 
    plt.plot(tau_l, g_l_fwd[:, 0, 0].real, '+', label='fwd iter')
    plt.plot(tau_l, g_l_root[:, 0, 0].real, 'x', label='root')
    plt.ylim(bottom=-1.1, top=0)
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$G(\tau)$')
    plt.legend(loc='best')

    plt.subplot(*subp); subp[-1] += 1
    plt.semilogy(tau_l, np.abs(g_l_fwd - g_l_root)[:, 0, 0].real, '+')
    plt.ylabel('Solver difference')
    plt.xlabel(r'$\tau$')
    
    
    plt.tight_layout()
    
    plt.show()
    
