"""Compare free Green's function solvers in imaginary time, 
DLR, and Matsubara space.

Copyright 2021 Hugo U.R. Strand

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
or implied. See the License for the specific language governing
permissions and limitations under the License."""


import itertools
import numpy as np

from pydlr import dlr, kernel


def test_free_greens_function_scalar_bosonic(verbose=False):

    xi = 1
    lamb = 30.
    beta = 10.234
    
    H_aa = np.array([[0.3]])
    S_aa = np.array([[1.]])

    N = 100
    w_n = 1.j*np.pi/beta * 2*np.arange(-N+1, N)
    
    # -- DLR
    
    d = dlr(lamb=lamb, python_impl=True, xi=xi)
    
    tau_l = d.get_tau(beta)

    G_laa_free = d.free_greens_function_tau(H_aa, beta, S_aa=S_aa)

    G_xaa = d.dlr_from_tau(G_laa_free)
    G_laa_dlr = d.tau_from_dlr(G_xaa)
    G_laa_eval = d.eval_dlr_tau(G_xaa, tau_l, beta)

    w_q = d.get_matsubara_frequencies(beta)
    G_qaa_dlr = d.matsubara_from_dlr(G_xaa, beta)

    G_qaa_exact = (1./(w_q - H_aa[0,0])).reshape((len(d), 1, 1))
    G_qaa_free = d.free_greens_function_matsubara(H_aa, beta, S_aa=S_aa)

    G_naa_exact = (1./(w_n - H_aa[0,0])).reshape((len(w_n), 1, 1))
    G_naa_eval = d.eval_dlr_freq(G_xaa, w_n, beta)

    G_xaa_matfree = d.dlr_from_matsubara(G_qaa_free, beta)
    G_laa_matfree = d.tau_from_dlr(G_xaa_matfree)

    G_xaa_free = d.free_greens_function_dlr(H_aa, beta)
    G_laa_dlrfree = d.tau_from_dlr(G_xaa_free)
    
    if verbose:

        triqs_ref = False
        
        if triqs_ref:
            # -- TRIQS ref calc
            from triqs.gf import Gf, MeshImFreq, inverse, iOmega_n, make_gf_from_fourier

            mesh = MeshImFreq(beta, 'Boson', N)
            g_n = Gf(mesh=mesh, target_shape=[])
            g_n << inverse(iOmega_n - H_aa[0,0])
            g_t = make_gf_from_fourier(g_n)
            #t_t = np.array([ complex(x) for x in g_t.mesh ])

            w_n_ref = np.array([ complex(x) for x in mesh ])
            np.testing.assert_array_almost_equal(w_n, w_n_ref)

            from triqs.plot.mpl_interface import oplot, oploti, oplotr, plt
        
        sidx = np.argsort(tau_l)
        for val in [tau_l, G_laa_free, G_laa_dlr, G_laa_eval]:
            val = val[sidx]

        # -- Viz

        import matplotlib.pyplot as plt

        plt.figure(figsize=(5, 8))

        subp = [2, 1, 1]    

        plt.subplot(*subp); subp[-1] += 1

        if triqs_ref: oplotr(g_t)

        plt.plot(tau_l, G_laa_free[:, 0, 0], '+', label='free')
        plt.plot(tau_l, G_laa_dlr[:, 0, 0], 'x', label='dlr')
        plt.plot(tau_l, G_laa_matfree[:, 0, 0], '.', label='mat free')
        plt.plot(tau_l, G_laa_dlrfree[:, 0, 0], '.', label='dlr free')
        plt.plot(tau_l, G_laa_eval[:, 0, 0], 'o', alpha=0.5, label='eval')

        plt.legend(loc='upper right')
        plt.xlabel(r'$\tau$')

        subp = [2, 2, 3]    

        plt.subplot(*subp); subp[-1] += 1

        if triqs_ref: oplotr(g_n)
        
        plt.plot(w_q.imag, G_qaa_free[:, 0, 0].real, '+', label='free')
        plt.plot(w_q.imag, G_qaa_dlr[:, 0, 0].real, 'x', label='dlr')
        plt.plot(w_n.imag, G_naa_eval[:, 0, 0].real, '.', label='eval', alpha=0.2)
        plt.legend(loc='upper right')
        plt.xlabel(r'$i\omega_n$')
        plt.ylabel(r'Re[$G(i\omega_n)$]')

        plt.subplot(*subp); subp[-1] += 1

        if triqs_ref: oploti(g_n)
        
        plt.plot(w_q.imag, G_qaa_free[:, 0, 0].imag, '+', label=f'free')        
        plt.plot(w_q.imag, G_qaa_dlr[:, 0, 0].imag, 'x', label='dlr')
        plt.plot(w_n.imag, G_naa_eval[:, 0, 0].imag, '.', label='eval', alpha=0.2)
        plt.legend(loc='upper right')
        plt.xlabel(r'$i\omega_n$')
        plt.ylabel(r'Im[$G(i\omega_n)$]')
        
        plt.tight_layout()
        plt.show()
        
    np.testing.assert_array_almost_equal(G_laa_free, G_laa_dlr)
    np.testing.assert_array_almost_equal(G_laa_free, G_laa_matfree)
    np.testing.assert_array_almost_equal(G_laa_free, G_laa_eval)
    np.testing.assert_array_almost_equal(G_laa_free, G_laa_dlrfree)
    
    np.testing.assert_array_almost_equal(G_qaa_free, G_qaa_exact)
    np.testing.assert_array_almost_equal(G_qaa_free, G_qaa_dlr)
    np.testing.assert_array_almost_equal(G_naa_eval, G_naa_exact)


def test_free_greens_function_scalar(verbose=False):

    lamb = 300.
    #lamb = 30.
    beta = 10.234
    
    H_aa = np.array([[0.3]])
    S_aa = np.array([[1.]])
        
    d = dlr(lamb=lamb, python_impl=True, verbose=True, eps=1e-13)

    tau_l = d.get_tau(beta)
    G_laa = d.free_greens_function_tau(H_aa, beta, S_aa=S_aa)
    G_xaa = d.dlr_from_tau(G_laa)
    G_qaa_ref = d.matsubara_from_dlr(G_xaa, beta)

    tau_i = np.linspace(0, beta, num=400)
    G_iaa = d.eval_dlr_tau(G_xaa, tau_i, beta)

    w_q = d.get_matsubara_frequencies(beta)

    G_qaa = d.free_greens_function_matsubara(H_aa, beta, S_aa=S_aa)
    G_xaa_ref = d.dlr_from_matsubara(G_qaa, beta)
    G_laa_ref = d.tau_from_dlr(G_xaa_ref).real

    print(f'diff tau = {np.max(np.abs(G_laa - G_laa_ref))}')
    print(f'diff Matsubara = {np.max(np.abs(G_qaa - G_qaa_ref))}')
    
    g_xaa = d.free_greens_function_dlr(H_aa, beta)
    g_laa = d.tau_from_dlr(g_xaa)
    
    print(f'diff DLR = {np.max(np.abs(G_laa - g_laa))}')
    
    if verbose:
        
        sidx = np.argsort(tau_l)
        tau_l = tau_l[sidx]
        g_laa = g_laa[sidx]
        G_laa = G_laa[sidx]
        G_laa_ref = G_laa_ref[sidx]
        # -- Viz

        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 9))

        subp = [2, 1, 1]    

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(tau_l, G_laa[:, 0, 0], '+')
        plt.plot(tau_l, G_laa_ref[:, 0, 0], 'x')
        plt.plot(tau_l, g_laa[:, 0, 0], '.')
        plt.plot(tau_i, G_iaa[:, 0, 0], '-')
        plt.xlabel(r'$\tau$')

        plt.subplot(*subp); subp[-1] += 1
        lr = plt.plot(w_q.imag, G_qaa[:, 0, 0].real, '+', label=f'Re')
        li = plt.plot(w_q.imag, G_qaa[:, 0, 0].imag, '+', label=f'Im')
        
        cr, ci = [x[0].get_color() for x in [lr, li]]
        plt.plot(w_q.imag, G_qaa_ref[:, 0, 0].real, 'x', color=cr)
        plt.plot(w_q.imag, G_qaa_ref[:, 0, 0].imag, 'x', color=ci)
        plt.legend(loc='upper right')
        plt.xlabel(r'$i\omega_n$')
        
        plt.tight_layout()
        plt.show()
        
    np.testing.assert_array_almost_equal(G_laa, G_laa_ref)
    np.testing.assert_array_almost_equal(G_qaa, G_qaa_ref)
    np.testing.assert_array_almost_equal(G_laa, g_laa)


def test_free_greens_function_matrix(verbose=False):

    lamb = 40.
    beta = 10.234
    
    H_aa = np.array([
        [1, 0.4, 0.1],
        [0.4, 2, 0.5],
        [0.1, 0.5, 3],
        ])

    S_aa = np.array([
        [1, 0.1, 0.2],
        [0.1, 1, 0.3],
        [0.2, 0.3, 1],
        ])
    
    d = dlr(lamb=lamb)

    tau_l = d.get_tau(beta)
    G_laa = d.free_greens_function_tau(H_aa, beta, S_aa=S_aa)
    G_xaa = d.dlr_from_tau(G_laa)
    G_qaa_ref = d.matsubara_from_dlr(G_xaa, beta)

    tau_i = np.linspace(0, beta, num=400)
    G_iaa = d.eval_dlr_tau(G_xaa, tau_i, beta)

    w_q = d.get_matsubara_frequencies(beta)

    G_qaa = d.free_greens_function_matsubara(H_aa, beta, S_aa=S_aa)
    G_xaa_ref = d.dlr_from_matsubara(G_qaa, beta)
    G_laa_ref = d.tau_from_dlr(G_xaa_ref).real

    print(f'diff tau = {np.max(np.abs(G_laa - G_laa_ref))}')
    print(f'diff Matsubara = {np.max(np.abs(G_qaa - G_qaa_ref))}')

    g_xaa = d.free_greens_function_dlr(H_aa, beta, S_aa=S_aa)
    g_laa = d.tau_from_dlr(g_xaa)
    
    print(f'diff tau DLR = {np.max(np.abs(G_laa - g_laa))}')
        
    if verbose:
        
        sidx = np.argsort(tau_l)
        tau_l = tau_l[sidx]
        g_laa = g_laa[sidx]
        G_laa = G_laa[sidx]
        G_laa_ref = G_laa_ref[sidx]
        # -- Viz

        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 9))

        subp = [3, 3, 1]    

        for i, j in itertools.product(range(3), repeat=2):
            plt.subplot(*subp); subp[-1] += 1
            lr = plt.plot(w_q.imag, G_qaa[:, i, j].real, '+', label=f'Re[G{i}{j}]')
            li = plt.plot(w_q.imag, G_qaa[:, i, j].imag, '+', label=f'Im[G{i}{j}]')
            cr, ci = [x[0].get_color() for x in [lr, li]]
            plt.plot(w_q.imag, G_qaa_ref[:, i, j].real, 'x', color=cr)
            plt.plot(w_q.imag, G_qaa_ref[:, i, j].imag, 'x', color=ci)
            plt.legend(loc='upper right')
            plt.xlabel(r'$i\omega_n$')

        plt.figure(figsize=(8, 9))

        subp = [3, 3, 1]    

        for i, j in itertools.product(range(3), repeat=2):
            plt.subplot(*subp); subp[-1] += 1
            plt.plot(tau_l, G_laa[:, i, j], '+', label=f'G{i}{j}')
            plt.plot(tau_l, G_laa_ref[:, i, j], 'x')
            plt.plot(tau_l, g_laa[:, i, j], '.')
            plt.plot(tau_i, G_iaa[:, i, j], '-')
            plt.legend(loc='upper right')
            plt.xlabel(r'$\tau$')

        plt.tight_layout()
        plt.show()
        
    np.testing.assert_array_almost_equal(G_laa, G_laa_ref)
    np.testing.assert_array_almost_equal(G_qaa, G_qaa_ref)
    np.testing.assert_array_almost_equal(g_laa, G_laa)


if __name__ == '__main__':

    test_free_greens_function_scalar_bosonic(verbose=True)
    test_free_greens_function_scalar(verbose=True)
    test_free_greens_function_matrix(verbose=True)
    
