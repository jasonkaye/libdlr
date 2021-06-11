
import itertools
import numpy as np

from pydlr import dlr, kernel


def free_gf(d, H_aa, beta, S_aa=None):
    

    return g_laa


def test_free_greens_function_scalar(verbose=False):

    lamb = 30.
    beta = 10.234
    
    H_aa = np.array([[0.3]])
    S_aa = np.array([[1.]])
        
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

    test_free_greens_function_scalar(verbose=False)
    test_free_greens_function_matrix(verbose=False)
