
""" Author: Hugo U.R. Strand (2021) """


import numpy as np

from pydlr import dlr


def test_convolution_scalar(verbose=False):

    beta = 3.337
    e1, e2 = 3., 3.3

    d = dlr(lamb=30.)

    g1_q = d.free_greens_function_matsubara(np.array([[e1]]), beta)
    g2_q = d.free_greens_function_matsubara(np.array([[e2]]), beta)

    # -- DLR Matsubara convolution

    gg_q = np.einsum('qab,qbc->qac', g1_q, g2_q)

    gg_x = d.dlr_from_matsubara(gg_q, beta)
    gg_l = d.tau_from_dlr(gg_x)

    # -- DLR coeff convolution

    g1_x = d.dlr_from_matsubara(g1_q, beta)
    g2_x = d.dlr_from_matsubara(g2_q, beta)

    gg_x_dlr = d.convolution(g1_x, g2_x, beta=beta)
    gg_l_dlr = d.tau_from_dlr(gg_x_dlr)

    # -- DLR convolution matrix

    n, na, _ = g1_x.shape
    
    C_xaxa = d.convolution_matrix(g1_x, beta=beta)
    C_AA = C_xaxa.reshape((n*na, n*na))
    
    B_Aa = g2_x.reshape((n*na, na))
    gg_x_mat = np.matmul(C_AA, B_Aa).reshape((n, na, na))
    gg_l_mat = d.tau_from_dlr(gg_x_mat)

    print(f'diff scalar = {np.max(np.abs(gg_l - gg_l_dlr))}')
    print(f'diff scalar (convmat) = {np.max(np.abs(gg_l - gg_l_mat))}')
    
    # --
    
    if verbose:

        # -- Viz

        tau_l = d.get_tau(beta)

        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 6))
        subp = [1, 1, 1]

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(tau_l, gg_l_dlr[:, 0, 0].real, 'o', label='DLR conv', alpha=0.5)
        plt.plot(tau_l, gg_l_mat[:, 0, 0].real, '+', label='DLR conv mat')
        plt.plot(tau_l, gg_l[:, 0, 0].real, 'x', label='Matsub conv')
        plt.ylabel(r'$(g_1 \ast g_2)(\tau)$')
        plt.xlabel(r'$\tau$')
        plt.legend(loc='best')
        plt.show()

    # -- Test

    np.testing.assert_array_almost_equal(gg_l, gg_l_dlr)
    np.testing.assert_array_almost_equal(gg_l, gg_l_mat)


def test_convolution_matrix():

    beta = 5.337

    E1_aa = np.array([
        [-1., 0.2, 0.4 + 1.j],
        [0.2, 0, 0.1j],
        [0.4 - 1.j, -0.1j, 1],
        ])

    np.testing.assert_array_almost_equal(E1_aa, E1_aa.conj().T)

    E2_aa = np.array([
        [-3., 0.3, 0.1 + 0.2j],
        [0.3, 1, 0.3j],
        [0.1 - 0.2j, -0.3j, 0],
        ])

    np.testing.assert_array_almost_equal(E2_aa, E2_aa.conj().T)
    
    d = dlr(lamb=40.)
    w_q = d.get_matsubara_frequencies(beta=beta)

    g1_qaa = d.free_greens_function_matsubara(E1_aa, beta)
    g2_qaa = d.free_greens_function_matsubara(E2_aa, beta)

    # -- DLR Matsubara convolution

    gg_qaa = np.einsum('qab,qbc->qac', g1_qaa, g2_qaa)

    gg_xaa = d.dlr_from_matsubara(gg_qaa, beta)
    gg_laa = d.tau_from_dlr(gg_xaa)

    # -- DLR coeff convolution
    
    g1_xaa = d.dlr_from_matsubara(g1_qaa, beta)
    g2_xaa = d.dlr_from_matsubara(g2_qaa, beta)

    gg_xaa_ref = d.convolution(g1_xaa, g2_xaa, beta)
    gg_laa_ref = d.tau_from_dlr(gg_xaa_ref)

    # -- DLR convolution matrix

    n, na, _ = g1_xaa.shape
    
    C_xaxa = d.convolution_matrix(g1_xaa, beta=beta)
    C_AA = C_xaxa.reshape((n*na, n*na))
    
    B_Aa = g2_xaa.reshape((n*na, na))
    gg_xaa_mat = np.matmul(C_AA, B_Aa).reshape((n, na, na))
    gg_laa_mat = d.tau_from_dlr(gg_xaa_mat)


    # -- Test

    print(f'diff matrix = {np.max(np.abs(gg_laa - gg_laa_ref))}')
    print(f'diff matrix (convmat) = {np.max(np.abs(gg_laa - gg_laa_mat))}')

    np.testing.assert_array_almost_equal(gg_laa, gg_laa_ref)
    np.testing.assert_array_almost_equal(gg_laa, gg_laa_mat)
    

if __name__ == '__main__':

    test_convolution_scalar(verbose=True)
    test_convolution_matrix()
