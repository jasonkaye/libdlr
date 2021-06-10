
""" Author: Hugo U.R. Strand (2021) """


import numpy as np

from pydlr import dlr


def test_convolution_scalar():

    beta = 1.337
    e1, e2 = 3., 3.3

    d = dlr(lamb=10)

    g1_q = d.free_greens_function_matsubara(np.array([[e1]]), beta)
    g2_q = d.free_greens_function_matsubara(np.array([[e2]]), beta)

    # -- DLR Matsubara convolution

    gg_q = np.einsum('qab,qbc->qac', g1_q, g2_q)

    gg_x = d.dlr_from_matsubara(gg_q)
    gg_l = d.tau_from_dlr(gg_x)

    # -- DLR coeff convolution

    g1_x = d.dlr_from_matsubara(g1_q)
    g2_x = d.dlr_from_matsubara(g2_q)

    gg_x_ref = d.convolution(g1_x, g2_x)
    gg_l_ref = d.tau_from_dlr(gg_x_ref)

    # -- Test

    print(f'diff scalar = {np.max(np.abs(gg_l - gg_l_ref))}')
    np.testing.assert_array_almost_equal(gg_l, gg_l_ref, decimal=9)


def test_convolution_matrix():

    beta = 1.337

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
    
    d = dlr(lamb=10)
    w_q = d.get_matsubara_frequencies(beta=beta)

    g1_qaa = d.free_greens_function_matsubara(E1_aa, beta)
    g2_qaa = d.free_greens_function_matsubara(E2_aa, beta)

    # -- DLR Matsubara convolution

    gg_qaa = np.einsum('qab,qbc->qac', g1_qaa, g2_qaa)

    gg_xaa = d.dlr_from_matsubara(gg_qaa)
    gg_laa = d.tau_from_dlr(gg_xaa)

    # -- DLR coeff convolution
    
    g1_xaa = d.dlr_from_matsubara(g1_qaa)
    g2_xaa = d.dlr_from_matsubara(g2_qaa)

    gg_xaa_ref = d.convolution(g1_xaa, g2_xaa)
    gg_laa_ref = d.tau_from_dlr(gg_xaa_ref)

    # -- Test

    print(f'diff matrix = {np.max(np.abs(gg_laa - gg_laa_ref))}')
    np.testing.assert_array_almost_equal(gg_laa, gg_laa_ref, decimal=9)
    

if __name__ == '__main__':

    test_convolution_scalar()
    test_convolution_matrix()
