
""" Author: Hugo U.R. Strand (2021) """


import numpy as np

from pydlr import dlr


def test_dyson_and_volterra_matsubara():

    beta = 3.421
    d = dlr(lamb=20.)

    e0, e1 = -0.55, 0.3
    V = 0.2

    E_aa = np.array([
        [e0, V],
        [V, e1],
        ])
    
    w_q = d.get_matsubara_frequencies(beta=beta)
 
    g0_qaa = d.free_greens_function_matsubara(np.array([[e0]]), beta)
    g1_qaa = d.free_greens_function_matsubara(np.array([[e1]]), beta)

    G_qaa = d.free_greens_function_matsubara(E_aa, beta)[:, 0, 0].reshape((len(w_q), 1, 1))

    G_qaa_dyson = d.dyson_matsubara(np.array([[e0]]), V**2 * g1_qaa, beta)
    G_qaa_volterra = d.volterra_matsubara(g0_qaa, V**2 * g1_qaa, beta)

    np.testing.assert_array_almost_equal(G_qaa, G_qaa_dyson)
    np.testing.assert_array_almost_equal(G_qaa, G_qaa_volterra)

    g1_xaa = d.dlr_from_matsubara(g1_qaa, beta)
    G_xaa_dyson_dlr = d.dyson_dlr(np.array([[e0]]), V**2 * g1_xaa, beta)
    G_qaa_dyson_dlr = d.matsubara_from_dlr(G_xaa_dyson_dlr, beta)

    np.testing.assert_array_almost_equal(G_qaa, G_qaa_dyson_dlr)
    

if __name__ == '__main__':

    test_dyson_and_volterra_matsubara()
