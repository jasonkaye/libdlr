
import numpy as np

from pydlr import dlr


def test_dlr_convolution(verbose=False):

    beta = 2.

    e1, e2 = 3., 3.3

    d = dlr(lamb=10)
    w_q = d.get_matsubara_frequencies(beta)

    g1_q = 1./ ( w_q - e1 )
    g2_q = 1./ ( w_q - e2 )

    g1_x = d.dlr_from_matsubara(g1_q, beta=beta)
    g2_x = d.dlr_from_matsubara(g2_q, beta=beta)

    g1_l = d.tau_from_dlr(g1_x)
    g2_l = d.tau_from_dlr(g2_x)

    # -- Analytic product in imaginary time
    
    gg_l_anal = g1_l/(e1 - e2) + g2_l/(e2 - e1)
    
    # -- DLR Matsubara convolution

    gg_q = g1_q * g2_q

    gg_x_mat = d.dlr_from_matsubara(gg_q, beta=beta)
    gg_l_mat = d.tau_from_dlr(gg_x_mat)

    # -- DLR coeff convolution

    g1_x = d.dlr_from_matsubara(g1_q, beta=beta)
    g2_x = d.dlr_from_matsubara(g2_q, beta=beta)

    gg_x = d.convolution(g1_x, g2_x, beta=beta)
    gg_l = d.tau_from_dlr(gg_x)

    # -- Viz

    if verbose:
        tau_l = d.get_tau(beta)
        import matplotlib.pyplot as plt
        plt.plot(tau_l, gg_l, '+', label='gg (DLR conv.)')
        plt.plot(tau_l, gg_l_mat, 'x', label='gg (Matsubara)')
        plt.plot(tau_l, gg_l_anal, '.', label='gg (analytic)')
        plt.legend()
        plt.show()

    # -- Test

    np.testing.assert_array_almost_equal(gg_l, gg_l_anal, decimal=9)
    np.testing.assert_array_almost_equal(gg_l, gg_l_mat, decimal=9)


if __name__ == '__main__':

    test_dlr_convolution(verbose=True)
