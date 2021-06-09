
import numpy as np

from pydlr import dlr


def test_convolution():

    beta = 1.337
    e1, e2 = 3., 3.3

    d = dlr(lamb=10)
    w_q = d.get_matsubara_frequencies(beta=beta)

    g1_q = 1./ ( w_q - e1 )
    g2_q = 1./ ( w_q - e2 )

    # -- DLR Matsubara convolution

    gg_q = g1_q * g2_q

    gg_x = d.dlr_from_matsubara(gg_q)
    gg_l = d.tau_from_dlr(gg_x)

    # -- DLR coeff convolution

    g1_x = d.dlr_from_matsubara(g1_q)
    g2_x = d.dlr_from_matsubara(g2_q)

    gg_x_ref = d.convolution(g1_x, g2_x)
    gg_l_ref = d.tau_from_dlr(gg_x_ref)

    # -- Test

    np.testing.assert_array_almost_equal(gg_l, gg_l_ref, decimal=9)


if __name__ == '__main__':

    test_convolution()
