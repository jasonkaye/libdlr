
""" Author: Hugo U.R. Strand (2021) """


import numpy as np

from pydlr.kernel import chebyschev_collocation_points_1st_kind
from pydlr.kernel import chebyschev_barycentric_weights_1st_kind
from pydlr.kernel import barycentric_interpolation


def test_barycentric_interp(verbose=False):

    N = 32
    x_i = chebyschev_collocation_points_1st_kind(N)
    w_i = chebyschev_barycentric_weights_1st_kind(N)

    f = lambda x : np.sin(4 * np.pi * x)
    
    f_i = f(x_i)

    x_j = np.linspace(-1, 1, num=400)
    f_j = f(x_j)

    f_j_interp = barycentric_interpolation(x_j, x_i, f_i, w_i)

    print(f'diff = {np.max(np.abs(f_j - f_j_interp))}')
    np.testing.assert_array_almost_equal(f_j, f_j_interp)
    
    if verbose:

        import matplotlib.pyplot as plt
        plt.plot(x_i, f_i, '.')
        plt.plot(x_j, f_j, '-')
        plt.plot(x_j, f_j_interp, '-', label='interp')
        plt.legend()
        plt.show()


if __name__ == '__main__':

    test_barycentric_interp(verbose=True)
    
