
""" Author: Hugo U.R. Strand (2021) """


import numpy as np

from pydlr.kernel import chebyshev_collocation_points_1st_kind
from pydlr.kernel import chebyshev_barycentric_weights_1st_kind
from pydlr.kernel import barycentric_chebyshev_interpolation


def test_barycentric_interp(verbose=False):

    N = 32
    x_i = chebyshev_collocation_points_1st_kind(N)
    w_i = chebyshev_barycentric_weights_1st_kind(N)

    f = lambda x : np.sin(4 * np.pi * x)
    
    f_i = f(x_i)

    x_j = np.linspace(-1, 1, num=1000)
    f_j = f(x_j)

    f_j_interp = barycentric_chebyshev_interpolation(x_j, x_i, f_i, w_i)

    print(f'diff = {np.max(np.abs(f_j - f_j_interp))}')
    np.testing.assert_array_almost_equal(f_j, f_j_interp)
    
    if verbose:

        import matplotlib.pyplot as plt
        plt.plot(x_i, f_i, '.')
        plt.plot(x_j, f_j, '-')
        plt.plot(x_j, f_j_interp, '-', label='interp')
        plt.legend()
        plt.show()


def test_barycentric_interp_on_grid(verbose=False):

    N = 32
    x_i = chebyshev_collocation_points_1st_kind(N)
    w_i = chebyshev_barycentric_weights_1st_kind(N)

    f = lambda x : np.sin(4 * np.pi * x)
    
    f_i = f(x_i)

    x_j = np.concatenate((x_i[::8]*0.99, [x_i[0], x_i[-1]], x_i[1::4]*0.99))
    f_j = f(x_j)

    f_j_interp = barycentric_chebyshev_interpolation(x_j, x_i, f_i, w_i)

    if verbose:
        print(f'f_j = \n{f_j}')
        print(f'f_j_interp = \n{f_j_interp}')
    
        print(f'diff = {np.max(np.abs(f_j - f_j_interp))}')
    
    np.testing.assert_array_almost_equal(f_j, f_j_interp)
            

if __name__ == '__main__':

    test_barycentric_interp_on_grid(verbose=True)
    test_barycentric_interp(verbose=True)
    
