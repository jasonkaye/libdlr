
""" Author: Hugo U.R. Strand (2021) """


import numpy as np

from pydlr import dlr
from pydlr.kernel import kernel_discretization


def test_gridding(verbose=False):

    lamb = 40.
    d = dlr(lamb=lamb, python_impl=False)
    kmat, t, w, err = kernel_discretization(lamb, error_est=True)
    print(f'err = {err:2.2E}')
    

    if verbose:

        import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 10))

        subp = [4, 1, 1]

        plt.subplot(*subp); subp[-1] += 1
        #plt.plot(t_panel_break_pt, 0.*t_panel_break_pt, 'or')
        plt.plot(d.t, 0.*d.t, '+')
        plt.plot(t, 0.*t, 'x')

        plt.subplot(*subp); subp[-1] += 1
        plt.semilogy(d.t, np.abs(t - d.t), '.-')

        plt.subplot(*subp); subp[-1] += 1
        #plt.plot(w_panel_break_pt, 0.*w_panel_break_pt, 'or')
        plt.plot(d.om, 0.*d.om, '+', alpha=0.5)
        plt.plot(w, 0.*w, 'x')

        plt.subplot(*subp); subp[-1] += 1
        plt.semilogy(d.om, np.abs(w - d.om), '.-')

        plt.show()

        
    np.testing.assert_array_almost_equal(t, d.t)
    np.testing.assert_array_almost_equal(w, d.om)
    np.testing.assert_array_almost_equal(d.kmat, kmat)
    assert(err < 1e-14)
        

if __name__ == '__main__':

    test_gridding(verbose=True)
    
