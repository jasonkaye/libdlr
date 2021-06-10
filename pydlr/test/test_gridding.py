
import numpy as np

from pydlr import chebyschev_collocation_points_1st_kind
from pydlr import chebyschev_barycentric_weights_1st_kind
from pydlr import barycentric_interpolation

from pydlr import dlr, kernel, gridparams


def kernel_discretization(lamb):

    order, npt, npo, nt, no = gridparams(lamb)

    #print(f'order = {order}, npt = {npt}, npo = {npo}, nt = {nt}, no = {no}')
    
    N = 24
    x_i = chebyschev_collocation_points_1st_kind(N)
    w_i = chebyschev_barycentric_weights_1st_kind(N)

    # -- Tau panel discretization
    
    i = np.arange(npt)
    t_panel_break_pt = np.zeros(npt + 1)
    t_panel_break_pt[1:] = 0.5 ** (npt - i)

    t = np.zeros(nt)
    for i in range(npt):
        a, b = t_panel_break_pt[i], t_panel_break_pt[i + 1]
        t[i*order:(i+1)*order] = a + (b - a)*0.5*(x_i+1)

    # -- Frequency panel discretization

    j = np.arange(npo)
    w_panel_break_pt = np.zeros(2*npo + 1)
    w_panel_break_pt[npo+1:] = lamb * 0.5 ** (npo - j - 1)
    w_panel_break_pt[:npo] = - w_panel_break_pt[npo+1:][::-1]

    w = np.zeros(no)
    for i in range(2*npo):
        a, b = w_panel_break_pt[i], w_panel_break_pt[i + 1]
        w[i*order:(i+1)*order] = a + (b - a)*0.5*(x_i+1)    

    kmat = kernel(t[:nt//2], w)
    kmat = np.vstack((kmat, kmat[::-1, ::-1]))

    # -- Error estimate
    
    x2_i = chebyschev_collocation_points_1st_kind(2*N)

    err = 0.

    for widx in range(no):
        for tp in range(npt):
            a, b = t_panel_break_pt[tp], t_panel_break_pt[tp + 1]
            X = a + (b - a)*0.5*(x2_i + 1)
            K = np.squeeze(kernel(X, np.array([w[widx]])))
            K_interp = barycentric_interpolation(x2_i, x_i, kmat[N*tp:N*(tp+1), widx], w_i)
            perr = np.max(np.abs(K - K_interp))
            err = np.max([err, perr])

    for tidx in range(nt//2):
        for wp in range(2*npo):
            a, b = w_panel_break_pt[wp], w_panel_break_pt[wp + 1]
            X = a + (b - a)*0.5*(x2_i + 1)
            K = np.squeeze(kernel(np.array([t[tidx]]), X))
            K_interp = barycentric_interpolation(x2_i, x_i, kmat[tidx, N*wp:N*(wp+1)], w_i)
            perr = np.max(np.abs(K - K_interp))
            err = np.max([err, perr])            
            
    return kmat, t, w, err


def test_gridding(verbose=False):

    lamb = 10.
    d = dlr(lamb=lamb, python_impl=False)
    kmat, t, w, err = kernel_discretization(lamb)
    
    if False:
        order, npt, npo, nt, no = gridparams(lamb)

        print(f'order = {order}, npt = {npt}, npo = {npo}, nt = {nt}, no = {no}')

        N = 24
        x_i = chebyschev_collocation_points_1st_kind(N)
        w_i = chebyschev_barycentric_weights_1st_kind(N)

        # -- Tau panel discretization

        i = np.arange(npt)
        t_panel_break_pt = np.zeros(npt + 1)
        t_panel_break_pt[1:] = 0.5 ** (npt - i)

        t = np.zeros(nt)
        for i in range(npt):
            a, b = t_panel_break_pt[i], t_panel_break_pt[i + 1]
            t[i*order:(i+1)*order] = a + (b - a)*0.5*(x_i+1)

        # -- Frequency panel discretization

        j = np.arange(npo)
        w_panel_break_pt = np.zeros(2*npo + 1)
        w_panel_break_pt[npo+1:] = lamb * 0.5 ** (npo - j - 1)
        w_panel_break_pt[:npo] = - w_panel_break_pt[npo+1:][::-1]

        w = np.zeros(no)
        for i in range(2*npo):
            a, b = w_panel_break_pt[i], w_panel_break_pt[i + 1]
            w[i*order:(i+1)*order] = a + (b - a)*0.5*(x_i+1)    

        kmat = kernel(t[:nt//2], w)
        kmat = np.vstack((kmat, kmat[::-1, ::-1]))
    

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


    if False:
        x2_i = chebyschev_collocation_points_1st_kind(2*N)

        err = 0.

        for widx in range(no):
            for tp in range(npt):
                a, b = t_panel_break_pt[tp], t_panel_break_pt[tp + 1]
                X = a + (b - a)*0.5*(x2_i + 1)
                K = np.squeeze(kernel(X, np.array([w[widx]])))
                K_interp = barycentric_interpolation(x2_i, x_i, kmat[N*tp:N*(tp+1), widx], w_i)
                perr = np.max(np.abs(K - K_interp))
                err = np.max([err, perr])

        for tidx in range(nt//2):
            for wp in range(2*npo):
                a, b = w_panel_break_pt[wp], w_panel_break_pt[wp + 1]
                X = a + (b - a)*0.5*(x2_i + 1)
                K = np.squeeze(kernel(np.array([t[tidx]]), X))
                K_interp = barycentric_interpolation(x2_i, x_i, kmat[tidx, N*wp:N*(wp+1)], w_i)
                perr = np.max(np.abs(K - K_interp))
                err = np.max([err, perr])            
            
    print(f'err = {err:2.2E}')

    assert(err < 1e-14)
        

if __name__ == '__main__':

    test_gridding(verbose=True)
    
