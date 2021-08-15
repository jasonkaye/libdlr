
""" Author: Hugo U.R. Strand (2021) """


import numpy as np


from scipy.linalg import qr as scipy_qr


def chebyschev_collocation_points_1st_kind(N):
    """
    Return the Chebyschev collocation points of the first kind.

    The Chebyshev collocation points of the first kind are defined as

    .. math:: x_j = \\cos \\left( \\pi \\frac{2j + 1}{2N} \\right)

    where :math:`N` is the order and :math:`j = 0, ..., N-1`.


    Parameters
    ----------

    N : int
        Order of the collocation point set.


    Returns
    -------

    x_j : ndarray
        Array with the collocation points :math:`x_j` sorted in
        increasing order.

    Examples
    --------

    >>> from pydlr.kernel import chebyschev_collocation_points_1st_kind
    >>> x_j = chebyschev_collocation_points_1st_kind(10)
    array([-0.98768834, -0.89100652, -0.70710678, -0.4539905 , -0.15643447,
            0.15643447,  0.4539905 ,  0.70710678,  0.89100652,  0.98768834])

    """
    j = np.arange(N)
    x_j = np.cos(np.pi * (2*j + 1)/(2*N))[::-1]
    return x_j


def chebyschev_barycentric_weights_1st_kind(N):
    """
    Return the Chebyschev barycentric interpolation weights of the first kind.

    The barycentric interpolation weights are defined as

    .. math:: w_j = (-1)^{j} \\sin \\left( \\pi \\frac{2j + 1}{2N} \\right)

    where :math:`N` is the order and :math:`j = 0, ..., N-1`.
    
    Parameters
    ----------

    N : int
        Order of the collocation point set.


    Returns
    -------

    w_j : ndarray
        Array with the barycentric weights :math:`w_j`.

    Examples
    --------

    >>> from pydlr.kernel import chebyschev_barycentric_weights_1st_kind
    >>> x_j = chebyschev_barycentric_weights_1st_kind(10)

    """
    j = np.arange(N)
    w_i = (-1)**j * np.sin(np.pi * (2*j + 1)/(2*N))
    return w_i


def barycentric_interpolation(x, x_i, f_i, w_i):

    # -- Return value if x is on the grid x_i
    idxs = np.argwhere(x_i == x)
    if len(idxs) > 0: return f_i[idxs[0]]

    # -- Barycentric interpolation off the grid
    q_xi = w_i[:, None] / (x[None, :] - x_i[:, None])
    val_x = np.sum(q_xi * f_i[:, None, ...], axis=0) / np.sum(q_xi, axis=0)

    return val_x


def fermi_function(E, beta):

    f = np.zeros_like(E)
    p, m = np.argwhere(E > 0), np.argwhere(E <= 0)

    f[p] = np.exp(-beta*E[p]) / (1. + np.exp(-beta*E[p]))
    f[m] = 1. / (np.exp(beta*E[m]) + 1.)

    return f


def kernel(tau, omega):

    kernel = np.empty((len(tau), len(omega)))

    p, = np.where(omega > 0.)
    m, = np.where(omega <= 0.)
    w_p, w_m = omega[p].T, omega[m].T

    tau = tau[:, None]

    kernel[:, p] = np.exp(-tau*w_p) / (1 + np.exp(-w_p))
    kernel[:, m] = np.exp((1. - tau)*w_m) / (1 + np.exp(w_m))

    return kernel


def gridparams(lamb, order=24):

    npt = int(np.max([np.ceil(np.log(lamb)/np.log(2.))-2, 1]))
    npo = int(np.max([np.ceil(np.log(lamb)/np.log(2.)), 1]))

    nt = 2 * order * npt
    no = 2 * order * npo

    return order, npt, npo, nt, no


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
