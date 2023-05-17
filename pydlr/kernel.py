"""Imagniary time analytic continuation kernel routines."""
"""
Copyright 2021 Hugo U.R. Strand

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
or implied. See the License for the specific language governing
permissions and limitations under the License."""


import numpy as np

from scipy.linalg import lu_solve, lu_factor
from scipy.linalg.interpolative import interp_decomp


def chebyshev_collocation_points_1st_kind(N):
    """
    Return the Chebyshev collocation points of the first kind.

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

    >>> from pydlr.kernel import chebyshev_collocation_points_1st_kind
    >>> x_j = chebyshev_collocation_points_1st_kind(10)
    array([-0.98768834, -0.89100652, -0.70710678, -0.4539905 , -0.15643447,
            0.15643447,  0.4539905 ,  0.70710678,  0.89100652,  0.98768834])

    """
    j = np.arange(N)
    x_j = np.cos(np.pi * (2*j + 1)/(2*N))[::-1]
    return x_j


def chebyshev_barycentric_weights_1st_kind(N):
    """
    Return the Chebyshev barycentric interpolation weights of the first kind.

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

    >>> from pydlr.kernel import chebyshev_barycentric_weights_1st_kind
    >>> w_j = chebyshev_barycentric_weights_1st_kind(10)
    array([ 0.15643447, -0.4539905 ,  0.70710678, -0.89100652,  0.98768834,
           -0.98768834,  0.89100652, -0.70710678,  0.4539905 , -0.15643447])

    """
    j = np.arange(N)
    w_j = (-1)**j * np.sin(np.pi * (2*j + 1)/(2*N))
    return w_j


def barycentric_chebyshev_interpolation(x, x_j, f_j, w_j):
    """
    Return the barycentric interpolation of a Chebyshev polynomial on arbitrary points.

    The Barycentric interpolation formula has the form

    .. math:: f(x) = \\left( \\sum_{j=0}^{N-1} \\frac{f_j w_j}{x - x_j} \\right) \Bigg/ \\left( \\sum_{j=0}^{N-1} \\frac{w_j}{x - x_j} \\right)

    Parameters
    ----------

    x : ndarray
        points :math:`x` to evaluate the Chebyshev polynomial at

    x_j : ndarray
        Chebyshev collocation points :math:`x_j` of the first kind

    f_j : ndarray
        Values  :math:`f_j` of the polynomial at the collocation points, :math:`f_j = f(x_j)`

    w_j : ndarray
        Chebyshev barycentric interpolation weights :math:`w_j`


    Returns
    -------

    f_x : ndarray
        Values :math:`f(x)` of the polynomial at the points :math:`x`

    
    Examples
    --------
    >>> import numpy as np
    >>> f = lambda x : np.cos(4*np.pi * x)
    >>> from pydlr.kernel import chebyshev_collocation_points_1st_kind
    >>> from pydlr.kernel import chebyshev_barycentric_weights_1st_kind
    >>> x_j = chebyshev_collocation_points_1st_kind(32)
    >>> w_j = chebyshev_barycentric_weights_1st_kind(32)
    >>> f_j = f(x_j)
    >>> from pydlr.kernel import barycentric_chebyshev_interpolation
    >>> x = np.linspace(-1, 1, num=1000)
    >>> f_x = barycentric_chebyshev_interpolation(x, x_j, f_j, w_j)
    >>> np.testing.assert_array_almost_equal(f_x, f(x))
    >>> print(np.max(np.abs(f_x - f(x))) < 1e-9)
    True

    """
    
    # -- Barycentric interpolation off the grid
    with np.testing.suppress_warnings() as sup:
        sup.filter(RuntimeWarning, "divide by zero encountered in true_divide")
        sup.filter(RuntimeWarning, "invalid value encountered in true_divide")
        
        q_xj = w_j[:, None] / (x[None, :] - x_j[:, None])
        f_x = np.sum(q_xj * f_j[:, None, ...], axis=0) / np.sum(q_xj, axis=0)

    # -- Direct value lookup if x is on the grid x_j
    idxs = np.argwhere(x[:, None] == x_j[None, :])
    if len(idxs) > 0:
        xidx, j = idxs.T
        f_x[xidx] = f_j[j]
    
    return f_x


def fermi_function(E, beta):
    """
    Evaluate the Fermi distribution function at energy :math:`E` and inverse temperature :math:`\beta`.

    The Fermi distribution function :math:`f_\\beta(E)` has the form

    .. math:: f_\\beta(E) = \\frac{1}{1 + e^{\\beta E}}

    the evaluation is stabilized using separate formulas for :math:`E \lessgtr 0`.

    Parameters
    ----------

    E : ndarray
        Energies

    beta : float
        Inverse temperature

    Returns
    -------

    f : ndarray
        The Fermi distribution function evaluated at :math:`E`
    
    """
    
    f = np.zeros_like(E)
    p, m = np.argwhere(E > 0), np.argwhere(E <= 0)

    f[p] = np.exp(-beta*E[p]) / (1. + np.exp(-beta*E[p]))
    f[m] = 1. / (np.exp(beta*E[m]) + 1.)

    return f


def kernel(tau, omega):
    """
    Evaluate the imaginary time and real frequency analytical continuation kernel :math:`K(\\tau, \\omega)`.

    The Fermionic analytical continuation kernel has the form
    
    .. math:: K(\\tau, \\omega) = \\frac{e^{-\\omega \\tau}}{1 + e^{-\\omega}}

    in normalized imaginary time :math:`\\tau \\in [0, 1]` and frequency :math:`\omega`.

    The evaluation is stabilized using separate formulas for :math:`E \lessgtr 0`.

    Parameters
    ----------

    tau : ndarray
        Points in imaginary time :math:`\\tau_j`

    omega : ndarray
        Points in real frequency :math:`\\omega_k`

    Returns
    -------

    kernel : ndarray
        Kernel :math:`K_{jk}` evaluated on the :math:`\\tau_j` and :math:`\\omega_k` grids,
        :math:`K_{jk} = K(\\tau_j, \\omega_k)`

    """

    kernel = np.empty((len(tau), len(omega)))

    p, = np.where(omega > 0.)
    m, = np.where(omega <= 0.)
    w_p, w_m = omega[p].T, omega[m].T

    tau = tau[:, None]

    kernel[:, p] = np.exp(-tau*w_p) / (1 + np.exp(-w_p))
    kernel[:, m] = np.exp((1. - tau)*w_m) / (1 + np.exp(w_m))

    return kernel


def gridparams(lamb, order=24, lambda_scale=1.):
    """
    Empirical discretization parameters of :math:`K(\\tau, \\omega)` for given :math:`\Lambda`

    Parameters
    ----------

    lamb : float
        Cutoff parameter :math:`\Lambda`

    order : int, optional
        Order of the polynomial expansion

    Returns
    -------

    order : int
        polynomial order
    
    npt : int
        number of panel refinements in imaginary time

    npo : int
        number of panel refinements in frequency

    nt : int
        total number of imaginary time points
    
    no : int
        total number of real frenquency points

    """

    lamb *= lambda_scale
    
    npt = int(np.max([np.ceil(np.log(lamb)/np.log(2.))-2, 1]))
    npo = int(np.max([np.ceil(np.log(lamb)/np.log(2.)), 1]))

    nt = 2 * order * npt
    no = 2 * order * npo

    return order, npt, npo, nt, no


def kernel_discretization(lamb, error_est=False):
    """
    Return kernel discretization correct to machine prescision for given :math:`\Lambda`

    Parameters
    ----------

    lamb : float
        Cut-off parameter :math:`\Lambda`

    error_est : bool, optional
        Perform convergence check of discretization (NB, slow)

    Returns
    -------

    kmat : ndarray
        Discretization of the analytical continuation kernel :math:`K_{jk}`

    t : ndarray
        Discretized imaginary time grid :math:`\\tau_j`

    w : ndarray
        Discretized real frequency grid :math:`\\omega_k`

    """
    
    order, npt, npo, nt, no = gridparams(lamb)

    #print(f'order = {order}, npt = {npt}, npo = {npo}, nt = {nt}, no = {no}')
    
    x_i = chebyshev_collocation_points_1st_kind(order)
    w_i = chebyshev_barycentric_weights_1st_kind(order)

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

    if not error_est:
        return kmat, t, w
    else:
        # -- Error estimate

        N = order
        x2_i = chebyshev_collocation_points_1st_kind(2*N)

        err = 0.

        for widx in range(no):
            for tp in range(npt):
                a, b = t_panel_break_pt[tp], t_panel_break_pt[tp + 1]
                X = a + (b - a)*0.5*(x2_i + 1)
                K = np.squeeze(kernel(X, np.array([w[widx]])))
                K_interp = barycentric_chebyshev_interpolation(
                    x2_i, x_i, kmat[N*tp:N*(tp+1), widx], w_i)
                perr = np.max(np.abs(K - K_interp))
                err = np.max([err, perr])

        for tidx in range(nt//2):
            for wp in range(2*npo):
                a, b = w_panel_break_pt[wp], w_panel_break_pt[wp + 1]
                X = a + (b - a)*0.5*(x2_i + 1)
                K = np.squeeze(kernel(np.array([t[tidx]]), X))
                K_interp = barycentric_chebyshev_interpolation(
                    x2_i, x_i, kmat[tidx, N*wp:N*(wp+1)], w_i)
                perr = np.max(np.abs(K - K_interp))
                err = np.max([err, perr])            
            
        return kmat, t, w, err


class KernelInterpolativeDecoposition:

    """
    Interpolative decomposition class for the imaginary-time analytical continuation kernel.

    This is the SciPy based driver class for the Discrete Lehmann Representation (DLR).

    Parameters
    ----------

    lamb : float
        DLR scale parameter :math:`\\Lambda`.
    eps : float
        Set accuracy of the DLR representation.
    xi : sign, optional,
        Statistical sign :math:`\\xi = \\pm 1` for bosons and fermions respectively.
    max_rank : int, optional
        Maximum rank of the DLR kernel decomposition. Default 500.
    nmax : int, optional
        Maxumum index of the Matsubara frequency grid. Default int(lamb).
    verbose : bool, optional
        Default `False`.

    """

    def __init__(self, lamb, eps=1e-15, xi=-1, max_rank=500, nmax=None, verbose=False):

        if verbose:
            import time        
            t_start = time.time()

        self.xi = xi # +1 for bosons, -1 for fermions
        self.lamb = lamb
        self.eps = eps

        if nmax is None: nmax = int(lamb)

        if verbose: t = time.time()
        #self.kmat, self.t, self.om, self.err = kernel_discretization(self.lamb, error_est=True)
        self.kmat, self.t, self.om = kernel_discretization(self.lamb, error_est=False)
        if verbose: print(f'kernel {time.time() - t} s')

        # -- Select real frequency points

        if verbose: t = time.time()

        #id_eps = self.eps * self.lamb
        id_eps = self.eps # Do not scale eps with lamb, to be consistent with Fortran impl.
        if id_eps >= 0.1: id_eps = 0.1
        self.rank, self.oidx, self.proj_w = interp_decomp(self.kmat, id_eps, rand=False)
        self.oidx = np.sort(self.oidx[:self.rank])
        self.dlrrf = self.om[self.oidx]
        if verbose: print(f'ID w {time.time() - t} s')
        
        # -- Select imaginary time points

        if verbose: t = time.time()
        self.tidx, self.proj_t = interp_decomp(self.kmat[:, self.oidx].T, self.rank, rand=False)
        self.tidx = np.sort(self.tidx[:self.rank])

        #self.dlrit = self.t[self.tidx]

        tt = self.t.copy()
        tt += (self.t[::-1] > 0) * (1 - self.t[::-1])
        self.dlrit = tt[self.tidx]
        
        if verbose: print(f'ID t {time.time() - t} s')
            
        # -- Matsubara frequency points

        if verbose: t = time.time()
        n = np.arange(-nmax, nmax+1)
        zeta = 0.5 * (1 - xi) # 0 for bosons, 1 for fermions
        iwn = 1.j * np.pi * (2*n + zeta)
        self.kmat_mf = -1./(iwn[:, None] - self.dlrrf[None, :])
        if verbose: print(f'kernel mats {time.time() - t} s')

        if verbose: t = time.time()
        self.mfidx, self.proj_mf = interp_decomp(self.kmat_mf.T, self.rank, rand=True)
        self.mfidx = np.sort(self.mfidx[:self.rank])
        if verbose: print(f'ID mats {time.time() - t} s')

        self.nmax = nmax
        self.dlrmf = n[self.mfidx]

        # -- Transform matrix DLR-tau (LU-decomposed)

        if verbose: t = time.time()
        self.T_lx = self.kmat[self.tidx][:, self.oidx]
        self.dlrit2cf, self.it2cfpiv = lu_factor(self.T_lx)
        if verbose: print(f'lu ix {time.time() - t} s')

        # -- Transform matrix DLR-Matsubara (LU-decomposed)

        if verbose: t = time.time()
        self.T_qx = self.kmat_mf[self.mfidx, :]
        self.dlrmf2cf, self.mf2cfpiv = lu_factor(self.T_qx)
        if verbose: print(f'lu mats {time.time() - t} s')

        if verbose: print(f'dlr init done {time.time() - t_start} s')
        
