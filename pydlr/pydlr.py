""" Discrete Lehman Representation (DLR) implementation
using Numpy and Scipy."""
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
import numpy.polynomial.legendre as leg

from scipy.linalg import eigh as scipy_eigh 
from scipy.linalg import lu_solve, lu_factor


from .kernel import kernel, KernelInterpolativeDecoposition


class dlr(object):
    
    """
    Discrete Lehmann Representation (DLR) class.

    Provides the DRL basis, transforms, and algorithms.

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
    python_impl : bool, optional
        Switch between the python and fortran library driver. Default `True`.

    """
    
    def __init__(self, lamb, eps=1e-15, xi=-1,
                 max_rank=500, nmax=None, verbose=False, python_impl=True, dense_imfreq=False):

        self.xi = xi
        self.lamb = lamb

        eps = np.abs(eps)
        min_eps = 10 * np.finfo(float).resolution

        self.eps = eps if eps > min_eps else min_eps

        opts = dict(eps=eps, xi=xi, max_rank=max_rank, nmax=nmax, verbose=verbose)
        
        if not python_impl:
            from .kernel_fortran import KernelInterpolativeDecopositionFortran
            KID = KernelInterpolativeDecopositionFortran
        else:
            KID = KernelInterpolativeDecoposition
            opts['dense_imfreq']= dense_imfreq
            
        kid = KID(lamb, **opts)

        members = [
            'rank', 'dlrit', 'dlrrf', 'dlrmf',
            'dlrit2cf', 'it2cfpiv', 'dlrmf2cf', 'mf2cfpiv', 'T_lx', 'T_qx',
            ]
        
        for member in members: setattr(self, member, getattr(kid, member))
        
        del kid

        # -- Split real-frequency nodes (assuming sorted)
        
        self.pm_idx = np.argwhere(self.dlrrf > 0)[0,0]
        self.dlrrf_p = self.dlrrf[self.pm_idx:]
        self.dlrrf_m = self.dlrrf[:self.pm_idx]
        self.kernel_nominator_p = 1 / (1 + np.exp(-self.dlrrf_p))
        self.kernel_nominator_m = 1 / (1 + np.exp(+self.dlrrf_m))

        # -- Auxilliary variables

        tau_l = self.get_tau(1.)
        w_x = self.dlrrf

        I = np.eye(len(w_x))
        self.W_xx = 1. / (I + w_x[:, None] - w_x[None, :]) - I

        self.k1_x = -np.squeeze(kernel(np.ones(1), w_x))
        self.TtT_xx = self.dlr_from_tau(tau_l[:, None] * self.T_lx)

        self.__bosonic_corr_freq = lambda omega: np.tanh(0.5 * omega)
        self.bosonic_corr_x = self.__bosonic_corr_freq(self.dlrrf)

        self.W_bc_xx = self.W_xx * self.bosonic_corr_x[:, None]
        self.TtT_bc_xx = self.TtT_xx * self.bosonic_corr_x[None, :]
        
        
    def __len__(self):
        """The length of the DLR expansion

        Returns
        -------

        rank : int
            Number of DLR expansion coefficients.
        """
        return self.rank

    
    def __xi_arg(self, xi):
        """Internal helper function to filter xi arguments."""

        if xi is None: xi = self.xi

        assert( np.abs(xi) == 1 and type(xi) == int )

        return xi
        

    def get_dlr_frequencies(self):
        """Get real frequency DLR grid :math:`\\omega_x \\in [-\\Lambda, \\Lambda]`

        Returns
        -------

        w_x : (n), ndarray
            Scaled real frequency points :math:`\\omega_x`
            for the DLR representation.
        """
        
        return self.dlrrf

    
    # -- Imaginary time

    def get_tau_over_beta(self):
        
        """Get grid in scaled imaginary time :math:`\\tau/\\beta \\in [0, 1]`

        Returns
        -------

        ttau_l : (n), ndarray
            Scaled imaginary time points :math:`\\tilde{\\tau}_l = \\tau_l / \\beta` 
            for the DLR representation, :math:`\\tilde{\\tau}_l \\in [0, 1]`.
        """
        
        ttau_l = self.dlrit
        return ttau_l

    
    def get_tau(self, beta):

        """Get :math:`\\tau`-grid in imaginary time :math:`\\tau \\in [0, \\beta]`.

        Parameters
        ----------

        beta : float
            Inverse temperature :math:`\\beta`

        Returns
        -------

        tau_l : (n), ndarray
            Imaginary time points :math:`\\tau_l` for the DLR representation, :math:`\\tau_l \\in [0, \\beta]`.
        """
        
        tau_l = self.get_tau_over_beta() * beta
        return tau_l


    def dlr_from_tau(self, G_laa):

        """Transform the rank-3 array_like Green's function `G_laa` from imaginary time to DLR space.

        Parameters
        ----------

        G_laa : (n,m,m), array_like
            Green's function in imaginary time with :math:`m \\times m` orbital indices.

        Returns
        -------

        G_xaa : (n,m,m), ndarray
            Green's function i DLR coefficient space with :math:`m \\times m` orbital indices.
        """

        G_xaa = lu_solve((self.dlrit2cf, self.it2cfpiv), G_laa)
        return G_xaa


    def tau_from_dlr(self, G_xaa):

        """Transform the rank-3 array_like Green's function `G_xaa` from DLR space to imaginary time.

        Parameters
        ----------

        G_xaa : (n,m,m), array_like
            Green's function in DLR coefficient space with :math:`m \\times m` orbital indices.

        Returns
        -------

        G_laa : (n,m,m), ndarray
            Green's function in imaginary time with :math:`m \\times m` orbital indices.
        """

        G_laa = np.tensordot(self.T_lx, G_xaa, axes=(1, 0))
        return G_laa


    def tau_from_legendre(self, G_naa):

        """Transform the rank-3 array_like Green's function `G_naa` from Legendre coefficient space to imaginary time.

        Parameters
        ----------

        G_naa : (n,m,m), array_like
            Green's function in Legendre coefficient space with :math:`m \\times m` orbital indices.

        Returns
        -------

        G_laa : (n,m,m), ndarray
            Green's function in imaginary time with :math:`m \\times m` orbital indices.
        """

        x = 2 * self.get_tau_over_beta() - 1
        G_laa = np.rollaxis(leg.legval(x, G_naa), -1)
        return G_laa


    def eval_dlr_tau(self, G_xaa, tau_k, beta):

        """Evaluate the DLR coefficient Green's function `G_xaa` at arbibrary points in imaginary time.

        Parameters
        ----------

        G_xaa : (n,m,m), array_like
            Green's function in DLR coefficient space with :math:`m \\times m` orbital indices.

        tau_k : (k), array_like
            Imaginary time points :math:`\\tau_k` where to evaluate the Green's function.

        beta : float
            Inverse temperature :math:`\\beta`

        Returns
        -------

        G_kaa : (k,m,m), ndarray
            Green's function at the imaginary time points :math:`\\tau_k` with :math:`m \\times m` orbital indices.
        """

        assert(len(G_xaa.shape) == 3 or len(G_xaa.shape) == 1)
        assert(len(tau_k.shape) == 1)

        tau_k = tau_k[:, None] / beta
        
        w_p = self.dlrrf_p[None, :]
        K_kp = np.exp(-tau_k*w_p) * self.kernel_nominator_p[None, :]
        G_kaa = np.einsum('kp,p...->k...', K_kp, G_xaa[self.pm_idx:])

        w_m = self.dlrrf_m[None, :]
        K_km = np.exp((1 - tau_k)*w_m) * self.kernel_nominator_m[None, :]
        G_kaa += np.einsum('km,m...->k...', K_km, G_xaa[:self.pm_idx])

        return G_kaa


    def lstsq_dlr_from_tau(self, tau_i, G_iaa, beta):
        """Return DLR coefficients by least squares fit to values on arbitrary imaginary time grid.

        Parameters
        ----------

        tau_i : (i), array_like
            Imaginary time points :math:`\\tau_i` where the Green's function is sampled.

        G_iaa : (i,m,m), array_like
            Green's function in imaginary time space :math:`G(\\tau_i)` with :math:`m \\times m` orbital indices.

        beta : float
            Inverse temperature :math:`\\beta`

        Returns
        -------

        G_xaa : (k,m,m), ndarray
            Green's function in DLR coefficient space with :math:`m \\times m` orbital indices.
        """

        shape_iaa = G_iaa.shape
        assert(len(shape_iaa) == 3)
        assert(len(tau_i.shape) == 1)
        
        shape_iA = (shape_iaa[0], shape_iaa[1]*shape_iaa[2])
        shape_xaa = (len(self), shape_iaa[1], shape_iaa[2])

        K_ix = kernel(tau_i/beta, self.dlrrf)

        G_xaa = np.linalg.lstsq(
            K_ix, G_iaa.reshape(shape_iA),
            rcond=None)[0].reshape(shape_xaa)

        return G_xaa
    
    
    # -- Matsubara Frequency

    def get_matsubara_frequencies(self, beta):

        """Get Matsubara frequency grid.

        Parameters
        ----------

        beta : float
            Inverse temperature :math:`\\beta`

        Returns
        -------

        w_q : (n), ndarray
            Matsubara frequency points :math:`i\\omega_q` for the DLR representation.
        """

        zeta = (1 - self.xi)/2
        w_q = 1.j * np.pi/beta * (2*self.dlrmf + zeta)
        return w_q


    def dlr_from_matsubara(self, G_qaa, beta, xi=None):

        """Transform the rank-3 array_like Green's function `G_qaa` from Matsbuara frequency to DLR space.

        Parameters
        ----------

        G_qaa : (n,m,m), array_like
            Green's function in Matsubara frequency with :math:`m \\times m` orbital indices.
        
        beta : float
            Inverse temperature :math:`\\beta`

        xi : int, optional
            Statistics sign, :math:`\\xi = +1` for bosons and :math:`\\xi = -1` for fermions.

        Returns
        -------

        G_xaa : (n,m,m), ndarray
            Green's function i DLR coefficient space with :math:`m \\times m` orbital indices.
        """
        xi = self.__xi_arg(xi)

        G_xaa = lu_solve((self.dlrmf2cf, self.mf2cfpiv), G_qaa / beta)

        if xi == 1: G_xaa /= self.bosonic_corr_x[:, None, None]

        return G_xaa


    def matsubara_from_dlr(self, G_xaa, beta, xi=None):

        """Transform the rank-3 array_like Green's function `G_xaa` from DLR space to Matsbuara frequency.

        Parameters
        ----------

        G_xaa : (n,m,m), array_like
            Green's function i DLR coefficient space with :math:`m \\times m` orbital indices.

        beta : float
            Inverse temperature :math:`\\beta`

        xi : int, optional
            Statistics sign, :math:`\\xi = +1` for bosons and :math:`\\xi = -1` for fermions.

        Returns
        -------

        G_qaa : (n,m,m), ndarray
            Green's function in Matsubara frequency with :math:`m \\times m` orbital indices.
        """
        xi = self.__xi_arg(xi)

        if xi == 1:
            G_qaa = beta * np.tensordot(
                self.T_qx * self.bosonic_corr_x[None, :], G_xaa, axes=(1, 0))
        else:
            G_qaa = beta * np.tensordot(self.T_qx, G_xaa, axes=(1, 0))

        if len(G_qaa.shape) == 3: G_qaa = np.transpose(G_qaa, axes=(0, 2, 1))

        return G_qaa
    
    
    def eval_dlr_freq(self, G_xaa, z, beta, xi=None):

        """Evaluate the DLR coefficient Green's function `G_xaa` at arbibrary points `z` in frequency space.

        Parameters
        ----------

        G_xaa : (n,m,m), array_like
            Green's function in DLR coefficient space with :math:`m \\times m` orbital indices.

        z : (k), array_like
            Frequency points :math:`z_k` where to evaluate the Green's function.

        beta : float
            Inverse temperature :math:`\\beta`

        xi : int, optional
            Statistics sign, :math:`\\xi = +1` for bosons and :math:`\\xi = -1` for fermions.

        Returns
        -------

        G_zaa : (k,m,m), ndarray
            Green's function at the frequency points :math:`z_k` with :math:`m \\times m` orbital indices.
        """
        xi = self.__xi_arg(xi)

        w_x = self.dlrrf / beta
        kernel_zx = 1./(z[:, None] + w_x[None, :])
        if xi == 1: kernel_zx *= self.bosonic_corr_x[None, :]
        G_zaa = np.einsum('x...,zx->z...', G_xaa, kernel_zx)
        if len(G_zaa.shape) == 3: G_zaa = np.transpose(G_zaa, axes=(0, 2, 1))
        G_zaa = G_zaa.conj()
        
        return G_zaa
    

    def lstsq_dlr_from_matsubara(self, w_q, G_qaa, beta):
        """Return DLR coefficients by least squares fit to values on arbitrary Matsubara frequency grid.

        Parameters
        ----------

        w_q : (q), array_like
            Imaginary frequency points :math:`i\\omega_q` where the Green's function is sampled.

        G_qaa : (q,m,m), array_like
            Green's function in imaginary frequency space :math:`G(i\\omega_q)` with :math:`m \\times m` orbital indices.

        beta : float
            Inverse temperature :math:`\\beta`

        Returns
        -------

        G_xaa : (k,m,m), ndarray
            Green's function in DLR coefficient space with :math:`m \\times m` orbital indices.
        """

        shape_qaa = G_qaa.shape
        assert(len(shape_qaa) == 3)
        
        shape_qA = (shape_qaa[0], shape_qaa[1]*shape_qaa[2])
        shape_xaa = (len(self), shape_qaa[1], shape_qaa[2])

        K_qx = -1./(w_q[:, None] - self.dlrrf[None, :]/beta)

        G_xaa = np.linalg.lstsq(
            K_qx, G_qaa.reshape(shape_qA),
            rcond=None)[0].reshape(shape_xaa)

        return G_xaa


    # -- Mathematical operations

    def quadratic_hamiltonian(self, G_xaa, beta, xi=None):

        """ Get quadratic Hamiltonian contribution of physical Green's function.

        Parameters
        ----------

        G_xaa : (k,m,m), ndarray
            Green's function in DLR coefficient space with :math:`m \\times m` orbital indices.        

        beta : float
            Inverse temperature :math:`\\beta`

        Returns
        -------

        H_aa : (m,m), ndarray
            Quadratic Hamiltonian contribution to the Green's function.

        """

        xi = self.__xi_arg(xi)

        w_x = self.dlrrf

        K0_x = kernel(np.array([0.]), w_x)
        K1_x = kernel(np.array([1.]), w_x)

        D_x = (-K0_x + xi * K1_x) * w_x / beta
        H_aa = np.tensordot(D_x, G_xaa, axes=([-1, 0]))[0]

        return H_aa
    

    def convolution(self, A_xaa, B_xaa, beta, xi=None):

        """ DLR convolution with :math:`\mathcal{O}(N^2)` scaling. Author: Hugo U.R. Strand (2021)

        Imaginary time convolution
        
        .. math:: C = A \\ast B

        reformulated in DLR coefficient space. The notation :math:`C = A \\ast B` is short hand for:

        .. math:: C_{ij}(\\tau) = \\sum_k \\int_{0}^\\beta d\\bar{\\tau} A_{ik}(\\tau - \\bar{\\tau}) B_{kj}(\\bar{\\tau})

        Parameters
        ----------

        A_xaa : (n,m,m)
            Green's function :math:`A` in DLR coefficient space with :math:`m \\times m` orbital indices.

        B_xaa : (n,m,m)
            Green's function :math:`B` in DLR coefficient space with :math:`m \\times m` orbital indices.

        beta : float
            Inverse temperature :math:`\\beta`

        xi : int, optional
            Statistics sign, :math:`\\xi = +1` for bosons and :math:`\\xi = -1` for fermions.

        Returns
        -------

        C_xaa : (n,m,m), ndarray
            Green's function :math:`C` in DLR coefficient space with :math:`m \\times m` orbital indices,
            given by the convolution :math:`C = A \\ast B`.
        """
        xi = self.__xi_arg(xi)
        
        W_xx = self.W_xx if xi == -1 else self.W_bc_xx
        TtT_xx = self.TtT_xx if xi == -1 else self.TtT_bc_xx

        n, na, _ = A_xaa.shape

        WA_xaa = np.matmul(W_xx.T, A_xaa.reshape((n, na*na))).reshape((n, na, na))
        C_xaa = np.matmul(WA_xaa, B_xaa)
        del WA_xaa

        WB_xaa = np.matmul(W_xx.T, B_xaa.reshape((n, na*na))).reshape((n, na, na))
        C_xaa += np.matmul(A_xaa, WB_xaa)
        del WB_xaa
            
        AB_xaa = np.matmul(A_xaa, B_xaa)
        C_xaa += -xi * self.k1_x[:, None, None] * AB_xaa

        C_xaa += np.matmul(TtT_xx, AB_xaa.reshape((n, na*na))).reshape((n, na, na))
        del AB_xaa
        
        C_xaa *= beta
        
        return C_xaa

    
    def convolution_matrix(self, A_xaa, beta, xi=None):

        """ DLR convolution matrix with :math:`\mathcal{O}(N^2)` scaling. Author: Hugo U.R. Strand (2021) 

        The imaginary time convolution matrix :math:`M` is given by
        
        .. math:: M = [A \\ast]

        i.e. the combination of the Green's function :math:`A` and the convolution operator :math:`\\ast`.

        Author: Hugo U.R. Strand (2021) 

        Parameters
        ----------

        A_xaa : (n,m,m)
            Green's function :math:`A` in DLR coefficient space with :math:`m \\times m` orbital indices.

        beta : float
            Inverse temperature :math:`\\beta`

        xi : int, optional
            Statistics sign, :math:`\\xi = +1` for bosons and :math:`\\xi = -1` for fermions.

        Returns
        -------

        M_xaxa : (n,m,n,m), ndarray
            Convolution matrix :math:`[A \\ast]` as a rank-4 tensor in DLR and orbital space.
        """
        xi = self.__xi_arg(xi)
        
        W_xx = self.W_xx if xi == -1 else self.W_bc_xx
        TtT_xx = self.TtT_xx if xi == -1 else self.TtT_bc_xx

        n, na, _ = A_xaa.shape

        Q_xaa = np.einsum('yx,yab->xab', W_xx, A_xaa)
        Q_xaa += -xi * self.k1_x[:,None,None] * A_xaa

        M_xxaa = np.einsum( 'xy,yab->yxab', W_xx,   A_xaa)
        M_xxaa += np.einsum('xy,yab->xyab', TtT_xx, A_xaa)
        M_xxaa += np.einsum('xy,yab->xyab', np.eye(n),   Q_xaa)
        M_xxaa *= beta
        
        M_xaxa = np.moveaxis(M_xxaa, 2, 1)
        
        return M_xaxa


    # -- Free Green's function solvers

    def free_greens_function_dlr(self, H_aa, beta, S_aa=None, xi=None):

        """ Return the free Green's function in DLR coefficent space.

        The free Green's function is the solution to the Dyson differential equation

        .. math:: (-\\partial_\\tau - H_{ij}) G_{jk}(\\tau) = \\delta(\\tau)

        Parameters
        ----------

        H_aa : (m,m), array_like
            Single-particle Hamiltonian matrix in :math:`m \\times m` orbital space.

        beta : float
            Inverse temperature :math:`\\beta`

        S_aa : (m,m), array_like, optional
            Overlap matrix for the generalized case of non-orthogonal basis functions,
            where the Dyson equation takes the form 
            :math:`(-S_{ij} \\partial_\\tau - H_{ij}) G_{jk}(\\tau) = \\delta(\\tau)`.
            Default is `S_aa = None`.

        xi : int, optional
            Statistics sign, :math:`\\xi = +1` for bosons and :math:`\\xi = -1` for fermions.

        Returns
        -------

        G_xaa : (n,m,m), ndarray
            Free Green's function :math:`G` in DLR coefficient space with :math:`m \\times m` orbital indices.

        Notes
        -----

        The fastest algorithm for calculation of free Green's functions is the `free_greens_function_tau` method.
        """
        xi = self.__xi_arg(xi)
        
        na = H_aa.shape[0]
        I_aa = np.eye(na)

        if S_aa is None: S_aa = I_aa

        w_x = self.dlrrf
        n = len(w_x)

        D_lx = self.T_lx * w_x[None, :] / beta
        D_AA = np.kron(D_lx, S_aa) - np.kron(self.T_lx, H_aa)
                
        bc_x = kernel(np.array([0.]), w_x) - xi * kernel(np.array([1.]), w_x)
        D_AA[(n-1)*na:, :] = np.kron(bc_x, S_aa)

        b_Aa = np.zeros((n*na, na))
        b_Aa[(n-1)*na:, :] = -I_aa

        g_xaa = np.linalg.solve(D_AA, b_Aa).reshape((n, na, na))

        return g_xaa
        
        
    def free_greens_function_tau(self, H_aa, beta, S_aa=None, xi=None):

        """ Return the free Green's function in imaginary time.

        The free Green's function is the solution to the Dyson differential equation

        .. math:: (-\\partial_\\tau - H_{ij}) G_{jk}(\\tau) = \\delta(\\tau)

        Parameters
        ----------

        H_aa : (m,m), array_like
            Single-particle Hamiltonian matrix in :math:`m \\times m` orbital space.

        beta : float
            Inverse temperature :math:`\\beta`

        S_aa : (m,m), array_like, optional
            Overlap matrix for the generalized case of non-orthogonal basis functions,
            where the Dyson equation takes the form 
            :math:`(-S_{ij} \\partial_\\tau - H_{ij}) G_{jk}(\\tau) = \\delta(\\tau)`.
            Default is `S_aa = None`.

        xi : int, optional
            Statistics sign, :math:`\\xi = +1` for bosons and :math:`\\xi = -1` for fermions.

        Returns
        -------

        G_laa : (n,m,m), ndarray
            Free Green's function :math:`G` in imaginary time with :math:`m \\times m` orbital indices.
        """
        xi = self.__xi_arg(xi)

        w_x = self.dlrrf

        if S_aa is None:
            E, U = np.linalg.eigh(H_aa)
        else:
            E, U = scipy_eigh(H_aa, S_aa)

        tau_l = self.get_tau(1.)
        g_lE = -kernel(tau_l, E*beta)
        
        if xi == 1: g_lE /= self.__bosonic_corr_freq(E*beta)[None, :]
        
        g_laa = np.einsum('lE,aE,Eb->lab', g_lE, U, U.T.conj())

        return g_laa    


    def free_greens_function_matsubara(self, H_aa, beta, S_aa=None):

        """ Return the free Green's function in Matsubara frequency.

        The free Green's function is the solution to the Dyson equation

        .. math:: G_{ij}(i \\omega_n ) = \\left[ i\\omega_n - H_{ij} \\right]^{-1}

        Parameters
        ----------

        H_aa : (m,m), array_like
            Single-particle Hamiltonian matrix in :math:`m \\times m` orbital space.

        beta : float
            Inverse temperature :math:`\\beta`

        S_aa : (m,m), array_like, optional
            Overlap matrix for the generalized case of non-orthogonal basis functions,
            where the Dyson equation takes the form 
            :math:`= G_{jk}(i\\omega_n) = (-S_{ij} i\\omega_n - H_{ij})^{-1}`.
            Default is `S_aa = None`.

        Returns
        -------

        G_qaa : (n,m,m), ndarray
            Free Green's function :math:`G` in Matsubara frequency with :math:`m \\times m` orbital indices.

        Notes
        -----

        The fastest algorithm for calculation of free Green's functions is the `free_greens_function_tau` method.
        """

        if S_aa is None: S_aa = np.eye(H_aa.shape[0])    
        w_q = self.get_matsubara_frequencies(beta)
        g_qaa = np.linalg.inv(w_q[:, None, None] * S_aa[None, ...] - H_aa[None, ...])

        return g_qaa    

    
    # -- Dyson equation solvers

    def dyson_matsubara(self, H_aa, Sigma_qaa, beta, S_aa=None):

        """ Solve the Dyson equation in Matsubara frequency.

        The Dyson equation gives the Green's function as

        .. math:: G_{ij}(i \\omega_n ) = \\left[ i\\omega_n - H_{ij} - \\Sigma(i\\omega_n) \\right]^{-1}

        Parameters
        ----------

        H_aa : (m,m), array_like
            Single-particle Hamiltonian matrix in :math:`m \\times m` orbital space.

        Sigma_qaa : (n,m,m), ndarray
            Self-energy :math:`\\Sigma` in Matsubara frequency with :math:`m \\times m` orbital indices.

        beta : float
            Inverse temperature :math:`\\beta`

        S_aa : (m,m), array_like, optional
            Overlap matrix for the generalized case of non-orthogonal basis functions. Default is `S_aa = None`.

        Returns
        -------

        G_qaa : (n,m,m), ndarray
            Green's function :math:`G` in Matsubara frequency with :math:`m \\times m` orbital indices.

        Notes
        -----

        The Matsubara frequency Dyson solver is the fastest Dyson solver, 
        albeit not as accurate as the DLR solver `dyson_dlr`.
        """

        if S_aa is None: S_aa = np.eye(H_aa.shape[0])
        w_q = self.get_matsubara_frequencies(beta)
        G_qaa = np.linalg.inv(w_q[:, None, None] * S_aa[None, ...] - H_aa[None, ...] - Sigma_qaa)

        return G_qaa


    def dyson_dlr(self, H_aa, Sigma_xaa, beta, S_aa=None,
                  iterative=False, lomem=False, verbose=False, tol=1e-12):

        """ Solve the Dyson equation in DLR coefficient space.

        The Dyson equation gives the imaginary time Green's function :math:`G_{ij}(\\tau)` as

        .. math:: (-\\partial_\\tau - H_{ij} - \\Sigma_{ij} \\ast \\, ) G_{jk}(\\tau) = \\delta(\\tau)

        Using the free Green's function :math:`g` this can be rewritten as the integral equation

        .. math:: (1 - g \\ast \\Sigma \\ast \\, ) G = g

        which here is solved directly in DLR coefficient space.

        Parameters
        ----------

        H_aa : (m,m), array_like
            Single-particle Hamiltonian matrix in :math:`m \\times m` orbital space.

        Sigma_xaa : (n,m,m), ndarray
            Self-energy :math:`\\Sigma` in DLR coefficient space with :math:`m \\times m` orbital indices.

        beta : float
            Inverse temperature :math:`\\beta`

        S_aa : (m,m), array_like, optional
            Overlap matrix for the generalized case of non-orthogonal basis functions. Default is `S_aa = None`.

        iterative : bool, optional
            Default `False`

        lomem : bool, optional
            Default `False`

        tol : float, optional
            Tolerance for iterative GMRES solver.

        Returns
        -------

        G_xaa : (n,m,m), ndarray
            Green's function :math:`G` in DLR coefficient space with :math:`m \\times m` orbital indices.

        Notes
        -----

        This DLR space Dyson solver is the most accurate algorithm for solving the Dyson equation. 

        By default it uses Lapack's direct solver on matrix problem in combined DLR and orbital space.
        This is fast and accurate for small problems. For larger problems the :math:`\\mathcal{O}(N^2M^2)`
        memory foot-print limits the performance.

        Hence, for large problems the solver should be run with `iterative=True` and `lomem=True` and will
        then use GMRES to solve the linear system using an implicit matrix formulation. 
        This formulation gives a memory foot print that is of the same order as storing a single Green's function.
        However, the requested GMRES tolerance `tol` has to be carefully tuned.
        """

        na = H_aa.shape[0]
        n = Sigma_xaa.shape[0]
        I_aa = np.eye(na)

        g0_iaa = self.free_greens_function_tau(H_aa, beta, S_aa=S_aa)
        g0_xaa = self.dlr_from_tau(g0_iaa)
        
        g0Sigma_xaa = self.convolution(g0_xaa, Sigma_xaa, beta)
        
        if iterative:

            from scipy.sparse.linalg import LinearOperator
            from scipy.sparse.linalg import gmres as scipy_gmres

            if lomem:
                def Amatvec(x_Aa):
                    x_xaa = x_Aa.reshape((n, na, na))
                    y_xaa = x_xaa.copy()
                    y_xaa -= self.convolution(g0Sigma_xaa, x_xaa, beta)
                    return self.tau_from_dlr(y_xaa).reshape((n*na, na))
            else:
                C_AA = self.convolution_matrix(g0Sigma_xaa, beta).reshape((n*na, n*na))
                def Amatvec(x_Aa):
                    y_Aa = x_Aa - C_AA @ x_Aa
                    y_xaa = y_Aa.reshape((n, na, na))
                    return self.tau_from_dlr(y_xaa).reshape((n*na, na))

            I_aa = np.eye(na)
            g0Sigma_qaa = self.matsubara_from_dlr(g0Sigma_xaa, beta)
            M_qaa = np.linalg.inv(I_aa[None, ...] - g0Sigma_qaa)

            def Mmatvec(x_Aa):
                x_iaa = x_Aa.reshape((n, na, na))
                x_qaa = self.matsubara_from_dlr(self.dlr_from_tau(x_iaa), beta)
                y_qaa = M_qaa @ x_qaa
                y_xaa = self.dlr_from_matsubara(y_qaa, beta)
                return y_xaa.reshape((n*na, na))

            D_AA = LinearOperator(matvec=Amatvec, shape=(n * na, n * na), dtype=complex)
            M_AA = LinearOperator(matvec=Mmatvec, shape=(n * na, n * na), dtype=complex)
            b_Aa = g0_iaa.reshape((n*na, na))

            x0_Aa = M_AA.matvec(b_Aa)
            
            g_Aa, info = solver_wrapper(scipy_gmres, verbose=verbose)(D_AA, b_Aa, x0=x0_Aa, M=M_AA, atol=tol)
            g_xaa = g_Aa.reshape((n, na, na))

            if verbose:
                diff = np.max(np.abs(self.tau_from_dlr(x0_Aa.reshape((n, na, na))) - self.tau_from_dlr(g_Aa.reshape((n,na,na)))))
                print(f'Corr on precond: {diff:2.2E}')

        else:

            if lomem: raise NotImplementedError

            #D_AA = np.kron(self.T_lx, I_aa) - self.tau_from_dlr(self.convolution_matrix(g0Sigma_xaa, beta)).reshape((n*na, n*na))        
            #b_Aa = g0_iaa.reshape((n*na, na))

            D_AA = np.kron(np.eye(n), I_aa) - self.convolution_matrix(g0Sigma_xaa, beta).reshape((n*na, n*na))
            b_Aa = g0_xaa.reshape((n*na, na))

            g_xaa = np.linalg.solve(D_AA, b_Aa).reshape((n, na, na))

        return g_xaa
    
    
    def dyson_dlr_integrodiff(self, H_aa, Sigma_xaa, beta, S_aa=None, xi=None):

        """ Solve the Dyson equation in DLR coefficient space.

        The Dyson equation gives the imaginary time Green's function :math:`G_{ij}(\\tau)` as

        .. math:: (-\\partial_\\tau - H_{ij} - \\Sigma_{ij} \\ast \\, ) G_{jk}(\\tau) = \\delta(\\tau)

        which here is solved directly in DLR coefficient space.

        Parameters
        ----------

        H_aa : (m,m), array_like
            Single-particle Hamiltonian matrix in :math:`m \\times m` orbital space.

        Sigma_xaa : (n,m,m), ndarray
            Self-energy :math:`\\Sigma` in DLR coefficient space with :math:`m \\times m` orbital indices.

        beta : float
            Inverse temperature :math:`\\beta`

        S_aa : (m,m), array_like, optional
            Overlap matrix for the generalized case of non-orthogonal basis functions. Default is `S_aa = None`.

        Returns
        -------

        G_xaa : (n,m,m), ndarray
            Green's function :math:`G` in DLR coefficient space with :math:`m \\times m` orbital indices.

        Notes
        -----

        This DLR space Dyson solver should not be used in production, use `dyson_dlr` instead. 
        The differential formulation gives a worse condition number, as compared to the integral 
        equation formulation. The method is kept here for educational purposes.
        """
        xi = self.__xi_arg(xi)

        na = H_aa.shape[0]
        I_aa = np.eye(na)

        if S_aa is None: S_aa = I_aa

        w_x = self.dlrrf
        n = len(w_x)

        D_lx = self.T_lx * w_x[None, :] / beta

        D_AA = -self.tau_from_dlr(self.convolution_matrix(Sigma_xaa, beta)).reshape((n*na, n*na))
        D_AA += np.kron(D_lx, S_aa) - np.kron(self.T_lx, H_aa)
        
        bc_x = kernel(np.array([0.]), w_x) - xi * kernel(np.array([1.]), w_x)
        D_AA[(n-1)*na:, :] = np.kron(bc_x, S_aa)

        b_Aa = np.zeros((n*na, na))
        b_Aa[(n-1)*na:, :] = -I_aa

        g_xaa = np.linalg.solve(D_AA, b_Aa).reshape((n, na, na))

        return g_xaa
    

class solver_wrapper:
    """
    Wrapper for scipy.sparse.linalg iterative linearsystems solvers
    """

    def __init__(self, solver, verbose=False):
        self.solver = solver
        self.verbose = verbose
        self.callback_kwargs = dict(callback=self.callback)
        
        # -- Hack to suppress missing paramter warning 'callback_type'
        # -- in new Scipy GMRES

        self.old_scipy = True
        from inspect import signature
        if 'callback_type' in signature(solver).parameters:
            self.callback_kwargs['callback_type'] = 'legacy'
            self.old_scipy = False

    def kwargs_filter(self, d):
        if self.old_scipy: d['tol'] = d.pop('atol')
        return d
                
    def callback(self, x):
        self.iter += 1

    def __call__(self, A, b, **kwargs):
        self.iter = 1
        ret = self.solver(A, b, **self.callback_kwargs,
                          **self.kwargs_filter(kwargs))
        if self.verbose: print('GMRES N iter:', self.iter)
        return ret
