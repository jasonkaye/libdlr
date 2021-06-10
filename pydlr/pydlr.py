""" 

Python wrapper module for dlrcode 

Author: Hugo U.R. Strand (2021) 

"""


import numpy as np
import numpy.polynomial.legendre as leg


from scipy.linalg import lu_solve, lu_factor
from scipy.linalg import eigh as scipy_eigh 
from scipy.linalg import qr as scipy_qr


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


def dlr_decomp(kmat, eps, lamb, eps_rank=None):

    if eps_rank is None:
        eps_rank = np.linalg.matrix_rank(kmat, tol=eps * lamb) 

    _, P_o = scipy_qr(kmat, pivoting=True, mode='r')
    P_o = P_o[:eps_rank]

    _, P_t = scipy_qr(kmat[:, P_o].T, pivoting=True, mode='r')
    P_t = P_t[:eps_rank]

    return P_o, P_t, eps_rank


def chebyschev_collocation_points_1st_kind(N):
    j = np.arange(N)
    x = np.cos(np.pi * (2*j + 1)/(2*N))[::-1]
    return x


def chebyschev_barycentric_weights_1st_kind(N):
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
    

def get_P(piv):
    """ Permutation matrix corresponding to Lapack piv index vector """
    P = np.eye(len(piv))
    for i, p in enumerate(piv):
        a = P[:, i].copy()
        b = P[:, p].copy()
        P[:, i], P[:, p] = b, a
    return P


def get_idx(piv):
    """ Numpy index vector corresponding to Lapack piv index vector """
    P = get_P(piv)
    idx = np.zeros_like(piv)
    for i in range(len(piv)):
        idx[i] = np.argwhere(P==1)

    return idx


def get_A(lu, piv):
    L, U = np.tril(lu, k=-1) + np.eye(lu.shape[0]), np.triu(lu)
    P = get_P(piv)
    A = P @ L @ U
    return A

import os
import glob
libname = glob.glob(os.path.dirname(__file__) + '/../lib/libdlr_c.*')[0]

# -- CFFI

from cffi import FFI

ffi = FFI()
ffi.cdef("void c_gridparams(double *lambda, int *p, int *npt, int *npo, int *nt, int *no);")
ffi.cdef("void c_kfine_cc(char *fb, double *lambda, int *p, int *npt, int *npo, double *t, double *om, double *kmat, double *err);")
ffi.cdef("void c_dlr_rf(double *lambda, double *eps, int *nt, int *no, double *om, double *kmat, int *rank, double *dlrrf, int *oidx);")
ffi.cdef("void c_dlr_it(double *lambda, int *nt, int *no, double *t, double *kmat, int *rank, int *oidx, double* dlrit, int *tidx);")
ffi.cdef("void c_dlr_it2cf(int *nt, int *no, double *kmat, int *rank, int *oidx, int *tidx, double *dlrit2cf, int *it2cfpiv);")
ffi.cdef("void c_dlr_mf(char *fb, int *nmax, int *rank, double *dlrrf, int *dlrmf);")
ffi.cdef("void c_dlr_mf2cf(char *fb, int *nmax, int *rank, double *dlrrf,int *dlrmf, double _Complex *dlrmf2cf, int *mf2cfpiv);")

lib = ffi.dlopen(libname)

class dlrBase(object):

    def __init__(self, lamb, eps=1e-15, fb=b'f',
                 max_rank=500, nmax=None, verbose=False, python_impl=False):

        self.xi = {b'f':-1., b'b':1.}[fb]
        
        self.lamb = lamb
        self.eps = eps

        if verbose:
            print(f'lambda = {self.lamb}')
            print(f'eps = {self.eps}')

        if nmax is None: nmax = int(lamb)
        
        # -- Determine kernel discretization from heuristics

        if False:
            lamb = ffi.new('double *', lamb)
            res = [ ffi.new('int [1]') for n in range(5) ]

            lib.c_gridparams(lamb, *res)

            p, npt, npo, nt, no = [ x for x in res ]
            self.p, self.npt, self.npo, self.nt, self.no = [ x[0] for x in res ]

            if python_impl:
                np.testing.assert_array_almost_equal(gridparams(self.lamb), [ x[0] for x in res ])


        self.p, self.npt, self.npo, self.nt, self.no = gridparams(self.lamb)
        
        if verbose:
            print(f'p = {self.p}')
            print(f'npt = {self.npt}')
            print(f'npo = {self.npo}')
            print(f'nt = {self.nt}')
            print(f'no = {self.no}')

        # -- Build analytical continuation kernel

        self.kmat, self.t, self.om, self.err = kernel_discretization(self.lamb)  

        if False:
            fb = ffi.new('char *', fb)
            t = ffi.new(f'double [{self.nt}]')
            om = ffi.new(f'double [{self.no}]')
            kmat = ffi.new(f'double [{self.nt*self.no}]')
            err = ffi.new('double [2]')

            lib.c_kfine_cc(fb, lamb, p, npt, npo, t, om, kmat, err)

            self.t_ref = np.frombuffer(ffi.buffer(t), dtype=np.float)
            self.om_ref = np.frombuffer(ffi.buffer(om), dtype=np.float)
            self.kmat_ref = np.frombuffer(ffi.buffer(kmat), dtype=np.float).reshape((self.no, self.nt)).T
            self.err_ref = np.frombuffer(ffi.buffer(err), dtype=np.float)

            np.testing.assert_array_almost_equal(self.t, self.t_ref)
            np.testing.assert_array_almost_equal(self.om, self.om_ref)
            np.testing.assert_array_almost_equal(self.kmat, self.kmat_ref)
        
        if verbose:
            print(f't.shape = {self.t.shape}')
            print(f'om.shape = {self.om.shape}')
            print(f'kmat.shape = {self.kmat.shape}')
            print(f'err.shape = {self.err.shape}')

        #print(f't = {self.t}')
        #print(f'om = {self.om}')

        # -- Select real frequency points

        self.oidx, self.tidx, self.rank = dlr_decomp(self.kmat, self.eps, self.lamb)

        self.dlrrf = self.om[self.oidx]
        self.dlrit = self.t[self.tidx]
        
        self.oidx += 1
        self.tidx += 1

        if False:
            eps = ffi.new('double *', eps)
            rank = ffi.new('int *', max_rank)
            oidx = ffi.new(f'int [{rank[0]}]')
            dlrrf = ffi.new(f'double [{rank[0]}]')

            lib.c_dlr_rf(lamb, eps, nt, no, om, kmat, rank, dlrrf, oidx)

            self.rank_ref = rank[0]
            self.oidx_ref = np.frombuffer(ffi.buffer(oidx), dtype=np.int32)[:self.rank]
            self.dlrrf_ref = np.frombuffer(ffi.buffer(dlrrf), dtype=np.float)[:self.rank]

            rank[0] = self.rank
            self.oidx_ref[:] = self.oidx
            self.dlrrf_ref[:] = self.dlrrf

        if verbose:
            print(f'rank = {self.rank}')
            print(f'oidx = {self.oidx}')
            print(f'dlrrf = {self.dlrrf}')

        # -- Select imaginary time points

        if False:
            tidx = ffi.new(f'int [{self.rank}]')
            dlrit = ffi.new(f'double [{self.rank}]')

            lib.c_dlr_it(lamb, nt, no, t, kmat, rank, oidx, dlrit, tidx)

            self.tidx_ref = np.frombuffer(ffi.buffer(tidx), dtype=np.int32)
            self.dlrit_ref = np.frombuffer(ffi.buffer(dlrit), dtype=np.float)

            self.tidx_ref[:] = self.tidx
            self.dlrit_ref[:] = self.dlrit

        if verbose:
            print(f'tidx = {self.tidx}')
            print(f'dlrit = {self.dlrit}')
        
        # -- Transform matrix (LU-decomposed)

        self.dlrit2cf, self.it2cfpiv = lu_factor(self.kmat[self.tidx - 1][:, self.oidx - 1])
        self.it2cfpiv += 1 # one based Fortran indexing

        if False:
            it2cfpiv = ffi.new(f'int [{self.rank}]')
            dlrit2cf = ffi.new(f'double [{self.rank**2}]')

            lib.c_dlr_it2cf(nt,no,kmat,rank,oidx,tidx,dlrit2cf,it2cfpiv)

            self.it2cfpiv_ref = np.frombuffer(ffi.buffer(it2cfpiv), dtype=np.int32)
            self.dlrit2cf_ref = np.frombuffer(
                ffi.buffer(dlrit2cf), dtype=np.float).reshape((self.rank, self.rank)).T

            A1 = get_A(self.dlrit2cf, self.it2cfpiv - 1)
            A2 = get_A(self.dlrit2cf_ref, self.it2cfpiv_ref - 1)
            np.testing.assert_array_almost_equal(A1, A2)
        
        if verbose:
            print(f'it2cfpiv = {self.it2cfpiv}')
            #print(f'dlrit2cf = \n{self.dlrit2cf}')

        # -- Matsubara frequency points

        n = np.arange(-nmax, nmax+1)
        iwn = 1.j * np.pi * (2*n + 1)

        self.kmat_mf = 1./(iwn[:, None] + self.dlrrf)

        _, P_mf = scipy_qr(self.kmat_mf.T, pivoting=True, mode='r')
        P_mf = P_mf[:self.rank]

        self.nmax = nmax
        self.dlrmf = n[P_mf]
        
        if False:
            nmax = ffi.new('int *', nmax)
            dlrmf = ffi.new(f'int [{self.rank}]')
            lib.c_dlr_mf(fb,nmax,rank,dlrrf,dlrmf)

            self.nmax_ref = nmax[0]
            self.dlrmf_ref = np.frombuffer(ffi.buffer(dlrmf), dtype=np.int32)

            np.testing.assert_array_almost_equal(np.sort(self.dlrmf), np.sort(self.dlrmf_ref))

            self.dlrmf_ref[:] = self.dlrmf
            
        if verbose:
            print(f'nmax = {self.nmax}')
            print(f'dlrmf = {self.dlrmf}')


        self.dlrmf2cf, self.mf2cfpiv = lu_factor(self.kmat_mf[P_mf, :])
        self.mf2cfpiv += 1 # one based Fortran indexing

        if False:
            mf2cfpiv = ffi.new(f'int [{self.rank}]')
            dlrmf2cf = ffi.new(f'double _Complex [{self.rank**2}]')
        
            lib.c_dlr_mf2cf(fb,nmax,rank,dlrrf,dlrmf,dlrmf2cf,mf2cfpiv)

            self.mf2cfpiv_ref = np.frombuffer(ffi.buffer(mf2cfpiv), dtype=np.int32)
            self.dlrmf2cf_ref = np.frombuffer(ffi.buffer(dlrmf2cf), dtype=np.complex).reshape((self.rank, self.rank)).T

            A1 = get_A(self.dlrmf2cf, self.mf2cfpiv - 1)
            A2 = get_A(self.dlrmf2cf_ref, self.mf2cfpiv_ref - 1)
            np.testing.assert_array_almost_equal(A1, A2)
            
        if verbose:
            print(f'mf2cfpiv = {self.mf2cfpiv}')
            #print(f'dlrmf2cf = \n{self.dlrmf2cf}')

        self.T_lx = self.kmat[self.tidx - 1][:, self.oidx - 1]
        self.T_qx = self.kmat_mf[P_mf, :]

        if False:
            T_lx_ref = get_A(self.dlrit2cf, self.it2cfpiv - 1)
            T_qx_ref = get_A(self.dlrmf2cf, self.mf2cfpiv - 1)
        
            np.testing.assert_array_almost_equal(self.T_lx, T_lx_ref)
            np.testing.assert_array_almost_equal(self.T_qx, T_qx_ref)
        

    # -- Imaginary time

    def get_tau_over_beta(self):
        tt = self.t.copy()
        tt += (self.t[::-1] > 0) * (1 - self.t[::-1])
        return tt[self.tidx - 1]

    
    def get_tau(self, beta=1.):
        tau_l = self.get_tau_over_beta() * beta
        return tau_l


    def dlr_from_tau(self, G_laa):
        G_xaa = lu_solve((self.dlrit2cf, self.it2cfpiv - 1), G_laa)
        return G_xaa


    def tau_from_dlr(self, G_xaa):
        G_laa = np.tensordot(self.T_lx, G_xaa, axes=(1, 0))
        return G_laa


    def tau_from_legendre(self, G_naa):
        x = 2 * self.get_tau_over_beta() - 1
        G_laa = np.rollaxis(leg.legval(x, G_naa), -1)
        return G_laa


    def eval_dlr_tau(self, G_xaa, tau, beta):
        assert( self.xi == -1. ) # Not implemented for bosons, yet
        
        w_x = self.om[self.oidx - 1] / beta
        
        p = np.argwhere(w_x > 0.)
        m = np.argwhere(w_x <= 0.)
        
        w_p, G_paa = w_x[p].T, G_xaa[p][:, 0]
        w_m, G_maa = w_x[m].T, G_xaa[m][:, 0]

        tau = tau[:, None]

        G_taa = np.einsum('p...,tp->t...', G_paa, np.exp(-tau*w_p) / (1 + np.exp(-beta*w_p))) + \
                np.einsum('m...,tm->t...', G_maa, np.exp((beta - tau)*w_m) / ( np.exp(beta*w_m) + 1 ))

        return G_taa

    # -- Matsubara Frequency

    def get_matsubara_frequencies(self, beta=1.):
        zeta = (1 - self.xi)/2
        w_q = 1.j * np.pi/beta * (2*self.dlrmf + zeta)
        return w_q


    def matsubara_from_dlr(self, G_xaa, beta=1.):
        G_qaa = beta * np.tensordot(self.T_qx, G_xaa, axes=(1, 0))
        if len(G_qaa.shape) == 3: G_qaa = np.transpose(G_qaa, axes=(0, 2, 1))
        G_qaa = G_qaa.conj()
        return G_qaa


    def dlr_from_matsubara(self, G_qaa, beta=1.):
        G_xaa = lu_solve((self.dlrmf2cf, self.mf2cfpiv - 1), G_qaa.conj() / beta)
        return G_xaa
    
    
    def eval_dlr_freq(self, G_xaa, iwn, beta):
        
        w_x = self.om[self.oidx - 1] / beta
        G_qaa = np.einsum('x...,qx->q...', G_xaa, 1./(iwn[:, None] + w_x[None, :]))
        if len(G_qaa.shape) == 3: G_qaa = np.transpose(G_qaa, axes=(0, 2, 1))
        G_qaa = G_qaa.conj()
        
        return G_qaa
    
    # -- Mathematical operations

    def convolution(self, A_xaa, B_xaa, beta=1.):

        tau_l = self.get_tau(1.)
        w_x = self.om[self.oidx - 1]

        I = np.eye(len(w_x))
        W_xx = 1. / (I + w_x[:, None] - w_x[None, :]) - I

        k1_x = -np.squeeze(kernel(np.ones(1), w_x))

        AB_xaa = np.matmul(A_xaa, B_xaa)

        C_xaa = k1_x[:, None, None] * AB_xaa + \
            self.dlr_from_tau(tau_l[:, None, None] * self.tau_from_dlr(AB_xaa)) + \
            np.matmul(np.tensordot(W_xx, A_xaa, axes=(0, 0)), B_xaa) + \
            np.matmul(A_xaa, np.tensordot(W_xx, B_xaa, axes=(0, 0)))        

        C_xaa *= beta
        
        return C_xaa
    

    def free_greens_function_tau(self, H_aa, beta, S_aa=None):

        w_x = self.om[self.oidx - 1]

        if S_aa is None:
            E, U = np.linalg.eigh(H_aa)
        else:
            E, U = scipy_eigh(H_aa, S_aa)

        tau_l = self.get_tau(1.)
        g_lE = -kernel(tau_l, E*beta)
        g_laa = np.einsum('lE,aE,Eb->lab', g_lE, U, U.T.conj())

        return g_laa    


    def free_greens_function_matsubara(self, H_aa, beta, S_aa=None):

        if S_aa is None: S_aa = np.eye(H_aa.shape[0])    
        w_q = self.get_matsubara_frequencies(beta)
        g_qaa = np.linalg.inv(w_q[:, None, None] * S_aa[None, ...] - H_aa[None, ...])

        return g_qaa    

    
    def dyson_matsubara(self, H_aa, Sigma_qaa, beta, S_aa=None):

        if S_aa is None: S_aa = np.eye(H_aa.shape[0])
        w_q = self.get_matsubara_frequencies(beta)
        G_qaa = np.linalg.inv(w_q[:, None, None] * S_aa[None, ...] - H_aa[None, ...] - Sigma_qaa)

        return G_qaa

    
    def volterra_matsubara(self, g_qaa, Sigma_qaa, beta):

        N = Sigma_qaa.shape[-1]    
        A_qaa = np.eye(N) - np.matmul(g_qaa, Sigma_qaa)        
        G_qaa = np.linalg.solve(A_qaa, g_qaa)

        return G_qaa
    

class dlrPython(dlrBase):

    def __init__(self, lamb, eps=1e-15, fb=b'f',
                 max_rank=500, nmax=None, verbose=False, python_impl=False):

        self.xi = {b'f':-1., b'b':1.}[fb]
        
        self.lamb = lamb
        self.eps = eps

        if nmax is None: nmax = int(lamb)
        
        self.p, self.npt, self.npo, self.nt, self.no = gridparams(self.lamb)
        self.kmat, self.t, self.om, self.err = kernel_discretization(self.lamb)  
        
        # -- Select real frequency points
        # -- Select imaginary time points

        self.oidx, self.tidx, self.rank = dlr_decomp(self.kmat, self.eps, self.lamb)

        self.dlrrf = self.om[self.oidx]
        self.dlrit = self.t[self.tidx]
        
        self.oidx += 1
        self.tidx += 1
        
        # -- Transform matrix (LU-decomposed)

        self.dlrit2cf, self.it2cfpiv = lu_factor(self.kmat[self.tidx - 1][:, self.oidx - 1])
        self.it2cfpiv += 1 # one based Fortran indexing

        # -- Matsubara frequency points

        n = np.arange(-nmax, nmax+1)
        iwn = 1.j * np.pi * (2*n + 1)

        self.kmat_mf = 1./(iwn[:, None] + self.dlrrf)

        _, P_mf = scipy_qr(self.kmat_mf.T, pivoting=True, mode='r')
        P_mf = P_mf[:self.rank]

        self.nmax = nmax
        self.dlrmf = n[P_mf]
                    
        # -- Transform matrix (LU-decomposed)

        self.dlrmf2cf, self.mf2cfpiv = lu_factor(self.kmat_mf[P_mf, :])
        self.mf2cfpiv += 1 # one based Fortran indexing
            
        self.T_lx = self.kmat[self.tidx - 1][:, self.oidx - 1]
        self.T_qx = self.kmat_mf[P_mf, :]

        
#from pydlr_fortran import dlrFortran
#dlr = dlrFortran
dlr = dlrPython
