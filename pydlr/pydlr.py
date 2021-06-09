""" 

Python wrapper module for dlrcode 

Author: Hugo U.R. Strand (2021) 

"""


import numpy as np
import numpy.polynomial.legendre as leg

from scipy.linalg import lu_solve
from scipy.linalg import eigh as scipy_eigh 


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

class dlr(object):

    def __init__(self, lamb, eps=1e-15, fb=b'f', max_rank=500, nmax=None, verbose=False):

        self.xi = {b'f':-1., b'b':1.}[fb]
        
        self.lamb = lamb
        self.eps = eps

        if verbose:
            print(f'lambda = {self.lamb}')
            print(f'eps = {self.eps}')

        if nmax is None: nmax = int(lamb)
        
        # -- Determine kernel discretization from heuristics

        lamb = ffi.new('double *', lamb)
        res = [ ffi.new('int [1]') for n in range(5) ]

        lib.c_gridparams(lamb, *res)

        p, npt, npo, nt, no = [ x for x in res ]
        self.p, self.npt, self.npo, self.nt, self.no = [ x[0] for x in res ]

        if verbose:
            print(f'p = {self.p}')
            print(f'npt = {self.npt}')
            print(f'npo = {self.npo}')
            print(f'nt = {self.nt}')
            print(f'no = {self.no}')

        # -- Build analytical continuation kernel

        fb = ffi.new('char *', fb)
        t = ffi.new(f'double [{self.nt}]')
        om = ffi.new(f'double [{self.no}]')
        kmat = ffi.new(f'double [{self.nt*self.no}]')
        err = ffi.new('double [2]')

        lib.c_kfine_cc(fb, lamb, p, npt, npo, t, om, kmat, err)

        self.t = np.frombuffer(ffi.buffer(t), dtype=np.float)
        self.om = np.frombuffer(ffi.buffer(om), dtype=np.float)
        self.kmat = np.frombuffer(ffi.buffer(kmat), dtype=np.float).reshape((self.no, self.nt)).T
        self.err = np.frombuffer(ffi.buffer(err), dtype=np.float)

        if verbose:
            print(f't.shape = {self.t.shape}')
            print(f'om.shape = {self.om.shape}')
            print(f'kmat.shape = {self.kmat.shape}')
            print(f'err.shape = {self.err.shape}')

        #print(f't = {self.t}')
        #print(f'om = {self.om}')

        # -- Select real frequency points

        eps = ffi.new('double *', eps)
        rank = ffi.new('int *', max_rank)
        oidx = ffi.new(f'int [{rank[0]}]')
        dlrrf = ffi.new(f'double [{rank[0]}]')

        lib.c_dlr_rf(lamb, eps, nt, no, om, kmat, rank, dlrrf, oidx)

        self.rank = rank[0]
        self.oidx = np.frombuffer(ffi.buffer(oidx), dtype=np.int32)[:self.rank]
        self.dlrrf = np.frombuffer(ffi.buffer(dlrrf), dtype=np.float)[:self.rank]
        
        if verbose:
            print(f'rank = {self.rank}')
            print(f'oidx = {self.oidx}')
            print(f'dlrrf = {self.dlrrf}')

        # -- Select imaginary time points

        tidx = ffi.new(f'int [{self.rank}]')
        dlrit = ffi.new(f'double [{self.rank}]')

        lib.c_dlr_it(lamb, nt, no, t, kmat, rank, oidx, dlrit, tidx)

        self.tidx = np.frombuffer(ffi.buffer(tidx), dtype=np.int32)
        self.dlrit = np.frombuffer(ffi.buffer(dlrit), dtype=np.float)

        if verbose:
            print(f'tidx = {self.tidx}')
            print(f'dlrit = {self.dlrit}')

        # -- Transform matrix (LU-decomposed)
        
        it2cfpiv = ffi.new(f'int [{self.rank}]')
        dlrit2cf = ffi.new(f'double [{self.rank**2}]')

        lib.c_dlr_it2cf(nt,no,kmat,rank,oidx,tidx,dlrit2cf,it2cfpiv)
        
        self.it2cfpiv = np.frombuffer(ffi.buffer(it2cfpiv), dtype=np.int32)
        self.dlrit2cf = np.frombuffer(ffi.buffer(dlrit2cf), dtype=np.float).reshape((self.rank, self.rank)).T

        if verbose:
            print(f'it2cfpiv = {self.it2cfpiv}')
            #print(f'dlrit2cf = \n{self.dlrit2cf}')

        # -- Matsubara frequency points
        
        nmax = ffi.new('int *', nmax)
        dlrmf = ffi.new(f'int [{self.rank}]')
        lib.c_dlr_mf(fb,nmax,rank,dlrrf,dlrmf)

        self.nmax = nmax[0]
        self.dlrmf = np.frombuffer(ffi.buffer(dlrmf), dtype=np.int32)

        if verbose:
            print(f'nmax = {self.nmax}')
            print(f'dlrmf = {self.dlrmf}')

        mf2cfpiv = ffi.new(f'int [{self.rank}]')
        dlrmf2cf = ffi.new(f'double _Complex [{self.rank**2}]')
        
        lib.c_dlr_mf2cf(fb,nmax,rank,dlrrf,dlrmf,dlrmf2cf,mf2cfpiv)

        self.mf2cfpiv = np.frombuffer(ffi.buffer(mf2cfpiv), dtype=np.int32)
        self.dlrmf2cf = np.frombuffer(ffi.buffer(dlrmf2cf), dtype=np.complex).reshape((self.rank, self.rank)).T

        if verbose:
            print(f'mf2cfpiv = {self.mf2cfpiv}')
            #print(f'dlrmf2cf = \n{self.dlrmf2cf}')

        self.T_lx = get_A(self.dlrit2cf, self.it2cfpiv - 1)
        self.T_qx = get_A(self.dlrmf2cf, self.mf2cfpiv - 1)

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

        G_taa = np.einsum('pab,tp->tab', G_paa, np.exp(-tau*w_p) / (1 + np.exp(-beta*w_p))) + \
                np.einsum('mab,tm->tab', G_maa, np.exp((beta - tau)*w_m) / ( np.exp(beta*w_m) + 1 ))

        return G_taa

    # -- Matsubara Frequency

    def get_matsubara_frequencies(self, beta=1.):
        zeta = (1 - self.xi)/2
        w_q = 1.j * np.pi/beta * (2*self.dlrmf + zeta)
        return w_q


    def matsubara_from_dlr(self, G_xaa, beta=1.):
        #G_qaa = beta * np.tensordot(self.T_qx, G_xaa, axes=(1, 0)).conj()        
        G_qaa = beta * np.tensordot(self.T_qx, G_xaa, axes=(1, 0))
        G_qaa = np.transpose(G_qaa, axes=(0, 2, 1)).conj()
        return G_qaa


    def dlr_from_matsubara(self, G_qaa, beta=1.):
        G_xaa = lu_solve((self.dlrmf2cf, self.mf2cfpiv - 1), G_qaa.conj() / beta)
        return G_xaa
    
    
    def eval_dlr_freq(self, G_xaa, iwn, beta):
        
        w_x = self.om[self.oidx - 1] / beta
        G_qaa = np.einsum('xab,qx->qab', G_xaa, 1./(iwn[:, None] + w_x[None, :]))
        G_qaa = np.transpose(G_qaa, axes=(0, 2, 1)).conj()
        
        return G_qaa
    
    # -- Mathematical operations

    def convolution(self, a_x, b_x):

        tau_l = self.get_tau(1.)
        w_x = self.om[self.oidx - 1]

        k_lx = self.kmat[self.tidx - 1][:, self.oidx - 1]

        I = np.eye(len(w_x))
        W_xx = 1. / (I + w_x[:, None] - w_x[None, :]) - I

        k1_x = -np.squeeze(kernel(np.ones(1), w_x))

        c_x = a_x * b_x * k1_x + \
            b_x * np.sum(a_x[:, None] * W_xx, axis=0) + \
            a_x * np.sum(b_x[:, None] * W_xx, axis=0) + \
            self.dlr_from_tau(tau_l * np.sum(k_lx * (a_x * b_x)[None, :], axis=-1))
        
        return c_x


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

    
    def dyson(self, H_aa, Sigma_xaa, beta):
        raise NotImplementedError
