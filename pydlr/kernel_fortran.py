""" 

Python wrapper module for dlrcode 

Author: Hugo U.R. Strand (2021) 

"""

import os
import glob

libname = glob.glob(os.path.dirname(__file__) + '/../lib/libdlr_c.*')[0]

# -- CFFI

from cffi import FFI

ffi = FFI()
ffi.cdef("void c_gridparams(double *lambda, int *p, int *npt, int *npo, int *nt, int *no);")
ffi.cdef("void c_kfine_cc(double *lambda, int *p, int *npt, int *npo, double *t, double *om, double *kmat, double *err);")
ffi.cdef("void c_dlr_rf(double *lambda, double *eps, int *nt, int *no, double *om, double *kmat, int *rank, double *dlrrf, int *oidx);")
ffi.cdef("void c_dlr_it(double *lambda, int *nt, int *no, double *t, double *kmat, int *rank, int *oidx, double* dlrit, int *tidx);")
ffi.cdef("void c_dlr_it2cf(int *rank, double *dlrrf, double *dlrit, double *dlrit2cf, int *it2cfpiv);")
ffi.cdef("void c_dlr_mf(int *nmax, int *rank, double *dlrrf, int *xi, int *dlrmf);")
ffi.cdef("void c_dlr_mf2cf(int *nmax, int *rank, double *dlrrf,int *dlrmf, int *xi, double _Complex *dlrmf2cf, int *mf2cfpiv);")

lib = ffi.dlopen(libname)


import numpy as np

from .pydlr import dlrBase


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


class KernelInterpolativeDecopositionFortran:

    def __init__(self, lamb, eps=1e-15, xi=-1,
                 max_rank=500, nmax=None, verbose=False):

        print('--> Fortran driver')
        
        self.xi = xi
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
        t = ffi.new(f'double [{self.nt}]')
        om = ffi.new(f'double [{self.no}]')
        kmat = ffi.new(f'double [{self.nt*self.no}]')
        err = ffi.new('double [2]')

        lib.c_kfine_cc(lamb, p, npt, npo, t, om, kmat, err)

        self.t = np.frombuffer(ffi.buffer(t), dtype=np.float)
        self.om = np.frombuffer(ffi.buffer(om), dtype=np.float)
        self.kmat = np.frombuffer(ffi.buffer(kmat), dtype=np.float).reshape((self.no, self.nt)).T
        self.err = np.frombuffer(ffi.buffer(err), dtype=np.float)
        
        if verbose:
            print(f't.shape = {self.t.shape}')
            print(f'om.shape = {self.om.shape}')
            print(f'kmat.shape = {self.kmat.shape}')
            print(f'err.shape = {self.err.shape}')

        # -- Select real frequency points

        eps = ffi.new('double *', eps)
        rank = ffi.new('int *', max_rank)
        oidx = ffi.new(f'int [{rank[0]}]')
        dlrrf = ffi.new(f'double [{rank[0]}]')

        lib.c_dlr_rf(lamb, eps, nt, no, om, kmat, rank, dlrrf, oidx)

        self.rank = rank[0]
        self.oidx = np.frombuffer(ffi.buffer(oidx), dtype=np.int32)[:self.rank] - 1
        self.dlrrf = np.frombuffer(ffi.buffer(dlrrf), dtype=np.float)[:self.rank]

        if verbose:
            print(f'rank = {self.rank}')
            print(f'oidx = {self.oidx}')
            print(f'dlrrf = {self.dlrrf}')

        # -- Select imaginary time points

        tidx = ffi.new(f'int [{self.rank}]')
        dlrit = ffi.new(f'double [{self.rank}]')

        lib.c_dlr_it(lamb, nt, no, t, kmat, rank, oidx, dlrit, tidx)

        self.tidx = np.frombuffer(ffi.buffer(tidx), dtype=np.int32) - 1
        self.dlrit = np.frombuffer(ffi.buffer(dlrit), dtype=np.float)
        self.dlrit = (self.dlrit > 0) * self.dlrit + (self.dlrit < 0) * (1 + self.dlrit)

        if verbose:
            print(f'tidx = {self.tidx}')
            print(f'dlrit = {self.dlrit}')
        
        # -- Transform matrix (LU-decomposed)

        it2cfpiv = ffi.new(f'int [{self.rank}]')
        dlrit2cf = ffi.new(f'double [{self.rank**2}]')

        lib.c_dlr_it2cf(rank,dlrrf,dlrit,dlrit2cf,it2cfpiv)

        self.it2cfpiv = np.frombuffer(ffi.buffer(it2cfpiv), dtype=np.int32) - 1
        self.dlrit2cf = np.frombuffer(
            ffi.buffer(dlrit2cf), dtype=np.float).reshape((self.rank, self.rank)).T
        
        if verbose:
            print(f'it2cfpiv = {self.it2cfpiv}')
            #print(f'dlrit2cf = \n{self.dlrit2cf}')

        # -- Matsubara frequency points

        nmax = ffi.new('int *', nmax)
        xi = ffi.new('int *', xi)
        dlrmf = ffi.new(f'int [{self.rank}]')
        lib.c_dlr_mf(nmax,rank,dlrrf,xi,dlrmf)

        self.nmax = nmax[0]
        self.dlrmf = np.frombuffer(ffi.buffer(dlrmf), dtype=np.int32)
            
        if verbose:
            print(f'nmax = {self.nmax}')
            print(f'dlrmf = {self.dlrmf}')

        mf2cfpiv = ffi.new(f'int [{self.rank}]')
        dlrmf2cf = ffi.new(f'double _Complex [{self.rank**2}]')

        lib.c_dlr_mf2cf(nmax,rank,dlrrf,dlrmf,xi,dlrmf2cf,mf2cfpiv)

        self.mf2cfpiv = np.frombuffer(ffi.buffer(mf2cfpiv), dtype=np.int32) - 1
        self.dlrmf2cf = np.frombuffer(ffi.buffer(dlrmf2cf), dtype=np.complex).reshape((self.rank, self.rank)).T
            
        if verbose:
            print(f'mf2cfpiv = {self.mf2cfpiv}')
            #print(f'dlrmf2cf = \n{self.dlrmf2cf}')

        self.T_lx = get_A(self.dlrit2cf, self.it2cfpiv)
        self.T_qx = get_A(self.dlrmf2cf, self.mf2cfpiv)
