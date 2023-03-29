"""Python wrapper module for dlrcode 

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


import os
import glob

#libname = glob.glob(os.path.dirname(__file__) + '/../lib/libdlr_c.*')[0]
libname = '/Users/hugstr/apps/libdlr/lib/libdlr_c.dylib'

# -- CFFI

from cffi import FFI

ffi = FFI()
ffi.cdef("void c_ccfine_init(double *lambda, int *p, int *npt, int *npo, int *nt, int *no);")
ffi.cdef("void c_ccfine(double *lambda, int *p, int *npt, int *npo, double *t, double *om);")
ffi.cdef("void c_dlr_kfine(double *lambda, int *p, int *npt, int *npo, double *t, double *om, double *kmat, double *err);")
#ffi.cdef("void c_dlr_rf(double *lambda, double *eps, int *nt, int *no, double *om, double *kmat, int *rank, double *dlrrf, int *oidx);")
#ffi.cdef("void c_dlr_it(double *lambda, int *nt, int *no, double *t, double *kmat, int *rank, int *oidx, double* dlrit, int *tidx);")
ffi.cdef("void c_dlr_it_build(double *lambda, double *eps, int *r, double *dlrrf, double *dlrit);")
ffi.cdef("void c_dlr_cf2it_init(int *rank, double *dlrrf, double *dlrit, double *cf2it);")
ffi.cdef("void c_dlr_it2cf_init(int *rank, double *dlrrf, double *dlrit, double *dlrit2cf, int *it2cfpiv);")
ffi.cdef("void c_dlr_mf(int *nmax, int *rank, double *dlrrf, int *xi, int *dlrmf);")
ffi.cdef("void c_dlr_cf2mf_init(int *rank, double *dlrrf,int *dlrmf, int *xi, double _Complex *cf2mf);")
ffi.cdef("void c_dlr_mf2cf_init(int *nmax, int *rank, double *dlrrf,int *dlrmf, int *xi, double _Complex *dlrmf2cf, int *mf2cfpiv);")

lib = ffi.dlopen(libname)


import numpy as np


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

        if verbose: print('--> Fortran driver')
        
        self.xi = xi
        self.lamb = lamb
        self.eps = eps

        if verbose:
            print(f'xi = {self.xi}')
            print(f'lambda = {self.lamb}')
            print(f'eps = {self.eps}')

        if nmax is None: nmax = int(lamb)

        
        r_r = ffi.new('int [1]', [max_rank])
        r_eps = ffi.new('double [1]', [eps])
        r_lamb = ffi.new('double [1]', [lamb])
        
        dlrit = ffi.new(f'double [{max_rank}]')
        dlrrf = ffi.new(f'double [{max_rank}]')

        lib.c_dlr_it_build(r_lamb, r_eps, r_r, dlrrf, dlrit)

        rank = r_r
        self.rank = int(r_r[0])

        self.dlrrf = np.frombuffer(ffi.buffer(dlrrf), dtype=float)[:self.rank].copy()
        self.dlrit = np.frombuffer(ffi.buffer(dlrit), dtype=float)[:self.rank].copy()
        del dlrrf, dlrit
        
        # -- Sort real frequencies

        sidx = np.argsort(self.dlrrf)
        self.dlrrf = self.dlrrf[sidx]
        dlrrf = ffi.cast('double *', self.dlrrf.ctypes.data)

        if verbose:
            print(f'rank = {self.rank}')
            print(f'dlrrf = {self.dlrrf}')
            print(f'dlrit = {self.dlrit}')

        # -- Sort imaginary times

        dlrit_pos = (self.dlrit > 0) * self.dlrit + (self.dlrit < 0) * (1 + self.dlrit)
        sidx = np.argsort(dlrit_pos)
        self.dlrit = self.dlrit[sidx]

        self.dlrit_org = self.dlrit.copy()
        dlrit = ffi.cast('double *', self.dlrit_org.ctypes.data)

        # -- Make times positive

        self.dlrit = (self.dlrit > 0) * self.dlrit + (self.dlrit < 0) * (1 + self.dlrit)

        if verbose:
            print(f'dlrit (trans) = {self.dlrit}')
        
        # -- Transform matrix (LU-decomposed)

        it2cfpiv = ffi.new(f'int [{self.rank}]')
        dlrit2cf = ffi.new(f'double [{self.rank**2}]')

        lib.c_dlr_it2cf_init(rank,dlrrf,dlrit,dlrit2cf,it2cfpiv)

        self.it2cfpiv = np.frombuffer(ffi.buffer(it2cfpiv), dtype=np.int32) - 1
        self.dlrit2cf = np.frombuffer(
            ffi.buffer(dlrit2cf), dtype=float).reshape((self.rank, self.rank)).T
        
        if verbose:
            print(f'it2cfpiv = {self.it2cfpiv}')
            #print(f'dlrit2cf = \n{self.dlrit2cf}')

        cf2it = ffi.new(f'double [{self.rank**2}]')

        lib.c_dlr_cf2it_init(rank,dlrrf,dlrit,cf2it)

        self.cf2it = np.frombuffer(
            ffi.buffer(cf2it), dtype=float).reshape((self.rank, self.rank)).T

        # -- Store raw ffi objects

        self.r_r = rank
        self.r_dlrrf = dlrrf
        self.r_dlrit = dlrit
        self.r_it2cf = dlrit2cf
        self.r_it2cfp = it2cfpiv
            
        # -- Matsubara frequency points

        if nmax < self.rank: nmax = self.rank
        
        if verbose:
            print(f'nmax = {nmax}')
        
        nmax = ffi.new('int *', nmax)
        xi = ffi.new('int *', int(xi))
        dlrmf = ffi.new(f'int [{self.rank}]')

        lib.c_dlr_mf(nmax,rank,dlrrf,xi,dlrmf)

        self.nmax = nmax[0]
        self.dlrmf = np.frombuffer(ffi.buffer(dlrmf), dtype=np.int32)
            
        if verbose:
            print(f'nmax = {self.nmax}')
            print(f'dlrmf = {self.dlrmf}')

        mf2cfpiv = ffi.new(f'int [{self.rank}]')
        dlrmf2cf = ffi.new(f'double _Complex [{self.rank**2}]')

        lib.c_dlr_mf2cf_init(nmax,rank,dlrrf,dlrmf,xi,dlrmf2cf,mf2cfpiv)

        self.mf2cfpiv = np.frombuffer(ffi.buffer(mf2cfpiv), dtype=np.int32) - 1
        self.dlrmf2cf = np.frombuffer(ffi.buffer(dlrmf2cf), dtype=np.complex).reshape((self.rank, self.rank)).T
            
        if verbose:
            print(f'mf2cfpiv = {self.mf2cfpiv}')
            #print(f'dlrmf2cf = \n{self.dlrmf2cf}')

        cf2mf = ffi.new(f'double _Complex [{self.rank**2}]')

        lib.c_dlr_cf2mf_init(rank,dlrrf,dlrmf,xi,cf2mf)

        self.cf2mf = np.frombuffer(ffi.buffer(cf2mf), dtype=np.complex).reshape((self.rank, self.rank)).T
            
        #self.T_lx = get_A(self.dlrit2cf, self.it2cfpiv)
        #self.T_qx = get_A(self.dlrmf2cf, self.mf2cfpiv)

        self.T_lx = self.cf2it
        self.T_qx = self.cf2mf
        
