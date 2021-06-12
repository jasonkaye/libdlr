
""" Discrete Lehman Representation (DLR) implementation
using Numpy and Scipy.

Author: Hugo U.R. Strand (2021) """


import numpy as np
import numpy.polynomial.legendre as leg


from scipy.linalg import qr as scipy_qr
from scipy.linalg import eigh as scipy_eigh 
from scipy.linalg import lu_solve, lu_factor


from .kernel import kernel, kernel_discretization, fermi_function


class dlrBase(object):

    def __init__(self, lamb, eps=1e-15, fb=b'f',
                 max_rank=500, nmax=None, verbose=False, python_impl=False):

        self.xi = {b'f':-1., b'b':1.}[fb]
        
        self.lamb = lamb
        self.eps = eps

        if nmax is None: nmax = int(lamb)
        
        self.kmat, self.t, self.om, self.err = kernel_discretization(self.lamb)  
        
        self.rank = np.linalg.matrix_rank(self.kmat, tol=eps * lamb) 

        # -- Select real frequency points

        _, oidx = scipy_qr(self.kmat, pivoting=True, mode='r')
        self.oidx = oidx[:self.rank]
        self.dlrrf = self.om[self.oidx]

        # -- Select imaginary time points

        _, tidx = scipy_qr(self.kmat[:, self.oidx].T, pivoting=True, mode='r')
        self.tidx = tidx[:self.rank]
        self.dlrit = self.t[self.tidx]
        
        # -- Transform matrix (LU-decomposed)

        self.T_lx = self.kmat[self.tidx][:, self.oidx]
        self.dlrit2cf, self.it2cfpiv = lu_factor(self.T_lx)

        # -- Matsubara frequency points

        n = np.arange(-nmax, nmax+1)
        iwn = 1.j * np.pi * (2*n + 1)

        self.kmat_mf = 1./(iwn[:, None] + self.dlrrf[None, :])

        _, mfidx = scipy_qr(self.kmat_mf.T, pivoting=True, mode='r')
        self.mfidx = mfidx[:self.rank]

        self.nmax = nmax
        self.dlrmf = n[self.mfidx]
                    
        # -- Transform matrix (LU-decomposed)

        self.T_qx = self.kmat_mf[self.mfidx, :]
        self.dlrmf2cf, self.mf2cfpiv = lu_factor(self.T_qx)


    # -- Imaginary time

    def get_tau_over_beta(self):
        tt = self.t.copy()
        tt += (self.t[::-1] > 0) * (1 - self.t[::-1])
        return tt[self.tidx]

    
    def get_tau(self, beta=1.):
        tau_l = self.get_tau_over_beta() * beta
        return tau_l


    def dlr_from_tau(self, G_laa):
        G_xaa = lu_solve((self.dlrit2cf, self.it2cfpiv), G_laa)
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
        
        w_x = self.dlrrf / beta
        
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
        G_xaa = lu_solve((self.dlrmf2cf, self.mf2cfpiv), G_qaa.conj() / beta)
        return G_xaa
    
    
    def eval_dlr_freq(self, G_xaa, iwn, beta):
        
        w_x = self.dlrrf / beta
        G_qaa = np.einsum('x...,qx->q...', G_xaa, 1./(iwn[:, None] + w_x[None, :]))
        if len(G_qaa.shape) == 3: G_qaa = np.transpose(G_qaa, axes=(0, 2, 1))
        G_qaa = G_qaa.conj()
        
        return G_qaa
    
    # -- Mathematical operations

    def convolution(self, A_xaa, B_xaa, beta=1.):

        tau_l = self.get_tau(1.)
        w_x = self.dlrrf

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


    def convolution_matrix(self, A_xaa, beta=1.):

        w_x = self.dlrrf
        tau_l = self.get_tau(1.)

        n = len(w_x)        
        I = np.eye(n)
        
        W_xx = 1. / (I + w_x[:, None] - w_x[None, :]) - I

        k1_x = -np.squeeze(kernel(np.ones(1), w_x))

        C_xxaa = A_xaa[:, None, ...] * W_xx.T[:, :, None, None] + \
            self.dlr_from_tau(np.einsum('l,lx,x...->lx...', tau_l, self.T_lx, A_xaa))
        C_xaa = k1_x[:, None, None] * A_xaa + np.tensordot(W_xx, A_xaa, axes=(0, 0))
        for i in range(n): C_xxaa[i, i, :, :] += C_xaa[i]    
        C_xxaa *= beta

        C_xaxa = np.moveaxis(C_xxaa, 2, 1)
        
        return C_xaxa
        

    def free_greens_function_dlr(self, H_aa, beta, S_aa=None):

        na = H_aa.shape[0]
        I_aa = np.eye(na)

        if S_aa is None: S_aa = I_aa

        w_x = self.dlrrf
        n = len(w_x)

        D_lx = self.T_lx * w_x[None, :] / beta
        D_AA = np.kron(D_lx, S_aa) - np.kron(self.T_lx, H_aa)
                
        bc_x = kernel(np.array([0.]), w_x) + kernel(np.array([1.]), w_x)
        D_AA[(n-1)*na:, :] = np.kron(bc_x, S_aa)

        b_Aa = np.zeros((n*na, na))
        b_Aa[(n-1)*na:, :] = -I_aa

        g_xaa = np.linalg.solve(D_AA, b_Aa).reshape((n, na, na))

        return g_xaa
        
        
    def free_greens_function_tau(self, H_aa, beta, S_aa=None):

        w_x = self.dlrrf

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
    

dlr = dlrBase

#from .pydlr_fortran import dlrFortran as dlr
        
