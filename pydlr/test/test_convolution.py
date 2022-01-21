"""Tests for the imaginary time convolution routines.

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


import unittest

import numpy as np

from pydlr import dlr


class TestConvolution(unittest.TestCase):

    
    def test_convolution_scalar(self, verbose=False):

        xi = +1
        beta = 2.337
        #beta = 3.337
        #beta = 100.

        lamb = np.max([20., 10. * beta])
        e1, e2 = 3., 3.3

        d = dlr(lamb=lamb, xi=xi)

        g1_l = d.free_greens_function_tau(np.array([[e1]]), beta)
        g2_l = d.free_greens_function_tau(np.array([[e2]]), beta)

        c = 1./(e1 - e2)
        gg_l_anal = c * g1_l - c * g2_l

        g1_q = d.free_greens_function_matsubara(np.array([[e1]]), beta)
        g2_q = d.free_greens_function_matsubara(np.array([[e2]]), beta)

        # -- DLR Matsubara convolution

        gg_q = np.einsum('qab,qbc->qac', g1_q, g2_q)
        gg_l_matsub = d.tau_from_dlr(d.dlr_from_matsubara(gg_q, beta))

        # -- DLR coeff convolution

        g1_x = d.dlr_from_matsubara(g1_q, beta)
        g2_x = d.dlr_from_matsubara(g2_q, beta)

        gg_x_dlr = d.convolution(g1_x, g2_x, beta=beta)
        gg_l_dlr = d.tau_from_dlr(gg_x_dlr)

        # -- DLR convolution matrix

        n, na, _ = g1_x.shape

        C_xaxa = d.convolution_matrix(g1_x, beta=beta)
        C_AA = C_xaxa.reshape((n*na, n*na))

        B_Aa = g2_x.reshape((n*na, na))
        gg_x_mat = np.matmul(C_AA, B_Aa).reshape((n, na, na))
        gg_l_mat = d.tau_from_dlr(gg_x_mat)

        if verbose:
            print(f'diff scalar = {np.max(np.abs(gg_l_anal - gg_l_matsub))} (matsubara)')
            print(f'diff scalar = {np.max(np.abs(gg_l_anal - gg_l_dlr))} (dlr)')
            print(f'diff scalar = {np.max(np.abs(gg_l_anal - gg_l_mat))} (dlr convmat)')

            # -- Viz

            tau_l = d.get_tau(beta)

            import matplotlib.pyplot as plt

            plt.figure(figsize=(6, 6))
            subp = [1, 1, 1]

            plt.subplot(*subp); subp[-1] += 1
            plt.plot(tau_l, gg_l_dlr[:, 0, 0].real, 'o', label='DLR conv', alpha=0.5)
            plt.plot(tau_l, gg_l_mat[:, 0, 0].real, '+', label='DLR conv mat')
            plt.plot(tau_l, gg_l_matsub[:, 0, 0].real, 's', label='Matsub conv')
            plt.plot(tau_l, gg_l_anal[:, 0, 0].real, 'x', label='Analytic')
            plt.ylabel(r'$(g_1 \ast g_2)(\tau)$')
            plt.xlabel(r'$\tau$')
            plt.legend(loc='best')
            plt.show()

        # -- Test

        self.assertTrue(np.allclose(gg_l_anal, gg_l_matsub))
        self.assertTrue(np.allclose(gg_l_anal, gg_l_dlr))
        self.assertTrue(np.allclose(gg_l_anal, gg_l_mat))


    def test_convolution_matrix(self, verbose=False):

        beta = 5.337

        E1_aa = np.array([
            [-1., 0.2, 0.4 + 1.j],
            [0.2, 0, 0.1j],
            [0.4 - 1.j, -0.1j, 1],
            ])

        self.assertTrue(np.allclose(E1_aa, E1_aa.conj().T))

        E2_aa = np.array([
            [-3., 0.3, 0.1 + 0.2j],
            [0.3, 1, 0.3j],
            [0.1 - 0.2j, -0.3j, 0],
            ])

        self.assertTrue(np.allclose(E2_aa, E2_aa.conj().T))

        d = dlr(lamb=40.)
        w_q = d.get_matsubara_frequencies(beta=beta)

        g1_qaa = d.free_greens_function_matsubara(E1_aa, beta)
        g2_qaa = d.free_greens_function_matsubara(E2_aa, beta)

        # -- DLR Matsubara convolution

        gg_qaa = np.einsum('qab,qbc->qac', g1_qaa, g2_qaa)

        gg_xaa = d.dlr_from_matsubara(gg_qaa, beta)
        gg_laa = d.tau_from_dlr(gg_xaa)

        # -- DLR coeff convolution

        g1_xaa = d.dlr_from_matsubara(g1_qaa, beta)
        g2_xaa = d.dlr_from_matsubara(g2_qaa, beta)

        gg_xaa_ref = d.convolution(g1_xaa, g2_xaa, beta)
        gg_laa_ref = d.tau_from_dlr(gg_xaa_ref)

        # -- DLR convolution matrix

        n, na, _ = g1_xaa.shape

        C_xaxa = d.convolution_matrix(g1_xaa, beta=beta)
        C_AA = C_xaxa.reshape((n*na, n*na))

        B_Aa = g2_xaa.reshape((n*na, na))
        gg_xaa_mat = np.matmul(C_AA, B_Aa).reshape((n, na, na))
        gg_laa_mat = d.tau_from_dlr(gg_xaa_mat)

        if verbose:
            print(f'diff matrix = {np.max(np.abs(gg_laa - gg_laa_ref))}')
            print(f'diff matrix (convmat) = {np.max(np.abs(gg_laa - gg_laa_mat))}')

        self.assertTrue(np.allclose(gg_laa, gg_laa_ref))
        self.assertTrue(np.allclose(gg_laa, gg_laa_mat))
    

if __name__ == '__main__':
    unittest.main()
