"""Test all Dyson equation solvers on exactly solvable dimer problem.

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


class TestDyson(unittest.TestCase):


    def setUp(self):
    
        xi = +1
        beta = 128.0
        d = dlr(lamb=40. * beta, xi=xi)

        d.tau_from_matsubara = lambda G_qaa : d.tau_from_dlr(d.dlr_from_matsubara(G_qaa, beta))
        d.matsubara_from_tau = lambda G_iaa : d.matsubara_from_dlr(d.dlr_from_tau(G_iaa), beta)
        
        #e0, e1 = -0.55, 0.3
        #V = 10.2 # NB! Large hybridization V gives kernel g0 * delta with non-trivial condition number

        e0, e1 = 1.0, 0.1
        V = 0.06

        E_aa = np.array([
            [e0, V],
            [V, e1],
            ])

        E0 = np.array([[e0]])
        E1 = np.array([[e1]])

        g0_iaa = d.free_greens_function_tau(E0, beta)
        g0_qaa = d.matsubara_from_tau(g0_iaa)

        g1_iaa = d.free_greens_function_tau(E1, beta)
        delta_iaa = V**2 * g1_iaa

        self.delta_xaa = d.dlr_from_tau(delta_iaa)
        self.delta_qaa = d.matsubara_from_dlr(self.delta_xaa, beta)

        self.G_iaa = d.free_greens_function_tau(E_aa, beta)[:, 0, 0].reshape((d.rank, 1, 1))

        self.d = d
        self.beta = beta
        self.E0 = E0


    
    def test_dyson_matsubara(self):
        G_iaa = self.d.tau_from_matsubara(self.d.dyson_matsubara(self.E0, self.delta_qaa, self.beta))
        print(f'err = {np.max(np.abs(G_iaa - self.G_iaa)):2.2E}')
        self.assertTrue(np.allclose(G_iaa, self.G_iaa))


    def test_dyson_dlr_integrodiff(self):
        G_iaa = self.d.tau_from_dlr(self.d.dyson_dlr_integrodiff(self.E0, self.delta_xaa, self.beta))
        print(f'err = {np.max(np.abs(G_iaa - self.G_iaa)):2.2E}')
        self.assertTrue(np.allclose(G_iaa, self.G_iaa))


    def test_dyson_dlr(self):
        G_iaa = self.d.tau_from_dlr(self.d.dyson_dlr(
            self.E0, self.delta_xaa, self.beta, iterative=False, lomem=False))
        print(f'err = {np.max(np.abs(G_iaa - self.G_iaa)):2.2E}')
        self.assertTrue(np.allclose(G_iaa, self.G_iaa))


    def test_dyson_dlr_iterative(self):
        G_iaa = self.d.tau_from_dlr(self.d.dyson_dlr(
            self.E0, self.delta_xaa, self.beta, iterative=True, lomem=False))
        print(f'err = {np.max(np.abs(G_iaa - self.G_iaa)):2.2E}')
        self.assertTrue(np.allclose(G_iaa, self.G_iaa))


    def test_dyson_dlr_iterative_lomem(self):
        G_iaa = self.d.tau_from_dlr(self.d.dyson_dlr(
            self.E0, self.delta_xaa, self.beta, iterative=True, lomem=True))
        print(f'err = {np.max(np.abs(G_iaa - self.G_iaa)):2.2E}')
        self.assertTrue(np.allclose(G_iaa, self.G_iaa))


if __name__ == '__main__':
    unittest.main()
