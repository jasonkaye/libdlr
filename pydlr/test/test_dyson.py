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


import numpy as np

from pydlr import dlr

xi = +1
#beta = 10.0
beta = 128.0
d = dlr(lamb=40. * beta, xi=xi)

tau_from_matsubara = lambda G_qaa : d.tau_from_dlr(d.dlr_from_matsubara(G_qaa, beta))
matsubara_from_tau = lambda G_iaa : d.matsubara_from_dlr(d.dlr_from_tau(G_iaa), beta)

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
g0_qaa = matsubara_from_tau(g0_iaa)

g1_iaa = d.free_greens_function_tau(E1, beta)
delta_iaa = V**2 * g1_iaa

delta_xaa = d.dlr_from_tau(delta_iaa)
delta_qaa = d.matsubara_from_dlr(delta_xaa, beta)

G_iaa_anal = d.free_greens_function_tau(E_aa, beta)[:, 0, 0].reshape((d.rank, 1, 1))


def check(msg, G_iaa, G_iaa_ref):
    err = np.max(np.abs(G_iaa - G_iaa_ref))
    print(f'err = {err:2.2E} -- {msg}')
    np.testing.assert_array_almost_equal(G_iaa, G_iaa_ref)


def test_dyson_and_volterra_matsubara():

    G_iaa_dyson = tau_from_matsubara(d.dyson_matsubara(E0, delta_qaa, beta))
    check('matsubara dyson', G_iaa_dyson, G_iaa_anal)
    

def test_dyson_dlr():
    
    G_iaa_dlr_dyson_intdiff = d.tau_from_dlr(d.dyson_dlr_integrodiff(E0, delta_xaa, beta))
    check('dlr intdiff', G_iaa_dlr_dyson_intdiff, G_iaa_anal)

    G_iaa_dlr_dyson = d.tau_from_dlr(
        d.dyson_dlr(E0, delta_xaa, beta, iterative=False, lomem=False))
    check('dlr int', G_iaa_dlr_dyson, G_iaa_anal)

    G_iaa_dlr_dyson_iter = d.tau_from_dlr(
        d.dyson_dlr(E0, delta_xaa, beta, iterative=True, lomem=False))
    check('dlr int iter', G_iaa_dlr_dyson_iter, G_iaa_anal)

    G_iaa_dlr_dyson_iter_lomem = d.tau_from_dlr(
        d.dyson_dlr(E0, delta_xaa, beta, iterative=True, lomem=True, verbose=True))
    check('dlr int iter lomem', G_iaa_dlr_dyson_iter_lomem, G_iaa_anal)
    
    
if __name__ == '__main__':

    test_dyson_and_volterra_matsubara()
    test_dyson_dlr()
