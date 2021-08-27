
""" Author: Hugo U.R. Strand (2021) """


import numpy as np

from pydlr import dlr


beta = 128.0
d = dlr(lamb=40. * beta)

tau_from_matsubara = lambda G_qaa : d.tau_from_dlr(d.dlr_from_matsubara(G_qaa, beta))
matsubara_from_tau = lambda G_iaa : d.matsubara_from_dlr(d.dlr_from_tau(G_iaa), beta)

e0, e1 = -0.55, 0.3
V = 10.2 # NB! Large hybridization V gives kernel g0 * delta with non-trivial condition number

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
