
""" Author: Hugo U.R. Strand (2021) """

import time
import numpy as np
import matplotlib.pyplot as plt

from pydlr import dlr


def test_dyson_and_volterra_matsubara(beta=10., verbose=False):

    t = time.time()
    lamb = 10. * beta
    print(f'lambd = {lamb}')
    
    d = dlr(lamb=lamb)
    print(f'init dlr {time.time() - t} s')

    e0, e1 = -0.55, 0.3
    V = 0.2
    #V = 1.2

    E_aa = np.array([
        [e0, V],
        [V, e1],
        ])
    
    w_q = d.get_matsubara_frequencies(beta=beta)

    tau_i = d.get_tau(beta)

    print(f'np = {len(tau_i)}')
    
    g0_iaa = d.free_greens_function_tau(np.array([[e0]]), beta)
    g1_iaa = d.free_greens_function_tau(np.array([[e1]]), beta)

    G_iaa = d.free_greens_function_tau(E_aa, beta)[:, 0, 0].reshape((len(tau_i), 1, 1))
        
    #g0_qaa = d.free_greens_function_matsubara(np.array([[e0]]), beta)
    #g1_qaa = d.free_greens_function_matsubara(np.array([[e1]]), beta)

    #G_qaa = d.free_greens_function_matsubara(E_aa, beta)[:, 0, 0].reshape((len(w_q), 1, 1))

    g0_xaa = d.dlr_from_tau(g0_iaa)
    g1_xaa = d.dlr_from_tau(g1_iaa)
    G_xaa = d.dlr_from_tau(G_iaa)

    g0_qaa = d.matsubara_from_dlr(g0_xaa, beta)
    g1_qaa = d.matsubara_from_dlr(g1_xaa, beta)

    G_qaa = d.matsubara_from_dlr(G_xaa, beta)
    
    G_qaa_mat_dyson = d.dyson_matsubara(np.array([[e0]]), V**2 * g1_qaa, beta)
    G_qaa_mat_volt = d.volterra_matsubara(g0_qaa, V**2 * g1_qaa, beta)

    G_iaa_mat_dyson = d.tau_from_dlr(d.dlr_from_matsubara(G_qaa_mat_dyson, beta))
    G_iaa_mat_volt = d.tau_from_dlr(d.dlr_from_matsubara(G_qaa_mat_volt, beta))
    
    #np.testing.assert_array_almost_equal(G_qaa, G_qaa_dyson)
    #np.testing.assert_array_almost_equal(G_qaa, G_qaa_volterra)

    G_xaa_dlr_dyson = d.dyson_dlr_integrodiff(np.array([[e0]]), V**2 * g1_xaa, beta)
    G_qaa_dlr_dyson = d.matsubara_from_dlr(G_xaa_dlr_dyson, beta)
    G_iaa_dlr_dyson = d.tau_from_dlr(G_xaa_dlr_dyson)

    t = time.time()
    G_xaa_dlr_dyson_int = d.dyson_dlr(np.array([[e0]]), V**2 * g1_xaa, beta)
    print(f'solve {time.time() - t} s')
    G_iaa_dlr_dyson_int = d.tau_from_dlr(G_xaa_dlr_dyson_int)

    t = time.time()
    G_xaa_dlr_dyson_int_iter = d.dyson_dlr(np.array([[e0]]), V**2 * g1_xaa, beta, iterative=True, lomem=True, verbose=True, tol=1e-14)
    print(f'gmres {time.time() - t} s')
    G_iaa_dlr_dyson_int_iter = d.tau_from_dlr(G_xaa_dlr_dyson_int_iter)

    #np.testing.assert_array_almost_equal(G_qaa, G_qaa_dyson_dlr)
    #np.testing.assert_array_almost_equal(G_iaa, G_iaa_dyson_dlr)

    err_mat_volt = np.max(np.abs(G_iaa - G_iaa_mat_volt))
    err_mat_dyson = np.max(np.abs(G_iaa - G_iaa_mat_dyson))
    err_dlr_dyson = np.max(np.abs(G_iaa - G_iaa_dlr_dyson))
    err_dlr_dyson_int = np.max(np.abs(G_iaa - G_iaa_dlr_dyson_int))
    err_dlr_dyson_int_iter = np.max(np.abs(G_iaa - G_iaa_dlr_dyson_int_iter))

    print(f'err_mat_volt  = {err_mat_volt}')
    print(f'err_mat_dyson = {err_mat_dyson}')
    print(f'err_dlr_dyson = {err_dlr_dyson}')
    print(f'err_dlr_dyson_int = {err_dlr_dyson_int}')
    print(f'err_dlr_dyson_int_iter = {err_dlr_dyson_int_iter}')

    if verbose:
        plt.figure(figsize=(10, 12))
        subp = [2, 1, 1]

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(tau_i, -np.squeeze(g0_iaa), '.-')
        plt.plot(tau_i, -np.squeeze(g1_iaa), '.-')
        plt.plot(tau_i, -np.squeeze(G_iaa), '.-')
        plt.semilogy([], [])
        
        plt.xlabel(r'$\tau$')
        plt.ylabel(r'$G(\tau)$')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(tau_i, -np.squeeze(G_iaa), '.-', label='ref')
        plt.plot(tau_i, -np.squeeze(G_iaa_mat_dyson).real, '.-', label='mat dyson')
        plt.plot(tau_i, -np.squeeze(G_iaa_mat_volt).real, '.-', label='mat volt')
        plt.plot(tau_i, -np.squeeze(G_iaa_dlr_dyson).real, '.-', label='dlr dyson')
        plt.plot(tau_i, -np.squeeze(G_iaa_dlr_dyson_int).real, '.-', label='dlr dyson int')
        plt.legend()
        plt.semilogy([], [])

        plt.show(); exit()

    return err_mat_volt, err_mat_dyson, err_dlr_dyson, err_dlr_dyson_int, err_dlr_dyson_int_iter
    

if __name__ == '__main__':

    betas = 8 * 2 ** np.arange(0, 10)
    #betas = [1024]
    print(f'betas = {betas}')

    errs = []
    for beta in betas:
        print('-'*72)
        print(f'beta = {beta}')
        mv, md, dd, di, dii = test_dyson_and_volterra_matsubara(beta=beta, verbose=False)
        errs.append((mv, md, dd, di, dii))

    mv, md, dd, di, dii = np.array(errs).T

    plt.plot(betas, mv, '.-', label='matsubara dyson integro')
    plt.plot(betas, md, '.-', label='matsubara dyson inverse')
    plt.plot(betas, dd, '.-', label='dlr dyson integro-diff')
    plt.plot(betas, di, '.-', label='dlr dyson integro')
    plt.plot(betas, dii, '.-', label='dlr dyson integro (gmres)')
    plt.legend(loc='best')
    plt.loglog([], [])
    plt.xlabel(r'$\beta$')

    plt.savefig('figure_dyson_alg_cf.pdf')
    plt.show()
