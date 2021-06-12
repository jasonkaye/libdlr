
""" Author: Hugo U.R. Strand (2021) """


import numpy as np

from pydlr import dlr


def test_convolution_scalar(verbose=False):

    beta = 3.337
    e1, e2 = 3., 3.3

    d = dlr(lamb=30.)

    g1_q = d.free_greens_function_matsubara(np.array([[e1]]), beta)
    g2_q = d.free_greens_function_matsubara(np.array([[e2]]), beta)

    # -- DLR Matsubara convolution

    gg_q = np.einsum('qab,qbc->qac', g1_q, g2_q)

    gg_x = d.dlr_from_matsubara(gg_q, beta)
    gg_l = d.tau_from_dlr(gg_x)

    # -- DLR coeff convolution

    g1_x = d.dlr_from_matsubara(g1_q, beta)
    g2_x = d.dlr_from_matsubara(g2_q, beta)

    gg_x_dlr = d.convolution(g1_x, g2_x, beta=beta)
    gg_l_dlr = d.tau_from_dlr(gg_x_dlr)

    # -- DLR convolution matrix

    a_x = np.squeeze(g1_x)
    b_x = np.squeeze(g2_x)

    # -- Scalar
    
    tau_l = d.get_tau(1.)
    w_x = d.dlrrf

    n = len(w_x)
    na = g1_x.shape[-1]

    I = np.eye(len(w_x))
    W_xx = 1. / (I + w_x[:, None] - w_x[None, :]) - I

    from pydlr import kernel
    
    k1_x = -np.squeeze(kernel(np.ones(1), w_x))

    #C_xx = np.diag(k1_x * a_x) + \
    #    d.dlr_from_tau(np.einsum('l,ln,n->ln', tau_l, d.T_lx, a_x)) + \
    #    np.diag(np.tensordot(a_x, W_xx, axes=(0, 0))) + a_x[:, None] * W_xx.T

    C_xx = a_x[:, None, ...] * W_xx.T + \
        d.dlr_from_tau(np.einsum('l,lx,x...->lx...', tau_l, d.T_lx, a_x))
    C_xx[np.diag_indices(n)] += k1_x * a_x + np.tensordot(a_x, W_xx, axes=(0, 0))
    C_xx *= beta

    gg_x_mat = np.matmul(C_xx, b_x)
    gg_l_mat = d.tau_from_dlr(gg_x_mat)
    
    diff = np.max(np.abs(np.squeeze(gg_l) - gg_l_mat))
    print(f'diff (mat) = {diff}')

    # -- Matrix

    A_xaa = g1_x
    B_xaa = g2_x

    C_xxaa = A_xaa[:, None, ...] * W_xx.T[:, :, None, None] + \
        d.dlr_from_tau(np.einsum('l,lx,x...->lx...', tau_l, d.T_lx, A_xaa))
    C_xaa = k1_x[:, None, None] * A_xaa + np.tensordot(W_xx, A_xaa, axes=(0, 0))
    for i in range(n): C_xxaa[i, i, :, :] += C_xaa[i]    
    C_xxaa *= beta

    print(C_xxaa.shape)
    
    C_xaxa = np.moveaxis(C_xxaa, 2, 1)
    C_AA = C_xaxa.reshape((n*na, n*na))
    
    B_Aa = B_xaa.reshape((n*na, na))
    gg_x_matmat = np.matmul(C_AA, B_Aa).reshape((n, na, na))
    gg_l_matmat = d.tau_from_dlr(gg_x_matmat)

    diff = np.max(np.abs(gg_l - gg_l_matmat))
    print(f'diff (mat mat) = {diff}')
    
    # -- matrix impl

    C_xaxa = d.convolution_matrix(A_xaa, beta=beta)
    C_AA = C_xaxa.reshape((n*na, n*na))
    
    B_Aa = B_xaa.reshape((n*na, na))
    gg_x_matmat = np.matmul(C_AA, B_Aa).reshape((n, na, na))
    gg_l_matmat = d.tau_from_dlr(gg_x_matmat)

    diff = np.max(np.abs(gg_l - gg_l_matmat))
    print(f'diff (mat mat) = {diff}')
    
    # -- Refs

    A_xaa = g1_x
    B_xaa = g2_x
    
    AB_xaa = np.matmul(A_xaa, B_xaa)

    C_xaa = k1_x[:, None, None] * AB_xaa + \
        d.dlr_from_tau(tau_l[:, None, None] * d.tau_from_dlr(AB_xaa)) + \
        np.matmul(np.tensordot(W_xx, A_xaa, axes=(0, 0)), B_xaa) + \
        np.matmul(A_xaa, np.tensordot(W_xx, B_xaa, axes=(0, 0)))        

    C_xaa *= beta
    C_laa = d.tau_from_dlr(C_xaa)

    diff = np.max(np.abs(gg_l - C_laa))
    print(f'diff (conv) = {diff}')
    
    # -- Term by term

    # 1
    
    C1_xx = beta * np.diag(k1_x * a_x)
    gg1_l_mat = d.tau_from_dlr(np.matmul(C1_xx, b_x))
    C1_laa = d.tau_from_dlr(beta * k1_x[:, None, None] * AB_xaa)
    diff = np.max(np.abs(gg1_l_mat - np.squeeze(C1_laa)))
    print(f'diff (1) = {diff}')

    # 2

    C2_xx = np.diag(np.tensordot(a_x, W_xx, axes=(0, 0))) + a_x[:, None] * W_xx.T
    gg2_l_mat = d.tau_from_dlr(np.matmul(C2_xx, b_x))
    C2_laa = d.tau_from_dlr(
        np.matmul(np.tensordot(W_xx, A_xaa, axes=(0, 0)), B_xaa) + \
        np.matmul(A_xaa, np.tensordot(W_xx, B_xaa, axes=(0, 0)))        
        )
    diff = np.max(np.abs(gg2_l_mat - np.squeeze(C2_laa)))
    print(f'diff (2) = {diff}')

    # 3

    #rank = np.linalg.matrix_rank(d.T_lx)
    #print(f'rank = {rank}, shape = {d.T_lx.shape}')
    #np.testing.assert_array_almost_equal(I, np.matmul(T_lx_inv, d.T_lx))
    
    C3_xx = d.dlr_from_tau(np.einsum('l,ln,n->ln', tau_l, d.T_lx, a_x))
    
    gg3_l_mat = d.tau_from_dlr(np.matmul(C3_xx, b_x))
    C3_laa = d.tau_from_dlr(
        d.dlr_from_tau(tau_l[:, None, None] * d.tau_from_dlr(AB_xaa))
        )
    diff = np.max(np.abs(gg3_l_mat - np.squeeze(C3_laa)))
    print(f'diff (3) = {diff}')
    
    print(gg_l_mat.shape)
    print(gg_l_dlr.shape)
    print(C_laa.shape)
    
    exit()
    
    # --
    
    if verbose:

        # -- Viz

        tau_l = d.get_tau(beta)

        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 6))
        subp = [1, 1, 1]

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(tau_l, gg_l_dlr[:, 0, 0].real, '+', label='DLR conv')
        plt.plot(tau_l, gg_l[:, 0, 0].real, 'x', label='Matsub conv')
        plt.ylabel(r'$(g_1 \ast g_2)(\tau)$')
        plt.xlabel(r'$\tau$')
        plt.legend(loc='best')
        plt.show()

    # -- Test

    print(f'diff scalar = {np.max(np.abs(gg_l - gg_l_dlr))}')
    np.testing.assert_array_almost_equal(gg_l, gg_l_dlr, decimal=9)


def test_convolution_matrix():

    beta = 5.337

    E1_aa = np.array([
        [-1., 0.2, 0.4 + 1.j],
        [0.2, 0, 0.1j],
        [0.4 - 1.j, -0.1j, 1],
        ])

    np.testing.assert_array_almost_equal(E1_aa, E1_aa.conj().T)

    E2_aa = np.array([
        [-3., 0.3, 0.1 + 0.2j],
        [0.3, 1, 0.3j],
        [0.1 - 0.2j, -0.3j, 0],
        ])

    np.testing.assert_array_almost_equal(E2_aa, E2_aa.conj().T)
    
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

    # -- Test

    print(f'diff matrix = {np.max(np.abs(gg_laa - gg_laa_ref))}')
    np.testing.assert_array_almost_equal(gg_laa, gg_laa_ref, decimal=9)
    

if __name__ == '__main__':

    test_convolution_scalar(verbose=True)
    test_convolution_matrix()
