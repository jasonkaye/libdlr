
import itertools
import numpy as np

from pydlr import dlr, kernel


def free_gf(d, H_aa, beta, S_aa=None):
    

    return g_laa


def test_free_greens_function_scalar(verbose=False):

    lamb = 9.
    beta = 1.234
    
    H_aa = np.array([[0.3]])
    S_aa = np.array([[1.]])
        
    d = dlr(lamb=lamb)

    tau_l = d.get_tau(beta)
    G_laa = d.free_greens_function_tau(H_aa, beta, S_aa=S_aa)
    G_xaa = d.dlr_from_tau(G_laa)
    G_qaa_ref = d.matsubara_from_dlr(G_xaa, beta)

    tau_i = np.linspace(0, beta, num=400)
    G_iaa = d.eval_dlr_tau(G_xaa, tau_i, beta)

    w_q = d.get_matsubara_frequencies(beta)

    G_qaa = d.free_greens_function_matsubara(H_aa, beta, S_aa=S_aa)
    G_xaa_ref = d.dlr_from_matsubara(G_qaa, beta)
    G_laa_ref = d.tau_from_dlr(G_xaa_ref).real

    print(f'diff tau = {np.max(np.abs(G_laa - G_laa_ref))}')
    print(f'diff Matsubara = {np.max(np.abs(G_qaa - G_qaa_ref))}')

    print(f'G_laa = {np.squeeze(G_laa)}')

    
    w_x = d.dlrrf
    n = len(w_x)
    zeros = np.zeros(n)

    D_x = w_x/beta - np.squeeze(H_aa)
    D_lx = d.T_lx * D_x[None, :]
    
    DG_x = D_x * np.squeeze(G_xaa)
    DG_l = d.tau_from_dlr(DG_x)
    np.testing.assert_array_almost_equal(DG_l, zeros)

    DG_l_ref = np.matmul(D_lx, np.squeeze(G_xaa))
    np.testing.assert_array_almost_equal(DG_l_ref, zeros)
    
    #print(f'DG_x = {DG_x}')
    #print(f'DG_l = {DG_l}')
    #print(f'DG_l_ref = {DG_l_ref}')

    k0_x = kernel(np.array([0.]), w_x)
    k1_x = kernel(np.array([1.]), w_x)
    bc_x = k0_x + k1_x
    bc = np.dot(bc_x, np.squeeze(G_xaa))
    print(f'bc = {bc}')

    D_lx[-1, :] = bc_x
    b_l = np.zeros(n)
    b_l[-1] = -1
    
    print(f'shape = {D_lx.shape}')
    print(f'rank = {np.linalg.matrix_rank(D_lx)}')
    
    g_x = np.linalg.solve(D_lx, b_l)
    g_l = d.tau_from_dlr(g_x)
    g_laa = g_l.reshape((n, 1, 1))

    print(f'g_x = {g_x}')
    print(f'g_l = {g_l}')

    print(f'diff tau DLR = {np.max(np.abs(G_laa - g_laa))}')

    # -----

    na = H_aa.shape[0]
    I_aa = np.eye(na)
    
    if S_aa is None: S_aa = I_aa
    
    w_x = d.dlrrf
    n = len(w_x)
    
    D_xaa = w_x[:, None, None]/beta * S_aa[None, ...] - H_aa[None, ...]
    D_lxaa = d.T_lx[:, :, None, None] * D_xaa[None, ...]
    
    D_laxa = np.moveaxis(D_lxaa, 2, 1)
    D_AA = D_lxaa.reshape(n*na, n*na)
    
    bc_x = kernel(np.array([0.]), w_x) + kernel(np.array([1.]), w_x)
    D_AA[(n-1)*na:, :] = np.kron(bc_x, S_aa)

    #print(f'D_lx = \n {D_lx}')
    #print(f'D_AA = \n {D_AA}')

    diff = np.max(np.abs(D_lx - D_AA))
    print(f'diff = {diff}')
    np.testing.assert_array_almost_equal(D_lx, D_AA)
    
    b_Aa = np.zeros((n*na, na))
    b_Aa[(n-1)*na:, :] = -I_aa
    
    diff = np.max(np.abs(b_l - np.squeeze(b_Aa)))
    print(f'diff = {diff}')
    np.testing.assert_array_almost_equal(b_l, np.squeeze(b_Aa))

    g_xaa_ref = np.linalg.solve(D_AA, b_Aa).reshape((n, na, na))
    g_laa_ref = d.tau_from_dlr(g_xaa_ref)
    
    np.testing.assert_array_almost_equal(g_laa, g_laa_ref)
    
    g_xaa_ref = d.free_greens_function_dlr(H_aa, beta)
    g_laa_ref = d.tau_from_dlr(g_xaa_ref)

    np.testing.assert_array_almost_equal(g_laa, g_laa_ref)
    
    exit()
    
    if verbose:
        
        sidx = np.argsort(tau_l)
        tau_l = tau_l[sidx]
        G_laa = G_laa[sidx]
        G_laa_ref = G_laa_ref[sidx]
        # -- Viz

        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 9))

        subp = [2, 1, 1]    

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(tau_l, G_laa[:, 0, 0], '+')
        plt.plot(tau_l, G_laa_ref[:, 0, 0], 'x')
        plt.plot(tau_l, g_laa[:, 0, 0], '.')
        plt.plot(tau_i, G_iaa[:, 0, 0], '-')
        plt.legend(loc='upper right')
        plt.xlabel(r'$\tau$')

        plt.subplot(*subp); subp[-1] += 1
        lr = plt.plot(w_q.imag, G_qaa[:, 0, 0].real, '+', label=f'Re')
        li = plt.plot(w_q.imag, G_qaa[:, 0, 0].imag, '+', label=f'Im')
        
        cr, ci = [x[0].get_color() for x in [lr, li]]
        plt.plot(w_q.imag, G_qaa_ref[:, 0, 0].real, 'x', color=cr)
        plt.plot(w_q.imag, G_qaa_ref[:, 0, 0].imag, 'x', color=ci)
        plt.legend(loc='upper right')
        plt.xlabel(r'$i\omega_n$')
        
        plt.show()
        
    np.testing.assert_array_almost_equal(G_laa, G_laa_ref)
    np.testing.assert_array_almost_equal(G_qaa, G_qaa_ref)


def test_free_greens_function_matrix(verbose=False):

    lamb = 30.
    beta = 10.234
    
    H_aa = np.array([
        [1, 0.4, 0.1],
        [0.4, 2, 0.5],
        [0.1, 0.5, 3],
        ])

    S_aa = np.array([
        [1, 0.1, 0.2],
        [0.1, 1, 0.3],
        [0.2, 0.3, 1],
        ])
        
    d = dlr(lamb=lamb)

    tau_l = d.get_tau(beta)
    G_laa = d.free_greens_function_tau(H_aa, beta, S_aa=S_aa)
    G_xaa = d.dlr_from_tau(G_laa)
    G_qaa_ref = d.matsubara_from_dlr(G_xaa, beta)

    tau_i = np.linspace(0, beta, num=400)
    G_iaa = d.eval_dlr_tau(G_xaa, tau_i, beta)

    w_q = d.get_matsubara_frequencies(beta)

    G_qaa = d.free_greens_function_matsubara(H_aa, beta, S_aa=S_aa)
    G_xaa_ref = d.dlr_from_matsubara(G_qaa, beta)
    G_laa_ref = d.tau_from_dlr(G_xaa_ref).real

    print(f'diff tau = {np.max(np.abs(G_laa - G_laa_ref))}')
    print(f'diff Matsubara = {np.max(np.abs(G_qaa - G_qaa_ref))}')

    w_x = d.dlrrf
    print(f'w_x = {w_x}')

    g_xaa = np.zeros_like(G_xaa)
    for idx, w in enumerate(w_x):
        g_xaa[idx] = np.linalg.inv(-w/beta * S_aa + H_aa)

    g_laa = d.tau_from_dlr(g_xaa)
    print(f'diff tau DLR = {np.max(np.abs(G_laa - g_laa))}')
        
    if verbose:
        
        sidx = np.argsort(tau_l)
        tau_l = tau_l[sidx]
        G_laa = G_laa[sidx]
        G_laa_ref = G_laa_ref[sidx]
        # -- Viz

        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 9))

        subp = [3, 3, 1]    

        #plt.plot(tau_i, G00_i, label=r'$G(\tau)$')
        for i, j in itertools.product(range(3), repeat=2):
            plt.subplot(*subp); subp[-1] += 1
            plt.plot(tau_l, G_laa[:, i, j], '+', label=f'G{i}{j}')
            plt.plot(tau_l, G_laa_ref[:, i, j], 'x')
            plt.plot(tau_l, g_laa[:, i, j], '.')
            plt.plot(tau_i, G_iaa[:, i, j], '-')
            plt.legend(loc='upper right')
        plt.xlabel(r'$\tau$')

        plt.figure(figsize=(8, 9))

        subp = [3, 3, 1]    

        #plt.plot(tau_i, G00_i, label=r'$G(\tau)$')
        for i, j in itertools.product(range(3), repeat=2):
            plt.subplot(*subp); subp[-1] += 1
            lr = plt.plot(w_q.imag, G_qaa[:, i, j].real, '+', label=f'Re[G{i}{j}]')
            li = plt.plot(w_q.imag, G_qaa[:, i, j].imag, '+', label=f'Im[G{i}{j}]')
            cr, ci = [x[0].get_color() for x in [lr, li]]
            plt.plot(w_q.imag, G_qaa_ref[:, i, j].real, 'x', color=cr)
            plt.plot(w_q.imag, G_qaa_ref[:, i, j].imag, 'x', color=ci)
            plt.legend(loc='upper right')
        plt.xlabel(r'$i\omega_n$')
        
        plt.show()
        
    np.testing.assert_array_almost_equal(G_laa, G_laa_ref)
    np.testing.assert_array_almost_equal(G_qaa, G_qaa_ref)


if __name__ == '__main__':

    test_free_greens_function_scalar(verbose=True)
    #test_free_greens_function_matrix(verbose=True)
