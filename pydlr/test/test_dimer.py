
import numpy as np

from pydlr import dlr, fermi_function

def test_dimer(verbose=False):
    
    # -- Setup test system

    N = 40
    xi = -1.
    beta = 1.3
    lamb = 10.
    Nw = int(lamb)

    n = np.arange(-Nw, Nw+1)
    zeta = (1 - xi)/2
    iwn = 1.j * np.pi/beta * (2*n + zeta) 

    e1, e2 = 3., 3.3
    V = 4

    E_aa = np.array([
        [e1, V],
        [V, e2],
        ])

    E, U = np.linalg.eigh(E_aa)
    rho = U.conj().T @ np.diag(fermi_function(E, beta)) @ U
    n00 = rho[0, 0]

    if verbose:
        print(f'rho =\n {rho}')
        print(f'n00 = {n00} (anal)')

    # -- DLR

    d = dlr(lamb=lamb)

    tau_l = d.get_tau(beta)
    shape = (len(tau_l), 1, 1)
    G00_l = d.free_greens_function_tau(E_aa, beta)[:, 0, 0].reshape(shape)
    
    G00_x = d.dlr_from_tau(G00_l)

    tau_i = np.linspace(0, beta, num=512)
    G00_i = d.eval_dlr_tau(G00_x, tau_i, beta)

    np.testing.assert_array_almost_equal(G00_l, d.tau_from_dlr(G00_x))
    np.testing.assert_array_almost_equal(G00_l, d.eval_dlr_tau(G00_x, tau_l, beta))

    G00_i_ref = d.eval_dlr_tau(G00_x, tau_i, beta)

    w_q = d.get_matsubara_frequencies(beta)
    G00_q = d.matsubara_from_dlr(G00_x, beta)

    np.testing.assert_array_almost_equal(G00_q, d.eval_dlr_freq(G00_x, w_q, beta))

    G00_q_anal = 1./ ( w_q - e1 - V**2/( w_q - e2 ) )
    G00_q_anal = G00_q_anal.reshape(shape)
    np.testing.assert_array_almost_equal(G00_q, G00_q_anal)

    G00_x_ref = d.dlr_from_matsubara(G00_q_anal, beta)

    G00_l_ref = d.tau_from_dlr(G00_x_ref)
    np.testing.assert_array_almost_equal(G00_l, G00_l_ref)

    G00_q_ref = d.matsubara_from_dlr(G00_x_ref, beta)
    np.testing.assert_array_almost_equal(G00_q, G00_q_ref)

    # -- Dyson solution
    
    g2_qaa = d.free_greens_function_matsubara(np.array([[e2]]), beta)
    G00_q_dyson_ref = d.dyson_matsubara(np.array([[e1]]), V**2 * g2_qaa, beta)
    np.testing.assert_array_almost_equal(G00_q, G00_q_dyson_ref)
    
    if verbose:
        
        # -- Viz

        G00_i = np.squeeze(G00_i)
        G00_l = np.squeeze(G00_l)
        G00_q = np.squeeze(G00_q)
        G00_x = np.squeeze(G00_x)
        G00_x_ref = np.squeeze(G00_x_ref)

        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 9))

        subp = [3, 1, 1]

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(tau_i, G00_i, label=r'$G(\tau)$')
        plt.plot(tau_l, G00_l, '.', label=r'$G(\tau_l)$ DLR points')
        plt.xlabel(r'$\tau$')
        plt.legend()

        plt.subplot(*subp); subp[-1] += 1
        s = np.argsort(w_q.imag)
        plt.plot(w_q[s].imag, G00_q[s].real, '.-')
        plt.plot(w_q[s].imag, G00_q[s].imag, '.-')
        plt.xlabel(r'$i\omega_n$')
        plt.ylabel(r'$G(i\omega_n)$')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(np.abs(G00_x.real), '.-', label='Re DLR (tau)')
        plt.plot(np.abs(G00_x.imag), '.-', label='Im DLR (tau)')
        plt.plot(np.abs(G00_x_ref.real), 'x-', label='Re DLR (tau)')
        plt.plot(np.abs(G00_x_ref.imag), 'x-', label='Im DLR (tau)')
        plt.semilogy([], [])
        plt.legend()
        plt.xlabel(r'$n$')
        plt.ylabel(r'$|G_n|$ expansion coeff')

        plt.tight_layout()
        plt.show()


if __name__ == '__main__':

    test_dimer(verbose=True)
