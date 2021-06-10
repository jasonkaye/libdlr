""" Evaluation of integral between a semi-circular 
spectral function and the imaginary time kernel, giving 
the imaginary time Green's function analytically.

The analytical result is compared with the DLR representation
on a dense tau mesh.

The result is also compared to the iterative solution of
the corresponding Bethe lattice self-consistency \Sigma = G / 4

Author: Hugo U.R. Strand (2021) """


import numpy as np
from scipy.integrate import quad

from pydlr import dlr, kernel


def eval_semi_circ_G_tau(t):
    I = lambda x: -2 / np.pi * kernel(np.array([t]), np.array([x]))[0,0]        
    g, res = quad(I, -1, 1, weight='alg', wvar=(0.5, 0.5))
    return g

eval_semi_circ_G_tau = np.vectorize(eval_semi_circ_G_tau)


def test_semi_cirular_G_tau(verbose=False):

    # -- Dense mesh analytic evaluation
    
    tau_i = np.linspace(0, 1, num=200)
    G_i = eval_semi_circ_G_tau(tau_i)

    # -- Evaluation on DLR-tau points
    # -- and reinterpolation on dense grid from DLR coefficients
    
    d = dlr(lamb=10.)
    tau_l = d.get_tau()
    shape = (len(tau_l), 1, 1)
    G_l = eval_semi_circ_G_tau(tau_l).reshape(shape)
    G_x = d.dlr_from_tau(G_l)
    G_i_ref = d.eval_dlr_tau(G_x, tau_i, 1.)[:, 0, 0]

    G_l = np.squeeze(G_l)
    
    print(f'diff = {np.max(np.abs(G_i - G_i_ref))}')
    np.testing.assert_array_almost_equal(G_i, G_i_ref)

    # -- Iterative determination using the Dyson equation

    max_iter = 20
    G_q = np.zeros(len(tau_l), dtype=complex)
    
    for iter in range(max_iter):
        G_q = d.dyson_matsubara(np.array([[0.]]), 0.25 * G_q.reshape(shape), 1.)[:, 0, 0]
        G_x_ref = d.dlr_from_matsubara(G_q)
        G_l_ref = d.tau_from_dlr(G_x_ref).real
        G_i_ref2 = d.eval_dlr_tau(G_x_ref, tau_i, 1.)

        diff = np.max(np.abs(G_i - G_i_ref2))
        print(f'diff = {diff}')

        if diff < 5e-14: break

    np.testing.assert_array_almost_equal(G_l, G_l_ref)
    np.testing.assert_array_almost_equal(G_i, G_i_ref2)
        
    if verbose:
        
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 6))

        subp = [2, 1, 1]

        plt.subplot(*subp); subp[-1] += 1
        plt.title('Semi-circular spectral function')
        plt.plot(tau_i, G_i, '-', label='analytic')
        plt.plot(tau_l, G_l, '.', label='tau DLR points')
        plt.plot(tau_l, G_l_ref, 'x', label='tau DLR points (iterative)')
        plt.xlabel(r'$\tau$')
        plt.ylabel(r'$G(\tau)$')
        plt.legend(loc='best')    

        plt.subplot(*subp); subp[-1] += 1
        plt.semilogy(tau_i, np.abs(G_i - G_i_ref), '.-', label='analytic')
        plt.semilogy(tau_i, np.abs(G_i - G_i_ref2), '.-', label='iterative')
        plt.xlabel(r'$\tau$')
        plt.ylabel(r'$|G(\tau) - G_{DLR}(\tau)|$')
        plt.legend(loc='best')    

        plt.tight_layout()
        plt.savefig('figure_test_semi_circ.pdf')

        plt.show()


if __name__ == '__main__':

    test_semi_cirular_G_tau(verbose=True)
    
