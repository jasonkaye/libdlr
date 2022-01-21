"""Compare Fortran and Python kernel implementations.

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

from pydlr import dlr, kernel


class TestKernel(unittest.TestCase):
    
    def test_kernel(self, verbose=False):

        d = dlr(lamb=10., python_impl=True, verbose=True)
        d.tt = (d.t > 0) * d.t + (d.t[::-1] > 0) * (1 - d.t[::-1])

        kmat = kernel(d.tt, d.om)

        self.assertTrue(np.allclose(kmat, d.kmat))

        if verbose:
            import matplotlib.pyplot as plt

            plt.figure(figsize=(6, 8))

            subp = [4, 1, 1]

            plt.subplot(*subp); subp[-1] += 1
            plt.title(r'Kernel $K_\Lambda(\tau, \omega)$, $\Lambda = %3.1f$, $\epsilon = %2.2E$' % (d.lamb, d.eps))
            plt.pcolormesh(d.tt, d.om, d.kmat.T, shading='nearest',vmin=0, vmax=1.0)
            plt.xlabel(r'$\tau$')
            plt.ylabel(r'$\omega$')
            plt.colorbar()

            plt.subplot(*subp); subp[-1] += 1
            plt.title(r'Kernel $K_\Lambda(\tau, \omega)$, $\Lambda = %3.1f$, $\epsilon = %2.2E$' % (d.lamb, d.eps))
            plt.pcolormesh(d.tt, d.om, kmat.T, shading='nearest',vmin=0, vmax=1.0)
            plt.xlabel(r'$\tau$')
            plt.ylabel(r'$\omega$')
            plt.colorbar()

            plt.subplot(*subp); subp[-1] += 1
            plt.title(r'$N_{dlr} = %i$' % d.rank)
            plt.plot(d.tt, 0*d.tt, '.-')
            plt.plot(0*d.om, d.om, '.-')
            plt.plot(d.tt[d.tidx], 0*d.tt[d.tidx], 'x', label=r'$\tau_i$ (DLR)')
            plt.plot(0.*d.om[d.oidx], d.om[d.oidx], 'x', label=r'$\omega_j$ (DLR)')
            plt.xlabel(r'$\tau$')
            plt.ylabel(r'$\omega$')
            plt.legend(loc='best')

            plt.subplot(*subp); subp[-1] += 1
            plt.plot(d.dlrmf, 0*d.dlrmf, 'o', label=r'$i\omega_n$ (DLR) subset')
            plt.xlabel(r'$i \omega_n$')
            plt.legend(loc='best')

            plt.tight_layout()
            plt.show()


if __name__ == '__main__':
    unittest.main()
