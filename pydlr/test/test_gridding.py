"""Test gridding of imaginary time analytical continuation kernel.

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
from pydlr.kernel import kernel_discretization


class TestBarycentricInterpolation(unittest.TestCase):


    def test_gridding(self, verbose=False):

        lamb = 40.
        d = dlr(lamb=lamb)
        kmat, t, w, err = kernel_discretization(lamb, error_est=True)
        print(f'err = {err:2.2E}')

        if verbose:

            import matplotlib.pyplot as plt
            plt.figure(figsize=(6, 10))

            subp = [4, 1, 1]

            plt.subplot(*subp); subp[-1] += 1
            #plt.plot(t_panel_break_pt, 0.*t_panel_break_pt, 'or')
            plt.plot(d.t, 0.*d.t, '+')
            plt.plot(t, 0.*t, 'x')

            plt.subplot(*subp); subp[-1] += 1
            plt.semilogy(d.t, np.abs(t - d.t), '.-')

            plt.subplot(*subp); subp[-1] += 1
            #plt.plot(w_panel_break_pt, 0.*w_panel_break_pt, 'or')
            plt.plot(d.om, 0.*d.om, '+', alpha=0.5)
            plt.plot(w, 0.*w, 'x')

            plt.subplot(*subp); subp[-1] += 1
            plt.semilogy(d.om, np.abs(w - d.om), '.-')

            plt.show()

        self.assertTrue(np.allclose(t, d.t))
        self.assertTrue(np.allclose(w, d.om))
        self.assertTrue(np.allclose(d.kmat, kmat))
        self.assertTrue(err < 1e-14)
        

if __name__ == '__main__':
    unittest.main()
    
