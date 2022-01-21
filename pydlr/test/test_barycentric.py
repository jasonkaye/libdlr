"""Tests for the barycentric interpolation routines.

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

from pydlr.kernel import chebyshev_collocation_points_1st_kind
from pydlr.kernel import chebyshev_barycentric_weights_1st_kind
from pydlr.kernel import barycentric_chebyshev_interpolation


class TestBarycentricInterpolation(unittest.TestCase):

    
    def test_barycentric_interp(self, verbose=False):
        """ Test Barycentric interpolation by evaluation of 
        known function on equidistant grid."""

        N = 32
        x_i = chebyshev_collocation_points_1st_kind(N)
        w_i = chebyshev_barycentric_weights_1st_kind(N)

        f = lambda x : np.sin(4 * np.pi * x)

        f_i = f(x_i)

        x_j = np.linspace(-1, 1, num=1000)
        f_j = f(x_j)

        f_j_interp = barycentric_chebyshev_interpolation(x_j, x_i, f_i, w_i)

        if verbose:
            print(f'f_j = \n{f_j}')
            print(f'f_j_interp = \n{f_j_interp}')
            print(f'diff = {np.max(np.abs(f_j - f_j_interp))}')

        self.assertTrue(np.allclose(f_j, f_j_interp))


    def test_barycentric_interp_on_grid(self, verbose=False):
        """The Barycentric interpolation diverges when the 
        interpolation point is equal to a Chebyshev collocation grid point.

        Here we check that this special case is correctly handled."""

        N = 32
        x_i = chebyshev_collocation_points_1st_kind(N)
        w_i = chebyshev_barycentric_weights_1st_kind(N)

        f = lambda x : np.sin(4 * np.pi * x)

        f_i = f(x_i)

        x_j = np.concatenate((x_i[::8]*0.99, [x_i[0], x_i[-1]], x_i[1::4]*0.99))
        f_j = f(x_j)

        f_j_interp = barycentric_chebyshev_interpolation(x_j, x_i, f_i, w_i)

        if verbose:
            print(f'f_j = \n{f_j}')
            print(f'f_j_interp = \n{f_j_interp}')
            print(f'diff = {np.max(np.abs(f_j - f_j_interp))}')

        self.assertTrue(np.allclose(f_j, f_j_interp))


if __name__ == '__main__':
    unittest.main()
