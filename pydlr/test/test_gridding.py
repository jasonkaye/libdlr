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

        for lamb in [10., 40., 80., 160., 320., 640.]:
            kmat, t, w, err = kernel_discretization(lamb, error_est=True)
            print(f'err = {err:2.2E}')
            self.assertTrue(err < 1e-14)
        

if __name__ == '__main__':
    unittest.main()
    
