"""Test least squares fit of DLR coefficients using arbitrary tau data.

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


import numpy as np

from pydlr import dlr
from pydlr import kernel

def test_lstsq_dlr_tau(verbose=False):

    d = dlr(lamb=100., eps=1e-10)

    beta = 10.

    tau_i = np.linspace(0, beta, num=400)
    shape = (len(tau_i), 1, 1)

    k = 100
    w_k = 1.0 + np.random.randn(k)
    c_k = np.random.rand(k)
    c_k /= np.sum(c_k)

    G_iaa = - (kernel(tau_i/beta, beta*w_k) @ c_k).reshape(shape)
    np.testing.assert_almost_equal(-1., G_iaa[0,0,0] + G_iaa[-1,0,0])

    G_xaa = d.lstsq_dlr_from_tau(tau_i, G_iaa, beta)
    
    G_laa = d.tau_from_dlr(G_xaa)
    tau_l = d.get_tau(beta)

    G_iaa_ref = d.eval_dlr_tau(G_xaa, tau_i, beta)
    
    if verbose:
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 12))
        subp = [1, 1, 1]

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(tau_i, G_iaa[:,0,0], 'x-')
        plt.plot(tau_l, G_laa[:,0,0], '.')
        plt.plot(tau_i, G_iaa_ref[:,0,0], '+')
        
        plt.xlabel(r'$\tau$')
        plt.ylabel(r'$G(\tau)$')

        plt.show(); exit()

if __name__ == '__main__':

    test_lstsq_dlr_tau(verbose=True)
