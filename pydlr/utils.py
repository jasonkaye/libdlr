""" Utilities functions used in testing and demnstration.

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

from scipy.integrate import quad

from pydlr import kernel

def analytic_bethe_G_tau(tau, beta):
    """Get Bethe graph imaginary time Green's function.

    Obtained by semi-analytic evaluation using adaptive integration of

    .. math:: G(\\tau) = -\\frac{2}{\\pi} \\int_{-1}^{1} K\\left(\\frac{\\tau}{\\beta}, \\beta \\omega\\right) \\sqrt{1 - \\omega^2} \, d\\omega  

    Parameters
    ----------

    tau : (n), ndarray
        Imaginary time :math:`\\tau` points to evaulate the Bethe Green's function on.

    beta : float
        Inverse temperature :math:`\\beta`

    Returns
    -------

    G_iaa : (n, 1, 1), ndarray
        Imaginary time Green's function :math:`G(\\tau)` evaluated on the given imaginary times.
    """
    
    I = lambda omega, tau: -2/np.pi * kernel(np.array([tau/beta]), np.array([beta*omega]))[0,0]
    integral = lambda tau: quad(I, -1, 1, weight='alg', wvar=(0.5, 0.5), args=tau)[0]
    return np.vectorize(integral)(tau).reshape((len(tau), 1, 1))
