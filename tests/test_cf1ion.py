import libMcPhase
import unittest
import numpy as np
import numpy.testing as nptest

class cf1ionTests(unittest.TestCase):

    ham_ref = np.array([[-6.9491e-01, -1.8754e+00, -6.8087e-02, -1.0278e+00,  9.7844e-01,  7.7351e-01,  7.5851e-01,           0,           0],
                        [-1.8754e+00,  1.3833e+00,  8.2511e-01,  3.0767e-01, -7.3185e-01, -9.5181e-01,  3.8676e-01,  1.0034e+00,           0],
                        [-6.8087e-02,  8.2511e-01, -1.1818e-01, -9.5223e-02,  1.3303e+00, -2.4313e-01, -2.0237e+00, -3.8676e-01,  7.5851e-01],
                        [-1.0278e+00,  3.0767e-01, -9.5223e-02, -4.3640e-01, -3.4219e-01,  1.8579e+00,  2.4313e-01, -9.5181e-01, -7.7351e-01],
                        [ 9.7844e-01, -7.3185e-01,  1.3303e+00, -3.4219e-01, -2.6766e-01,  3.4219e-01,  1.3303e+00,  7.3185e-01,  9.7844e-01],
                        [ 7.7351e-01, -9.5181e-01, -2.4313e-01,  1.8579e+00,  3.4219e-01, -4.3640e-01,  9.5223e-02,  3.0767e-01,  1.0278e+00],
                        [ 7.5851e-01,  3.8676e-01, -2.0237e+00,  2.4313e-01,  1.3303e+00,  9.5223e-02, -1.1818e-01, -8.2511e-01, -6.8087e-02],
                        [          0,  1.0034e+00, -3.8676e-01, -9.5181e-01,  7.3185e-01,  3.0767e-01, -8.2511e-01,  1.3833e+00,  1.8754e+00],
                        [          0,           0,  7.5851e-01, -7.7351e-01,  9.7844e-01,  1.0278e+00, -6.8087e-02,  1.8754e+00, -6.9491e-01]])

    def test_creation(self):
        with self.assertRaises(RuntimeError):   # J must be integer or half-integer
            cfp = libMcPhase.cf1ion(2.3)
        with self.assertRaises(RuntimeError):   # Unknown ion name (must be a rare earth 3+)
            cfp = libMcPhase.cf1ion('Pd3+')
        # Random parameters, compared to results from Saficf and cfield
        cfp = libMcPhase.cf1ion('Pr3+', B20=9.1381e-04, B21=1.0579e-01, B22=1.0386e-01,
                                B40=-7.2952e-04, B41=6.4828e-03, B42=-3.3059e-03, B43=1.5821e-02, B44=-4.5648e-03,
                                B60=-2.1369e-05, B61=4.3678e-04, B62=8.8995e-05, B63=1.0380e-04, B64=4.7701e-04, B65=-1.1485e-03, B66=3.9818e-04)
        nptest.assert_allclose(cfp.hamiltonian(), self.ham_ref, rtol=1e-4)
