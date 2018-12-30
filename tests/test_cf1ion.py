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

    ham_refi = np.array([[         0, -1.1372e+00,  5.7374e-01, -5.6652e-01,  2.2417e-01, -5.4185e-01, -1.2659e+00,           0,           0],
                        [ 1.1372e+00,           0,  3.4309e-01,  2.0210e-01,  6.7990e-01, -1.4629e+00, -2.7092e-01, -1.6747e+00,           0],
                        [-5.7374e-01, -3.4309e-01,           0,  7.2829e-01, -1.1407e+00,  5.1338e-01, -2.3456e+00,  2.7092e-01, -1.2659e+00],
                        [ 5.6652e-01, -2.0210e-01, -7.2829e-01,           0,  3.2722e-01, -1.8480e+00, -5.1338e-01, -1.4629e+00,  5.4185e-01],
                        [-2.2417e-01, -6.7990e-01,  1.1407e+00, -3.2722e-01,           0, -3.2722e-01, -1.1407e+00, -6.7990e-01,  2.2417e-01],
                        [ 5.4185e-01,  1.4629e+00, -5.1338e-01,  1.8480e+00,  3.2722e-01,           0, -7.2829e-01,  2.0210e-01,  5.6652e-01],
                        [ 1.2659e+00,  2.7092e-01,  2.3456e+00,  5.1338e-01,  1.1407e+00,  7.2829e-01,           0, -3.4309e-01,  5.7374e-01],
                        [          0,  1.6747e+00, -2.7092e-01,  1.4629e+00,  6.7990e-01, -2.0210e-01,  3.4309e-01,           0,  1.1372e+00],
                        [          0,           0,  1.2659e+00, -5.4185e-01, -2.2417e-01, -5.6652e-01, -5.7374e-01, -1.1372e+00,           0]])

    # Eigenvalues from Saficf
    eval_refr = np.array([-4.0496e+00, -3.4501e+00, -1.8251e+00, -1.7653e+00,  7.7521e-01, 1.3800e+00,  1.6579e+00,  3.5307e+00,  3.7464e+00])
    eval_refi = np.array([-4.1784e+00, -3.2256e+00, -1.2661e+00, -3.2787e-01, -5.7147e-16, 3.2787e-01,  1.2661e+00,  3.2256e+00,  4.1784e+00])
    eval_reff = np.array([-5.8293e+00, -3.9243e+00, -2.8666e+00, -1.3448e+00, -9.6973e-01, 1.6797e+00,  2.8187e+00,  5.0176e+00,  5.4188e+00])

    # Eigenvalues from cfield
    #results-full/levels.cef:#! Eigenvalues= 0 1.90492 2.96261 4.48442 4.85954 7.50893 8.64798 10.8469 11.2481  meV   
    #results-imag/levels.cef:#! Eigenvalues= 0 0.952819 2.9123 3.85057 4.17843 4.5063 5.44457 7.40405 8.35687  meV   
    #results-real/levels.cef:#! Eigenvalues= 0 0.599541 2.22457 2.28429 4.82484 5.42965 5.70748 7.58031 7.796  meV
    en_reff = np.array([0, 1.90492, 2.96261, 4.48442, 4.85954, 7.50893, 8.64798, 10.8469, 11.2481])
    en_refi = np.array([0, 0.952819, 2.9123, 3.85057, 4.17843, 4.5063, 5.44457, 7.40405, 8.35687])
    en_refr = np.array([0, 0.599541, 2.22457, 2.28429, 4.82484, 5.42965, 5.70748, 7.58031, 7.796])

    # For Wybourne normalisation
    en_refw = np.array([0, 0, 0.0464263, 0.0495632, 0.0742188, 0.0875256, 0.0949677, 0.11971, 0.11971])

    def test_creation_and_diagonalisation(self):
        with self.assertRaises(RuntimeError):   # J must be integer or half-integer
            cfp = libMcPhase.cf1ion(2.3)
        with self.assertRaises(RuntimeError):   # Unknown ion name (must be a rare earth 3+)
            cfp = libMcPhase.cf1ion('Pd3+')
        # Random parameters, compared to results from Saficf and cfield
        cfp = libMcPhase.cf1ion('Pr3+', B20=9.1381e-04, B21=1.0579e-01, B22=1.0386e-01,
                                B40=-7.2952e-04, B41=6.4828e-03, B42=-3.3059e-03, B43=1.5821e-02, B44=-4.5648e-03,
                                B60=-2.1369e-05, B61=4.3678e-04, B62=8.8995e-05, B63=1.0380e-04, B64=4.7701e-04, B65=-1.1485e-03, B66=3.9818e-04)
        nptest.assert_allclose(cfp.hamiltonian(), self.ham_ref, rtol=1e-4)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E, self.eval_refr, atol=1e-4)
        nptest.assert_allclose(E - E[0], self.en_refr, atol=1e-4)
        nptest.assert_allclose(V @ np.diag(E) @ V.conj().T, cfp.hamiltonian(), atol=1e-4)

    def test_imaginary_parameters(self):
        cfp = libMcPhase.cf1ion('Pr3+', B22S=-6.6761e-02, B21S=2.3224e-02, 
                                B44S=-8.1745e-03, B43S=-2.3769e-03, B42S=4.8620e-03, B41S=6.6729e-03, 
                                B66S=-6.6456e-04, B65S=8.0453e-04, B64S=3.4691e-04, B63S=5.1962e-04, B62S=-1.2112e-04, B61S=3.5266e-05)
        nptest.assert_allclose(cfp.hamiltonian(), 1j*self.ham_refi, rtol=1e-4)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E, self.eval_refi, atol=1e-4)
        nptest.assert_allclose(E - E[0], self.en_refi, atol=1e-4)
        nptest.assert_allclose(V @ np.diag(E) @ V.conj().T, cfp.hamiltonian(), atol=1e-4)

    def test_all_parameters(self):
        cfp = libMcPhase.cf1ion('Pr3+', B20=9.1381e-04, B21=1.0579e-01, B22=1.0386e-01,
                                B40=-7.2952e-04, B41=6.4828e-03, B42=-3.3059e-03, B43=1.5821e-02, B44=-4.5648e-03,
                                B60=-2.1369e-05, B61=4.3678e-04, B62=8.8995e-05, B63=1.0380e-04, B64=4.7701e-04, B65=-1.1485e-03, B66=3.9818e-04,
                                B22S=-6.6761e-02, B21S=2.3224e-02, 
                                B44S=-8.1745e-03, B43S=-2.3769e-03, B42S=4.8620e-03, B41S=6.6729e-03, 
                                B66S=-6.6456e-04, B65S=8.0453e-04, B64S=3.4691e-04, B63S=5.1962e-04, B62S=-1.2112e-04, B61S=3.5266e-05)
        nptest.assert_allclose(cfp.hamiltonian(), self.ham_ref + 1j*self.ham_refi, rtol=1e-4)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E, self.eval_reff, atol=1e-4)
        nptest.assert_allclose(E - E[0], self.en_reff, atol=1e-4)
        nptest.assert_allclose(V @ np.diag(E) @ V.conj().T, cfp.hamiltonian(), atol=1e-4)

    def test_wybourne(self):
        cfp = libMcPhase.cf1ion('Pr3+', type='Llm', B20=9.1381e-04, B21=1.0579e-01, B22=1.0386e-01,
                                B40=-7.2952e-04, B41=6.4828e-03, B42=-3.3059e-03, B43=1.5821e-02, B44=-4.5648e-03,
                                B60=-2.1369e-05, B61=4.3678e-04, B62=8.8995e-05, B63=1.0380e-04, B64=4.7701e-04, B65=-1.1485e-03, B66=3.9818e-04,
                                B22S=-6.6761e-02, B21S=2.3224e-02, 
                                B44S=-8.1745e-03, B43S=-2.3769e-03, B42S=4.8620e-03, B41S=6.6729e-03, 
                                B66S=-6.6456e-04, B65S=8.0453e-04, B64S=3.4691e-04, B63S=5.1962e-04, B62S=-1.2112e-04, B61S=3.5266e-05)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E - E[0], self.en_refw, atol=1e-3)
