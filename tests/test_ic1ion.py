import libmcphase
import unittest
import numpy as np
import numpy.testing as nptest

class ic1ionTests(unittest.TestCase):
    # Eigenvalues from McPhase 5.4
    en_refr = np.array([-1400.18, -1399.68, -1398.16, -1397.99, -1395.58, -1394.96, -1394.66, -1393, -1392.79, -1138.31, -1138.08, -1136.94, -1136.62,
                        -1134.23, -1134.15, -1133.97, -1133.61, -1133.59, -1132.48, -1132.38, -866.993, -866.41, -865.065, -864.368, -861.593, -861.056,
                        -861.049, -860.474, -859.914, -859.474, -858.926, -858.058, -856.903, -790.393, -789.612, -789.31, -789.161, -788.855, -616.178,
                        -615.661, -615.332, -614.225, -614.118, -614.007, -613.522, -558.719, -557.832, -557.34, -556.635, -556.411, -555.567, -555.043,
                        -554.843, -554.372, -192.949, -189.563, -188.162, -187.462, -186.38, -186.378, -184.65, -184.329, -181.205, 688.838, 692.736,
                        693.247, 694.493, 695.46, 1161.97, 1239.94, 1241.31, 1242.12, 1242.19, 1242.47, 1242.77, 1243.14, 1247.49, 1247.63, 1248.64,
                        1249.32, 1249.98, 1250.73, 1251.11, 1252.69, 1252.71, 1392.6, 1393.44, 1394.54, 1395.23, 1395.74, 4395.67])
    en_refi = np.array([-1401.95, -1400, -1399.35, -1397.6, -1397.35, -1394.78, -1393.62, -1391.6, -1391.2, -1138.68, -1138.51, -1137.76, -1137.14,
                        -1136.26, -1135.27, -1133.49, -1133.13, -1132.8, -1130.93, -1130.77, -867.828, -866.929, -865.631, -865.456, -864.899, -863.497,
                        -861.647, -859.534, -859.274, -858.897, -857.456, -855.493, -854.43, -790.395, -789.715, -789.491, -788.752, -788.509, -616.586,
                        -615.789, -615.547, -614.621, -613.821, -613.495, -612.842, -559.522, -558.528, -558.2, -557.096, -555.609, -555.275, -554.533,
                        -554.178, -553.457, -194.468, -191.231, -190.087, -188.514, -186.898, -184.581, -183.756, -181.359, -180.12, 687.351, 691.495,
                        694.211, 695.196, 696.31, 1161.97, 1238.78, 1238.85, 1239.82, 1241.21, 1242.38, 1243.71, 1244, 1246.92, 1247.8, 1248.35,
                        1249.18, 1250.33, 1251.35, 1252.72, 1254.41, 1254.74, 1392.24, 1393.39, 1394.75, 1395.45, 1395.9, 4395.7])

    pars_orth = {'B20':9.1381e-04, 'B22':1.0386e-01, 
                 'B40':-7.2952e-04, 'B42':-3.3059e-03, 'B44':-4.5648e-03, 
                 'B60':-2.1369e-05, 'B62':8.8995e-05, 'B64':4.7701e-04, 'B66':3.9818e-04}
    pars_real = dict(pars_orth, **{'B21':1.0579e-01, 'B41':6.4828e-03, 'B43':1.5821e-02, 'B61':4.3678e-04, 'B63':1.0380e-04, 'B65':-1.1485e-03})
    pars_imag = {'B22S':-6.6761e-02, 'B21S':2.3224e-02, 
                 'B44S':-8.1745e-03, 'B43S':-2.3769e-03, 'B42S':4.8620e-03, 'B41S':6.6729e-03, 
                 'B66S':-6.6456e-04, 'B65S':8.0453e-04, 'B64S':3.4691e-04, 'B63S':5.1962e-04, 'B62S':-1.2112e-04, 'B61S':3.5266e-05}
    all_pars = dict(pars_real, **pars_imag)

    def test_creation_and_diagonalisation(self):
        with self.assertRaises(RuntimeError):   # Can only be constructed from ion name
            cfp = libmcphase.ic1ion(2.5)
        with self.assertRaises(RuntimeError):   # Unknown ion name
            cfp = libmcphase.ic1ion('O2-')
        # Random parameters, compared to results from McPhase normal
        cfp = libmcphase.ic1ion('Pr3+', **self.pars_real)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E, self.en_refr, atol=1e-2)
        nptest.assert_allclose(V @ np.diag(E) @ V.conj().T, cfp.hamiltonian(), atol=1e-4)

    def test_all_parameters(self):
        cfp = libmcphase.ic1ion('Pr3+', **self.all_pars)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E, self.en_refi, atol=1e-2)
        nptest.assert_allclose(V @ np.diag(E) @ V.conj().T, cfp.hamiltonian(), atol=1e-4)

    def test_magnetisation(self):
        # Test compared to McPhase 5.4
        ref_mag = np.array([3.35996e-13, 0.299607, 0.585826, 0.848715, 1.0833, 1.28896, 1.46781, 1.62318, 1.75859, 1.87726, 1.98196])
        cfp = libmcphase.ic1ion('Pr3+', **self.all_pars)
        mag = cfp.magnetisation(np.arange(0, 21, 2), [1, 0, 0], 1, 'bohr')
        nptest.assert_allclose(mag, ref_mag, atol=1e-4)

    def test_susceptibility(self):
        # Test compared to McPhase 5.4 magnetisation vs T at 1T
        ref_sus = np.array([0.0377984, 0.0581256, 0.0392299, 0.0284567, 0.0221967, 0.0181698, 0.0153757, 0.0133276, 0.0117635, 0.0105307, 0.00953459])
        cfp = libmcphase.ic1ion('Pr3+', **self.pars_orth)
        sus = cfp.susceptibility(np.arange(1, 302, 30), [1, 0, 0], 'bohr')
        nptest.assert_allclose(sus, ref_sus, atol=1e-4)
        

