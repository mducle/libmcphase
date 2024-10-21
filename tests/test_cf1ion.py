import libmcphase
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

    # For different normalisation (cfield)
    #results-akq/levels.cef:#! Eigenvalues= 0 5.70236e-06 0.0218494 0.0226182 0.034851 0.0453383 0.047894 0.0758639 0.0760302  meV   
    #results-lkq/levels.cef:#! Eigenvalues= 0 0.000225895 0.0464263 0.0495632 0.0742188 0.0875256 0.0949677 0.11971 0.120501  meV   
    en_ref_a = np.array([0, 5.70236e-06, 0.0218494, 0.0226182, 0.034851, 0.0453383, 0.047894, 0.0758639, 0.0760302])   # Stevens Alm (orth only)
    en_ref_l = np.array([0, 0.000225895, 0.0464263, 0.0495632, 0.0742188, 0.0875256, 0.0949677, 0.11971, 0.120501])    # Wybourne Llm

    pars_orth = {'B20':9.1381e-04, 'B22':1.0386e-01, 
                 'B40':-7.2952e-04, 'B42':-3.3059e-03, 'B44':-4.5648e-03, 
                 'B60':-2.1369e-05, 'B62':8.8995e-05, 'B64':4.7701e-04, 'B66':3.9818e-04}
    pars_real = dict(pars_orth, **{'B21':1.0579e-01, 'B41':6.4828e-03, 'B43':1.5821e-02, 'B61':4.3678e-04, 'B63':1.0380e-04, 'B65':-1.1485e-03})
    pars_imag = {'B22S':-6.6761e-02, 'B21S':2.3224e-02, 
                 'B44S':-8.1745e-03, 'B43S':-2.3769e-03, 'B42S':4.8620e-03, 'B41S':6.6729e-03, 
                 'B66S':-6.6456e-04, 'B65S':8.0453e-04, 'B64S':3.4691e-04, 'B63S':5.1962e-04, 'B62S':-1.2112e-04, 'B61S':3.5266e-05}
    all_pars = dict(pars_real, **pars_imag)

    # Parameters for physical properties check
    pp_cfpars = {'B20':0.37737, 'B22':3.9770, 'B40':-0.031787, 'B42':-0.11611, 'B44':-0.12544}

    def test_creation_and_diagonalisation(self):
        with self.assertRaises(RuntimeError):   # J must be integer or half-integer
            cfp = libmcphase.cf1ion(2.3)
        with self.assertRaises(RuntimeError):   # Unknown ion name (must be a rare earth 3+)
            cfp = libmcphase.cf1ion('Pd3+')
        # Random parameters, compared to results from Saficf and cfield
        cfp = libmcphase.cf1ion('Pr3+', **self.pars_real)
        nptest.assert_allclose(cfp.hamiltonian(), self.ham_ref, rtol=1e-4)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E, self.eval_refr, atol=1e-4)
        nptest.assert_allclose(E - E[0], self.en_refr, atol=1e-4)
        nptest.assert_allclose(V @ np.diag(E) @ V.conj().T, cfp.hamiltonian(), atol=1e-4)

    def test_imaginary_parameters(self):
        cfp = libmcphase.cf1ion('Pr3+', **self.pars_imag)
        nptest.assert_allclose(cfp.hamiltonian(), 1j*self.ham_refi, rtol=1e-4)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E, self.eval_refi, atol=1e-4)
        nptest.assert_allclose(E - E[0], self.en_refi, atol=1e-4)
        nptest.assert_allclose(V @ np.diag(E) @ V.conj().T, cfp.hamiltonian(), atol=1e-4)

    def test_all_parameters(self):
        cfp = libmcphase.cf1ion('Pr3+', **self.all_pars)
        nptest.assert_allclose(cfp.hamiltonian(), self.ham_ref + 1j*self.ham_refi, rtol=1e-4)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E, self.eval_reff, atol=1e-4)
        nptest.assert_allclose(E - E[0], self.en_reff, atol=1e-4)
        nptest.assert_allclose(V @ np.diag(E) @ V.conj().T, cfp.hamiltonian(), atol=1e-4)

    def test_types(self):
        cfp = libmcphase.cf1ion('Pr3+', type='Alm', **self.pars_orth)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E - E[0], self.en_ref_a, atol=1e-6)
        cfp = libmcphase.cf1ion('Pr3+', type='Llm', **self.all_pars)
        V, E = cfp.eigensystem()
        nptest.assert_allclose(E - E[0], self.en_ref_l, atol=1e-4)

    def test_heatcapacity(self):
        # Compared to values from Mantid
        cf = libmcphase.cf1ion('Ce3+', **self.pp_cfpars)
        # Test Heat capacity calculations
        Cv = cf.heatcapacity(np.linspace(1,300,300))
        #self.assertAlmostEqual(TCv[150], 151, 4)
        self.assertAlmostEqual(Cv[100], 4.2264, 3)
        self.assertAlmostEqual(Cv[150], 5.9218, 3)
        self.assertAlmostEqual(Cv[200], 5.4599, 3)
        cf = libmcphase.cf1ion('Ce3+', **{k:v*11.6045 for k,v in self.pp_cfpars.items()}, unit='K')
        Cv = cf.heatcapacity(np.linspace(1,300,300))
        self.assertAlmostEqual(Cv[100], 4.2264, 3)
        self.assertAlmostEqual(Cv[150], 5.9218, 3)
        self.assertAlmostEqual(Cv[200], 5.4599, 3)

    def test_susceptibility(self):
        # Test susceptibility calculations
        cf = libmcphase.cf1ion('Ce3+', **self.pp_cfpars)
        #Tchi_powder, chi_powder = cf.getSusceptibility(np.linspace(1, 300, 50), Hdir="powder")
        tt = np.linspace(1, 300, 50)
        cc = [np.array(cf.susceptibility(tt, hdir, 'cgs')) for hdir in [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]]
        chi_powder = (cc[0] + cc[1] + cc[2]) / 3
        #self.assertAlmostEqual(Tchi_powder[10], 62.02, 2)
        self.assertAlmostEqual(chi_powder[5], 1.92026e-2, 6)
        self.assertAlmostEqual(chi_powder[10], 1.03471e-2, 6)
        self.assertAlmostEqual(chi_powder[15], 0.73004e-2, 6)

    def test_magnetisation_vs_T(self):
        # Test M(T) calculations
        cf = libmcphase.cf1ion('Ce3+', **self.pp_cfpars)
        tt = np.linspace(1, 300, 50)
        #Tmt_powder, mt_powder = cf.getMagneticMoment(1.0, Temperature=np.linspace(1, 300, 50), Hdir="powder", Unit="cgs")
        cc = [np.array(cf.susceptibility(tt, hdir, 'cgs')) for hdir in [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]]
        mm = [np.array(cf.magnetisation([1.0], hdir, tt, 'cgs')) for hdir in [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]]
        chi_powder = (cc[0] + cc[1] + cc[2]) / 3
        mt_powder = np.squeeze(mm[0] + mm[1] + mm[2]) / 3
        self.assertAlmostEqual(chi_powder[5], mt_powder[5], 6)
        self.assertAlmostEqual(chi_powder[10], mt_powder[10], 6)
        self.assertAlmostEqual(chi_powder[15], mt_powder[15], 6)
        #_, invmt_powder_SI = cf.getMagneticMoment(1.0, Temperature=np.linspace(1, 300, 50), Hdir="powder", Unit="SI", Inverse=True)
        #self.assertAlmostEqual(chi_powder[5] * 10, 1 / invmt_powder_SI[5], 2)
        #self.assertAlmostEqual(chi_powder[10] * 10, 1 / invmt_powder_SI[10], 2)
        #self.assertAlmostEqual(chi_powder[15] * 10, 1 / invmt_powder_SI[15], 2)

        # Test different Hmag
        #_, h_mag_10 = cf.getMagneticMoment(Hmag=10, Temperature=np.linspace(1, 300, 50), Hdir="powder", Unit="bohr")
        mm = [np.array(cf.magnetisation([10.], hdir, tt, 'bohr')) for hdir in [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]]
        h_mag_10 = np.squeeze(mm[0] + mm[1] + mm[2]) / 3
        self.assertAlmostEqual(h_mag_10[5], 0.323607, 5)
        self.assertAlmostEqual(h_mag_10[10], 0.182484, 5)
        self.assertAlmostEqual(h_mag_10[15], 0.129909, 5)
        #_, h_mag_5 = cf.getMagneticMoment(Hmag=5, Temperature=np.linspace(1, 300, 50), Hdir="powder", Unit="bohr")
        mm = [np.array(cf.magnetisation([5.], hdir, tt, 'bohr')) for hdir in [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]]
        h_mag_5 = np.squeeze(mm[0] + mm[1] + mm[2]) / 3
        self.assertAlmostEqual(h_mag_5[5], 0.16923426, 6)
        self.assertAlmostEqual(h_mag_5[10], 0.09228022, 6)
        self.assertAlmostEqual(h_mag_5[15], 0.06525625, 6)

    def test_magnetisation(self):
        # Test M(H) calculations
        cf = libmcphase.cf1ion('Ce3+', **self.pp_cfpars)
        #Hmag_SI, mag_SI = cf.getMagneticMoment(np.linspace(0, 30, 15), Temperature=10, Hdir=[0, 1, -1], Unit="SI")
        mag_SI = np.squeeze(cf.magnetisation(np.linspace(0, 30, 15), [0, 1, -1], [10], 'SI'))
        self.assertAlmostEqual(mag_SI[1], 1.8139, 3)
        self.assertAlmostEqual(mag_SI[5], 6.7859, 3)
        self.assertAlmostEqual(mag_SI[9], 8.2705, 3)
        #_, mag_bohr = cf.getMagneticMoment(np.linspace(0, 30, 15), Temperature=10, Hdir=[0, 1, -1], Unit="bohr")
        mag_bohr = np.squeeze(cf.magnetisation(np.linspace(0, 30, 15), [0, 1, -1], [10], 'bohr'))
        self.assertAlmostEqual(mag_SI[1] / 5.5849, mag_bohr[1], 3)
        self.assertAlmostEqual(mag_SI[5] / 5.5849, mag_bohr[5], 3)
        self.assertAlmostEqual(mag_SI[9] / 5.5849, mag_bohr[9], 3)
