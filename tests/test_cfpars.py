import libMcPhase
import unittest

class cfparsTests(unittest.TestCase):

    def test_creation(self):
        with self.assertRaises(RuntimeError):   # J must be integer or half-integer
            cfp = libMcPhase.cfpars(2.3)
        with self.assertRaises(RuntimeError):   # Unknown ion name (must be a rare earth 3+)
            cfp = libMcPhase.cfpars('Pd3+')
        cfp = libMcPhase.cfpars('Pr3+', B20=0.1)
        self.assertEqual(cfp.B20, 0.1)
        self.assertEqual(cfp.ion.lower(), 'pr3+')
        self.assertEqual(cfp.J, 4)
        self.assertEqual(cfp.alpha, -1.0 * 2*2*13/3/3/5/5/11)
        self.assertEqual(cfp.beta, -1.0 * 2*2/3/3/5/11/11)
        self.assertEqual(cfp.gamma, 1.0 * 2*2*2*2*17/3/3/3/3/5/7/11/11/13)
        self.assertEqual(cfp.B22, 0.)
        cfp = libMcPhase.cfpars('Pr3+', type='Llm', B20=0.1)
        self.assertEqual(cfp.type, cfp.Types.Llm)
        with self.assertRaises(RuntimeError):   # cannot set type without ionname
            cfp = libMcPhase.cfpars(2.5, type='Llm', B20=0.1)
        cfp = libMcPhase.cfpars(2.5, unit='cm', B20=0.1)
        self.assertEqual(cfp.ion, '')
        self.assertEqual(cfp.J, 2.5)
        self.assertEqual(cfp.alpha, 1.)
        self.assertEqual(cfp.unit, cfp.Units.cm)

    def test_units(self):
        # conversion factors from: https://www.physics.nist.gov/cuu/Constants/energy.html
        cfp = libMcPhase.cfpars(2.5, unit='meV', B20=1)
        self.assertEqual(cfp.B20, 1)
        cfp.unit = 'cm'
        self.assertAlmostEqual(cfp.B20, 8.065544005)
        cfp.unit = 'K'
        self.assertAlmostEqual(cfp.B20, 11.6045221)
        cfp.unit = 'meV'
        self.assertAlmostEqual(cfp.B20, 1, places=15)
        cfp = libMcPhase.cfpars(2.5, unit='cm', B20=1)
        self.assertEqual(cfp.B20, 1)
        cfp.unit = 'meV'
        self.assertAlmostEqual(cfp.B20, 0.12398419739)
        cfp.unit = 'K'
        self.assertAlmostEqual(cfp.B20, 1.43877736)
        cfp.unit = 'cm'
        self.assertAlmostEqual(cfp.B20, 1, places=15)
        cfp = libMcPhase.cfpars(2.5, unit='K', B20=1)
        self.assertEqual(cfp.B20, 1)
        cfp.unit = 'meV'
        self.assertAlmostEqual(cfp.B20, 0.086173303)
        cfp.unit = 'cm'
        self.assertAlmostEqual(cfp.B20, 0.69503457)
        cfp.unit = 'K'
        self.assertAlmostEqual(cfp.B20, 1, places=15)

    def test_types(self):
        par0 = [0.1, -0.01, 0.001]
        cfp = libMcPhase.cfpars('Nd3+', type='Blm', B20=par0[0], B40=par0[1], B60=par0[2])
        self.assertEqual(cfp.B20, par0[0])
        self.assertEqual(cfp.B40, par0[1])
        self.assertEqual(cfp.B60, par0[2])
        # Conversion constants for Nd3+
        theta_k = [-1.0 * 7/3/3/11/11, -1.0 * 2*2*2*17/3/3/3/11/11/11/13, -1.0 * 5*17*19/3/3/3/7/11/11/11/13/13]
        rk = [1.114, 2.910, 15.03]
        lamb = [0.5, 0.125, 0.0625]
        cfp.type = 'ARlm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii] / theta_k[ii])
        cfp.type = 'Alm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii] / theta_k[ii] / rk[ii])
        cfp.type = 'Vlm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii])
        cfp.type = 'Wlm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii] / rk[ii])
        cfp.type = 'Llm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii] / lamb[ii] / theta_k[ii])
        cfp.type = 'Blm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertEqual(par, par0[ii])
        # Try starting from Wybourne parameters
        cfp = libMcPhase.cfpars('Nd3+', type='Llm', B20=par0[0], B40=par0[1], B60=par0[2])
        cfp.type = 'Blm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii] * lamb[ii] * theta_k[ii])
        cfp.type = 'ARlm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii] * lamb[ii])
        cfp.type = 'Alm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii] * lamb[ii] / rk[ii])
        cfp.type = 'Vlm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii] * lamb[ii] * theta_k[ii])
        cfp.type = 'Wlm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii] * lamb[ii] * theta_k[ii] / rk[ii])
        cfp.type = 'Llm'
        for ii, par in enumerate([cfp.B20, cfp.B40, cfp.B60]):
            self.assertAlmostEqual(par, par0[ii], places=15)

    def test_ion(self):
        cfp = libMcPhase.cfpars('Nd3+', B20 = 1)
        self.assertEqual(cfp.ion.lower(), 'nd3+')
        self.assertEqual(cfp.B20, 1)
        cfp.ion = 'Pr3+'
        self.assertEqual(cfp.ion.lower(), 'pr3+')
        self.assertAlmostEqual(cfp.B20, (-1.0 * 2*2*13/3/3/5/5/11) / (-1.0 * 7/3/3/11/11))
        cfp.J = 4
        self.assertEqual(cfp.ion, '')
        self.assertEqual(cfp.type, cfp.Types.Blm)
        with self.assertRaises(RuntimeError):   # lost all ion info, cannot know ion dependent factors
            cfp.type = 'Llm'
        self.assertAlmostEqual(cfp.B20, (-1.0 * 2*2*13/3/3/5/5/11) / (-1.0 * 7/3/3/11/11))

    def test_conversions(self):
        cfp = libMcPhase.cfpars('Nd3+', type='Llm', unit='K', B20=10)
        self.assertEqual(cfp.ion.lower(), 'nd3+')
        self.assertEqual(cfp.unit, cfp.Units.K)
        self.assertEqual(cfp.type, cfp.Types.Llm)
        self.assertEqual(cfp.B20, 10)
        cfp.unit = 'cm'
        self.assertEqual(cfp.unit, cfp.Units.cm)
        cfp.type = 'Alm'
        self.assertEqual(cfp.type, cfp.Types.Alm)
        self.assertAlmostEqual(cfp.B20, 10 * 0.69503457 * 0.5 / 1.114)
        cfp.unit = 'meV'
        self.assertEqual(cfp.unit, cfp.Units.meV)
        self.assertAlmostEqual(cfp.B20, 10 * 0.086173303 * 0.5 / 1.114)
        cfp.ion = 'Pr3+'
        self.assertEqual(cfp.ion.lower(), 'pr3+')
        self.assertAlmostEqual(cfp.B20, 10 * 0.086173303 * 0.5 / 1.1963)
        cfp.type = 'Blm'
        self.assertEqual(cfp.type, cfp.Types.Blm)
        self.assertAlmostEqual(cfp.B20, 10 * 0.086173303 * 0.5 * (-1.0 * 2*2*13/3/3/5/5/11))
