/* physprop.hpp
 *
 * Header file for the physical properties base class
 *
 * (c) 2024 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the COPYING file
 */

#pragma once

namespace libMcPhase {

// Basic physical constants (needs cm constants as internal energy units in ic1ion is cm)
static const double K_B = 0.08617343183;       // meV/K - Boltzmann constant
static const double MU_B = 0.0578838263;       // meV/T - Bohr magneton
static const double K_Bc = 0.6950348004;       // cm/K - Boltzmann constant
static const double MU_Bc = 0.46686447783;     // cm/T - Bohr magneton

// Conversion factors for different magnetic units for magnetisation. Order: [bohr, cgs, SI].
// NAMUB is N_A * MU_B in J/T/mol == Am^2/mol is the SI unit. 
// The cgs unit is N_A * MU_B in erg/G/mol == emu/mol is different only by a factor of 1000 larger
static const double NAMUB = 5.5849397;  // N_A*mu_B - J/T/mol - product of Bohr magneton and Avogadro's number
static const std::array<double, 3> MAGCONV = {1., NAMUB*1000, NAMUB};

// Note these constants are strange because the default energy unit in this module is cm-1
static const double NAMUBSQ_ERG = 0.26074098;     // N_A * muB[erg/G] * muB[cm/T] * [Tesla_to_Gauss=1e-4]
static const double NAMUBSQ_JOULE = 3.276568e-06; // N_A * muB[J/T] * muB[cm/T] * mu0
// Factor of mu0 is needed in SI due to different definition of the magnetisation, B and H fields

// Conversion factors for different magnetic units for magnetic susceptibility. Order: [bohr, cgs, SI].
// The susceptibility prefactor is (in principle) N_A * MU_B^2 but we need to account for various units...
// The atomic (bohr) susceptibility is in uB/T/ion; cgs is in erg/G^2/mol==cm^3/mol; SI in J/T^2/mol==m^3/mol
// Note that chi_SI = (4pi*10^-6)chi_cgs [*not* 4pi*10-7!]
static const std::array<double, 3> SUSCCONV = {MU_B, NAMUBSQ_ERG, NAMUBSQ_JOULE};

// Conversion factor for heat capacity calculations
static const double NAMEV = 96.48533212;          // J/mol = N_A * meV

// EPSILON to determine if energy levels are degenerate or not
static const double DELTA_EPS = 1e-6;

} // namespace libMcPhase
