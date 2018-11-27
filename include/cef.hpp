/* cef.hpp
 * 
 * A rewrite of the crystal field module of McPhase (cfield / so1ion) originally by Peter Fabi (nee Hoffmann).
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "cfpars.hpp"
#include "eigen.hpp"

RowMatrixXd cef_hamiltonian(cfpars &pars);
