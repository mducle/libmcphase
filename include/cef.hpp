/* cef.hpp
 * 
 * A rewrite of the crystal field module of McPhase (cfield / so1ion) originally by Peter Fabi (nee Hoffmann).
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#ifndef CEF_H
#define CEF_H

#include "cfpars.hpp"
#include "eigen.hpp"

using namespace libMcPhase;

//namespace libMcPhase {

RowMatrixXd cef_hamiltonian(cfpars &pars);

//} // namespace libMcPhase

#endif
