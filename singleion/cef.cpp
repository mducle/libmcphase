/* cef.cpp
 * 
 * A rewrite of the crystal field module of McPhase (cfield / so1ion) originally by Peter Fabi (nee Hoffmann).
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "cef.hpp"
#include <iostream>

RowMatrixXd cef_hamiltonian(cfpars &pars)
{
}


int main()
{
  cfpars pars;
  pars.set(cfpars::Blm::B20, 0.2);
  pars.set(cfpars::Blm::B22, -0.05);
  pars.set(cfpars::Blm::B40, 0.01);
  pars.set(cfpars::Blm::B42, 0.002);
  pars.set(cfpars::Blm::B44, -0.001);
  RowMatrixXd Hcf = cef_hamiltonian(pars);
  std::cout << Hcf << std::endl;
}
