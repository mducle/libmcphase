/* cef.cpp
 * 
 * A rewrite of the crystal field module of McPhase (cfield / so1ion) originally by Peter Fabi (nee Hoffmann).
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "cef.hpp"

RowMatrixXd cef_hamiltonian(cfpars &pars)
{
}


int main()
{
  cfpars pars;
  pars.B20 = 0.2;
  pars.B22 = -0.05;
  pars.B40 = 0.01;
  pars.B42 = 0.002;
  pars.B44 = -0.001;
  RowMatrixXd Hcf = cef_hamiltonian(pars);
  std::cout << Hcf << std::endl;
}
