/* cf1ion.hpp
 * 
 * A class for calculating the crystal field Hamiltonian in Russell-Saunders (LS-) coupling.
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#ifndef CF1ION_H
#define CF1ION_H

#include "eigen.hpp"
#include "cfpars.hpp"

namespace libMcPhase {

class cf1ion: public cfpars {

    public:
        // Constructors
        cf1ion() : cfpars() {};
        cf1ion(const int J2) : cfpars(J2) {};
        cf1ion(const double J) : cfpars(J) {};
        cf1ion(const std::string &ionname) : cfpars(ionname) {};
        // Methods
        RowMatrixXcd hamiltonian(bool upper=true);
        std::tuple<RowMatrixXcd, VectorXd> eigensystem();

}; // class cf1ion

} // namespace libMcPhase

#endif
