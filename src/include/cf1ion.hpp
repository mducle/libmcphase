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

    protected:
        bool m_ham_calc = false;                              // Flag to indicate if Hamiltonian calculated
        bool m_ev_calc = false;                               // Flag to indicate if eigenvectors/values calculated
        RowMatrixXcd m_hamiltonian;                           // Cached Hamiltonian
        RowMatrixXcd m_eigenvectors;                          // Cached eigenvectors
        VectorXd m_eigenvalues;                               // Cached eigenvalues

    public:
        // Setters
        virtual void set_unit(const Units newunit);
        virtual void set_type(const Type newtype);
        virtual void set_name(const std::string &ionname);
        virtual void set_J(const double J);
        virtual void set(const Blm blm, double val);
        virtual void set(int l, int m, double val);
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
