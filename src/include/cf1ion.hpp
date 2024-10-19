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
#include "physprop.hpp"

namespace libMcPhase {

class cf1ion: public cfpars, public physprop {

    protected:
        bool m_ham_calc = false;                              // Flag to indicate if Hamiltonian calculated
        bool m_upper = true;                                  // Flag to indicate if upper triangle of Ham is calc
        bool m_ev_calc = false;                               // Flag to indicate if eigenvectors/values calculated
        RowMatrixXcd m_hamiltonian;                           // Cached Hamiltonian
        RowMatrixXcd m_eigenvectors;                          // Cached eigenvectors
        VectorXd m_eigenvalues;                               // Cached eigenvalues
        RowMatrixXcd _hamiltonian(bool upper=true);
        void fill_upper();

    public:
        // Setters
        void set_unit(const Units newunit);
        void set_type(const Type newtype);
        void set_name(const std::string &ionname);
        void set_J(const double J);
        void set(const Blm blm, double val);
        void set(int l, int m, double val);
        // Constructors
        cf1ion() : cfpars() {};
        cf1ion(const int J2) : cfpars(J2) {};
        cf1ion(const double J) : cfpars(J) {};
        cf1ion(const std::string &ionname) : cfpars(ionname) {};
        // Methods
        RowMatrixXcd hamiltonian();
        std::tuple<RowMatrixXcd, VectorXd> eigensystem();
        RowMatrixXcd zeeman_hamiltonian(double H, std::vector<double> Hdir);
        std::vector<RowMatrixXcd> calculate_moments_matrix(RowMatrixXcd ev);

}; // class cf1ion

} // namespace libMcPhase

#endif
