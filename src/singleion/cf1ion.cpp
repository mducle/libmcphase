/* cf1ion.cpp
 * 
 * A class for calculating the crystal field Hamiltonian in Russell-Saunders (LS-) coupling.
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "cf1ion.hpp"

namespace libMcPhase {

static const std::array<std::array<int, 4>, 12> idq = { {{2,2,0,4}, {2,1,1,3}, {4,4,5,13}, {4,3,6,12}, {4,2,7,11}, {4,1,8,10}, 
                                                        {6,6,14,26}, {6,5,15,25}, {6,4,16,24}, {6,3,17,23}, {6,2,18,22}, {6,1,19,21}} };
static const std::array<std::array<int, 2>, 3> idq0 = { {{2,2}, {4,9}, {6,20}} };

// --------------------------------------------------------------------------------------------------------------- //
// Setter/getter methods for cfpars class
// --------------------------------------------------------------------------------------------------------------- //
void cf1ion::set(const Blm blm, double val) {
    cfpars::set(blm, val);
    m_ham_calc = false;
    m_ev_calc = false;
}

void cf1ion::set(int l, int m, double val) {
    cfpars::set(l, m, val);
    m_ham_calc = false;
    m_ev_calc = false;
}

void cf1ion::set_unit(cfpars::Units const newunit) {
    cfpars::set_unit(newunit);
    physprop::m_meVconv = m_econv;
    m_ham_calc = false;
    m_ev_calc = false;
}

void cf1ion::set_type(const cfpars::Type newtype) {
    cfpars::set_type(newtype);
    m_ham_calc = false;
    m_ev_calc = false;
}

void cf1ion::set_name(const std::string &ionname) {
    cfpars::set_name(ionname);
    m_ham_calc = false;
    m_ev_calc = false;
}

void cf1ion::set_J(const double J) {
    cfpars::set_J(J);
    m_ham_calc = false;
    m_ev_calc = false;
}

// --------------------------------------------------------------------------------------------------------------- //
// General methods for the cf1ion class
// --------------------------------------------------------------------------------------------------------------- //
void cf1ion::fill_upper() {
    int dimj = m_J2 + 1;
    for (int i=0; i<dimj; i++) {
        for (int j=i+1; j<dimj; j++) {
            m_hamiltonian(i,j) = std::conj(m_hamiltonian(j,i));
        }
    }
}

RowMatrixXcd cf1ion::_hamiltonian(bool upper) {
    if (m_ham_calc) {
        if (upper && !m_upper) {
            fill_upper();
        }
        return m_hamiltonian;
    }
    if (m_J2 <= 0) {
        throw std::runtime_error("Invalid value of J - must be strictly greater than zero");
    }

    int dimj = m_J2 + 1;
    m_hamiltonian = RowMatrixXcd::Zero(dimj, dimj);
    // For J=1/2, all Stevens operators are zero
    if (dimj == 2) {
        m_ham_calc = true;
        return m_hamiltonian;
    }

    // The crystal field potential operator is given by:
    //
    //         ---   k    [  k        q  k    ]     ---  k  k     ---   k [  k       q  k  ]
    // V   = i >    B     | O   - (-1)  O     |  +  >   B  O   +  >    B  | O  + (-1)  O   |
    //  cf     ---   -|q| [  |q|         -|q| ]     ---  0  0     ---   q [  q          -q ]
    //        k,q<0                                  k           k,q>0  
    //
    // where O^k_q are spherical harmonic tensor operators of rank k.
    // Reference: C. Rudowicz, J. Phys. C: Solid State Phys., vol 18, pp1415-1430 (1985).
    //
    // The terms in the square brackets above are the "tesseral harmonics" which means that
    // B^k_q parameters are strictly real.
    //
    // The energy matrix elements <LSJM_j|V_CF|L'SJ'M_j'> are given by summing
    // over the matrix elements of the tensor operator given by:
    //
    //          k                 J-M_j                           k
    // <LSJM | O  |L'SJ'M'> = (-1)      ( J  k  J' )   *   (LSJ||O ||L'SJ')
    //      j   q        j              (-Mj q  Mj')
    //
    //   where this is the 3-j symbol---^  and reduced matrix element--^
    //
    // The reduced matrix elements for the tensor operators are given by:
    //                   _______________
    //      k        -k | (2J + k + 1)!
    // <j||O ||j> = 2   | -------------
    //                 \|   (2J - k)!
    //
    // Reference: D.Smith and J.H.M. Thornley, Proc. Phys. Soc., 1966, vol 89, pp779.

    // Calculate the reduced matrix elements
    std::array<double, 7> rme{}; 
    for (int k=2; k<=6; k+=2) {
        rme[k] = (m_J2 > k) ? pow(2., -k) * sqrt( m_racah.f(m_J2 + k + 1) / m_racah.f(m_J2 - k) ) : 0.;
    }

    // The Hamiltonian must be hermitian, so we only need to compute values of the upper 
    // (or lower) triangle. In this case we choose to calculate the lower triangle as this is
    // what is needed by Eigen's SelfAdjointEigenSolver.

    // Each order q operator has only nonzero matrix elements along a diagonal of the matrix
    // Eg. Order 0 terms has it along the diagonal. Order 1 terms has it one element left,
    // and so on. We can use this to calculate only for the orders with nonzero parameters.

    // NB. the matrix elements calculated here using the threej symbols are actually matrix 
    // elements of the Wybourne normalised operator equivalent to spherical harmonics rather
    // than the Steven's operator equivalent to tesseral harmonics. 
    // Thus the internal parameters are in Wybourne normalisation but are divided by the
    // Stevens operator equivalent factors because they ignore the nature of ground multiplet.
    // See the Rudowicz paper for clarification.

    // First the diagonal elements (q=0 terms)
    for (auto iq: idq0) {
        int k = iq[0], m = iq[1];
        if (std::fabs(m_Bi[m]) > 1e-12) {
            for (int i=0; i<dimj; i++) {
                int mj = 2*i - m_J2;
                m_hamiltonian(i,i) += pow(-1., (m_J2-mj)/2.) * m_racah.threej(m_J2, 2*k, m_J2, -mj, 0, mj) * rme[k] * m_Bi[m] * m_econv;
            }
        }
    }

    // Now the off-diagonal terms - using a helper vector defined globally above to index
    for (auto iq: idq) {
        int k = iq[0], q = iq[1], m = iq[2], p = iq[3];
        if (std::fabs(m_Bi[m]) > 1e-12) {
            for (int i=0; i<(dimj-q); i++) {
                int mj = 2*i - m_J2, mjp = 2*(i+q) - m_J2;
                double tjp = m_racah.threej(m_J2, 2*k, m_J2, -mj, 2*q, mjp) - m_racah.threej(m_J2, 2*k, m_J2, -mj, -2*q, mjp);
                m_hamiltonian(i+q,i) += std::complex<double>(0., pow(-1., (m_J2-mj)/2.) * tjp * rme[k] * m_Bi[m] * m_econv);
            }
        }
        if (std::fabs(m_Bi[p]) > 1e-12) {
            for (int i=0; i<(dimj-q); i++) {
                int mj = 2*i - m_J2, mjp = 2*(i+q) - m_J2;
                double tjp = m_racah.threej(m_J2, 2*k, m_J2, -mj, 2*q, mjp) + m_racah.threej(m_J2, 2*k, m_J2, -mj, -2*q, mjp);
                m_hamiltonian(i+q,i) += pow(-1., (m_J2-mj)/2.) * tjp * rme[k] * m_Bi[p] * m_econv;
            }
        }
    }

    // Fill in upper triangle
    if (upper) {
        fill_upper();
    }

    m_ham_calc = true;
    m_upper = upper;
    return m_hamiltonian; 
}

RowMatrixXcd cf1ion::hamiltonian() {
    return _hamiltonian(true);
}

std::tuple<RowMatrixXcd, VectorXd> cf1ion::eigensystem() {
    if (!m_ev_calc) {
        SelfAdjointEigenSolver<RowMatrixXcd> es(_hamiltonian(false));
        m_eigenvectors = es.eigenvectors();
        m_eigenvalues = es.eigenvalues();
        m_ev_calc = true;
    }
    return std::tuple<RowMatrixXcd, VectorXd>(m_eigenvectors, m_eigenvalues);
}

RowMatrixXcd cf1ion::zeeman_hamiltonian(double H, std::vector<double> Hdir) {
    RowMatrixXcd zeeman = RowMatrixXcd::Zero(m_J2+1, m_J2+1);
    return zeeman;
}

std::vector<RowMatrixXcd> cf1ion::calculate_moments_matrix(RowMatrixXcd ev) {
    std::vector<RowMatrixXcd> moments;
    for (size_t ii=0; ii<3; ii++) {
        RowMatrixXcd Jmat = RowMatrixXcd::Zero(ev.rows(), ev.cols());
        moments.push_back((ev.adjoint()) * (Jmat * ev));
    }
    return moments;
}


} // namespace libMcPhase
