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
RowMatrixXcd cf1ion::hamiltonian(bool upper) {
    if (m_ham_calc) {
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
        for (int i=0; i<dimj; i++) {
            for (int j=i+1; j<dimj; j++) {
                m_hamiltonian(i,j) = std::conj(m_hamiltonian(j,i));
            }
        }
    }

    m_ham_calc = true;
    return m_hamiltonian; 
}

std::tuple<RowMatrixXcd, VectorXd> cf1ion::eigensystem() {
    if (!m_ev_calc) {
        SelfAdjointEigenSolver<RowMatrixXcd> es(hamiltonian(false));
        m_eigenvectors = es.eigenvectors();
        m_eigenvalues = es.eigenvalues();
        m_ev_calc = true;
    }
    return std::tuple<RowMatrixXcd, VectorXd>(m_eigenvectors, m_eigenvalues);
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates bulk properties (magnetisation, susceptibility)
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> cf1ion::calculate_boltzmann(VectorXd en, double T)
{
    std::vector<double> boltzmann, en_meV;
    // Need kBT in external energy units. K_B is in meV/K
    double beta = 1 / (K_B * T * m_econv);
    double Emin = std::numeric_limits<double>::max();
    for (size_t i=0; i < (size_t)en.size(); i++) {
        Emin = (en(i) < Emin) ? en(i) : Emin;
    }
    for (size_t i=0; i < (size_t)en.size(); i++) {
        const double expi = exp(-(en(i) - Emin) * beta);
        boltzmann.push_back((fabs(expi) > DELTA_EPS) ? expi : 0.);
    }
    return boltzmann;
}

std::vector<double> cf1ion::heatcapacity(std::vector<double> Tvec) {
    if(!m_ev_calc) {
        auto evsystem = eigensystem(); }
    std::vector<double> out;
    out.reserve(Tvec.size());
    std::vector<double> en;
    en.reserve(m_eigenvalues.size());
    for (size_t i=0; i < (size_t)m_eigenvalues.size(); i++) {
        en.push_back(m_eigenvalues(i) / m_econv);
    }
    for (auto T: Tvec) {
        double Z = 0., U = 0., U2 = 0.;
        std::vector<double> expfact = calculate_boltzmann(m_eigenvalues, T);
        for (size_t i=0; i < (size_t)m_eigenvalues.size(); i++) {
            Z += expfact[i];
            U += en[i] * expfact[i];
            U2 += en[i] * en[i] * expfact[i];
        }
        U /= Z;
        U2 /= Z;
        out.push_back( ((U2 - U * U) / (K_B * T * T)) * NAMEV );
    }
    return out;
}
/*
std::vector<double> cf1ion::magnetisation(std::vector<double> Hvec, std::vector<double> Hdir, double T, MagUnits unit_type)
{
    // Normalise the field direction vector
    double Hnorm = sqrt(Hdir[0] * Hdir[0] + Hdir[1] * Hdir[1] + Hdir[2] * Hdir[2]);
    if (fabs(Hnorm) < 1.e-6) {
        throw std::runtime_error("cf1ion::magnetisation(): Direction vector cannot be zero");
    }
    std::vector<double> nHdir;
    std::transform(Hdir.begin(), Hdir.end(), std::back_inserter(nHdir), [Hnorm](double Hd){ return Hd / Hnorm; });
    // Calculates Magnetisation M(H) at specified T
    if (!m_ham_calc)
        calculate_hamiltonian();
    std::vector<double> M;
    M.reserve(Hvec.size());
    // Loops through all the input field magnitudes and calculates the magnetisation
    for (auto H: Hvec) {
        if (unit_type == MagUnits::cgs) {
            H /= 1e4;   // For cgs, input field is in Gauss, need to convert to Tesla for Zeeman calculation
        }
        RowMatrixXcd ham = m_hamiltonian - zeeman_hamiltonian(H, Hdir);
        SelfAdjointEigenSolver<RowMatrixXcd> es(ham);
        // calculate_moments returns a vector of 3 moments *squared* vectors, in the x, y, z directions
        std::vector< std::vector<double> > moments_vec = calculate_moments(es.eigenvectors());
        std::vector<double> boltzmann = calculate_boltzmann(es.eigenvalues(), T);
        std::vector<double> Mdir;
        for (auto moments: moments_vec) {
            double Mexp = 0., Z = 0.;
            //std::inner_product(moments.begin(), moments.end(), boltzmann.begin(), Mexp);
            //std::accumulate(boltzmann.begin(), boltzmann.end(), Z);
            for (int ii=0; ii<ham.cols(); ii++) {
                Mexp += moments[ii] * boltzmann[ii];
                Z += boltzmann[ii];
            }
            Mdir.push_back(Mexp / Z);
        }
        M.push_back(sqrt(Mdir[0] * Mdir[0] + Mdir[1] * Mdir[1] + Mdir[2] * Mdir[2]) * MAGCONV[(int)unit_type]);
    }
    return M;
}

std::vector<double> cf1ion::susceptibility(std::vector<double> Tvec, std::vector<double> Hdir, MagUnits unit_type)
{
    // Normalise the field direction vector
    double Hnorm = sqrt(Hdir[0] * Hdir[0] + Hdir[1] * Hdir[1] + Hdir[2] * Hdir[2]);
    if (fabs(Hnorm) < 1.e-6) {
        throw std::runtime_error("cf1ion::magnetisation(): Direction vector cannot be zero");
    }
    std::vector<double> nHdir;
    std::transform(Hdir.begin(), Hdir.end(), std::back_inserter(nHdir), [Hnorm](double Hd){ return Hd / Hnorm; });
    // Calculates the susceptibility chi(T)
    if (!m_ev_calc)
        calculate_eigensystem();
    std::vector<double> chi;
    chi.reserve(Tvec.size());
    // Calculates the moments matrices in the x, y, z directions, and get the resultant
    std::vector<RowMatrixXcd> moments_mat_vec = calculate_moments_matrix(m_eigenvectors);
    RowMatrixXcd moments_mat = moments_mat_vec[0] * nHdir[0]
                               + moments_mat_vec[1] * nHdir[1]
                               + moments_mat_vec[2] * nHdir[2];
    // Now calculate the first and second order terms in the Van Vleck equation
    size_t nlev = m_eigenvectors.cols();
    std::vector<double> mu(nlev, 0.);
    std::vector<double> mu2(nlev, 0.);
    for (size_t ii=0; ii<nlev; ii++) {
        for (size_t jj=0; jj<nlev; jj++) {
            const double delta = m_eigenvalues[ii] - m_eigenvalues[jj];
            const double matel = (moments_mat(ii, jj) * std::conj(moments_mat(ii, jj))).real();
            if (fabs(delta) < DELTA_EPS) {
                mu[ii] += matel;           // First order term
            } else {
                mu2[ii] += matel / delta;  // Second order term
            }
        }
    }

    // Loops through all the input temperatures and calculates the susceptibility using:
    //                                 2                     2
    //           N_A --- [ <V_n|mu|V_n>      --- <V_n|mu|V_m>  ]
    // chi(T) =  --- >   [ ------------  - 2 >   ------------  ] exp(-E/k_BT)
    //            Z  --- [    k_B T          ---   En - Em     ]
    //                n                     m!=n

    for (auto T: Tvec) {
        std::vector<double> boltzmann = calculate_boltzmann(m_eigenvalues, T);
        const double beta = 1. / (K_B * T);
        double U = 0., Z = 0.;
        for (size_t ii=0; ii<nlev; ii++) {
            U += ((mu[ii] * beta) - (2 * mu2[ii])) * boltzmann[ii];
            Z += boltzmann[ii];
        }
        chi.push_back(SUSCCONV[(int)unit_type] * U / Z);
    }
    return chi;
}
*/

} // namespace libMcPhase
