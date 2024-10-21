/* physprop.cpp
 * 
 * A class for calculating the physical properties (magnetisation, susceptibility, heat capacity) associated with
 * a crystal field type single-ion Hamiltonian.
 *
 * (C) 2024 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "physprop.hpp"
#include "eigen.hpp"

namespace libMcPhase {

// --------------------------------------------------------------------------------------------------------------- //
// Calculates bulk properties (heat capacity, magnetisation, susceptibility)
// --------------------------------------------------------------------------------------------------------------- //

std::vector<double> physprop::calculate_boltzmann(VectorXd en, double T) // Note that energy must be in meV in this routine
{
    std::vector<double> boltzmann;
    // Need kBT in external energy units. K_B is in cm-1/K
    double beta = 1. / (K_B * T * m_meVconv);
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

std::vector<double> physprop::heatcapacity(std::vector<double> Tvec) {
    auto es = eigensystem();
    std::vector<double> out;
    out.reserve(Tvec.size());
    std::vector<double> en;
    double Emin = std::numeric_limits<double>::max();
    size_t sz = std::get<1>(es).size();
    en.reserve(sz);
    for (size_t i=0; i < sz; i++) {
        en.push_back(std::get<1>(es)(i) / m_meVconv);
        Emin = (en[i] < Emin) ? en[i] : Emin;
    }
    for (size_t i=0; i < sz; i++) {
        en[i] -= Emin;
    }
    for (auto T: Tvec) {
        double Z = 0., U = 0., U2 = 0.;
        std::vector<double> expfact = calculate_boltzmann(std::get<1>(es), T);
        for (size_t i=0; i < sz; i++) {
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

std::vector< std::vector<double> > physprop::magnetisation(std::vector<double> Hvec, std::vector<double> Hdir, std::vector<double> Tvec, MagUnits unit_type)
{
    // Normalise the field direction vector
    double Hnorm = sqrt(Hdir[0] * Hdir[0] + Hdir[1] * Hdir[1] + Hdir[2] * Hdir[2]);
    if (fabs(Hnorm) < 1.e-6) {
        throw std::runtime_error("physprop::magnetisation(): Direction vector cannot be zero");
    }
    std::vector<double> nHdir;
    std::transform(Hdir.begin(), Hdir.end(), std::back_inserter(nHdir), [Hnorm](double Hd){ return Hd / Hnorm; });
    // Calculates Magnetisation M(H) at specified T
    RowMatrixXcd ham0 = hamiltonian();
    std::vector< std::vector<double> > M;
    M.reserve(Hvec.size());
    // Loops through all the input field magnitudes and calculates the magnetisation
    std::vector<RowMatrixXcd> mag_ops = calculate_moments_matrix(RowMatrixXcd::Identity(ham0.rows(), ham0.cols()));
    RowMatrixXcd Jmat = nHdir[0] * mag_ops[0] + nHdir[1] * mag_ops[1] + nHdir[2] * mag_ops[2];
    for (auto H: Hvec) {
        if (unit_type == MagUnits::cgs) {
            H /= 1e4;   // For cgs, input field is in Gauss, need to convert to Tesla for Zeeman calculation
        }
        RowMatrixXcd ham = ham0 - zeeman_hamiltonian(H, Hdir);
        SelfAdjointEigenSolver<RowMatrixXcd> es(ham);
        std::vector<double> Mt;
        Mt.reserve(Tvec.size());
        for (auto T: Tvec) {
            std::vector<double> boltzmann = calculate_boltzmann(es.eigenvalues(), T);
            RowMatrixXcd me = (es.eigenvectors().adjoint()) * (Jmat * es.eigenvectors());
            double Mexp = 0., Z = 0.;
            for (int ii=0; ii<ham.cols(); ii++) {
                Mexp += me(ii,ii).real() * boltzmann[ii];
                Z += boltzmann[ii];
            }
            Mt.push_back((Mexp / Z) * MAGCONV[(int)unit_type]);
        }
        M.push_back(Mt);
    }
    return M;
}

std::vector<double> physprop::susceptibility(std::vector<double> Tvec, std::vector<double> Hdir, MagUnits unit_type)
{
    // Normalise the field direction vector
    double Hnorm = sqrt(Hdir[0] * Hdir[0] + Hdir[1] * Hdir[1] + Hdir[2] * Hdir[2]);
    if (fabs(Hnorm) < 1.e-6) {
        throw std::runtime_error("physprop::magnetisation(): Direction vector cannot be zero");
    }
    std::vector<double> nHdir;
    std::transform(Hdir.begin(), Hdir.end(), std::back_inserter(nHdir), [Hnorm](double Hd){ return Hd / Hnorm; });
    // Calculates the susceptibility chi(T)
    std::tuple<RowMatrixXcd, VectorXd> es = eigensystem();
    std::vector<double> chi;
    chi.reserve(Tvec.size());
    // Calculates the moments matrices in the x, y, z directions, and get the resultant
    std::vector<RowMatrixXcd> moments_mat_vec = calculate_moments_matrix(std::get<0>(es));
    RowMatrixXcd moments_mat = moments_mat_vec[0] * nHdir[0]
                             + moments_mat_vec[1] * nHdir[1]
                             + moments_mat_vec[2] * nHdir[2];
    // Now calculate the first and second order terms in the Van Vleck equation
    size_t nlev = std::get<0>(es).cols();
    std::vector<double> mu(nlev, 0.);
    std::vector<double> mu2(nlev, 0.);
    std::vector<double> eigenvalues(nlev, 0.);
    for (size_t ii=0; ii<nlev; ii++) {
        eigenvalues[ii] = std::get<1>(es)[ii] / m_meVconv; }  // Convert energies to meV
    for (size_t ii=0; ii<nlev; ii++) {
        for (size_t jj=0; jj<nlev; jj++) {
            const double delta = eigenvalues[ii] - eigenvalues[jj];
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
        std::vector<double> boltzmann = calculate_boltzmann(std::get<1>(es), T);
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

} // namespace libMcPhase
