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

VectorXd physprop::calculate_boltzmann(VectorXd en, double T) // Note that energy must be in meV in this routine
{
    VectorXd boltzmann = VectorXd::Zero(en.size());
    // Need kBT in external energy units. K_B is in cm-1/K
    double beta = 1. / (K_B * T * m_meVconv);
    double Emin = std::numeric_limits<double>::max();
    for (size_t i=0; i < (size_t)en.size(); i++) {
        Emin = (en(i) < Emin) ? en(i) : Emin;
    }
    for (size_t i=0; i < (size_t)en.size(); i++) {
        const double expi = exp(-(en(i) - Emin) * beta);
        boltzmann(i) = (fabs(expi) > DELTA_EPS) ? expi : 0.;
    }
    return boltzmann;
}

VectorXd physprop::heatcapacity(std::vector<double> Tvec) {
    auto es = eigensystem();
    VectorXd out = VectorXd::Zero(Tvec.size());
    size_t sz = std::get<1>(es).size();
    VectorXd en = VectorXd::Zero(sz);
    double Emin = std::numeric_limits<double>::max();
    for (size_t i=0; i < sz; i++) {
        en(i) = std::get<1>(es)(i) / m_meVconv;
        Emin = (en[i] < Emin) ? en[i] : Emin;
    }
    for (size_t i=0; i < sz; i++) {
        en[i] -= Emin;
    }
    for (size_t tt=0; tt<Tvec.size(); tt++) {
        double Z = 0., U = 0., U2 = 0., T = Tvec[tt];
        VectorXd expfact = calculate_boltzmann(std::get<1>(es), T);
        for (size_t i=0; i < sz; i++) {
            Z += expfact[i];
            U += en[i] * expfact[i];
            U2 += en[i] * en[i] * expfact[i];
        }
        U /= Z;
        U2 /= Z;
        out(tt) = ((U2 - U * U) / (K_B * T * T)) * NAMEV;
    }
    return out;
}

RowMatrixXd physprop::magnetisation(std::vector<double> Tvec, std::vector<double> Hvec, std::vector<double> Hdir, MagUnits unit_type)
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
    RowMatrixXd M = RowMatrixXd::Zero(Hvec.size(), Tvec.size());
    // Loops through all the input field magnitudes and calculates the magnetisation
    std::vector<RowMatrixXcd> mag_ops = calculate_moments_matrix(RowMatrixXcd::Identity(ham0.rows(), ham0.cols()));
    RowMatrixXcd Jmat = nHdir[0] * mag_ops[0] + nHdir[1] * mag_ops[1] + nHdir[2] * mag_ops[2];
    for (size_t hh=0; hh<Hvec.size(); hh++) {
        double H = Hvec[hh];
        if (unit_type == MagUnits::cgs) {
            H /= 1e4;   // For cgs, input field is in Gauss, need to convert to Tesla for Zeeman calculation
        }
        RowMatrixXcd ham = ham0 - zeeman_hamiltonian(H, Hdir);
        SelfAdjointEigenSolver<RowMatrixXcd> es(ham);
        for (size_t tt=0; tt<Tvec.size(); tt++) {
            VectorXd boltzmann = calculate_boltzmann(es.eigenvalues(), Tvec[tt]);
            RowMatrixXcd me = (es.eigenvectors().adjoint()) * (Jmat * es.eigenvectors());
            double Mexp = 0., Z = 0.;
            for (int ii=0; ii<ham.cols(); ii++) {
                Mexp += me(ii,ii).real() * boltzmann[ii];
                Z += boltzmann[ii];
            }
            M(hh, tt) = (Mexp / Z) * MAGCONV[(int)unit_type];
        }
    }
    return M;
}

VectorXd physprop::susceptibility(std::vector<double> Tvec, std::vector<double> Hdir, MagUnits unit_type)
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
    VectorXd chi = VectorXd::Zero(Tvec.size());
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

    for (size_t tt=0; tt<Tvec.size(); tt++) {
        VectorXd boltzmann = calculate_boltzmann(std::get<1>(es), Tvec[tt]);
        const double beta = 1. / (K_B * Tvec[tt]);
        double U = 0., Z = 0.;
        for (size_t ii=0; ii<nlev; ii++) {
            U += ((mu[ii] * beta) - (2 * mu2[ii])) * boltzmann[ii];
            Z += boltzmann[ii];
        }
        chi(tt) = SUSCCONV[(int)unit_type] * U / Z;
    }
    return chi;
}

RowMatrixXd physprop::peaks(double T)
{
    auto es = eigensystem();
    std::vector<RowMatrixXcd> moments_mat_vec = calculate_moments_matrix(std::get<0>(es));
    RowMatrixXcd trans = moments_mat_vec[0].cwiseProduct(moments_mat_vec[0].conjugate()) +
                         moments_mat_vec[1].cwiseProduct(moments_mat_vec[1].conjugate()) +
                         moments_mat_vec[2].cwiseProduct(moments_mat_vec[2].conjugate());
    size_t sz = std::get<1>(es).size();
    double Z = 0.;
    VectorXd expfact = calculate_boltzmann(std::get<1>(es), T);
    for (size_t i=0; i < sz; i++) {
        Z += expfact[i];
        trans.col(i) *= expfact[i];
    }
    trans *= MAGXSEC_MBSR / Z;  // The magnetic cross-section in milibarn/sr
    // Match transition matrix elements to energies
    std::vector< std::pair<double, double> > pkl;
    pkl.reserve(sz * sz);
    for (int i=0; i<sz; i++) {
        for (int j=0; j<sz; j++) {
            if (trans(i, j).real() > 1e-6) {
                pkl.push_back(std::make_pair(std::get<1>(es)[i] - std::get<1>(es)[j], trans(i,j).real()));
            }
        }
    }
    // Sums degenerate transitions
    std::sort(pkl.begin(), pkl.end(), 
        [](std::pair<double, double> a, std::pair<double, double> b) { return a.first > b.first; });
    std::vector<int> ndegen;
    int n_ex = 0, j = 0;
    ndegen.reserve(pkl.size());
    ndegen.push_back(j);
    for (int i=1; i<pkl.size(); i++) {
        double dE = pkl[j].first - pkl[i].first;
        if (dE < 1e-3 || (dE/pkl[i].first) < 1e-2) {
            pkl[j].second += pkl[i].second;
        } else {
            j = i;
            ndegen.push_back(j);
        }
    }
    RowMatrixXd rv = RowMatrixXd::Zero(ndegen.size(), 2);
    // Sort by intensity
    std::sort(ndegen.begin(), ndegen.end(), [&pkl](int a, int b) { return pkl[a].second > pkl[b].second; });
    for (int i=0; i<ndegen.size(); i++) {
        rv.row(i) << pkl[ndegen[i]].first, pkl[ndegen[i]].second;
    }
    return rv;
}

} // namespace libMcPhase
