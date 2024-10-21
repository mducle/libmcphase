/* pyic1ion.cpp
 *
 * Python bindings for the ic1ion class of libMcPhase
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "ic1ion.hpp"
#include <unordered_map>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "pycfpars.hpp"

namespace py = pybind11;
using namespace libMcPhase;

static const std::unordered_map<std::string, ic1ion::CoulombType> coulomb_names = {
    {"Slater", ic1ion::CoulombType::Slater}, {"CondonShortley", ic1ion::CoulombType::CondonShortley}, {"Racah", ic1ion::CoulombType::Racah} };

static const std::unordered_map<std::string, ic1ion::SpinOrbType> spinorb_names = {
    {"Zeta", ic1ion::SpinOrbType::Zeta}, {"Lambda", ic1ion::SpinOrbType::Lambda} };

ic1ion *ic1ion_init(py::args args, py::kwargs kwargs) {
    ic1ion *cls = new ic1ion;
    bool is_ic1ion = true;
    cf_parse(static_cast<cfpars*>(cls), args, kwargs, is_ic1ion);
    if (kwargs.contains("zeta")) {
        cls->set_spinorbit(kwargs["zeta"].cast<double>(), ic1ion::SpinOrbType::Zeta);
    }
    if (kwargs.contains("slater")) {
        cls->set_coulomb(kwargs["slater"].cast<std::vector<double>>(), ic1ion::CoulombType::Slater);
    }
    return cls;
}

void wrap_ic1ion(py::module &m) {

    py::class_<ic1ion, cfpars> pyic1ion(m, "ic1ion");

    pyic1ion.def(py::init<>())
        .def(py::init<const std::string &>(), py::arg("ionname"))
        .def(py::init(&ic1ion_init), cfpars_init_str)
        .def("set_coulomb", [](ic1ion &self, std::vector<double> val, std::string type) { self.set_coulomb(val,
             set_enum(type, coulomb_names, "Invalid normalisation, must be one of: Slater, CondonShortley, Racah")); })
        .def("set_spinorbit", [](ic1ion &self, double val, std::string type) { self.set_spinorbit(val,
             set_enum(type, spinorb_names, "Invalid normalisation, must be one of: Zeta, Lambda")); })
        .def("set_ci", &ic1ion::set_ci)
        .def("get_coulomb", &ic1ion::get_coulomb)
        .def("get_spinorbit", &ic1ion::get_spinorbit)
        .def("get_ci", &ic1ion::get_ci)
        .def("hamiltonian", &ic1ion::hamiltonian, "the crystal field Hamiltonian")
        .def("eigensystem", &ic1ion::eigensystem, "the eigenvectors and eigenvalues of the crystal field Hamiltonian")
        .def_property("GJ", [](ic1ion const &self) { return self.get_GJ(); }, [](ic1ion &self, double v) { self.set_GJ(v); })
        .def_property("zeta", [](ic1ion const &self) { return self.get_spinorbit(); }, [](ic1ion &self, double v) { self.set_spinorbit(v, ic1ion::SpinOrbType::Zeta); })
        .def_property("slater", [](ic1ion const &self) { return self.get_coulomb(); }, [](ic1ion &self, std::vector<double> v) { self.set_coulomb(v, ic1ion::CoulombType::Slater); })
        .def("zeeman_hamiltonian", &ic1ion::zeeman_hamiltonian, "the Zeeman Hamiltonian")
        .def("calculate_boltzmann", &ic1ion::calculate_boltzmann, "")
        .def("heatcapacity", &ic1ion::heatcapacity, "the heat capacity of the crystal field Hamiltonian in J/mol/K")
        .def("magnetisation", [](ic1ion &self, std::vector<double> H, std::vector<double> Hdir, std::vector<double> T, std::string unit) { return self.magnetisation(H, Hdir, T,
             set_enum(unit, mag_unit_names, "Invalid magnetic unit, must be one of: 'bohr', 'cgs', or 'SI'")); })
        .def("susceptibility", [](ic1ion &self, std::vector<double> T, std::vector<double> Hdir, std::string unit) { return self.susceptibility(T, Hdir,
             set_enum(unit, mag_unit_names, "Invalid magnetic unit, must be one of: 'bohr', 'cgs', or 'SI'")); })
        .def("get_states", &ic1ion::get_states, "Gets the list of states for this ion configuration");
}


