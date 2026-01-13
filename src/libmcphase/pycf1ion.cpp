/* pycf1ion.cpp
 *
 * Python bindings for the cf1ion class of libMcPhase
 *
 * (C) 2020 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "cf1ion.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "pycfpars.hpp"

namespace py = pybind11;
using namespace libMcPhase;
using namespace pybind11::literals;

cf1ion *cf1ion_init(py::args args, py::kwargs kwargs) {
    cf1ion *cls = new cf1ion;
    cf_parse(static_cast<cfpars*>(cls), args, kwargs);
    return cls;
}

void wrap_cf1ion(py::module &m) {

    py::class_<cf1ion, cfpars> pycf1ion(m, "cf1ion");

    pycf1ion.def(py::init<>())
        .def(py::init<const std::string &>(), py::arg("ionname"))
        .def(py::init<const double &>(), py::arg("J"))
        .def(py::init(&cf1ion_init), cfpars_init_str)
        .def_property("GJ", [](cf1ion const &self) { return self.get_GJ(); }, [](cf1ion &self, double v) { self.set_GJ(v); })
        .def("hamiltonian", &cf1ion::hamiltonian, "the crystal field Hamiltonian")
        .def("eigensystem", &cf1ion::eigensystem, "the eigenvectors and eigenvalues of the crystal field Hamiltonian")
        .def("zeeman_hamiltonian", &cf1ion::zeeman_hamiltonian, "the Zeeman Hamiltonian")
        .def("calculate_boltzmann", &cf1ion::calculate_boltzmann, "")
        .def("heatcapacity", &cf1ion::heatcapacity, "the heat capacity of the crystal field Hamiltonian in J/mol/K")
        .def("magnetisation", [](cf1ion &self, std::vector<double> T, std::vector<double> H, std::vector<double> Hdir, std::string unit) { return self.magnetisation(T, H, Hdir,
             set_enum(unit, mag_unit_names, "Invalid magnetic unit, must be one of: 'bohr', 'cgs', or 'SI'")); })
        .def("susceptibility", [](cf1ion &self, std::vector<double> T, std::vector<double> Hdir, std::string unit) { return self.susceptibility(T, Hdir,
             set_enum(unit, mag_unit_names, "Invalid magnetic unit, must be one of: 'bohr', 'cgs', or 'SI'")); })
        .def("peaks", &cf1ion::peaks, "list of peaks with intensity in mb/sr")
        .def("split2range", [](cf1ion &self, double E, bool s) { return self.split2range(E, s); }, py::arg("Energy_splitting"), py::arg("use_sym")=false,
            "returns maximum absolute values of Blm which when sampled will give on average E0 splitting")
        .def("fitengy", &cf1ion::fitengy, py::arg("E"), py::arg("use_sym"), "update parameters to fit input energies using Newman-Ng algorithm");
}


