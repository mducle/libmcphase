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
#include "pycfpars.hpp"

namespace py = pybind11;
using namespace libMcPhase;
using namespace pybind11::literals;

void cf_parse(cfpars *cls, py::args args, py::kwargs kwargs);

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
        .def("hamiltonian", &cf1ion::hamiltonian, "the crystal field Hamiltonian", "upper"_a=true)
        .def("eigensystem", &cf1ion::eigensystem, "the eigenvectors and eigenvalues of the crystal field Hamiltonian");
}


