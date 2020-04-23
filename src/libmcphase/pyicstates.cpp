/* pyicstates.cpp
 *
 * Python bindings for the fstates class of libMcPhase
 *
 * (C) 2020 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "ic_states.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace libMcPhase;

void wrap_icstates(py::module &m) {

    py::class_<fstates_t> pyicstate(m, "ic_state");

    py::enum_<orbital>(pyicstate, "orbital")
        .value("S", orbital::S)
        .value("P", orbital::P)
        .value("D", orbital::D)
        .value("F", orbital::F)
        .value("Fp", orbital::Fp)
        .value("G", orbital::G)
        .value("Gp", orbital::Gp)
        .value("H", orbital::H)
        .value("Hp", orbital::Hp)
        .value("I", orbital::I)
        .value("Ip", orbital::Ip)
        .value("K", orbital::K)
        .value("Kp", orbital::Kp)
        .value("L", orbital::L)
        .value("Lp", orbital::Lp)
        .value("M", orbital::M)
        .value("N", orbital::N)
        .value("O", orbital::O)
        .value("Q", orbital::Q)
        .export_values();

    pyicstate
        .def_property_readonly("S2", [](const fstates_t &self){ return self.S2; }, "Twice the total spin momentum")
        .def_property_readonly("L", [](const fstates_t &self){ return (int)abs(self.L); }, "The orbital momentum")
        .def_property_readonly("id", [](const fstates_t &self){ return self.id; }, "id string")
        .def_property_readonly("J2", [](const fstates_t &self){ return self.J2; }, "Twice the total momentum")
        .def_property_readonly("mJ2", [](const fstates_t &self){ return self.mJ2; }, "Twice the total azimuth momentum");
}


