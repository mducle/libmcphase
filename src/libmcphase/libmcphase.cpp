/* libmcphase.cpp
 *
 * Python bindings for libMcPhase
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include <pybind11/pybind11.h>

namespace py = pybind11;

void wrap_cfpars(py::module &);
void wrap_cf1ion(py::module &);
void wrap_ic1ion(py::module &);
void wrap_icstates(py::module &);

PYBIND11_MODULE(libmcphase, m) {
    m.doc() = "Python bindings for libMcPhase";

    wrap_cfpars(m);
    wrap_cf1ion(m);
    wrap_ic1ion(m);
    wrap_icstates(m);
}


