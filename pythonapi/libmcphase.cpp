/* libmcphase.cpp
 * 
 * Python bindings for libMcPhase
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "cfpars.hpp"
#include "cf1ion.hpp"
#include "ic1ion.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace libMcPhase;
using namespace pybind11::literals;

namespace { // Anonymous namespace for local functions

static const std::unordered_map<std::string, cfpars::Blm> Blm_names = {
    {"B22S", cfpars::Blm::B22S}, {"B21S", cfpars::Blm::B21S}, {"B20", cfpars::Blm::B20}, {"B21", cfpars::Blm::B21}, {"B22", cfpars::Blm::B22},
    {"B44S", cfpars::Blm::B44S}, {"B43S", cfpars::Blm::B43S}, {"B42S", cfpars::Blm::B42S}, {"B41S", cfpars::Blm::B41S}, 
    {"B40", cfpars::Blm::B40}, {"B41", cfpars::Blm::B41}, {"B42", cfpars::Blm::B42}, {"B43", cfpars::Blm::B43}, {"B44", cfpars::Blm::B44},
    {"B66S", cfpars::Blm::B66S}, {"B65S", cfpars::Blm::B65S}, {"B64S", cfpars::Blm::B64S}, {"B63S", cfpars::Blm::B63S}, 
    {"B62S", cfpars::Blm::B62S}, {"B61S", cfpars::Blm::B61S}, {"B60", cfpars::Blm::B60}, {"B61", cfpars::Blm::B61}, {"B62", cfpars::Blm::B62}, 
    {"B63", cfpars::Blm::B63}, {"B64", cfpars::Blm::B64}, {"B65", cfpars::Blm::B65}, {"B66", cfpars::Blm::B66} };

static const std::unordered_map<std::string, cfpars::Type> type_names = {
    {"Alm", cfpars::Type::Alm}, {"Blm", cfpars::Type::Blm}, {"Llm", cfpars::Type::Llm}, {"ARlm", cfpars::Type::ARlm} };

static const std::unordered_map<std::string, cfpars::Units> unit_names = { 
    {"meV", cfpars::Units::meV}, {"cm", cfpars::Units::cm}, {"K", cfpars::Units::K} };

static const char* cfpars_init_str = "Construct a cfpars object\n"
                                     "    args (one of):\n"
                                     "        J - total angular momentum of multiplet\n"
                                     "        ionname - name of ion, e.g. 'Pr3+'\n"
                                     "    kwargs:\n"
                                     "        unit - energy units of parameters. one of 'meV' (default), 'cm' or 'K'\n"
                                     "        type - type of parameters one of 'Alm', 'ARlm', 'Blm' (default), 'Llm' 'Vlm', or 'Wlm'\n"
                                     "               see online documentation or McPhase webpage for definitions\n"
                                     "        Blm - value of parameters, e.g. cfp = cfpars('Pr3+', B20=0.1, B22=-0.01, B40=0.001)\n";

void cfpars_set_type(cfpars *cls, std::string typestr) {
    std::unordered_map<std::string, cfpars::Type>::const_iterator type_it = type_names.find(typestr);
    if (type_it == type_names.end()) {
        throw std::runtime_error("Invalid type name, must be one of 'Alm', 'ARlm', 'Blm', 'Llm', 'Vlm' or 'Wlm'");
    } else {
        cls->set_type(type_it->second);
    }
}

void cfpars_set_unit(cfpars *cls, std::string unitstr) {
    std::unordered_map<std::string, cfpars::Units>::const_iterator unit_it = unit_names.find(unitstr);
    if (unit_it == unit_names.end()) {
        throw std::runtime_error("Invalid unit, must be one of 'meV', 'cm', or 'K'");
    } else {
        cls->set_unit(unit_it->second);
    }
}

void cf_parse(cfpars *cls, py::args args, py::kwargs kwargs) {
    if (!args && !kwargs) {
        return;
    }
    if (args && args.size() > 0) {
        try {
            double J = args[0].cast<double>();
            cls->set_J(J);
        } catch (py::cast_error) {
            try {
                std::string ionname = args[0].cast<std::string>();
                cls->set_name(ionname);
            } catch (py::cast_error) {
                throw std::runtime_error("Invalid first argument: must be the ion name as a string "
                                         "or the total angular momentum quantum number J");
            }
        }
    }
    else if (kwargs.contains("J")) {
        cls->set_J(kwargs["J"].cast<double>());
    } 
    else if (kwargs.contains("ionname")) {
        cls->set_name(kwargs["ionname"].cast<std::string>());
    } 
    if (kwargs.contains("type")) {
        cfpars_set_type(cls, kwargs["type"].cast<std::string>());
    }
    if (kwargs.contains("unit")) {
        cfpars_set_unit(cls, kwargs["unit"].cast<std::string>());
    }
    if (kwargs) {
        for (auto const &bname : Blm_names) {
            if (kwargs.contains(bname.first.c_str())) {
                cls->set(bname.second, kwargs[bname.first.c_str()].cast<double>());
            }
        }
    }
}

cfpars *cfpars_init(py::args args, py::kwargs kwargs) {
    cfpars *cls = new cfpars;
    cf_parse(cls, args, kwargs);
    return cls;
}

cf1ion *cf1ion_init(py::args args, py::kwargs kwargs) {
    cf1ion *cls = new cf1ion;
    cf_parse(static_cast<cfpars*>(cls), args, kwargs);
    return cls;
}

ic1ion *ic1ion_init(py::args args, py::kwargs kwargs) {
    ic1ion *cls = new ic1ion;
    cf_parse(static_cast<cfpars*>(cls), args, kwargs);
    return cls;
}

} // Anonymous namespace

PYBIND11_MODULE(libMcPhase, m) {
    m.doc() = "Python bindings for libMcPhase";

    py::class_<cfpars> pycfpars(m, "cfpars");

    py::enum_<cfpars::Units>(pycfpars, "Units")
        .value("meV", cfpars::Units::meV)
        .value("cm", cfpars::Units::cm)
        .value("K", cfpars::Units::K);

    py::enum_<cfpars::Normalisation>(pycfpars, "Normalisation")
        .value("Stevens", cfpars::Normalisation::Stevens)
        .value("Wybourne", cfpars::Normalisation::Wybourne);

    py::enum_<cfpars::Type>(pycfpars, "Types")
        .value("Alm", cfpars::Type::Alm)
        .value("Blm", cfpars::Type::Blm)
        .value("Llm", cfpars::Type::Llm)
        .value("ARlm", cfpars::Type::ARlm);

    pycfpars.def(py::init<>())
        .def(py::init<const std::string &>(), py::arg("ionname"))
        .def(py::init<const double &>(), py::arg("J"))
        .def(py::init(&cfpars_init), cfpars_init_str)
        .def_property("unit", &cfpars::get_unit, &cfpars_set_unit, "energy unit of the parameters")
        .def_property("type", &cfpars::get_type, &cfpars_set_type, "type of CF parameters. When changed, parameters values will be automatically converted.")
        .def_property("ion", &cfpars::get_name, &cfpars::set_name, "ion type. If reset, parameters will be updated by scaling by Stevens factors.")
        .def_property("J", &cfpars::get_J, &cfpars::set_J, "the total angular momentum quantum number of this multiplet")
        .def_property_readonly("normalisation", &cfpars::get_normalisation, "normalisation of CF parameters")
        .def_property_readonly("alpha", &cfpars::alpha, "the 2nd order Stevens operator equivalent factor")
        .def_property_readonly("beta", &cfpars::beta, "the 4th order Stevens operator equivalent factor")
        .def_property_readonly("gamma", &cfpars::gamma, "the 6th order Stevens operator equivalent factor")
        .def_property("B22S", [](cfpars const &self) { return self.get(cfpars::Blm::B22S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B22S, v); })
        .def_property("B21S", [](cfpars const &self) { return self.get(cfpars::Blm::B21S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B21S, v); })
        .def_property("B20", [](cfpars const &self) { return self.get(cfpars::Blm::B20); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B20, v); })
        .def_property("B21", [](cfpars const &self) { return self.get(cfpars::Blm::B21); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B21, v); })
        .def_property("B22", [](cfpars const &self) { return self.get(cfpars::Blm::B22); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B22, v); })
        .def_property("B44S", [](cfpars const &self) { return self.get(cfpars::Blm::B44S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B44S, v); })
        .def_property("B43S", [](cfpars const &self) { return self.get(cfpars::Blm::B43S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B43S, v); })
        .def_property("B42S", [](cfpars const &self) { return self.get(cfpars::Blm::B42S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B42S, v); })
        .def_property("B41S", [](cfpars const &self) { return self.get(cfpars::Blm::B41S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B41S, v); })
        .def_property("B40", [](cfpars const &self) { return self.get(cfpars::Blm::B40); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B40, v); })
        .def_property("B41", [](cfpars const &self) { return self.get(cfpars::Blm::B41); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B41, v); })
        .def_property("B42", [](cfpars const &self) { return self.get(cfpars::Blm::B42); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B42, v); })
        .def_property("B43", [](cfpars const &self) { return self.get(cfpars::Blm::B43); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B43, v); })
        .def_property("B44", [](cfpars const &self) { return self.get(cfpars::Blm::B44); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B44, v); })
        .def_property("B66S", [](cfpars const &self) { return self.get(cfpars::Blm::B66S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B66S, v); })
        .def_property("B65S", [](cfpars const &self) { return self.get(cfpars::Blm::B65S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B65S, v); })
        .def_property("B64S", [](cfpars const &self) { return self.get(cfpars::Blm::B64S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B64S, v); })
        .def_property("B63S", [](cfpars const &self) { return self.get(cfpars::Blm::B63S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B63S, v); })
        .def_property("B62S", [](cfpars const &self) { return self.get(cfpars::Blm::B62S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B62S, v); })
        .def_property("B61S", [](cfpars const &self) { return self.get(cfpars::Blm::B61S); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B61S, v); })
        .def_property("B60", [](cfpars const &self) { return self.get(cfpars::Blm::B60); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B60, v); })
        .def_property("B61", [](cfpars const &self) { return self.get(cfpars::Blm::B61); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B61, v); })
        .def_property("B62", [](cfpars const &self) { return self.get(cfpars::Blm::B62); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B62, v); })
        .def_property("B63", [](cfpars const &self) { return self.get(cfpars::Blm::B63); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B63, v); })
        .def_property("B64", [](cfpars const &self) { return self.get(cfpars::Blm::B64); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B64, v); })
        .def_property("B65", [](cfpars const &self) { return self.get(cfpars::Blm::B65); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B65, v); })
        .def_property("B66", [](cfpars const &self) { return self.get(cfpars::Blm::B66); }, [](cfpars &self, double v) { self.set(cfpars::Blm::B66, v); });

    py::class_<cf1ion, cfpars> pycf1ion(m, "cf1ion");
    
    pycf1ion.def(py::init<>())
        .def(py::init<const std::string &>(), py::arg("ionname"))
        .def(py::init<const double &>(), py::arg("J"))
        .def(py::init(&cf1ion_init), cfpars_init_str)
        .def("hamiltonian", &cf1ion::hamiltonian, "the crystal field Hamiltonian", "upper"_a=true)
        .def("eigensystem", &cf1ion::eigensystem, "the eigenvectors and eigenvalues of the crystal field Hamiltonian");

    py::class_<ic1ion, cfpars> pyic1ion(m, "ic1ion");
    
    pyic1ion.def(py::init<>())
        .def(py::init<const std::string &>(), py::arg("ionname"))
        .def(py::init(&ic1ion_init), cfpars_init_str)
        .def("racah_so", &ic1ion::racah_so)
        .def("racah_Umat", &ic1ion::racah_Umat, py::arg("k"))
        .def("racah_ukq", &ic1ion::racah_ukq, py::arg("k"), py::arg("q"))
        .def("racah_emat", &ic1ion::emat)
        .def("racah_ci", &ic1ion::ci)
        .def("racahW", &ic1ion::racahW)
        .def("hamiltonian", &ic1ion::hamiltonian, "the crystal field Hamiltonian")
        .def("eigensystem", &ic1ion::eigensystem, "the eigenvectors and eigenvalues of the crystal field Hamiltonian");

}


