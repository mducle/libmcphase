/* pycfpars.cpp
 *
 * Python bindings for the cfpars class of libMcPhase
 *
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#include "pycfpars.hpp"

namespace py = pybind11;
using namespace libMcPhase;

static const std::unordered_map<std::string, cfpars::Blm> Blm_names = {
    {"B22S", cfpars::Blm::B22S}, {"B21S", cfpars::Blm::B21S}, {"B20", cfpars::Blm::B20}, {"B21", cfpars::Blm::B21}, {"B22", cfpars::Blm::B22},
    {"B44S", cfpars::Blm::B44S}, {"B43S", cfpars::Blm::B43S}, {"B42S", cfpars::Blm::B42S}, {"B41S", cfpars::Blm::B41S},
    {"B40", cfpars::Blm::B40}, {"B41", cfpars::Blm::B41}, {"B42", cfpars::Blm::B42}, {"B43", cfpars::Blm::B43}, {"B44", cfpars::Blm::B44},
    {"B66S", cfpars::Blm::B66S}, {"B65S", cfpars::Blm::B65S}, {"B64S", cfpars::Blm::B64S}, {"B63S", cfpars::Blm::B63S},
    {"B62S", cfpars::Blm::B62S}, {"B61S", cfpars::Blm::B61S}, {"B60", cfpars::Blm::B60}, {"B61", cfpars::Blm::B61}, {"B62", cfpars::Blm::B62},
    {"B63", cfpars::Blm::B63}, {"B64", cfpars::Blm::B64}, {"B65", cfpars::Blm::B65}, {"B66", cfpars::Blm::B66} };

static const std::array<std::string, 27> Blmvec = {{"B22S", "B21S", "B20", "B21", "B22", "B44S", "B43S", "B42S", "B41S", 
    "B40", "B41", "B42", "B43", "B44", "B66S", "B65S", "B64S", "B63S", "B62S", "B61S", "B60", "B61", "B62", "B63", "B64", "B65", "B66" }};

static const std::unordered_map<std::string, cfpars::Type> type_names = {
    {"Alm", cfpars::Type::Alm}, {"Blm", cfpars::Type::Blm}, {"Llm", cfpars::Type::Llm}, {"ARlm", cfpars::Type::ARlm}, {"Nlm", cfpars::Type::Nlm} };
static const std::string type_err = "Invalid type name, must be one of 'Alm', 'ARlm', 'Blm', 'Llm', 'Nlm', 'Vlm' or 'Wlm'";

static const std::unordered_map<std::string, cfpars::Units> unit_names = {
    {"meV", cfpars::Units::meV}, {"cm", cfpars::Units::cm}, {"K", cfpars::Units::K} };
static const std::string unit_err = "Invalid unit, must be one of 'meV', 'cm', or 'K'";

static const std::unordered_map<std::string, cfpars::Sym> sym_names = {
    {"Ci", cfpars::Sym::Ci}, {"C1", cfpars::Sym::C1}, {"C2", cfpars::Sym::C2}, {"Cs", cfpars::Sym::Cs}, {"C2h", cfpars::Sym::C2h},
    {"C2v", cfpars::Sym::C2v}, {"D2", cfpars::Sym::D2}, {"D2h", cfpars::Sym::D2h}, {"C4", cfpars::Sym::C4}, {"S4", cfpars::Sym::S4},
    {"C4h", cfpars::Sym::C4h}, {"D4", cfpars::Sym::D4}, {"C4v", cfpars::Sym::C4v}, {"D2d", cfpars::Sym::D2d}, {"D4h", cfpars::Sym::D4h},
    {"C3", cfpars::Sym::C3}, {"S6", cfpars::Sym::S6}, {"D3", cfpars::Sym::D3}, {"C3v", cfpars::Sym::C3v}, {"D3d", cfpars::Sym::D3d},
    {"C6", cfpars::Sym::C6}, {"C3h", cfpars::Sym::C3h}, {"C6h", cfpars::Sym::C6h}, {"D6", cfpars::Sym::D6}, {"C6v", cfpars::Sym::C6v},
    {"D3h", cfpars::Sym::D3h}, {"D6h", cfpars::Sym::D6h}, {"T", cfpars::Sym::T}, {"Th", cfpars::Sym::Th}, {"Td", cfpars::Sym::Td},
    {"O", cfpars::Sym::O}, {"Oh", cfpars::Sym::Oh} };
static const std::string sym_err = "Invalid point symmetry, must be a Schoenflies symbol";

void cf_parse(cfpars *cls, py::args args, py::kwargs kwargs, bool is_ic1ion) {
    if (!args && !kwargs) {
        return;
    }
    if (args && args.size() > 0) {
        try {
            std::string ionname = args[0].cast<std::string>();
            cls->set_name(ionname);
        } catch (py::cast_error) {
            if (is_ic1ion) {
                throw std::runtime_error("ic1ion module can only be initialised by an ion name");
            }
            try {
                double J = args[0].cast<double>();
                cls->set_J(J);
            } catch (py::cast_error) {
                throw std::runtime_error("Invalid first argument: must be the ion name as a string "
                                         "or the total angular momentum quantum number J");
            }
        }
    }
    else if (!is_ic1ion && kwargs.contains("J")) {
        cls->set_J(kwargs["J"].cast<double>());
    }
    else if (kwargs.contains("ionname")) {
        cls->set_name(kwargs["ionname"].cast<std::string>());
    }
    if (kwargs.contains("type")) {
        cls->set_type(set_enum(kwargs["type"].cast<std::string>(), type_names, type_err));
    }
    if (kwargs.contains("unit")) {
        cls->set_unit(set_enum(kwargs["unit"].cast<std::string>(), unit_names, unit_err));
    }
    if (kwargs.contains("sym")) {
        cls->set_sym(set_enum(kwargs["sym"].cast<std::string>(), sym_names, sym_err));
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

void wrap_cfpars(py::module &m) {

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
        .value("ARlm", cfpars::Type::ARlm)
        .value("Nlm", cfpars::Type::Nlm);

    py::enum_<cfpars::Sym>(pycfpars, "Sym")
        .value("Ci", cfpars::Sym::Ci).value("C1", cfpars::Sym::C1)
        .value("C2", cfpars::Sym::C2).value("Cs", cfpars::Sym::Cs).value("C2h", cfpars::Sym::C2h)
        .value("C2v", cfpars::Sym::C2v).value("D2", cfpars::Sym::D2).value("D2h", cfpars::Sym::D2h)
        .value("C4", cfpars::Sym::C4).value("S4", cfpars::Sym::S4).value("C4h", cfpars::Sym::C4h)
        .value("D4", cfpars::Sym::D4).value("C4v", cfpars::Sym::C4v).value("D2d", cfpars::Sym::D2d) .value("D4h", cfpars::Sym::D4h)
        .value("C3", cfpars::Sym::C3).value("S6", cfpars::Sym::S6)
        .value("D3", cfpars::Sym::D3).value("C3v", cfpars::Sym::C3v).value("D3d", cfpars::Sym::D3d)
        .value("C6", cfpars::Sym::C6).value("C3h", cfpars::Sym::C3h).value("C6h", cfpars::Sym::C6h)
        .value("D6", cfpars::Sym::D6).value("C6v", cfpars::Sym::C6v).value("D3h", cfpars::Sym::D3h) .value("D6h", cfpars::Sym::D6h)
        .value("T", cfpars::Sym::T).value("Th", cfpars::Sym::Th).value("Td", cfpars::Sym::Td).value("O", cfpars::Sym::O).value("Oh", cfpars::Sym::Oh);

    pycfpars.def(py::init<>())
        .def(py::init<const std::string &>(), py::arg("ionname"))
        .def(py::init<const double &>(), py::arg("J"))
        .def(py::init(&cfpars_init), cfpars_init_str)
        .def("get_sym_allowed_pars", [](cfpars &self) { std::vector<cfpars::Blm> r0 = self.get_sym_allowed_Blm(); std::vector<std::string> rv(r0.size());
            std::transform(r0.begin(), r0.end(), rv.begin(), [](cfpars::Blm b) { return Blmvec[(int)b]; }); return rv; })
        .def_property("unit", &cfpars::get_unit, [](cfpars &self, std::string unit) { self.set_unit(set_enum(unit, unit_names, unit_err)); },
            "energy unit of the parameters")
        .def_property("type", &cfpars::get_type, [](cfpars &self, std::string type) { self.set_type(set_enum(type, type_names, type_err)); },
            "type of CF parameters. When changed, parameters values will be automatically converted.")
        .def_property("sym", &cfpars::get_sym, [](cfpars &self, std::string sym) { self.set_sym(set_enum(sym, sym_names, sym_err)); },
            "point symmetry of site in Schoenflies notation")
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
}


