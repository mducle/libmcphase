/* pycfpars.hpp
 *
 * Headers for python bindings for the cfpars class of libMcPhase
 *
 * (C) 2020 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#pragma once

template <typename T> T set_enum(std::string key, std::unordered_map<std::string, T> enum_map, std::string errmsg) {
    auto it = enum_map.find(key);
    if (it == enum_map.end()) {
        throw std::runtime_error(errmsg);
    } else {
        return it->second;
    }
}

static const char* cfpars_init_str = "Construct a cfpars object\n"
                                     "    args (one of):\n"
                                     "        J - total angular momentum of multiplet\n"
                                     "        ionname - name of ion, e.g. 'Pr3+'\n"
                                     "    kwargs:\n"
                                     "        unit - energy units of parameters. one of 'meV' (default), 'cm' or 'K'\n"
                                     "        type - type of parameters one of 'Alm', 'ARlm', 'Blm' (default), 'Llm' 'Vlm', or 'Wlm'\n"
                                     "               see online documentation or McPhase webpage for definitions\n"
                                     "        Blm - value of parameters, e.g. cfp = cfpars('Pr3+', B20=0.1, B22=-0.01, B40=0.001)\n";

