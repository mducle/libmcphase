### [v0.1.3](https://github.com/mducle/libmcphase/compare/v0.1.2...v0.1.3)

Library / dependencies update

* Update build system - now supports MacOS (Intel and ARM).
* Cleans up some compiler warnings

Example scripts (after installing with `pip install libmcphase`)

```
import libmcphase

cfp = libmcphase.cf1ion('Pr3+', B20=0.1, B40=0.01, B60=0.001, B66=-0.02)
V, E = cfp.eigensystem()
print(E)
```


### [v0.1.2](https://github.com/mducle/libmcphase/compare/v0.1.0...v0.1.2)

Pre-alpha release with only single-ion calculations of energy levels and some physical properties (magnetisation, susceptibility).

Very little documentation.


### [v0.1.0](https://github.com/mducle/libmcphase/compare/08c4a2e22b61c145b841ca531c7a11a3959d10fa...v0.1.0)

Initial commit of libmcphase.

* Reorganised McPhase code and changed directory structure.
* Rewrite of `cf1ion` module to remove old Fabi code and use 3j/6j symbols.
* Python interface using Pybind11.
