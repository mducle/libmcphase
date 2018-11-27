# libMcPhase

libMcPhase is a fork of the [McPhase](http://mcphase.de) 
[project](https://github.com/mducle/mcphase) intended to be a partial rewrite of 
the code as a Python library.

McPhase was originally developed as a series of C, C++, Perl and Basic programs.
The C/C++04 codes used the [MATPACK](http://www.matpack.de/) and 
[LAPACK](http://www.netlib.org/lapack/) linear algebra libraries.

libMcPhase will update the code to C++11 and use the 
[Eigen](http://eigen.tuxfamily.org) header library instead of MATPACK/LAPACK.
[PyBind11](https://github.com/pybind/pybind11) is used as the python binding.

Initially only the single-ion calculation modules `so1ion`, `ic1ion` and `cluster`
will be ported to C+11.

