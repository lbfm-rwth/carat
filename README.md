[![Travis build Status](https://travis-ci.org/lbfm-rwth/carat.svg?branch=master)](https://travis-ci.org/lbfm-rwth/carat)
[![AppVeyor build Status](https://ci.appveyor.com/api/projects/status/github/lbfm-rwth/carat?branch=master&svg=true)](https://ci.appveyor.com/project/lbfm-rwth/carat)
[![Code Coverage](https://codecov.io/github/lbfm-rwth/carat/coverage.svg?branch=master&token=)](https://codecov.io/gh/lbfm-rwth/carat)

# CARAT 

CARAT is a package for solving certain problems in mathematical
crystallography.

It is distributed via

     Lehrstuhl B fuer Mathematik
     RWTH-Aachen
     Prof. i.R. Plesken
     Pontdriesch 10-16
     52064 Aachen
     Germany
     email: carat@momo.math.rwth-aachen.de

NOTE: CARAT was developed for crystallographic groups in dimensions up to 6.
Most algorithms also work in higher dimensions. However, integer overflow
is not trapped in general.

## Website

You can find the official CARAT homepage and documentation on
<https://lbfm-rwth.github.io/carat/>.

## Installation

If you want to compile CARAT directly from GitHub, change
into the base directory (the one containing this README.md)
and issue the command

    ./autogen.sh

Please note that this requires you to have autoconf/automake
installed. 
Afterwards, or if you are using a released version of CARAT,
enter the following commands in the CARAT base directory:

     ./configure && make

CARAT relies on the GMP library (along with its header files) being
installed at a standard location, so that the compiler and linker will
find it. Otherwise, you will have to indicate their location with
appropriate compiler/linker flags.

In order to clean up from a previous installation, you may want to to a

     make clean

before building CARAT again.


## Bug reports and feature requests

Please submit bug reports and support requests via our GitHub issue tracker:

  <https://github.com/lbfm-rwth/carat/issues>


## License

CARAT is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

For details see the file LICENSE.
