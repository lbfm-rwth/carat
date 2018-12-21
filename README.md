[![Travis build Status](https://travis-ci.org/lbfm-rwth/carat.svg?branch=master)](https://travis-ci.org/lbfm-rwth/carat)
[![AppVeyor build Status](https://ci.appveyor.com/api/projects/status/github/lbfm-rwth/carat?branch=master&svg=true)](https://ci.appveyor.com/project/lbfm-rwth/carat)
[![Code Coverage](https://codecov.io/github/lbfm-rwth/carat/coverage.svg?branch=master&token=)](https://codecov.io/gh/lbfm-rwth/carat)

# CARAT 

Carat is a package for solving certain problems in mathematical
crystallography.

It is distributed via

     Lehrstuhl B fuer Mathematik
     RWTH-Aachen
     Prof. Plesken
     Pontdriesch 10-16
     52064 Aachen
     Germany
     email: carat@momo.math.rwth-aachen.de

NOTE: CARAT was developed for crystallographic groups in dimensions up to 6.
      Most algorithms also work in higher dimensions. However, integer overflow
      is not trapped in general.

INSTALLATION:

The easiest way to compile CARAT is to go to the CARAT base directory
(the one containing this REAMDME), and issue the make command with some
variables set on the command line, to order override their defaults in
the Makefile:

     make TOPDIR=`pwd` CFLAGS="-O -g" CC="gcc"

You may of course choose your own compiler and compiler flags. The value
of TOPDIR must be the CARAT base directory.

CARAT relies on the GMP library (along with its header files) being
installed at a standard location, so that the compiler and linker will
find it. Otherwise, you will have to indicate their location with
appropriate compiler/linker flags.

In order to clean up from a previous installation, you may want to to a

     make clean

before building CARAT again.

If you have any problems installing CARAT, please feel free to contact
me at carat@momo.math.rwth-aachen.de

Faithfully yours

Dominik Bernhardt


