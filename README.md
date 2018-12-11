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

Edit the Makefile in the directory you installed CARAT in, and
change the variables
 
     TOPDIR (which is the output of pwd when in the directory this makefile is stored in)
     CFLAGS

according to your needs. Some additional explanations/examples are given
there.

For Mac users only: You need to delete the line 

    #include <malloc.h>
    
from `include/typedef.h`. This will be fixed once there is 
a `configure` script.

Afterwards do

    make

or to be sure everything is made from scratch:

    make clean ; make 

If you have any problems installing CARAT, please feel free to contact
me at carat@momo.math.rwth-aachen.de

Faithfully yours

Dominik Bernhardt


