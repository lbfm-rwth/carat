############################################################################
# Makefile for the CARAT package, distributed by
#    Lehrstuhl B fuer Mathematik
#    RWTH Aachen
#    Templergraben 64
#    52064 Aachen
#    Germany
#    email: carat@momo.math.rwth-aachen.de
#
# To install the CARAT package, just edit the variables
# TOPDIR (which is the present directory), CC and CFLAGS
# below.
# If you need additional include pathes, please add them
# to the compiler, ie CC= gcc -I/YOUR_INCUDE_PATH
#
# If you encounter any problems, please contact
#    carat@momo.math.rwth-aachen.de
############################################################################


TOPDIR= /usb/carat
CC = gcc

# There are some special preprocessor flags which set some
# memory diagnostics:
# Use -DDIAG1 to check for not properly used memory only
# when calling malloc(...) and free(..), and -DDIAG2 to
# have a general control whats going on. THIS IS VERY SLOW.
# For the normal user we recommend neither to use -DDIAG1 nor -DDIAG2!

CFLAGS = -g -Wall -DDIAG1 -fwritable-strings # for a HP-UX-machine using gcc (momo)
                                       # the flag -fwritable-strings is
                                       # required for the use with gcc

# CFLAGS = -g -Aa                        # for a HP-UX-machine using cc

# CFLAGS = -m486 -O2                     # on a Linux machine (i486)


# The part below doesn't (better: shouldn't) need any editing
#============================================================================

# Some SUN-OS will need ranlib to be run on the library.
# This part is stolen from the gmp-Makefile
RANLIB_TEST = [ -f /usr/bin/ranlib -o -f /bin/ranlib ]
RANLIB = ranlib


ALL: Makefile\
     Autgrp\
     Base\
     Bravais\
     Contrib\
     Datei\
     Getput\
     Hyperbolic\
     Idem\
     Links\
     Longtools\
     Name\
     Matrix\
     M_alloc\
     Orbit\
     Polyeder\
     Presentation\
     Qcatalog\
     Reduction\
     Sort\
     Symm\
     Tools\
     Voronoi\
     Zassen\
     ZZ\
     Gmp\
     Executables

Autgrp: Makefile functions/Autgrp/Makefile
	cd functions/Autgrp; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Base: Makefile functions/Base/Makefile
	cd functions/Base; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Bravais: Makefile functions/Bravais/Makefile
	cd functions/Bravais;make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Contrib: Makefile functions/Contrib/Makefile
	cd functions/Contrib;make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Datei: Makefile functions/Datei/Makefile
	cd functions/Datei;make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Getput: Makefile functions/Getput/Makefile
	cd functions/Getput;make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Hyperbolic: Makefile functions/Hyperbolic/Makefile
	cd functions/Hyperbolic; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Idem: Makefile functions/Idem/Makefile
	cd functions/Idem; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Links:
	cd tables/symbol; make
	cd lib ; rm -f libfunctions.a ; ln -f -s functions.a libfunctions.a
	rm -f $(TOPDIR)/functions/Gmp/m_alloc.h
	ln -f -s $(TOPDIR)/include/m_alloc.h $(TOPDIR)/functions/Gmp/m_alloc.h

Longtools: Makefile functions/Longtools/Makefile
	cd functions/Longtools; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Name: Makefile functions/Name/Makefile
	cd functions/Name;make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Matrix: Makefile functions/Matrix/Makefile
	cd functions/Matrix;make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

M_alloc: Makefile functions/M_alloc/Makefile
	cd functions/M_alloc;make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Orbit: Makefile functions/Orbit/Makefile
	cd functions/Orbit; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Polyeder: Makefile functions/Polyeder/Makefile
	cd functions/Polyeder; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Presentation: Makefile functions/Presentation/Makefile
	cd functions/Presentation; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Qcatalog: Makefile tables/qcatalog.tar.gz
	cd tables; if [ !  -d qcatalog ] ; then tar xvzf qcatalog.tar.gz ; fi

Reduction: Makefile functions/Reduction/Makefile
	cd functions/Reduction; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Sort: Makefile functions/Sort/Makefile
	cd functions/Sort; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Symm: Makefile functions/Symm/Makefile
	cd functions/Symm; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Tools: Makefile functions/Tools/Makefile
	cd functions/Tools; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Voronoi: Makefile functions/Voronoi/Makefile
	cd functions/Voronoi; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Zassen: Makefile functions/Zassen/Makefile
	cd functions/Zassen; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

ZZ: Makefile functions/ZZ/Makefile
	cd functions/ZZ; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

Gmp: functions/Gmp/Makefile
	cd functions/Gmp ; ./configure --prefix=../.. ; make CFLAGS="$(CFLAGS)" CC="$(CC)" install

Executables: bin/Makefile
	if $(RANLIB_TEST) ; then $(RANLIB) lib/functions.a; else true; fi
	if $(RANLIB_TEST) ; then $(RANLIB) lib/libpresentation.a; else true; fi
	if $(RANLIB_TEST) ; then $(RANLIB) lib/libm_alloc.a; else true; fi
	cd bin; make CC="$(CC)" CFLAGS="$(CFLAGS)" TOPDIR=$(TOPDIR)

clean:
#	cd bin; make clean
	cd functions/Autgrp; make clean
	cd functions/Base; make clean
	cd functions/Bravais; make clean
	cd functions/Contrib; make clean
	cd functions/Datei; make clean
	cd functions/Getput; make clean
	cd functions/Hyperbolic; make clean
	cd functions/Idem; make clean
	cd functions/Longtools; make clean
	cd functions/Name; make clean
	cd functions/Matrix; make clean
	cd functions/M_alloc; make clean
	cd functions/Orbit; make clean
	cd functions/Polyeder; make clean
	cd functions/Presentation; make clean
	cd functions/Reduction; make clean
	cd functions/Sort; make clean
	cd functions/Symm; make clean
	cd functions/Tools; make clean
	cd functions/Voronoi; make clean
	cd functions/Zassen; make clean
	cd functions/ZZ; make clean
	cd functions/Gmp; make clean
	rm -f lib/libgmp.a
	rm -f lib/libm_alloc.a
	rm -f lib/libpresentation.a
	rm -f lib/functions.a
