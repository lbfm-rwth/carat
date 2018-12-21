CC = gcc
CFLAGS = -g -fwritable-strings -DDIAG1
#TOPDIR = /usb/carat

INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) -c -I$(INCL)

OBJS = almost_decomposable_lattice.o\
       bravais_catalog.o\
       centr.o\
       min_pol.o\
       symbol.o\
       v4_catalog.o\
       z_equivalent.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
