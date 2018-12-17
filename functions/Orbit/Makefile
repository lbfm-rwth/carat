CC = gcc
#TOPDIR = /usb/carat
CFLAGS = -g -fwritable-strings -DDIAG1

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) -I$(INCL)

OBJS = orb_division.o\
       orbit_alg.o\
       orbit_subdivision.o\
       row_spin.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
