TOPDIR = /usb/carat
CC = gcc
CFLAGS = -fwritable-strings -DDIAG1 -g

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rvuc $(TOPDIR)/lib/functions.a

GLOBAL = -DTOPDIR=\"$(TOPDIR)\" -DTABLES=\"$(TOPDIR)/tables/\" \
         -DATOMS=\"$(TOPDIR)/tables/atoms/\" \
         -DTABLEDIM=\"$(TOPDIR)/tables/dim\"

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) $(GLOBAL) -I$(INCL)

OBJS = brav_from_datei.o\
       free_bravais.o\
       get_symbol.o\
       get_zentr.o\
       gittstab.o\
       gittstabneu.o\
       lattice.o\
       lattice_tools.o\
       read_symbol.o\
       right_order.o\
       super_lattice.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o

