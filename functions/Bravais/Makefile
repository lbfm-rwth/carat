#TOPDIR = /usb/carat
CFLAGS = -g -fwritable-strings -DDIAG1
CC = gcc

INCL = $(TOPDIR)/include

AR = ar rvuc $(TOPDIR)/lib/functions.a

GLOBAL = -DTOPDIR=\"$(TOPDIR)\" -DTABLES=\"$(TOPDIR)/tables/\" \
         -DATOMS=\"$(TOPDIR)/tables/atoms/\" \
         -DTABLEDIM=\"$(TOPDIR)/tables/dim\"

COMP = $(CC) $(CFLAGS) -c $(GLOBAL) -I$(INCL)

OBJS = bravais_tools.o\
       formspace.o\
       init_bravais.o\
       invar_space.o\
       konj_bravais.o\
       lincomb.o\
       normlin.o\
       p_formspace.o\
       rform.o\
       tr_bravais.o\
       trace_bifo.o\
       normalisator.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o

