TOPDIR =/usb/carat
CFLAGS = -g -fwritable-strings -DDIAG1
CC = gcc

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rvuc $(TOPDIR)/lib/libpresentation.a

GLOBAL = -DTOPDIR=\"$(TOPDIR)\" -DTABLES=\"$(TOPDIR)/tables/\" \
         -DATOMS=\"$(TOPDIR)/tables/atoms/\" \
         -DTABLEDIM=\"$(TOPDIR)/tables/dim\"

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) $(GLOBAL) -I$(INCL)

OBJS = Back.o\
       barzen.o\
       KKGV.o\
       Pres.o\
       T_polyeder.o\
       T_ungleich_wand.o\
       barzen.o\
       easy_torus_orbit.o\
       fub.o\
       fub_tools.o\
       modshort.o\
       neugram.o\
       neukandidat.o\
       presentation_point_grp.o\
       prog.o\
       sscal.o\
       stab_orbit.o\
       tietze_tools.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
