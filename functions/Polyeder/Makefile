CFLAGS = -DDIAG1 -fwritable-strings -g
CC = gcc
#TOPDIR = /usb/carat

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rvuc $(TOPDIR)/lib/functions.a

GLOBAL = -DTOPDIR=\"$(TOPDIR)\" -DTABLES=\"$(TOPDIR)/tables/\" \
         -DATOMS=\"$(TOPDIR)/tables/atoms/\" \
         -DTABLEDIM=\"$(TOPDIR)/tables/dim\"

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) $(GLOBAL) -I$(INCL)

OBJS = first_fuber.o\
       first_polyeder.o\
       fuber_tools.o\
       polyeder_to_vecs.o\
       polyeder_tools.o\
       refine_fuber.o\
       refine_polyeder.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
