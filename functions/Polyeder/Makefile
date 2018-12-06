CFLAGS = -DDIAG1 -fwritable-strings -g
CC = gcc
#TOPDIR = /usb/carat

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) -I$(INCL)

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
