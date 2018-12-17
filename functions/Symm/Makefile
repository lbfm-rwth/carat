CC = gcc
CFLAGS = -fwritable-strings -DDIAG1 -g
#TOPDIR = /usb/carat

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) -I$(INCL)

OBJS = dsylv.o\
       rest_short.o\
       short.o\
       shortest.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
