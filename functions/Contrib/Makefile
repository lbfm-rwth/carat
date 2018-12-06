CC = gcc
#TOPDIR = /usb/carat
CFLAGS = -fwritable-strings -DDIAG1 -g


INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) -c -I$(INCL)

OBJS = suche_kand.o\
       torsionfree.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
