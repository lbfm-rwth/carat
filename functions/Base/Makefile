#TOPDIR = /usb/carat
CFLAGS = -g -fwritable-strings -DDIAG1
CC = gcc

INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) -c -I$(INCL)

OBJS = base.o\
       base2.o\
       conjugate.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o

