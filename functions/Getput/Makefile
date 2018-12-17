CFLAGS = -g -DDIAG1 -fwritable-strings
#TOPDIR = /usb/carat
CC = gcc

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) -I$(INCL)

OBJS = get_bravais.o\
       get_mat.o\
       put_bravais.o\
       put_mat.o\
       put_order.o\
       read_header.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o

