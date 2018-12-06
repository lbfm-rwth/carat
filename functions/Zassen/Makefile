CC = gcc
#TOPDIR = /usb/carat
CFLAGS = -fwritable-strings -DDIAG1 -g


INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) -c -I$(INCL)

OBJS = coboundary.o\
       cong_solve.o\
       convert_cocycle_to_column.o\
       convert_to_cozycle.o\
       put_cocycle.o\
       reget_gen.o\
       reverse_valuation.o\
       normalop.o\
       zass.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
