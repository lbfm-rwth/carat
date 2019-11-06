CC = gcc
#TOPDIR = /usb/carat
CFLAGS = -fwritable-strings -g

INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) -c -I$(INCL)

OBJS = ZZ.o\
       ZZ_cen_fun.o\
       ZZ_gen_vs.o\
       ZZ_irr_const.o\
       ZZ_lll.o\
       ZZ_zclass.o\
       get_groups.o\
       q2z.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
