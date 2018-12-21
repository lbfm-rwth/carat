CC = gcc
#TOPDIR = /usb/carat
CFLAGS = -fwritable-strings -DDIAG1 -g

INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) -c -I$(INCL)

OBJS = PoaN.o\
	name_fct.o\
	tsubgroups_fct.o\
	tsubgroups_fct_db.o\
	tsubgroups_tools.o\
	tSUPERgroups_fct_db.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
