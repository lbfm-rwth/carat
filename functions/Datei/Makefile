#TOPDIR = /usb/carat
CC = gcc
CFLAGS = -fwritable-strings -DDIAG1 -g

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) -I$(INCL)

OBJS = brav_from_datei.o\
       free_bravais.o\
       get_symbol.o\
       get_zentr.o\
       gittstab.o\
       gittstabneu.o\
       lattice.o\
       lattice_tools.o\
       read_symbol.o\
       right_order.o\
       super_lattice.o\
       get_data_dir.o

ALL:  $(OBJS)
	$(AR) *.o

get_data_dir.o:  get_data_dir.c
	$(COMP) -DTOPDIR=\"$(TOPDIR)\" $<

%.o: %.c
	$(COMP) $< -o $@

clean:
	rm -f *.o

strip:
	strip *.o

