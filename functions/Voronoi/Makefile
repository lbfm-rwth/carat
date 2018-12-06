#TOPDIR = /usb/carat
CFLAGS = -g -fwritable-strings -DDIAG1
CC = gcc

INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) -c -I$(INCL)

OBJS = all_vor_neighbours.o\
       bravais_flok.o\
       bravais_flok_datei.o\
       calc_vor_data.o\
       first_perfect.o\
       init_voronoi.o\
       normalizer.o\
       pair_red_inv.o\
       red_normal.o\
       vor_neighbour.o\
       vor_vertices.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
