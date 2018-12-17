CC = gcc
#TOPDIR = /usb/carat/
CFLAGS = -fwritable-strings -DDIAG1 -g

INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) -c -I$(INCL)

OBJS = Q_data.o\
	graph_tools.o\
	normalop_o.o\
	subgroupgraph.o\
	stabs.o\
	orbit_fkts.o\
	coho_size.o\
	lattices.o\
	phi.o\
	graph_mapped_word.o\
	aff_normal.o\
	super-k-groups-fcts.o\
	sub-super-tools.o\
	sub-k-groups-fcts.o\
	is_k_subgroup.o\
	calculate_representatives.o\
	sub-t-groups-fcts.o\
	obergruppenzahl.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
