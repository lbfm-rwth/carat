CC = gcc
TOPDIR = /usb/carat
CFLAGS = -fwritable-strings -DDIAG1 -g

INCL = $(TOPDIR)/include

AR = ar rvuc $(TOPDIR)/lib/libgraph.a

GLOBAL = -DTOPDIR=\"$(TOPDIR)\" -DTABLES=\"$(TOPDIR)/tables/\" \
         -DATOMS=\"$(TOPDIR)/tables/atoms/\" \
         -DTABLEDIM=\"$(TOPDIR)/tables/dim\"

COMP = $(CC) $(CFLAGS) -c $(GLOBAL) -I$(INCL)

OBJS =  Q_data.o\
	graph_tools.o\
	normalop_o.o\
	subgroupgraph.o\
	stabs.o\
	orbit_fkts.o\
	coho_size.o\
	phi.o\
	graph_mapped_word.o\
	aff_normal.o\
	super-k-groups-fcts.o\
	sub-super-tools.o\
	sub-k-groups-fcts.o\
	is_k_subgroup.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
