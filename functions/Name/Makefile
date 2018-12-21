CC = gcc
CFLAGS = -g -fwritable-strings -DDIAG1
#TOPDIR = /usb/carat

INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) -c -I$(INCL)

OBJS = HM_symbol.o\
       Q_catalog.o\
       aff_class_inf.o\
       compute_q_matrix.o\
       point_group.o\
       q_class_inf.o\
       z_class_inf.o\
	reverse_name_fct.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
