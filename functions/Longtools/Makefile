#TOPDIR = /usb/carat
CC = gcc
CFLAGS = -DDIAG1 -fwritable-strings -g

INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) -c -I$(INCL)

OBJS = MP_conv_mat.o\
       MP_gauss.o\
       MP_hnf.o\
       MP_pair_red.o\
       MP_red_sort.o\
       MP_solve.o\
       dump_MP_mat.o\
       long_elt.o\
       long_gauss.o\
       long_hnf.o\
       long_kernel_mat.o\
       long_mat_inv.o\
       long_qbase.o\
       long_rein_mat.o\
       long_solve_mat.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
