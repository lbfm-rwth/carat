#TOPDIR=/usb/carat
CFLAGS= -g -fwritable-strings -DDIAG1
CC = gcc

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) -I$(INCL)

OBJS = add_mat.o\
       col_row_ops_mat.o\
       comp_mat.o\
       construct_mat.o\
       divide_by_gcd.o\
       elt_div_mat.o\
       find_max_entry_mat.o\
       gauss_mat.o\
       inv_mat.o\
       kgv2rat_mat.o\
       kron_mat.o\
       modp_mat.o\
       mul_mat.o\
       null_mat.o\
       p_gauss_mat.o\
       p_lse_solve.o\
       p_mat_det.o\
       p_solve_mat.o\
       real_mat.o\
       red_mat.o\
       scal_pr_mat.o\
       solve_mat.o\
       tools_mat.o\
       tr_pose_mat.o\
       trace_mat.o\
       unity_mat.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
