INCL = $(TOPDIR)/include

AR = ar rvuc $(TOPDIR)/lib/libpresentation.a

GLOBAL = -DTOPDIR=\"$(TOPDIR)\" -DTABLES=\"$(TOPDIR)/tables/\" \
         -DATOMS=\"$(TOPDIR)/tables/atoms/\" \
         -DTABLEDIM=\"$(TOPDIR)/tables/dim\"

COMP = $(CC) $(CFLAGS) -c $(GLOBAL) -I$(INCL)

OBJS = calc_orders.o\
       carat_out.o\
       cayley_out.o\
       delta.o\
       derived_rc.o\
       dumprelators.o\
       dumprelators_file.o\
       exist.o\
       fundi.o\
       gap_out.o\
       generate_pres.o\
       inverse_prods.o\
       is_equal.o\
       is_identity.o\
       magma_out.o\
       normalize.o\
       old_as_new.o\
       shorter.o\
       sub_relator.o\
       substit.o\
       tietzetrans_rc.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
