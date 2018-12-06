
OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) -I$(INCL)

OBJS = autgrp.o\
       isometry.o\
       pr_aut.o\
       pr_isom.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
