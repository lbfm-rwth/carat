TOPDIR =/usb/carat
CFLAGS = -g -fwritable-strings -DDIAG1
CC = gcc

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rvuc $(TOPDIR)/lib/libpresentation.a

GLOBAL = -DTOPDIR=\"$(TOPDIR)\" -DTABLES=\"$(TOPDIR)/tables/\" \
         -DATOMS=\"$(TOPDIR)/tables/atoms/\" \
         -DTABLEDIM=\"$(TOPDIR)/tables/dim\"

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) $(GLOBAL) -I$(INCL)

OBJS = check_base.o\
       mapped_word.o\
       presentation.o\
       put_word.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
