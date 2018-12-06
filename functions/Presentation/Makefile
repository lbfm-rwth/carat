#TOPDIR =/usb/carat
CFLAGS = -g -fwritable-strings -DDIAG1
CC = gcc

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/libpresentation.a

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) -I$(INCL)

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
