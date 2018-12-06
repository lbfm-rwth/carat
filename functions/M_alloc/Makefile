#TOPDIR =/usb/carat
CFLAGS = -g -fwritable-strings -DDIAG1
CC = gcc

OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/libm_alloc.a

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) -I$(INCL)

OBJS = m_alloc.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
