
OBJFLAGS = -c
INCL = $(TOPDIR)/include

AR = ar rDc $(TOPDIR)/lib/functions.a

COMP = $(CC) $(CFLAGS) $(OBJFLAGS) -I$(INCL)

OBJS = chin_remainder.o\
       intpow.o\
       itoa.o\
       malloc2dim.o\
       min_div.o\
       ovfl_mul.o\
       prime_tools.o\
       random.o\
       tools.o

.c.o:
	$(COMP) $< -o $@

ALL:  $(OBJS)
	$(AR) *.o

clean:
	rm -f *.o

strip:
	strip *.o
