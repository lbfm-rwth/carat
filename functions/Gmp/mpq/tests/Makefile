# This file was generated automatically by configure.  Do not edit.
VPATH = .
links = 
host_alias = hppa1.1-hp-hpux10.20
host_cpu = hppa1.1
host_vendor = hp
host_os = hpux10.20
host_canonical = hppa1.1-hp-hpux10.20
target_alias = hppa1.1-hp-hpux10.20
target_cpu = hppa1.1
target_vendor = hp
target_os = hpux10.20
target_canonical = hppa1.1-hp-hpux10.20
build_alias = hppa1.1-hp-hpux10.20
build_cpu = hppa1.1
build_vendor = hp
build_os = hpux10.20
build_canonical = hppa1.1-hp-hpux10.20
# Makefile for mpq/tests for GNU MP

srcdir = .

CC = gcc -O2

TEST_LIBS = ../../libgmp.a
INCLUDES = -I../../mpn -I$(srcdir)/../..
CFLAGS = -g -O

.c.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) $(XCFLAGS) $<

TEST_SRCS = t-cmp.c t-cmp_ui.c t-get_d.c
TEST_OBJS = t-cmp.o t-cmp_ui.o t-get_d.o
TESTS = t-cmp t-cmp_ui t-get_d

check: Makefile st-cmp st-cmp_ui st-get_d
	@echo "The tests passed."

st-cmp: t-cmp
	./t-cmp
	touch $@
st-cmp_ui: t-cmp_ui
	./t-cmp_ui
	touch $@
st-get_d: t-get_d
	./t-get_d
	touch $@

t-cmp: t-cmp.o $(TEST_LIBS)
	$(CC) -o $@ $@.o $(TEST_LIBS) $(CFLAGS)
t-cmp_ui: t-cmp_ui.o $(TEST_LIBS)
	$(CC) -o $@ $@.o $(TEST_LIBS) $(CFLAGS)
t-get_d: t-get_d.o $(TEST_LIBS)
	$(CC) -o $@ $@.o $(TEST_LIBS) $(CFLAGS)

t-cmp.o: $(srcdir)/t-cmp.c
t-cmp_ui.o: $(srcdir)/t-cmp_ui.c
t-get_d.o: $(srcdir)/t-get_d.c

clean mostlyclean:
	rm -f *.o st-* $(TESTS)
distclean maintainer-clean: clean
	rm -f Makefile config.status

Makefile: $(srcdir)/Makefile.in
	$(SHELL) ./config.status
