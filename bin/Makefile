CC = gcc
CFLAGS = -fwritable-strings -DDIAG1 -pg
# TOPDIR = /usb/carat

# LIB = -L$(TOPDIR)/lib -lcarat -lm
LIB = -L$(TOPDIR)/lib -lpresentation -lfunctions -lpresentation -lgmp -lm_alloc -lm
SRC = $(TOPDIR)/src
INCL = $(TOPDIR)/include
OBJ = $(TOPDIR)/obj
LIBS = $(TOPDIR)/lib/libpresentation.a $(TOPDIR)/lib/functions.a  $(TOPDIR)/lib/libm_alloc.a

COMP = $(CC) $(CFLAGS) -I$(INCL)

ALL: DIR

#
#   PROGRAMS
#

PROGRAMS: Add Aut_grp Bravais_equiv\
          Bravais_grp Bravais_inclusions\
          Conj Conjugated Conj_bravais\
          Conv Datei\
          Elt Extract Extensions\
          First_perfect Form_space Form_elt\
          Formtovec Full Gauss Graph\
          Idem Inv Invar_space Isometry Is_finite\
          KSubgroups KSupergroups\
          Kron Long_solve Ltm\
          Modp Mtl Mul Mink_red Minpol\
          Name Normalizer Normalizer_in_N Normlin\
          Orbit Order\
          Pair_red Pdet Perfect_neighbours Polyeder P_lse_solve\
          Presentation Q_catalog QtoZ Red_gen Rest_short\
          Rein Rform Reverse_name\
          Same_generators Scalarmul Scpr Short\
          Shortest Short_reduce Signature\
          Simplify_mat Symbol\
	  TSubgroups TSupergroups\
          Torsionfree Tr Trace Trbifo Tr_bravais\
          Vor_vertices Vectoform ZZprog Zass_main Z_equiv


DIR: ./Makefile ./config.guess
	if [ ! -d `./config.guess`-`basename $(CC)` ] ; then \
	   mkdir `./config.guess`-`basename $(CC)` ; \
	fi
	cp Makefile `./config.guess`-`basename $(CC)`
	cd `./config.guess`-`basename $(CC)` ;\
         make TOPDIR="$(TOPDIR)" CC="$(CC)" CFLAGS="$(CFLAGS)" PROGRAMS

#-------------------------------------------------------------

Add:  $(SRC)/add.c $(LIBS)
	$(COMP) -o Add \
        $(SRC)/add.c  $(LIB)

Aut_grp:  $(SRC)/autgrp.c $(LIBS)
	$(COMP) -o Aut_grp \
        $(SRC)/autgrp.c  $(LIB)

Bravais_equiv:  $(SRC)/bravais_equiv.c $(LIBS)
	$(COMP) -o Bravais_equiv \
        $(SRC)/bravais_equiv.c  $(LIB)

Bravais_grp:  $(SRC)/bravaisgroup.c $(LIBS)
	$(COMP) -o Bravais_grp \
        $(SRC)/bravaisgroup.c  $(LIB)

Bravais_inclusions:  $(SRC)/bravais_inclusions.c $(LIBS)
	$(COMP) -o Bravais_inclusions \
        $(SRC)/bravais_inclusions.c  $(LIB)

Conj:  $(SRC)/conj.c $(LIBS)
	$(COMP) -o Conj \
        $(SRC)/conj.c  $(LIB)

Conjugated:  $(SRC)/conjugated.c $(LIBS)
	$(COMP) -o Conjugated \
        $(SRC)/conjugated.c  $(LIB)

Conj_bravais:  $(SRC)/conj_bravais.c $(LIBS)
	$(COMP) -o Conj_bravais \
        $(SRC)/conj_bravais.c  $(LIB)

Conv:  $(SRC)/conv.c $(LIBS)
	$(COMP) -o Conv \
        $(SRC)/conv.c  $(LIB)

Datei:  $(SRC)/datei.c $(LIBS)
	$(COMP) -o Datei \
        $(SRC)/datei.c  $(LIB)
	rm -f Bravais_catalog
	ln -f -s Datei Bravais_catalog

Elt:  $(SRC)/elt.c $(LIBS)
	$(COMP) -o Elt \
        $(SRC)/elt.c  $(LIB)

Extract:  $(SRC)/extract.c $(LIBS)
	$(COMP) -o Extract \
        $(SRC)/extract.c  $(LIB)
	rm -f Standard_affine_form
	ln -f -s Extract Standard_affine_form

Extensions:  $(SRC)/extensions.c $(LIBS)
	$(COMP) -o Extensions \
        $(SRC)/extensions.c  $(LIB)
	rm -f Vector_systems
	ln -f -s Extensions Vector_systems

First_perfect:  $(SRC)/first_perfect.c $(LIBS)
	$(COMP) -o First_perfect \
        $(SRC)/first_perfect.c  $(LIB)

Form_space:  $(SRC)/form_space.c $(LIBS)
	$(COMP) -o Form_space \
        $(SRC)/form_space.c  $(LIB)

Form_elt: $(SRC)/form_elt.c $(LIBS)
	$(COMP) -o Form_elt \
        $(SRC)/form_elt.c  $(LIB)

Formtovec:  $(SRC)/formtovec.c $(LIBS)
	$(COMP) -o Formtovec \
        $(SRC)/formtovec.c  $(LIB)

Full:  $(SRC)/full.c $(LIBS)
	$(COMP) -o Full \
        $(SRC)/full.c  $(LIB)

#Fundamental_domain:  $(SRC)/fundamental_domain.c $(LIBS)
#	$(COMP) -o Fundamental_domain \
#        $(SRC)/fundamental_domain.c  $(LIB)

Gauss:  $(SRC)/gauss.c $(LIBS)
	$(COMP) -o Gauss \
        $(SRC)/gauss.c  $(LIB)

Graph:  $(SRC)/graph.c $(LIBS)
	$(COMP) -o Graph \
        $(SRC)/graph.c  $(LIB)

#Gittstab:  $(SRC)/gittstab.c $(LIBS)
#	$(COMP) -o Gittstab \
#        $(SRC)/gittstab.c  $(LIB)

#Hypisom:  $(SRC)/hypisom.c $(LIBS)
#	$(COMP) -o Hypisom \
#        $(SRC)/hypisom.c  $(LIB)
#
#Hypstab:  $(SRC)/hypstab.c $(LIBS)
#	$(COMP) -o Hypstab \
#        $(SRC)/hypstab.c  $(LIB)

Idem:  $(SRC)/idem.c $(LIBS)
	$(COMP) -o Idem \
        $(SRC)/idem.c  $(LIB)

Inv:  $(SRC)/inv.c $(LIBS)
	$(COMP) -o Inv \
        $(SRC)/inv.c  $(LIB)

Invar_space:  $(SRC)/invar_space.c $(LIBS)
	$(COMP) -o Invar_space \
        $(SRC)/invar_space.c  $(LIB)

Isometry:  $(SRC)/isometry.c $(LIBS)
	$(COMP) -o Isometry \
        $(SRC)/isometry.c  $(LIB)

Is_finite:  $(SRC)/is_finite.c $(LIBS)
	$(COMP) -o Is_finite \
        $(SRC)/is_finite.c  $(LIB)

Kron:  $(SRC)/kron.c $(LIBS)
	$(COMP) -o Kron \
        $(SRC)/kron.c  $(LIB)

KSubgroups:  $(SRC)/sub-k-groups.c $(LIBS)
	$(COMP) -o KSubgroups \
        $(SRC)/sub-k-groups.c  $(LIB)

KSupergroups:  $(SRC)/super-k-groups.c $(LIBS)
	$(COMP) -o KSupergroups \
        $(SRC)/super-k-groups.c  $(LIB)

Long_solve:  $(SRC)/long_solve.c $(LIBS)
	$(COMP) -o Long_solve \
        $(SRC)/long_solve.c  $(LIB)

Ltm:  $(SRC)/ltm.c $(LIBS)
	$(COMP) -o Ltm \
        $(SRC)/ltm.c  $(LIB)

Modp:  $(SRC)/modp.c $(LIBS)
	$(COMP) -o Modp \
        $(SRC)/modp.c  $(LIB)

Mtl:  $(SRC)/mtl.c $(LIBS)
	$(COMP) -o Mtl \
        $(SRC)/mtl.c  $(LIB)

Mink_red:  $(SRC)/mink_red.c $(LIBS)
	$(COMP) -o Mink_red \
        $(SRC)/mink_red.c  $(LIB)

Minpol:  $(SRC)/minpol.c $(LIBS)
	$(COMP) -o Minpol \
        $(SRC)/minpol.c  $(LIB)

Mul:  $(SRC)/mul.c $(LIBS)
	$(COMP) -o Mul \
        $(SRC)/mul.c  $(LIB)

Name:  $(SRC)/name.c $(LIBS)
	$(COMP) -o Name \
        $(SRC)/name.c  $(LIB)

Normalizer:  $(SRC)/normalizer.c $(LIBS)
	$(COMP) -o Normalizer \
        $(SRC)/normalizer.c  $(LIB)

Normalizer_in_N:  $(SRC)/normalizer_in_n.c $(LIBS)
	$(COMP) -o Normalizer_in_N \
        $(SRC)/normalizer_in_n.c  $(LIB)

Normlin:  $(SRC)/normlin.c $(LIBS)
	$(COMP) -o Normlin \
        $(SRC)/normlin.c  $(LIB)

Orbit:  $(SRC)/orbit.c $(LIBS)
	$(COMP) -o Orbit \
        $(SRC)/orbit.c  $(LIB)

#Orbit_representatives:  $(SRC)/orbit_representatives.c $(LIBS)
#	$(COMP) -o Orbit_representatives \
#        $(SRC)/orbit_representatives.c  $(LIB)

#Orbrep:  $(SRC)/orbrep.c $(LIBS)
#	$(COMP) -o Orbrep \
#        $(SRC)/orbrep.c  $(LIB)

Order:  $(SRC)/order.c $(LIBS)
	$(COMP) -o Order \
        $(SRC)/order.c  $(LIB)

Pair_red:  $(SRC)/pair_red.c $(LIBS)
	$(COMP) -o Pair_red \
        $(SRC)/pair_red.c  $(LIB)

Pdet:  $(SRC)/pdet.c $(LIBS)
	$(COMP) -o Pdet \
        $(SRC)/pdet.c  $(LIB)

Perfect_neighbours:  $(SRC)/perfect_neighbours.c $(LIBS)
	$(COMP) -o Perfect_neighbours \
        $(SRC)/perfect_neighbours.c  $(LIB)

#Poincare:  $(SRC)/poincare.c $(LIBS)
#	$(COMP) -o Poincare \
#        $(SRC)/poincare.c  $(LIB)

Polyeder:  $(SRC)/polyeder.c $(LIBS)
	$(COMP) -o Polyeder \
        $(SRC)/polyeder.c  $(LIB)

P_lse_solve:  $(SRC)/p_lse_solve.c $(LIBS)
	$(COMP) -o P_lse_solve \
        $(SRC)/p_lse_solve.c  $(LIB)

Q_catalog:  $(SRC)/Q_catalog.c $(LIBS)
	$(COMP) -o Q_catalog \
        $(SRC)/Q_catalog.c  $(LIB)

QtoZ:  $(SRC)/qtoz.c $(LIBS)
	$(COMP) -o QtoZ \
        $(SRC)/qtoz.c  $(LIB)

Presentation:  $(SRC)/presentation.c $(LIBS)
	$(COMP) -o Presentation \
        $(SRC)/presentation.c $(LIB)

Red_gen:  $(SRC)/red_gen.c $(LIBS)
	$(COMP) -o Red_gen \
        $(SRC)/red_gen.c  $(LIB)

Rein:  $(SRC)/rein.c $(LIBS)
	$(COMP) -o Rein \
        $(SRC)/rein.c  $(LIB)

Rest_short:  $(SRC)/rest_short.c $(LIBS)
	$(COMP) -o Rest_short \
        $(SRC)/rest_short.c  $(LIB)

Rform:  $(SRC)/rform.c $(LIBS)
	$(COMP) -o Rform \
        $(SRC)/rform.c  $(LIB)

Reverse_name:  $(SRC)/reverse_name.c $(LIBS)
	$(COMP) -o Reverse_name \
        $(SRC)/reverse_name.c  $(LIB)

Same_generators:  $(SRC)/same_generators.c $(LIBS)
	$(COMP) -o Same_generators \
        $(SRC)/same_generators.c  $(LIB)

Scalarmul:  $(SRC)/scalarmul.c $(LIBS)
	$(COMP) -o Scalarmul \
        $(SRC)/scalarmul.c  $(LIB)

Scpr:  $(SRC)/scpr.c $(LIBS)
	$(COMP) -o Scpr \
        $(SRC)/scpr.c  $(LIB)

Short:  $(SRC)/short.c $(LIBS)
	$(COMP) -o Short \
        $(SRC)/short.c  $(LIB)

Shortest:  $(SRC)/shortest.c $(LIBS)
	$(COMP) -o Shortest \
        $(SRC)/shortest.c  $(LIB)

Short_reduce:  $(SRC)/short_reduce.c $(LIBS)
	$(COMP) -o Short_reduce \
        $(SRC)/short_reduce.c  $(LIB)

Signature:  $(SRC)/signature.c $(LIBS)
	$(COMP) -o Signature \
        $(SRC)/signature.c  $(LIB)

Simplify_mat:  $(SRC)/simplify_mat.c $(LIBS)
	$(COMP) -o Simplify_mat \
        $(SRC)/simplify_mat.c  $(LIB)

#Solve:  $(SRC)/solve.c $(LIBS)
#	$(COMP) -o Solve\
#        $(SRC)/solve.c  $(LIB)

Symbol:  $(SRC)/symbol.c $(LIBS)
	$(COMP) -o Symbol \
        $(SRC)/symbol.c  $(LIB)
	rm -f Bravais_type
	ln -f -s Symbol Bravais_type

TSubgroups:  $(SRC)/tsubgroups.c $(LIBS)
	$(COMP) -o TSubgroups \
        $(SRC)/tsubgroups.c  $(LIB)

TSupergroups:  $(SRC)/tsupergroups.c $(LIBS)
	$(COMP) -o TSupergroups \
        $(SRC)/tsupergroups.c  $(LIB)

Torsionfree:  $(SRC)/torsionfree.c $(LIBS)
	$(COMP) -o Torsionfree \
        $(SRC)/torsionfree.c  $(LIB)

Tr:  $(SRC)/tr.c $(LIBS)
	$(COMP) -o Tr \
        $(SRC)/tr.c  $(LIB)

Trace:  $(SRC)/trace.c $(LIBS)
	$(COMP) -o Trace \
        $(SRC)/trace.c  $(LIB)

Trbifo:  $(SRC)/trbifo.c $(LIBS)
	$(COMP) -o Trbifo \
        $(SRC)/trbifo.c  $(LIB)

Tr_bravais:  $(SRC)/tr_bravais.c $(LIBS)
	$(COMP) -o Tr_bravais \
        $(SRC)/tr_bravais.c  $(LIB)

Vor_vertices:  $(SRC)/vor_vertices.c $(LIBS)
	$(COMP) -o Vor_vertices \
        $(SRC)/vor_vertices.c  $(LIB)

Vectoform:  $(SRC)/vectoform.c $(LIBS)
	$(COMP) -o Vectoform \
        $(SRC)/vectoform.c  $(LIB)

ZZprog:  $(SRC)/ZZprog.c $(LIBS)
	$(COMP) -o ZZprog \
        $(SRC)/ZZprog.c  $(LIB)
	rm -f Sublattices
	ln -f -s ZZprog Sublattices

Zass_main:  $(SRC)/zass_main.c $(LIBS)
	$(COMP) -o Zass_main \
        $(SRC)/zass_main.c  $(LIB)

Z_equiv:  $(SRC)/z_equiv.c $(LIBS)
	$(COMP) -o Z_equiv \
        $(SRC)/z_equiv.c  $(LIB)

test:  $(SRC)/test.c $(LIBS)
	$(COMP) -o test \
        $(SRC)/test.c  $(LIB)

minset:  $(SRC)/minset.c $(LIBS)
	$(COMP) -o minset \
        $(SRC)/minset.c  $(LIB)

min_in:  $(SRC)/min_in.c $(LIBS)
	$(COMP) -o min_in \
        $(SRC)/min_in.c  $(LIB)

clean:
	rm -f */Add
	rm -f */Aut_grp
	rm -f */Bravais_flok
	rm -f */Bravais_catalog
	rm -f */Bravais_inclusions
	rm -f */Con
	rm -f */Conv
	rm -f */Datei
	rm -f */Dsylv
	rm -f */Elt
	rm -f */First_perfect
	rm -f */Form_elt
	rm -f */Form_space
	rm -f */Formtovec
	rm -f */Full
	rm -f */Fundamental_domain
	rm -f */Gauss
	rm -f */Graph
	rm -f */Gittstab
	rm -f */Hypisom
	rm -f */Hypstab
	rm -f */Idem
	rm -f */Inv
	rm -f */Invar_space
	rm -f */Isometry
	rm -f */Is_finite
	rm -f */Conj_bravais
	rm -f */KSubgroups
	rm -f */KSupergroups
	rm -f */Kron
	rm -f */Long_solve
	rm -f */Ltm
	rm -f */Malloc_liste
	rm -f */Mink_red
	rm -f */Modp
	rm -f */Mtl
	rm -f */Mul
	rm -f */Normalizer
	rm -f */Normalizer_in_N
	rm -f */Normlin
	rm -f */Orbit
	rm -f */Orbit_representatives
	rm -f */Orbrep
	rm -f */Order
	rm -f */P_lse_solve
	rm -f */Pair_red
	rm -f */Pdet
	rm -f */Perfect_neighbours
	rm -f */Polyeder
	rm -f */Poincare
	rm -f */Rein
	rm -f */Reverse_name
	rm -f */Rest_short
	rm -f */Rform
	rm -f */Same_generators
	rm -f */Scalarmul
	rm -f */Scpr
	rm -f */Short
	rm -f */Short_reduce
	rm -f */Shortest
	rm -f */Signature
	rm -f */Simplify_mat
	rm -f */Solve
	rm -f */Symbol ; rm -f */Bravais_type
	rm -f */TSubgroups
	rm -f */TSupergroups
	rm -f */Torsionfree
	rm -f */Tr
	rm -f */Tr_bravais
	rm -f */Trace
	rm -f */Trbifo
	rm -f */Vor_vertices
	rm -f */ZZprog ; rm -f */Sublattices
	rm -f */Zass
	rm -f */Zass_main
	rm -f */Bravais_flok
	rm -f */Extensions ; rm -f */Vector_systems
	rm -f */Extract ; rm -f */Standard_affine_form
	rm -f */Name
	rm -f */Minpol
	rm -f */Red_gen
	rm -f */Conjugated
	rm -f */Order
	rm -f */Malloc_liste
	rm -f */Presentation
	rm -f */Q_catalog
	rm -f */QtoZ
	rm -f */Z_equiv
	rm -f */Bravais_equiv
	rm -f */Bravais_grp
	rm -f */Vectoform
	rm -f */Malloc_liste
	rm -f */core

strip:
	strip *


