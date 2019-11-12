AM_CFLAGS = -Wall -Wno-unused-result
AM_CPPFLAGS = -I$(top_srcdir)/include $(GMP_CPPFLAGS) -DDATADIR="\"$(datadir)\""
AM_LDFLAGS = $(GMP_LDFLAGS)
LIBS += $(GMP_LIBS)

# internal library
noinst_LIBRARIES = libfunctions.a

libfunctions_a_SOURCES = \
    functions/Autgrp/autgrp.c \
    functions/Autgrp/auttools.c \
    functions/Autgrp/bachtools.c \
    functions/Autgrp/iotools.c \
    functions/Autgrp/isometry.c \
    functions/Autgrp/isotools.c \
    functions/Autgrp/lattools.c \
    functions/Autgrp/mattools.c \
    functions/Autgrp/orbtools.c \
    functions/Autgrp/perfecttools.c \
    functions/Autgrp/pr_aut.c \
    functions/Autgrp/pr_isom.c \
    functions/Autgrp/preproc.c \
    functions/Autgrp/sorttools.c \
    functions/Base/base.c \
    functions/Base/base2.c \
    functions/Base/conjugate.c \
    functions/Bravais/bravais_tools.c \
    functions/Bravais/formspace.c \
    functions/Bravais/init_bravais.c \
    functions/Bravais/invar_space.c \
    functions/Bravais/konj_bravais.c \
    functions/Bravais/lincomb.c \
    functions/Bravais/normalisator.c \
    functions/Bravais/normlin.c \
    functions/Bravais/p_formspace.c \
    functions/Bravais/rform.c \
    functions/Bravais/tr_bravais.c \
    functions/Bravais/trace_bifo.c \
    functions/Contrib/suche_kand.c \
    functions/Contrib/torsionfree.c \
    functions/Datei/brav_from_datei.c \
    functions/Datei/free_bravais.c \
    functions/Datei/get_data_dir.c \
    functions/Datei/get_symbol.c \
    functions/Datei/get_zentr.c \
    functions/Datei/gittstab.c \
    functions/Datei/gittstabneu.c \
    functions/Datei/lattice_tools.c \
    functions/Datei/lattice.c \
    functions/Datei/read_symbol.c \
    functions/Datei/right_order.c \
    functions/Datei/super_lattice.c \
    functions/Getput/get_bravais.c \
    functions/Getput/get_mat.c \
    functions/Getput/put_bravais.c \
    functions/Getput/put_mat.c \
    functions/Getput/put_order.c \
    functions/Getput/read_header.c \
    functions/Graph/aff_normal.c \
    functions/Graph/calculate_representatives.c \
    functions/Graph/coho_size.c \
    functions/Graph/graph_mapped_word.c \
    functions/Graph/graph_tools.c \
    functions/Graph/is_k_subgroup.c \
    functions/Graph/lattices.c \
    functions/Graph/normalop_o.c \
    functions/Graph/obergruppenzahl.c \
    functions/Graph/orbit_fkts.c \
    functions/Graph/phi.c \
    functions/Graph/Q_data.c \
    functions/Graph/stabs.c \
    functions/Graph/sub-k-groups-fcts.c \
    functions/Graph/sub-super-tools.c \
    functions/Graph/subgroupgraph.c \
    functions/Graph/super-k-groups-fcts.c \
    functions/Hyperbolic/hyp_isom.c \
    functions/Hyperbolic/hyp_stabilizer.c \
    functions/Idem/almost_decomposable_lattice.c \
    functions/Idem/bravais_catalog.c \
    functions/Idem/centr.c \
    functions/Idem/min_pol.c \
    functions/Idem/symbol.c \
    functions/Idem/v4_catalog.c \
    functions/Idem/z_equivalent.c \
    functions/Longtools/dump_MP_mat.c \
    functions/Longtools/long_elt.c \
    functions/Longtools/long_gauss.c \
    functions/Longtools/long_hnf.c \
    functions/Longtools/long_kernel_mat.c \
    functions/Longtools/long_mat_inv.c \
    functions/Longtools/long_qbase.c \
    functions/Longtools/long_rein_mat.c \
    functions/Longtools/long_solve_mat.c \
    functions/Longtools/MP_conv_mat.c \
    functions/Longtools/MP_gauss.c \
    functions/Longtools/MP_hnf.c \
    functions/Longtools/MP_pair_red.c \
    functions/Longtools/MP_red_sort.c \
    functions/Longtools/MP_solve.c \
    functions/Matrix/add_mat.c \
    functions/Matrix/col_row_ops_mat.c \
    functions/Matrix/comp_mat.c \
    functions/Matrix/construct_mat.c \
    functions/Matrix/divide_by_gcd.c \
    functions/Matrix/elt_div_mat.c \
    functions/Matrix/find_max_entry_mat.c \
    functions/Matrix/gauss_mat.c \
    functions/Matrix/inv_mat.c \
    functions/Matrix/kgv2rat_mat.c \
    functions/Matrix/kron_mat.c \
    functions/Matrix/modp_mat.c \
    functions/Matrix/mul_mat.c \
    functions/Matrix/null_mat.c \
    functions/Matrix/p_gauss_mat.c \
    functions/Matrix/p_lse_solve.c \
    functions/Matrix/p_mat_det.c \
    functions/Matrix/p_solve_mat.c \
    functions/Matrix/real_mat.c \
    functions/Matrix/red_mat.c \
    functions/Matrix/scal_pr_mat.c \
    functions/Matrix/solve_mat.c \
    functions/Matrix/tools_mat.c \
    functions/Matrix/tr_pose_mat.c \
    functions/Matrix/trace_mat.c \
    functions/Matrix/unity_mat.c \
    functions/Name/aff_class_inf.c \
    functions/Name/compute_q_matrix.c \
    functions/Name/HM_symbol.c \
    functions/Name/point_group.c \
    functions/Name/Q_catalog.c \
    functions/Name/q_class_inf.c \
    functions/Name/reverse_name_fct.c \
    functions/Name/z_class_inf.c \
    functions/Orbit/orb_division.c \
    functions/Orbit/orbit_alg.c \
    functions/Orbit/orbit_subdivision.c \
    functions/Orbit/row_spin.c \
    functions/Polyeder/first_fuber.c \
    functions/Polyeder/first_polyeder.c \
    functions/Polyeder/fuber_tools.c \
    functions/Polyeder/polyeder_to_vecs.c \
    functions/Polyeder/polyeder_tools.c \
    functions/Polyeder/refine_fuber.c \
    functions/Polyeder/refine_polyeder.c \
    functions/Presentation/check_base.c \
    functions/Presentation/mapped_word.c \
    functions/Presentation/presentation.c \
    functions/Presentation/put_word.c \
    functions/Reduction/mink_red.c \
    functions/Reduction/pair_red.c \
    functions/Reduction/reduction_sort.c \
    functions/Reduction/short_reduce.c \
    functions/Sort/compare.c \
    functions/Sort/quicksort.c \
    functions/Sort/search.c \
    functions/Symm/dsylv.c \
    functions/Symm/rest_short.c \
    functions/Symm/short.c \
    functions/Symm/shortest.c \
    functions/Tools/chin_remainder.c \
    functions/Tools/intpow.c \
    functions/Tools/itoa.c \
    functions/Tools/malloc2dim.c \
    functions/Tools/min_div.c \
    functions/Tools/ovfl_mul.c \
    functions/Tools/prime_tools.c \
    functions/Tools/random.c \
    functions/Tools/tools.c \
    functions/TSubgroups/name_fct.c \
    functions/TSubgroups/PoaN.c \
    functions/TSubgroups/tsubgroups_fct_db.c \
    functions/TSubgroups/tsubgroups_fct.c \
    functions/TSubgroups/tsubgroups_tools.c \
    functions/TSubgroups/tSUPERgroups_fct_db.c \
    functions/Voronoi/all_vor_neighbours.c \
    functions/Voronoi/bravais_flok_datei.c \
    functions/Voronoi/bravais_flok.c \
    functions/Voronoi/calc_vor_data.c \
    functions/Voronoi/first_perfect.c \
    functions/Voronoi/init_voronoi.c \
    functions/Voronoi/normalizer.c \
    functions/Voronoi/pair_red_inv.c \
    functions/Voronoi/red_normal.c \
    functions/Voronoi/vor_neighbour.c \
    functions/Voronoi/vor_vertices.c \
    functions/Zassen/coboundary.c \
    functions/Zassen/cong_solve.c \
    functions/Zassen/convert_cocycle_to_column.c \
    functions/Zassen/convert_to_cocycle.c \
    functions/Zassen/normalop.c \
    functions/Zassen/put_cocycle.c \
    functions/Zassen/reget_gen.c \
    functions/Zassen/reverse_valuation.c \
    functions/Zassen/zass.c \
    functions/ZZ/get_groups.c \
    functions/ZZ/q2z.c \
    functions/ZZ/ZZ_cen_fun.c \
    functions/ZZ/ZZ_gen_vs.c \
    functions/ZZ/ZZ_irr_const.c \
    functions/ZZ/ZZ_lll.c \
    functions/ZZ/ZZ_zclass.c \
    functions/ZZ/ZZ.c

# unused files:
#  functions/ZZ/ZZtest.c


# binaries
bin_PROGRAMS = \
    bin/Add \
    bin/Aut_grp \
    bin/Bravais_equiv \
    bin/Bravais_grp \
    bin/Bravais_inclusions \
    bin/Conj \
    bin/Conj_bravais \
    bin/Conjugated \
    bin/Conv \
    bin/Datei \
    bin/Elt \
    bin/Extensions \
    bin/Extract \
    bin/First_perfect \
    bin/Form_elt \
    bin/Form_space \
    bin/Formtovec \
    bin/Full \
    bin/Gauss \
    bin/Graph \
    bin/Idem \
    bin/Inv \
    bin/Invar_space \
    bin/Is_finite \
    bin/Isometry \
    bin/Kron \
    bin/KSubgroups \
    bin/KSupergroups \
    bin/Long_solve \
    bin/Ltm \
    bin/Mink_red \
    bin/Minpol \
    bin/Modp \
    bin/Mtl \
    bin/Mul \
    bin/Name \
    bin/Normalizer \
    bin/Normalizer_in_N \
    bin/Normlin \
    bin/Orbit \
    bin/Order \
    bin/P_lse_solve \
    bin/Pair_red \
    bin/Pdet \
    bin/Perfect_neighbours \
    bin/Polyeder \
    bin/Presentation \
    bin/Q_catalog \
    bin/QtoZ \
    bin/Red_gen \
    bin/Rein \
    bin/Rest_short \
    bin/Reverse_name \
    bin/Rform \
    bin/Same_generators \
    bin/Scalarmul \
    bin/Scpr \
    bin/Short \
    bin/Short_reduce \
    bin/Shortest \
    bin/Signature \
    bin/Simplify_mat \
    bin/Symbol \
    bin/Torsionfree \
    bin/Tr \
    bin/Tr_bravais \
    bin/Trace \
    bin/Trbifo \
    bin/TSubgroups \
    bin/TSupergroups \
    bin/Vectoform \
    bin/Vor_vertices \
    bin/Z_equiv \
    bin/Zass_main \
    bin/ZZprog


# currently unused:
#Fundamental_domain_SOURCES = src/fundamental_domain.c
#Gittstab_SOURCES = src/gittstab.c
#Hypisom_SOURCES = src/hypisom.c
#Hypstab_SOURCES = src/hypstab.c
#Orbit_representatives_SOURCES = src/orbit_representatives.c
#Orbrep_SOURCES = src/orbrep.c
#Poincare_SOURCES = src/poincare.c
#Solve_SOURCES = src/solve.c

bin_Add_SOURCES = src/add.c
bin_Add_LDADD = libfunctions.a

bin_Aut_grp_SOURCES = src/autgrp.c
bin_Aut_grp_LDADD = libfunctions.a

bin_Bravais_equiv_SOURCES = src/bravais_equiv.c
bin_Bravais_equiv_LDADD = libfunctions.a

bin_Bravais_grp_SOURCES = src/bravaisgroup.c
bin_Bravais_grp_LDADD = libfunctions.a

bin_Bravais_inclusions_SOURCES = src/bravais_inclusions.c
bin_Bravais_inclusions_LDADD = libfunctions.a

bin_Conj_bravais_SOURCES = src/conj_bravais.c
bin_Conj_bravais_LDADD = libfunctions.a

bin_Conj_SOURCES = src/conj.c
bin_Conj_LDADD = libfunctions.a

bin_Conjugated_SOURCES = src/conjugated.c
bin_Conjugated_LDADD = libfunctions.a

bin_Conv_SOURCES = src/conv.c
bin_Conv_LDADD = libfunctions.a

bin_Datei_SOURCES = src/datei.c
bin_Datei_LDADD = libfunctions.a

bin_Elt_SOURCES = src/elt.c
bin_Elt_LDADD = libfunctions.a

bin_Extensions_SOURCES = src/extensions.c
bin_Extensions_LDADD = libfunctions.a

bin_Extract_SOURCES = src/extract.c
bin_Extract_LDADD = libfunctions.a

bin_First_perfect_SOURCES = src/first_perfect.c
bin_First_perfect_LDADD = libfunctions.a

bin_Form_elt_SOURCES = src/form_elt.c
bin_Form_elt_LDADD = libfunctions.a

bin_Form_space_SOURCES = src/form_space.c
bin_Form_space_LDADD = libfunctions.a

bin_Formtovec_SOURCES = src/formtovec.c
bin_Formtovec_LDADD = libfunctions.a

bin_Full_SOURCES = src/full.c
bin_Full_LDADD = libfunctions.a

bin_Gauss_SOURCES = src/gauss.c
bin_Gauss_LDADD = libfunctions.a

bin_Graph_SOURCES = src/graph.c
bin_Graph_LDADD = libfunctions.a

bin_Idem_SOURCES = src/idem.c
bin_Idem_LDADD = libfunctions.a

bin_Inv_SOURCES = src/inv.c
bin_Inv_LDADD = libfunctions.a

bin_Invar_space_SOURCES = src/invar_space.c
bin_Invar_space_LDADD = libfunctions.a

bin_Is_finite_SOURCES = src/is_finite.c
bin_Is_finite_LDADD = libfunctions.a

bin_Isometry_SOURCES = src/isometry.c
bin_Isometry_LDADD = libfunctions.a

bin_Kron_SOURCES = src/kron.c
bin_Kron_LDADD = libfunctions.a

bin_KSubgroups_SOURCES = src/sub-k-groups.c
bin_KSubgroups_LDADD = libfunctions.a

bin_KSupergroups_SOURCES = src/super-k-groups.c
bin_KSupergroups_LDADD = libfunctions.a

bin_Long_solve_SOURCES = src/long_solve.c
bin_Long_solve_LDADD = libfunctions.a

bin_Ltm_SOURCES = src/ltm.c
bin_Ltm_LDADD = libfunctions.a

bin_Mink_red_SOURCES = src/mink_red.c
bin_Mink_red_LDADD = libfunctions.a

bin_Minpol_SOURCES = src/minpol.c
bin_Minpol_LDADD = libfunctions.a

bin_Modp_SOURCES = src/modp.c
bin_Modp_LDADD = libfunctions.a

bin_Mtl_SOURCES = src/mtl.c
bin_Mtl_LDADD = libfunctions.a

bin_Mul_SOURCES = src/mul.c
bin_Mul_LDADD = libfunctions.a

bin_Name_SOURCES = src/name.c
bin_Name_LDADD = libfunctions.a

bin_Normalizer_in_N_SOURCES = src/normalizer_in_n.c
bin_Normalizer_in_N_LDADD = libfunctions.a

bin_Normalizer_SOURCES = src/normalizer.c
bin_Normalizer_LDADD = libfunctions.a

bin_Normlin_SOURCES = src/normlin.c
bin_Normlin_LDADD = libfunctions.a

bin_Orbit_SOURCES = src/orbit.c
bin_Orbit_LDADD = libfunctions.a

bin_Order_SOURCES = src/order.c
bin_Order_LDADD = libfunctions.a

bin_P_lse_solve_SOURCES = src/p_lse_solve.c
bin_P_lse_solve_LDADD = libfunctions.a

bin_Pair_red_SOURCES = src/pair_red.c
bin_Pair_red_LDADD = libfunctions.a

bin_Pdet_SOURCES = src/pdet.c
bin_Pdet_LDADD = libfunctions.a

bin_Perfect_neighbours_SOURCES = src/perfect_neighbours.c
bin_Perfect_neighbours_LDADD = libfunctions.a

bin_Polyeder_SOURCES = src/polyeder.c
bin_Polyeder_LDADD = libfunctions.a

bin_Presentation_SOURCES = src/presentation.c
bin_Presentation_LDADD = libfunctions.a

bin_Q_catalog_SOURCES = src/Q_catalog.c
bin_Q_catalog_LDADD = libfunctions.a

bin_QtoZ_SOURCES = src/qtoz.c
bin_QtoZ_LDADD = libfunctions.a

bin_Red_gen_SOURCES = src/red_gen.c
bin_Red_gen_LDADD = libfunctions.a

bin_Rein_SOURCES = src/rein.c
bin_Rein_LDADD = libfunctions.a

bin_Rest_short_SOURCES = src/rest_short.c
bin_Rest_short_LDADD = libfunctions.a

bin_Reverse_name_SOURCES = src/reverse_name.c
bin_Reverse_name_LDADD = libfunctions.a

bin_Rform_SOURCES = src/rform.c
bin_Rform_LDADD = libfunctions.a

bin_Same_generators_SOURCES = src/same_generators.c
bin_Same_generators_LDADD = libfunctions.a

bin_Scalarmul_SOURCES = src/scalarmul.c
bin_Scalarmul_LDADD = libfunctions.a

bin_Scpr_SOURCES = src/scpr.c
bin_Scpr_LDADD = libfunctions.a

bin_Short_reduce_SOURCES = src/short_reduce.c
bin_Short_reduce_LDADD = libfunctions.a

bin_Short_SOURCES = src/short.c
bin_Short_LDADD = libfunctions.a

bin_Shortest_SOURCES = src/shortest.c
bin_Shortest_LDADD = libfunctions.a

bin_Signature_SOURCES = src/signature.c
bin_Signature_LDADD = libfunctions.a

bin_Simplify_mat_SOURCES = src/simplify_mat.c
bin_Simplify_mat_LDADD = libfunctions.a

bin_Symbol_SOURCES = src/symbol.c
bin_Symbol_LDADD = libfunctions.a

bin_Torsionfree_SOURCES = src/torsionfree.c
bin_Torsionfree_LDADD = libfunctions.a

bin_Tr_bravais_SOURCES = src/tr_bravais.c
bin_Tr_bravais_LDADD = libfunctions.a

bin_Tr_SOURCES = src/tr.c
bin_Tr_LDADD = libfunctions.a

bin_Trace_SOURCES = src/trace.c
bin_Trace_LDADD = libfunctions.a

bin_Trbifo_SOURCES = src/trbifo.c
bin_Trbifo_LDADD = libfunctions.a

bin_TSubgroups_SOURCES = src/tsubgroups.c
bin_TSubgroups_LDADD = libfunctions.a

bin_TSupergroups_SOURCES = src/tsupergroups.c
bin_TSupergroups_LDADD = libfunctions.a

bin_Vectoform_SOURCES = src/vectoform.c
bin_Vectoform_LDADD = libfunctions.a

bin_Vor_vertices_SOURCES = src/vor_vertices.c
bin_Vor_vertices_LDADD = libfunctions.a

bin_Z_equiv_SOURCES = src/z_equiv.c
bin_Z_equiv_LDADD = libfunctions.a

bin_Zass_main_SOURCES = src/zass_main.c
bin_Zass_main_LDADD = libfunctions.a

bin_ZZprog_SOURCES = src/ZZprog.c
bin_ZZprog_LDADD = libfunctions.a


# FIXME: the following is a bit of a hack. Not sure how to do this
# better, though.
all-local:
	@( cd $(top_srcdir)/tables; if [ ! -d qcatalog ] ; then tar xf qcatalog.tar.gz ; fi )
	@( cd $(top_srcdir)/tables/symbol && make )
	@mkdir -p bin
	@ln -f -s Datei bin/Bravais_catalog
	@ln -f -s Extract bin/Standard_affine_form
	@ln -f -s Extensions bin/Vector_systems
	@ln -f -s Symbol bin/Bravais_type
	@ln -f -s ZZprog bin/Sublattices
