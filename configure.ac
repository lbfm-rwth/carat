dnl
dnl Process this file with autoconf to produce a configure script.
dnl 

dnl
dnl Setup autoconf
dnl
AC_PREREQ([2.68])
AC_INIT([carat], [TODO_VERSION], [https://github.com/lbfm-rwth/carat/issues], [carat], [https://lbfm-rwth.github.io/carat/])
AC_CONFIG_SRCDIR([src/autgrp.c])
AC_CONFIG_HEADER([include/config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([cnf])

dnl
dnl Get canonical host info
dnl
AC_CANONICAL_HOST

dnl
dnl C is the language
dnl
AC_PROG_CC
AC_LANG([C])

dnl
dnl Setup automake
dnl
AC_PROG_RANLIB
AM_INIT_AUTOMAKE([1.11 -Wall foreign subdir-objects])
AM_SILENT_RULES([yes])
AM_PROG_AR

dnl
dnl Detect GMP
dnl
AC_ARG_WITH([gmp],
  [AS_HELP_STRING([--with-gmp@<:@=builtin|PREFIX@:>@],
    [prefix of GMP installation. e.g. /usr/local])],
    [],[with_gmp=yes])

GMP_CPPFLAGS=
GMP_LDFLAGS=
GMP_LIBS="-lgmp"
AS_CASE([$with_gmp],
  [no],[
    AC_MSG_ERROR([Building without GMP is not supported])
  ],
  [yes],[],
  [*],[
    GMP_CPPFLAGS="-I${with_gmp}/include"
    GMP_LDFLAGS="-L${with_gmp}/lib"
  ]
)

AX_CHECK_LIBRARY([GMP], [gmp.h], [gmp], [__gmpz_init],
    [],
    [AC_MSG_ERROR([GMP not found at prefix $with_gmp])])

AC_SUBST([GMP_CPPFLAGS])
AC_SUBST([GMP_LDFLAGS])
AC_SUBST([GMP_LIBS])

dnl
dnl Other checks
dnl

dnl check whether libm / -lm is available and necessary
AC_SEARCH_LIBS([floor], [m], [], [
  AC_MSG_ERROR([unable to find the floor() function])
])


dnl
dnl Define TOPDIR
dnl
dnl TODO: get rid of this? hardcoding paths in binaries seems problematic
dnl TODO: at the very least, when installing this should be recompiled to
dnl       a value relative to the $prefix
dnl
AC_DEFINE_UNQUOTED([TOPDIR], ["$ac_pwd"], [for functions/Datei/get_data_dir.c])

dnl
dnl Finally, generate the Makefiles and output everything
dnl
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
