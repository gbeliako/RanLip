#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([libranlip],[1.0],[gleb@deakin.edu.au esteban@v7w.com])
AM_INIT_AUTOMAKE([libranlip],[1.0],[])

AC_CONFIG_SRCDIR([src/ranlipdist.h])
AM_CONFIG_HEADER([config.h])

#get libtool going
AM_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_LANG_CPLUSPLUS

# Checks for libraries.
AC_CHECK_LIB([m], [cos])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
