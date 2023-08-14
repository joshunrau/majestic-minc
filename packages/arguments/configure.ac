dnl Process this file with autoconf to produce a configure script.

AC_INIT([arguments],[0.2.2],[Jonathan Harlap <jharlap@bic.mni.mcgill.ca>])
AC_CONFIG_SRCDIR([src/arguments.cpp])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE

#------------------------------------------------------------------------
# Checks for compiler toolset.
#------------------------------------------------------------------------
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET

AC_DISABLE_SHARED
AC_PROG_LIBTOOL


#------------------------------------------------------------------------
# Check that we have some required header files (and macros in them)
#------------------------------------------------------------------------
AC_HEADER_STDC
#AC_CHECK_HEADERS(float.h limits.h malloc.h)


#------------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#------------------------------------------------------------------------
#AC_C_INLINE
#AC_C_CONST
#AC_TYPE_SIZE_T
#AC_CHECK_TYPE(bool,int)


#------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------
#AC_FUNC_VPRINTF
#AC_CHECK_FUNCS(strtod strtol)


AC_CHECK_PROGS(HAVE_DOT, [dot], NO)
if test "$HAVE_DOT" = dot; then
    HAVE_DOT=YES
fi

AC_OUTPUT([Makefile src/Makefile doc/Makefile doc/doxygen_config epm-header])

