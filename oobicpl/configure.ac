dnl Process this file with autoconf to produce a configure script.

# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

AC_INIT([oobicpl], [0.4.9], [claude@bic.mni.mcgill.ca])
AC_CONFIG_SRCDIR([src/mniVolume.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([config.h])

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
# Look for perl (and make sure it's the right version)
#------------------------------------------------------------------------

AC_PATH_PROGS(PERL, [perl], no)

if test "$PERL" = no; then
    AC_MSG_WARN(No perl program found!)
    PERL=
dnl else check perl version ....
fi


#------------------------------------------------------------------------
# Check for third-party libraries.
#------------------------------------------------------------------------

smr_WITH_BUILD_PATH

mni_REQUIRE_BICPL

# the regular expression C library
smr_REQUIRED_LIB(pcre, pcre_compile, pcre.h)

# the C++ wrapper for the pcre library
AC_LANG_PUSH(C++)
AC_MSG_CHECKING([for -lpcre++])
LIBS="-lpcre++ $LIBS"
AC_TRY_LINK([#include <pcre++.h>],
            [using namespace pcrepp; Pcre EE("[a-z]", "i");],,
            [AC_MSG_ERROR(cannot find pcre++ library)])
AC_MSG_RESULT([yes])
AC_LANG_POP

# the arguments library
AC_LANG_PUSH(C++)
AC_MSG_CHECKING([for -larguments])
LIBS="-larguments $LIBS"
AC_TRY_LINK([#include <arguments.h>],
            [Arguments cArg("test", "description");],,
            [AC_MSG_ERROR(cannot find arguments library)])
AC_MSG_RESULT([yes])
AC_LANG_POP

#------------------------------------------------------------------------
# Check that we have some required header files (and macros in them)
#------------------------------------------------------------------------
AC_HEADER_STDC
AC_CHECK_HEADERS(float.h limits.h malloc.h)


#------------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#------------------------------------------------------------------------
AC_C_INLINE
AC_C_CONST
AC_TYPE_SIZE_T
#AC_CHECK_TYPE(bool,int)


#------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------
AC_FUNC_VPRINTF
AC_CHECK_FUNCS(strtod strtol)

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT(epm-header)


