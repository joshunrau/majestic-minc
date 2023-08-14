# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

AC_INIT([glim_image],[1.2],[Jonathan Taylor])
AC_CONFIG_SRCDIR([glim_image.c])

AM_INIT_AUTOMAKE

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADERS([config.h])

# This macro eases building the package; see m4/README.
smr_WITH_BUILD_PATH

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

# check for MINC
mni_REQUIRE_VOLUMEIO

# Finish up by writing output files.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT(epm-header)
