dnl Process this file with autoconf to produce a configure script.

AC_INIT(conglomerate, 1.6.8, bert@bic.mni.mcgill.ca)
AC_CONFIG_SRCDIR([surface_mask2.c])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

smr_WITH_BUILD_PATH

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_PROG_CC
AC_PROG_CXX

AC_CHECK_HEADERS(float.h)

AC_PROG_LIBTOOL

mni_REQUIRE_BICPL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT(epm-header
nlfit_smr
mritotal_suppress
smooth_mask
compute_icbm_vols
multispectral_stx_registration
sphere_resample_obj
suppress_fat
regional_thickness
remap_to_lobes
)
