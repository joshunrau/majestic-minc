ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = glim_image

glim_image_SOURCES = glim.c \
	glim_image.c \
	matrix.c \
	smoothness.c

noinst_HEADERS = glim.h \
	matrix.h \
	smoothness.h

EXTRA_DIST = $(m4_files)

m4_files = m4/acx_pthread.m4 \
           m4/ax_check_gl.m4 \
           m4/ax_check_glu.m4 \
           m4/ax_check_glut.m4 \
           m4/ax_lang_compiler_ms.m4 \
           m4/lf_link_headers.m4 \
           m4/mni_REQUIRE_LIB.m4 \
           m4/mni_REQUIRE_MNILIBS.m4 \
           m4/mni_REQUIRE_OPENINVENTOR.m4 \
           m4/mni_cxx_have_koenig_lookup.m4 \
           m4/smr_CGAL_MAKEFILE.m4 \
           m4/smr_OPTIONAL_LIB.m4 \
           m4/smr_REQUIRED_LIB.m4 \
           m4/smr_WITH_BUILD_PATH.m4

