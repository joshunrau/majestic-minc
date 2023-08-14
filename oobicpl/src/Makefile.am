AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = \
	create_lines \
	obj_colour_to_texture \
	object_hemispheres \
	surface_area_roi \
	surface_probability_map \
	surface_volume_coordinates \
	vertstat_to_volume \
	vertstats_average \
	vertstats_colour_object \
	vertstats_extract \
	vertstats_find_peaks \
	vertstats_info \
	vertstats_math \
	vertstats_stats \
	white_cortex_validity

include_HEADERS = mniVolume.h \
	mniLabelVolume.h \
	mniBaseVolume.h \
	mniVertstatsFile.h \
	mniVertstatsMath.h \
	tree.hh 

LDADD = liboobicpl.la

lib_LTLIBRARIES = liboobicpl.la

TESTS = vertstats_tests \
	single_column_file_test \
	multiple_columns_file_test

check_PROGRAMS = $(TESTS)

liboobicpl_la_SOURCES = mniBaseVolume.cc \
	mniVolume.cc \
	mniLabelVolume.cc \
	mniVertstatsFile.cc \
	mniVertstatsMath.cc

# Programmes
create_lines_SOURCES = create_lines.cc
#mniVertstatsMath.cc
obj_colour_to_texture_SOURCES = obj_colour_to_texture.cc
object_hemispheres_SOURCES = object_hemispheres.cc
surface_area_roi_SOURCES = surface_area_roi.cc
surface_probability_map_SOURCES = surface_probability_map.cc
surface_volume_coordinates_SOURCES = surface_volume_coordinates.cc
vertstat_to_volume_SOURCES = vertstat_to_volume.cc
vertstats_average_SOURCES = vertstats_average.cc
vertstats_colour_object_SOURCES = vertstats_colour_object.cc
vertstats_extract_SOURCES = vertstatsExtract.cc
vertstats_find_peaks_SOURCES = vertstats_find_peaks.cc
vertstats_info_SOURCES = vertstatsinfo.cc
vertstats_math_SOURCES = vertstats_math.cc
vertstats_stats_SOURCES = vertstats_stats.cc 
white_cortex_validity_SOURCES = white_cortex_validity.cc

# Tests
#negative_values_test_SOURCES = negativeValuesTest.cc
#volume_tests_SOURCES = volumeTests.cc
vertstats_tests_SOURCES = vertstatsTests.cc
single_column_file_test_SOURCES = singleColumnFileTest.cc
multiple_columns_file_test_SOURCES = multipleColumnsFileTest.cc


