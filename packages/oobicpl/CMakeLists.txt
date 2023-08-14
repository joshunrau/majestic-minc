# CMakeFiles.txt for the arguments library
#
# Vladimir S. FONOV - vladimir.fonov@gmail.com

CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

PROJECT(oobicpl)

SET(OOBICPL_PACKAGE_VERSION_MAJOR 1)
SET(OOBICPL_PACKAGE_VERSION_MINOR 4)
SET(OOBICPL_PACKAGE_VERSION_PATCH 60)

SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")
# ENABLE_TESTING()
# INCLUDE(CTest)

#ADD_SUBDIRECTORY( doc )

# Packaging defines
IF(NOT MINC_TOOLKIT_BUILD)
 
  SET(CPACK_GENERATOR TGZ)
  SET(CPACK_PACKAGE_VERSION_MAJOR ${OOBICPL_PACKAGE_VERSION_MAJOR})
  SET(CPACK_PACKAGE_VERSION_MINOR ${OOBICPL_PACKAGE_VERSION_MINOR})
  SET(CPACK_PACKAGE_VERSION_PATCH ${OOBICPL_PACKAGE_VERSION_PATCH})
  INCLUDE(CPack)
  
  FIND_PACKAGE( LIBMINC REQUIRED )
  FIND_PACKAGE( BICPL REQUIRED )
  FIND_PACKAGE( ARGUMENTS REQUIRED )
  
  FIND_PACKAGE( PCRE REQUIRED )
#   FIND_PACKAGE( PCRE++ REQUIRED )
  
  OPTION(OOBICPL_BUILD_SHARED_LIBS "Build shared libraries." OFF)
ENDIF(NOT MINC_TOOLKIT_BUILD)

SET(PACKAGE "oobicpl")
SET(PACKAGE_BUGREPORT "jharlap@bic.mni.mcgill.ca")

SET(PACKAGE_NAME "oobicpl")
SET(PACKAGE_VERSION "${OOBICPL_PACKAGE_VERSION_MAJOR}.${OOBICPL_PACKAGE_VERSION_MINOR}.${OOBICPL_PACKAGE_VERSION_PATCH}")
SET(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")


INCLUDE( ${LIBMINC_USE_FILE} )
INCLUDE( ${BICPL_USE_FILE} )
INCLUDE( ${ARGUMENTS_USE_FILE} )

ADD_DEFINITIONS(-DHAVE_CONFIG_H)

# add for building relocatable library
IF(UNIX)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF(UNIX)

set(CMAKE_MACOSX_RPATH ON)

# aliases
SET(VERSION "${PACKAGE_VERSION}")

# check for prereqs
INCLUDE(CheckFunctionExists)
INCLUDE(CheckIncludeFiles)

CHECK_FUNCTION_EXISTS(mkstemp  HAVE_MKSTEMP)
CHECK_FUNCTION_EXISTS(tmpnam   HAVE_TMPNAM)
CHECK_FUNCTION_EXISTS(tempnam  HAVE_TEMPNAM)
CHECK_FUNCTION_EXISTS(strerror HAVE_STRERROR)
CHECK_FUNCTION_EXISTS(popen    HAVE_POPEN)
CHECK_FUNCTION_EXISTS(fork     HAVE_FORK)
CHECK_FUNCTION_EXISTS(vfork    HAVE_VFORK)
CHECK_FUNCTION_EXISTS(fdopen   HAVE_FDOPEN)
CHECK_FUNCTION_EXISTS(strdup   HAVE_STRDUP)
CHECK_FUNCTION_EXISTS(getpwnam HAVE_GETPWNAM) 
CHECK_FUNCTION_EXISTS(select   HAVE_SELECT)
CHECK_FUNCTION_EXISTS(strerror HAVE_STRERROR) 
CHECK_FUNCTION_EXISTS(strtol   HAVE_STRTOL) 
CHECK_FUNCTION_EXISTS(strtod   HAVE_STRTOD) 
CHECK_FUNCTION_EXISTS(sysconf  HAVE_SYSCONF)
CHECK_FUNCTION_EXISTS(system   HAVE_SYSTEM)
CHECK_FUNCTION_EXISTS(vprintf  HAVE_VPRINTF)

CHECK_INCLUDE_FILES(float.h     HAVE_FLOAT_H)
CHECK_INCLUDE_FILES(limit.h     HAVE_LIMITS_H)
CHECK_INCLUDE_FILES(sys/stat.h  HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(values.h    HAVE_VALUES_H)
CHECK_INCLUDE_FILES(unistd.h    HAVE_UNISTD_H)
CHECK_INCLUDE_FILES(dirent.h    HAVE_DIRENT_H)
CHECK_INCLUDE_FILES(memory.h    HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(malloc.h    HAVE_MALLOC_H)
CHECK_INCLUDE_FILES(stdlib.h    HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(fcntl.h     HAVE_FCNTL_H)
CHECK_INCLUDE_FILES(vfork.h     HAVE_VFORK_H)
CHECK_INCLUDE_FILES(inttypes.h  HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(string.h    HAVE_STRING_H)
CHECK_INCLUDE_FILES(strings.h   HAVE_STRINGS_H)
CHECK_INCLUDE_FILES(dlfcn.h     HAVE_DLFCN_H)


# # config files for build
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY)

# others
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/epm-header.in ${CMAKE_CURRENT_BINARY_DIR}/epm-header)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} ${PCRE_INCLUDE_DIR}  src )#${PCREPP_INCLUDE_DIR}

SET(oobicpl_HEADERS
  src/mniVolume.h 
  src/mniLabelVolume.h 
  src/mniBaseVolume.h 
  src/mniVertstatsFile.h 
  src/mniVertstatsMath.h 
  src/tree.hh
)

# some variables
SET(oobicpl_LIB_SRCS
  ${oobicpl_HEADERS}
  src/mniBaseVolume.cc 
  src/mniVolume.cc 
  src/mniLabelVolume.cc 
  src/mniVertstatsFile.cc 
  src/mniVertstatsMath.cc
)

IF(OOBICPL_BUILD_SHARED_LIBS)
  SET(LIBRARY_TYPE SHARED)
  SET(LIBRARY_INSTALL LIBRARY)
ELSE(OOBICPL_BUILD_SHARED_LIBS)
  SET(LIBRARY_TYPE STATIC)
  SET(LIBRARY_INSTALL ARCHIVE)
ENDIF(OOBICPL_BUILD_SHARED_LIBS)

ADD_LIBRARY(oobicpl ${LIBRARY_TYPE} ${oobicpl_LIB_SRCS})
TARGET_LINK_LIBRARIES(oobicpl ${PCRECPP_LIBRARY} ${PCRE_LIBRARY} ${BICPL_LIBRARIES} ${VOLUME_IO_LIBRARIES} ${LIBMINC_LIBRARIES} )

LINK_LIBRARIES(oobicpl ${ARGUMENTS_LIBRARIES})

ADD_EXECUTABLE(create_lines              src/create_lines.cc)
ADD_EXECUTABLE(obj_colour_to_texture     src/obj_colour_to_texture.cc)
ADD_EXECUTABLE(object_hemispheres        src/object_hemispheres.cc)
ADD_EXECUTABLE(surface_area_roi          src/surface_area_roi.cc)
ADD_EXECUTABLE(surface_probability_map   src/surface_probability_map.cc)
ADD_EXECUTABLE(surface_volume_coordinates src/surface_volume_coordinates.cc)
ADD_EXECUTABLE(vertstat_to_volume        src/vertstat_to_volume.cc)
ADD_EXECUTABLE(vertstats_average         src/vertstats_average.cc)
ADD_EXECUTABLE(vertstats_colour_object   src/vertstats_colour_object.cc)
ADD_EXECUTABLE(vertstats_extract         src/vertstatsExtract.cc)
ADD_EXECUTABLE(vertstats_find_peaks      src/vertstats_find_peaks.cc)
ADD_EXECUTABLE(vertstats_info            src/vertstatsinfo.cc)
ADD_EXECUTABLE(vertstats_math            src/vertstats_math.cc)
ADD_EXECUTABLE(vertstats_stats           src/vertstats_stats.cc )
ADD_EXECUTABLE(white_cortex_validity     src/white_cortex_validity.cc)

IF(BUILD_TESTING)
  #ADD_EXECUTABLE(negative_values_test  negativeValuesTest.cc)
  #ADD_EXECUTABLE(volume_tests  volumeTests.cc)
  ADD_EXECUTABLE(vertstats_tests  src/vertstatsTests.cc)
  ADD_EXECUTABLE(single_column_file_test  src/singleColumnFileTest.cc)
  ADD_EXECUTABLE(multiple_columns_file_test  src/multipleColumnsFileTest.cc)
ENDIF(BUILD_TESTING)

INSTALL(TARGETS 
    create_lines 
    obj_colour_to_texture 
    object_hemispheres 
    surface_area_roi 
    surface_probability_map 
    surface_volume_coordinates 
    vertstat_to_volume 
    vertstats_average 
    vertstats_colour_object 
    vertstats_extract 
    vertstats_find_peaks 
    vertstats_info 
    vertstats_math 
    vertstats_stats 
    white_cortex_validity    
  DESTINATION bin)

INSTALL(FILES ${oobicpl_HEADERS}
   DESTINATION include)

INSTALL(TARGETS oobicpl ${LIBRARY_INSTALL} DESTINATION lib${LIB_SUFFIX})

set(OOBICPL_USE_FILE_CONFIG      "${CMAKE_CURRENT_BINARY_DIR}/UseOOBICPL.cmake")
set(OOBICPL_INCLUDE_DIRS_CONFIG   ${CMAKE_CURRENT_SOURCE_DIR}/src )
set(OOBICPL_LIBRARY_DIRS_CONFIG   ${CMAKE_CURRENT_BINARY_DIR} )
set(OOBICPL_LIBRARIES_CONFIG      oobicpl ${PCRECPP_LIBRARY} ${PCRE_LIBRARY} )#${PCREPP_LIBRARY}

configure_file(OOBICPLConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/OOBICPLConfig.cmake @ONLY )

configure_file(UseOOBICPL.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/UseOOBICPL.cmake @ONLY)

# config for install dir
set(OOBICPL_USE_FILE_CONFIG      "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}/UseOOBICPL.cmake")
set(OOBICPL_INCLUDE_DIRS_CONFIG  ${CMAKE_INSTALL_PREFIX}/include )
set(OOBICPL_LIBRARY_DIRS_CONFIG  ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX} )
set(OOBICPL_LIBRARIES_CONFIG     oobicpl ${PCRECPP_LIBRARY} ${PCRE_LIBRARY})# ${PCREPP_LIBRARY}

IF(SUPERBUILD_STAGING_PREFIX)
  STRING(REPLACE "${SUPERBUILD_STAGING_PREFIX}/" "" OOBICPL_LIBRARIES_CONFIG    "${OOBICPL_LIBRARIES_CONFIG}")
  STRING(REPLACE "${SUPERBUILD_STAGING_PREFIX}/" "" OOBICPL_INCLUDE_DIRS_CONFIG "${OOBICPL_INCLUDE_DIRS_CONFIG}")
  STRING(REPLACE "${SUPERBUILD_STAGING_PREFIX}/" "" OOBICPL_LIBRARY_DIRS_CONFIG "${OOBICPL_LIBRARY_DIRS_CONFIG}")
ENDIF(SUPERBUILD_STAGING_PREFIX)

configure_file(OOBICPLConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/OOBICPLConfig.cmake @ONLY )

configure_file(UseOOBICPL.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/UseOOBICPL.cmake @ONLY)

INSTALL(
  FILES 
   ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/UseOOBICPL.cmake 
   ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/OOBICPLConfig.cmake
  DESTINATION lib${LIB_SUFFIX}
  COMPONENT Development)


#TODO
IF(BUILD_TESTING)
#   ADD_SUBDIRECTORY( Testing )
ENDIF(BUILD_TESTING)
