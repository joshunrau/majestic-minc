# FindNetCDF.cmake module


FIND_PATH(NETCDF_INCLUDE_DIR netcdf.h /usr/include /usr/local/include /usr/local/bic/include)

FIND_LIBRARY(NETCDF_LIBRARY NAMES netcdf PATHS /usr/lib /usr/local/lib /usr/local/bic/lib)


IF (NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)
   SET(NETCDF_FOUND TRUE)
ENDIF (NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)


IF (NETCDF_FOUND)
   IF (NOT NETCDF_FIND_QUIETLY)
      MESSAGE(STATUS "Found NetCDF headers: ${NETCDF_INCLUDE_DIR}")
      MESSAGE(STATUS "Found NetCDF library: ${NETCDF_LIBRARY}")
   ENDIF (NOT NETCDF_FIND_QUIETLY)
ELSE (NETCDF_FOUND)
   IF (NETCDF_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Cound not find NetCDF")
   ENDIF (NETCDF_FIND_REQUIRED)
ENDIF (NETCDF_FOUND)
