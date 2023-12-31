# CMakeFiles.txt for the postf
#
# Vladimir S. FONOV - vladimir.fonov@gmail.com

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(POSTF_PACKAGE_VERSION_MAJOR 1)
SET(POSTF_PACKAGE_VERSION_MINOR 0)
SET(POSTF_PACKAGE_VERSION_PATCH 03)

IF(NOT MINC_TOOLKIT_BUILD)

  PROJECT(postf )

  # Packaging defines
  SET(CPACK_GENERATOR TGZ)
  SET(CPACK_PACKAGE_VERSION_MAJOR ${POSTF_PACKAGE_VERSION_MAJOR})
  SET(CPACK_PACKAGE_VERSION_MINOR ${POSTF_PACKAGE_VERSION_MINOR})
  SET(CPACK_PACKAGE_VERSION_PATCH ${POSTF_PACKAGE_VERSION_PATCH})
  INCLUDE(CPack)

  FIND_PACKAGE( LIBMINC REQUIRED )
  FIND_PACKAGE( X11 REQUIRED )
ENDIF(NOT MINC_TOOLKIT_BUILD)

SET(PACKAGE "postf")
SET(PACKAGE_BUGREPORT "a.janke@gmail.com")

SET(PACKAGE_NAME "postf")
SET(PACKAGE_VERSION "${POSTF_PACKAGE_VERSION_MAJOR}.${POSTF_PACKAGE_VERSION_MINOR}.${POSTF_PACKAGE_VERSION_PATCH}")
SET(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")


INCLUDE( ${LIBMINC_USE_FILE} )

SET(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${LIBMINC_INCLUDE_DIRS} )

ADD_DEFINITIONS(-DHAVE_CONFIG_H)

# aliases
SET(VERSION "${PACKAGE_VERSION}")

# others
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/epm-header.in  ${CMAKE_CURRENT_BINARY_DIR}/epm-header)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${X11_INCLUDE_DIR})

ADD_EXECUTABLE(postf postf.c mnifile.c conversions.c )
TARGET_LINK_LIBRARIES(postf ${LIBMINC_LIBRARIES} ${X11_LIBRARIES})

INSTALL(TARGETS postf DESTINATION bin)
