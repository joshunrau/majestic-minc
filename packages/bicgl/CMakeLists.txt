# CMakeFiles.txt for the BICGL
#
# Vladimir S. FONOV - vladimir.fonov@gmail.com
PROJECT(bicgl "C")

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(BICGL_PACKAGE_VERSION_MAJOR 2)
SET(BICGL_PACKAGE_VERSION_MINOR 0)
SET(BICGL_PACKAGE_VERSION_PATCH 1)
SET(BICGL_USE_GLUT ON)


IF(NOT MINC_TOOLKIT_BUILD)
  # Packaging defines
  SET(CPACK_GENERATOR TGZ)
  SET(CPACK_PACKAGE_VERSION_MAJOR ${BICGL_PACKAGE_VERSION_MAJOR})
  SET(CPACK_PACKAGE_VERSION_MINOR ${BICGL_PACKAGE_VERSION_MINOR})
  SET(CPACK_PACKAGE_VERSION_PATCH ${BICGL_PACKAGE_VERSION_PATCH})
  INCLUDE(CPack)

  FIND_PACKAGE( LIBMINC )
  FIND_PACKAGE( BICPL )
  FIND_PACKAGE( Threads )
  FIND_PACKAGE( OpenGL REQUIRED )
  FIND_PACKAGE( GLUT REQUIRED )
  FIND_PACKAGE( X11 QUIET )

  OPTION(BICGL_BUILD_SHARED_LIBS "Build shared libraries." OFF)

  INCLUDE( ${LIBMINC_USE_FILE} )
  INCLUDE( ${BICPL_USE_FILE} )
  ENABLE_TESTING()
  INCLUDE(CTest)

ENDIF(NOT MINC_TOOLKIT_BUILD)

IF(APPLE)
  SET(HAVE_APPLE_OPENGL_FRAMEWORK ON)
  SET(X_DISPLAY_MISSING ON)
ENDIF(APPLE)

SET(PACKAGE "bicgl")
SET(PACKAGE_BUGREPORT "robert.d.vincent@mcgill.ca")

SET(PACKAGE_NAME "bicgl")
SET(PACKAGE_VERSION "${BICGL_PACKAGE_VERSION_MAJOR}.${BICGL_PACKAGE_VERSION_MINOR}.${BICGL_PACKAGE_VERSION_PATCH}")
SET(PACKAGE_STRING  "${PACKAGE_NAME} ${PACKAGE_VERSION}")


#SET(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${LIBMINC_INCLUDE_DIRS} ${BICPL_INCLUDE_DIRS} )

ADD_DEFINITIONS(-DHAVE_CONFIG_H)
IF(NOT APPLE)
  ADD_DEFINITIONS(-D_XOPEN_SOURCE=500)
ENDIF(NOT APPLE)

# aliases
SET(VERSION "${PACKAGE_VERSION}")


# check for prereqs
INCLUDE(CheckIncludeFiles)

CHECK_INCLUDE_FILES(float.h     HAVE_FLOAT_H)
CHECK_INCLUDE_FILES(dlfcn.h     HAVE_FLFCN_H)
CHECK_INCLUDE_FILES(sys/stat.h  HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(values.h    HAVE_VALUES_H)
CHECK_INCLUDE_FILES(unistd.h    HAVE_UNISTD_H)
CHECK_INCLUDE_FILES(memory.h    HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(malloc.h    HAVE_MALLOC_H)
CHECK_INCLUDE_FILES(stdlib.h    HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(fcntl.h     HAVE_FCNTL_H)
CHECK_INCLUDE_FILES(vfork.h     HAVE_VFORK_H)
CHECK_INCLUDE_FILES(inttypes.h  HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(string.h    HAVE_STRING_H)
CHECK_INCLUDE_FILES(strings.h   HAVE_STRINGS_H)


IF(FREEGLUT_FOUND)
  SET(HAVE_FREEGLUT_H ON)
ELSE()
  CHECK_INCLUDE_FILES(GL/freeglut.h HAVE_FREEGLUT_H)
ENDIF()

IF(HAVE_FREEGLUT_H)
  message("Found freeglut.h")
ELSE()
  message("NOT Found freeglut.h")
ENDIF()

# CHECK_INCLUDE_FILES(bicgl/deform.h HAVE_DEFORM_H)

IF(CMAKE_USE_PTHREADS_INIT)
  SET(HAVE_PTHREAD ON)
ELSE(CMAKE_USE_PTHREADS_INIT)
  SET(HAVE_PTHREAD OFF)
ENDIF(CMAKE_USE_PTHREADS_INIT)

IF(WINDOWS)
  CHECK_INCLUDE_FILES(windows.h   HAVE_WINDOWS_H)
ELSE(WINDOWS)
  SET(HAVE_WINDOWS_H OFF)
ENDIF(WINDOWS)

# others
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/epm-header.in  ${CMAKE_CURRENT_BINARY_DIR}/epm-header)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${GLUT_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( OpenGL_graphics/Include Include )

IF(BICGL_USE_GLUT)
  INCLUDE_DIRECTORIES( GLUT_windows/Include )
ELSE()
  INCLUDE_DIRECTORIES(  X_windows/Include GLX_windows/Include )
ENDIF()

INCLUDE_DIRECTORIES( User_interface/Include Functionality/Include UI_calls_F/Include  )

IF(X11_FOUND)
  INCLUDE_DIRECTORIES( ${X11_X11_INCLUDE_PATH} )
ENDIF(X11_FOUND)

SET(PRINT_X_COLOURS print_X_colours)

SET(GLX_SRC
    GLX_windows/glx_windows.c
    GLX_windows/stored_font.c
  )

SET(GLX_headers
  GLX_windows/Include/glx_window_prototypes.h
  GLX_windows/Include/WS_windows.h
)

SET(GLUT_SRC
    GLUT_windows/glut_windows.c
#    GLUT_windows/copy_x_colours.c
  )

SET(GLUT_headers
  GLUT_windows/Include/WS_graphics.h
  GLUT_windows/Include/glut_window_prototypes.h
)

SET(X_SRC
    X_windows/x_windows.c
   )

# SET(GL_SRC
#     IrisGL_graphics/colour_def.c
#     IrisGL_graphics/draw.c
#     IrisGL_graphics/event_loop.c
#     IrisGL_graphics/lights.c
#     IrisGL_graphics/render.c
#     IrisGL_graphics/view.c
#     IrisGL_graphics/windows.c
#    )

SET(OPENGL_SRC
    OpenGL_graphics/colour_def.c
    OpenGL_graphics/draw.c
    OpenGL_graphics/event_loop.c
    OpenGL_graphics/lights.c
    OpenGL_graphics/render.c
    OpenGL_graphics/view.c
    OpenGL_graphics/windows.c
  )

SET(G_SRC
    G_graphics/draw.c
    G_graphics/draw_objects.c
    G_graphics/event_loop.c
    G_graphics/graphics_structs.c
    G_graphics/lights.c
    G_graphics/random_order.c
    G_graphics/render.c
    G_graphics/view.c
    G_graphics/windows.c
  )

IF(BICGL_USE_GLUT)
SET(bicgl_HEADERS
  ${GLUT_headers}
)
ELSE()
SET(bicgl_HEADERS
  ${GLX_headers}
)
ENDIF()

SET(bicgl_HEADERS
  ${bicgl_HEADERS}
  Include/graphics_base.h
  Include/GS_graphics.h
  Include/graphics.h
#  Include/GL/glut.h
#  Include/GL/freeglut_std.h
  Include/save.h
  Include/graphics_prototypes.h
  Include/random_order.h
  Include/graphics_struct.h
  # IrisGL_graphics/Include/GS_graphics.h
  # IrisGL_graphics/Include/def_gl.h
  # IrisGL_graphics/Include/irisgl_graphics_prototypes.h
  GLX_windows/Include/glx_window_prototypes.h
  GLX_windows/Include/WS_windows.h
  OpenGL_graphics/Include/opengl_graphics_prototypes.h
  OpenGL_graphics/Include/GS_graphics.h
  X_windows/Include/x_windows.h
  X_windows/Include/x_window_prototypes.h
 )


IF(BICGL_USE_GLUT)
  SET(OPENGL_X_SRC ${G_SRC} ${OPENGL_SRC} ${GLUT_SRC} ${bicgl_HEADERS})
ELSE()
  SET(OPENGL_X_SRC ${G_SRC} ${OPENGL_SRC} ${GLX_SRC} ${bicgl_HEADERS})
ENDIF()

IF(BICGL_BUILD_SHARED_LIBS)
  SET(LIBRARY_TYPE    SHARED)
  SET(LIBRARY_INSTALL LIBRARY)
ELSE(BICGL_BUILD_SHARED_LIBS)
  SET(LIBRARY_TYPE    STATIC)
  SET(LIBRARY_INSTALL ARCHIVE)
ENDIF(BICGL_BUILD_SHARED_LIBS)


SET(BICGL_SRC
  ${OPENGL_X_SRC}
)

ADD_LIBRARY(bicgl ${LIBRARY_TYPE}  ${BICGL_SRC})

IF(LIBMINC_STATIC_LIBRARIES AND BICPL_STATIC_LIBRARIES)

  IF(BICGL_USE_GLUT)
    TARGET_LINK_LIBRARIES(bicgl ${GLUT_LIBRARIES} ${OPENGL_gl_LIBRARY}
        ${BICPL_STATIC_LIBRARIES}  ${LIBMINC_STATIC_LIBRARIES} )
  ELSE()
    TARGET_LINK_LIBRARIES(bicgl ${OPENGL_gl_LIBRARY}
    ${BICPL_STATIC_LIBRARIES}  ${LIBMINC_STATIC_LIBRARIES} )
  ENDIF()

ELSE(LIBMINC_STATIC_LIBRARIES AND BICPL_STATIC_LIBRARIES)
  IF(BICGL_USE_GLUT)
    TARGET_LINK_LIBRARIES(bicgl ${GLUT_LIBRARIES} ${OPENGL_gl_LIBRARY}
        ${BICPL_LIBRARIES} ${LIBMINC_LIBRARIES} )
  ELSE()
    TARGET_LINK_LIBRARIES(bicgl ${OPENGL_gl_LIBRARY}
    ${BICPL_LIBRARIES} ${LIBMINC_LIBRARIES} )
  ENDIF()
ENDIF(LIBMINC_STATIC_LIBRARIES AND BICPL_STATIC_LIBRARIES)

# should we really install the library if it's static?
INSTALL(TARGETS bicgl
  ${LIBRARY_INSTALL} DESTINATION lib${LIB_SUFFIX})


INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/BICGLConfig.cmake
  DESTINATION lib${LIB_SUFFIX}
  COMPONENT Development)


INSTALL(FILES ${bicgl_HEADERS}
   DESTINATION include/bicgl)


set(BICGL_INCLUDE_DIRS_CONFIG
    ${CMAKE_CURRENT_SOURCE_DIR}/OpenGL_graphics/Include
    ${CMAKE_CURRENT_SOURCE_DIR}/Include
    ${CMAKE_CURRENT_SOURCE_DIR}/GLUT_windows/Include )

set(BICGL_LIBRARY_DIRS_CONFIG   ${CMAKE_CURRENT_BINARY_DIR} )
set(BICGL_LIBRARIES_CONFIG      bicgl )

configure_file(BICGLConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/BICGLConfig.cmake @ONLY )

# config for install dir
set(BICGL_INCLUDE_DIRS_CONFIG
    ${CMAKE_INSTALL_PREFIX}/include/bicgl/OpenGL_graphics/Include
    ${CMAKE_INSTALL_PREFIX}/include/bicgl/Include
    ${CMAKE_INSTALL_PREFIX}/include/bicgl/GLUT_windows/Include )

set(BICGL_LIBRARY_DIRS_CONFIG  ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX} )
set(BICGL_LIBRARIES_CONFIG     bicgl )

configure_file(BICGLConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/BICGLConfig.cmake @ONLY )

INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/BICGLConfig.cmake
  DESTINATION lib${LIB_SUFFIX}
  COMPONENT Development)

# VF: TODO
#IF(BUILD_TESTING)
#  ADD_SUBDIRECTORY( Testing )
#ENDIF(BUILD_TESTING)
