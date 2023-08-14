#-- #NOTES from http: // www.cmake.org/Wiki/CMake_Testing_With_CTest
#set(CTEST_CUSTOM_MEMCHECK_IGNORE
#    $ {CTEST_CUSTOM_MEMCHECK_IGNORE}
#    DummyExcludeMemcheckIgnoreTestSetGet
#    )

#-- #set(CTEST_CUSTOM_WARNING_MATCH
#-- #${CTEST_CUSTOM_WARNING_MATCH}
#-- #"{standard input}:[0-9][0-9]*: Warning: "
#-- #)

#
# For further details regarding this file,
# see http://www.cmake.org/Wiki/CMake_Testing_With_CTest#Customizing_CTest
#
# and
# http://www.kitware.com/blog/home/post/27
#
#----------------------------------------------------------------------

#-- #Reset maximum number of warnings so that they all show up.
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 1000)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS   1000)

# The following tests should not be run under valgrind
set(CTEST_CUSTOM_MEMCHECK_IGNORE
  @CTEST_CUSTOM_MEMCHECK_IGNORE@
  )

set(CTEST_EXTRA_COVERAGE_GLOB
  Source/.*/*.h
  Source/.*/*.hxx
  Source/.*/*.cxx
  )

set(CTEST_CUSTOM_COVERAGE_EXCLUDE
  ${CTEST_CUSTOM_COVERAGE_EXCLUDE}
  ".*Temporary/.*"
  ".*boost.*"
  # Exclude try_compile sources from coverage results:
  "/CMakeFiles/CMakeTmp/"

  # Exclude files generated by the moc pre-compiler
  ".*/moc_.*"

  # Exclude files generated by the uic pre-compiler
  ".*/ui_.*"

  # Exclude files generated by the resource pre-compiler
  ".*/qrc_.*"

  # Exclude files from the Testing directories
  ".*/Testing/.*"

  # Exclude generated python files
  ".*Python.cxx"
  ".*PythonInit.cxx"

  # Exclude Qt Designer plugins
  ".*/DesignerPlugins/.*"

  # Exclude generated cpp files
  ".*/generated_cpp/.*"
  )


set(CTEST_CUSTOM_WARNING_EXCEPTION
  ${CTEST_CUSTOM_WARNING_EXCEPTION}
  #"vtkparse.tab.c"
  "Microsoft SDKs."
  "VC.include."

  # Ignore 'detached HEAD' warnings
  ".*Note: checking out.*"
  ".*detached.*HEAD.*"
  ".*warning generated..*"

  # NIPYPE warnings:
  ".*NIPYPE-prefix.*warning.*"

  # PCRE warnings:
  ".*PCRE-prefix.*warning.*"
  ".*PCRE.*warning.*"
  ".*pcre.*has no symbols.*"
  ".*pcre.*warning*"

  ".*checking for C compiler warning.*"
  ".*gmake.*: warning.*"

  # JPEG warnings:
  ".*/JPEG.*"
  ".*warning: unused parameter .*cinfo.*"
  ".*/TIFF.*"

  # Swig warnings:
  ".*Swig-prefix.*warning.*"
  ".*Note.*SWIG.*"
  ".*CParse.*warning.*"
  ".*parser.y.*"
  ".*maximum warning verbosity.*"
  ".*Swig.*note.*"
  ".*Swig/Source.*"
  ".*warning.*argument unused during compilation.*"

  # Open-CV warnings:
  ".*OpenCV-build.*warning.*"
  ".*OpenCV-install.*warning.*"
  ".*OpenCV.*"

  # FFTW warnings:
  ".*fftw.*"

  # Numpy warnings
  "_configtest.*warning"
  # C4244: 'conversion' conversion from 'type1' to 'type2', possible loss of data
  # C4028: formal parameter 'number' different from declaration
  # C4996: 'function': was declared deprecated
  # C4133: '=' : incompatible types - from 'PyArray_Descr *' to 'PyObject *'
  # C4146: unary minus operator applied to unsigned type, result still unsigned
  # C4716: function '...' : must return a value
  # C4723: Potential divide by zero
  "numpy.(core|numarray|linalg|random).*warning (C4244|C4028|C4996|C4133|C4146|C4716|C4723)"
  # warning: assignment from incompatible pointer type
  # warning: ‘...’ defined but not used
  # warning: ‘...’ may be used uninitialized in this function
  "numpy.(core).*warning.*(assignment|defined but not used|uninitialized)"
  "NUMPY.*Warning"
  # Mac
  "numpy.core.src.multiarray.descriptor.*was declared here"

  # Tcl
  "tcl.tcl.unix.*warning: cast"
  #  warning: '...' is deprecated (declared at ...)
  "tcl.tcl.unix.*warning:.*is deprecated"
  "tcl.tcl.unix.*warning:.*ignoring return value of"

  # Tk
  "tcl.tk.unix.*warning: cast"
  "System.Library.Frameworks.Tk.framework.Headers.X11.Xlib.h.*warning: function declaration isn't a prototype"
  "tcl.tk.unix.*warning:.*ignoring return value of"

  # incrTcl
  "generic.itk_(option|archetype).*warning"
  "generic.itcl_.*warning"

  # qt suppressions from vtk...
  # Some Slicer dashboards include building bits of Qt which produce lots of
  # the following warnings when built with the MS compilers. Qt guys should
  # fix their code. Until they do, keep the Qt chatter off the Slicer dashboard
  # results:
  "include.[Qq]t([Cc]ore|[Gg]ui).*warning C4127: conditional expression is constant"
  "[Qq]t.*h.*warning.*declaration of .* shadows a member of .this"
  "[Qq]t.*h.*warning.*(copy constructor|assignment operator) could not be generated"
  # Tiger / 4.6.2 warning
  "include.[Qq]t[Cc]ore.qtextcodec.h.*warning.*is already a friend of"
  "include.QtGui.(qtextformat|qtablewidget).*warning"
  # Snowleopard / 4.6.2 warning
  "QtGui.framework.Headers.(qtextformat|qtablewidget).*warning"

  # Suppress warning caused when QT 'foreach' loops are combined
  ".*warning: declaration of '_container_' shadows a previous local"

  # STL - Tiger
  "include.c.*bits.stl_algo.h.*warning: comparison between signed and unsigned integer expressions"

  # Make
  "warning: jobserver unavailable"

  # Suppressing warnings about GL_GLEXT_LEGACY, the link reported below
  # report a similar problem with GL_GLEXT_PROTOTYPE.
  # http://lists.apple.com/archives/mac-opengl/2009/Dec/msg00081.html
  # That problem could be solved installing a newer version of X11 SDK
  # See http://xquartz.macosforge.org/trac/changeset/343
  ".*warning.*GL_GLEXT_LEGACY.*redefined"

  # ITK suppressions
  "[Uu]tilities.gdcm"
  "[Uu]tilities.vxl"
  "[Uu]tilities.itktiff"
  "([Ii]nsight|ITKv3).[Cc]ode.[Cc]ommon"
  "([Ii]nsight|ITKv3).[Cc]ode.[Nn]umerics"
  "([Ii]nsight|ITKv3).[Cc]ode.(IO|io)"
  "([Ii]nsight|ITKv3).[Cc]ode.[Ss]patial[Oo]bject"
  "([Ii]nsight|ITKv3).[Uu]tilities.[Nn]rrd(IO|io)"
  "([Ii]nsight|ITKv3).[Uu]tilities.(openjpeg|nifti)"
  "([Ii]nsight|ITKv3|BRAINSFit).*Informational: catch(...) semantics changed since Visual C\\+\\+ 7.1; structured exceptions (SEH) are no longer caught"

  # VTK suppressions
  "vtkfreetype"
  "Utilities.vtktiff"
  "VTK.*IO.vtkMySQLQuery.cxx"
  "VTK.*Utilities.vtkexodus2"
  "VTK.*Utilities.vtklibproj"
  "VTK.*Utilities.vtksqlite"
  "VTK.*Utilities.VPIC.*cxx"
  "VTK.*warn_unused_result"
  "VTK.*Filtering.*cxx"
  "VTK.*IO.*cxx"
  "VTK.*Infovis.*cxx"
  "VTK.*vtk.*warning"
  # exception specific to Mac/Carbon
  "VTK.Rendering.vtkCarbonRenderWindow.*warning.*(NewRgn|DiffRgn|EqualRgn|DisposeRgn).*is deprecated"
  # exception specific to Mac/X11
  "VTK.Rendering.vtkOpenGL.*warning: this is the location of the previous definition"

  # CTK - log4qt
  "logobjectptr.obj : warning LNK4221: no public symbols found; archive member will be inaccessible"
  "/usr/bin/ranlib: .*libLog4Qt.a.*has no symbols"
  "log4qt.rollingfileappender.h.*Warning: Property declaration maxFileSize has no READ accessor function."
  "ld.*has different visibility.*libLog4Qt.a"

  # CTK - dcmtk
  ".*dcmdata/dcovlay.h:93: warning: use of old-style cast.*"
  ".*/usr/bin/ranlib: .*lib(dcmtls|oflog).a.*has no symbols.*"
  ".*ld.*has different visibility.*(libdcmdata|libdcmnet|libdcmimgle|liboflog|libofstd).a.*"
  ".*DCMTK.(ofstd|dcmdata|dcmjpls|dcmnet|dcmimage|dcmimgle|dcmpstat|dcmqrdb).(lib|include|apps).*conversion from '(size_t|SOCKET)' to '.*', possible loss of data.*"
  ".*DCMTK.*warning.*"

  # Libs/OpenIGTLink
  "(OpenIGTLink|openigtlink).[Ss]ource.igtl"

  # Batchmake
  "BatchMake.Utilities.Zip.(zip|unzipcmd|zipcmd).*warning"

  # Libs/tclap
  "tclap.include.tclap.*Arg.h.*warning C4512"

  # Teem
  # Mac - teem/src/nrrd/superset.c:433: warning: format '%d' expects type 'int', but argument 6 has type 'ptrdiff_t'
  "teem.src.nrrd.superset.c.*warning"

  # Python - Windows
  "Modules.zlib.gzio"
  "Modules._ctypes.libffi_msvc.ffi.*warning C4018"
  "Modules.audioop.c.*warning C4018"
  "Modules._multiprocessing.*warning"
  "(Python|Objects|Modules|modules|PC).*conversion from '(Py_uintptr_t|Py_ssize_t|INT_PTR|size_t|__int64)' to '.*', possible loss of data"
  # Python - Linux
  "dist.py.*UserWarning.*licence.*distribution option is deprecated"
  "Objects.unicodeobject.c.*warning:.*differ in signedness"
  "[Ii]nclude.(string|unicodeobject).h.*note: expected .const char *. but argument is of type .unsigned char *."
  "Modules.(getpath|signalmodule).c.*warning: ignoring return value of.*declared with attribute warn_unused_result"
  "Modules.expat.xmlparse.*warning.*discards qualifiers from pointer target type"
  # Python - Mac
  "ranlib: file:.*libpython2.6.a.*has no symbols"
  "python.Modules._cursesmodule.c.*warning.*may be used uninitialized in this function"
  "python.Mac.Modules.(cf|Nav).*warning: (cast|unused)"
  "python.Modules._ctypes.*warning: function declaration isn't a prototype"
  "QuickTime.framework.QuickTime, missing required architecture x86_64 in file"
  "python.Modules._ssl.*incompatible pointer type"
  "python.Mac.Modules.carbonevt.*defined but not used"
  "python.Mac.Modules.qt._Qtmodule.*used uninitialized in this function"
  "Modules.main.c.*warning: format not a string literal and no format arguments"
  # About redefinition of symbols
  "pyconfig.h.*warning:.*redefined"
  "features.h.*"

  # curl suppressions
  "cmcurl.*warning.*conditional expression is constant"
  "cmcurl.*warning.*conversion from.*possible loss of data"
  # C4701: potentially uninitialized local variable '...' used
  # C4057: 'function' : '...' differs in indirection to slightly different base types from '...'
  # C4245: '=' : conversion from '...' to '...', signed/unsigned mismatch
  # C4706: assignment within conditional expression
  # C4232: nonstandard extension used : '...' : address of dllimport '...' is not static, identity not guaranteed
  "cmcurl.(transfer|ftp|file|cookie|url|telnet|multi|hostip4|formdata|easy).c.*warning (C4244|C4701|C4057|C4245|C4706|C4232)"
  # C4131: uses old-style declarator
  # C4244: conversion from '...' to '...', possible loss of data
  # C4127: conditional expression is constant
  "y.tab.c.*warning (C4131|C4244|C4127|C4701)"
  # C4100: unreferenced formal parameter
  "getdate.y.*warning(C4100|C4127|C4244)"
  # Mac
  "curl.mprintf.*warning.*redefined"
  "usr.include.secure._stdio.h.*warning: this is the location of the previous definition"
  "ranlib: file:.*bin.libslicerlibcurl.a.*has no symbols"

  #PythonQt
  "PythonQt.src.*conversion from 'size_t' to '.*', possible loss of data"

  #Libarchive
  "LibArchive.libArchive.*signed/unsigned mismatch"
  "LibArchive.libArchive.*conversion from 'size_t' to '.*', possible loss of data"

  # Visual studio spurious warnings...
  "The following environment variables were not found"

  # Since NUMPY has test that return build errors, let's add the following exception
  "WARNING non-zero return value in ctest from"


  ## HACK: THIS SHOULD NOT BE SUPPRESSED, NEED TO FIX OPENCV!!
  "BRAINSCutApplyModel.cxx.*warning.*increases required alignment"

  "ModuleFactory.*warning.*SymbolPointer"

  ## HACK: THIS SHOULD NOT BE SUPPRESSED, NEED TO FIX IN ANTS
  "ANTS.*warning"
  ## External Packages
  "Note.*SWIG"
  ".*parser.y.*"
  "maximum warning verbosity"
  ".*OpenCV.*"
  ".*VTK.*"
  ".*has no symbols."
  "ITKv5"
  "SlicerExecutionModel"
  "SimpleITK"

  ".*has no symbols"
  "CMake Warning:"
  "note: expanded from macro"
  ": note:"
  "vrscanl.c.* warning: unused parameter"
  "/usr/bin/libtool: warning same member name"
  )

set(CTEST_CUSTOM_WARNING_MATCH
  ${CTEST_CUSTOM_WARNING_MATCH}
  #"CMake Warning[ :]"
  )

if(APPLE)
set(CTEST_CUSTOM_WARNING_EXCEPTION
  ${CTEST_CUSTOM_WARNING_EXCEPTION}
  "warning -.: directory name .* does not exist"

  # Suppressing warnings about duplicate libraries in Darwin
  # At some point this may be addressed by CMake feature request:
  # http://public.kitware.com/Bug/view.php?id=10179
  "ld.*warning.*duplicate dylib.*"
  )
endif()

set(CTEST_CUSTOM_ERROR_MATCH
  ${CTEST_CUSTOM_ERROR_MATCH}
  "CMake Error[ :]"
  )

set(CTEST_CUSTOM_ERROR_EXCEPTION
  ${CTEST_CUSTOM_ERROR_EXCEPTION}
  # Numpy errors
  "NUMPY.*Warning"
  "NUMPY._configtest.*undefined reference"
  "_configtest.*error"
  "collect2: ld returned 1 exit status"

  #SWIG
  "Note.*SWIG"
  ".*parser.y.*"
  "maximum warning verbosity"
  ".*OpenCV.*"
  ".*VTK.*"
  ".*has no symbols."
  "ITKv5"
  "SlicerExecutionModel"
  "SimpleITK"

  )
