project(mni-perllib)
CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

IF(NOT MINC_TOOLKIT_BUILD)
  FIND_PACKAGE( Perl REQUIRED )
ENDIF(NOT MINC_TOOLKIT_BUILD)


# quick and dirty hack to install mni-perllib using CMake
# TODO: make it proper
SET(DEFAULT_DATA_DIR "${CMAKE_INSTALL_PREFIX}/share/mni-models")

CONFIGURE_FILE(MNI/DataDir.pm.in ${CMAKE_CURRENT_BINARY_DIR}/DataDir.pm @ONLY)

INSTALL(FILES 
   MNI.pm
   mni_perllib.pod
  DESTINATION 
   perl)

SET( MNI_PERLLIB_MODULES 
#    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/Compress.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/FileUtilities.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/MincUtilities.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/MiscUtilities.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/NumericUtilities.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/PathUtilities.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/Spawn.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/Startup.pm
#    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/Tag.pm
#    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/TagPoint.pm
#    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/TagSet.pm
#    ${CMAKE_CURRENT_BINARY_DIR}/DataDir.pm
)

INSTALL(FILES 
    ${MNI_PERLLIB_MODULES}
    MNI/Spawn.pod
    ${CMAKE_CURRENT_BINARY_DIR}/DataDir.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/Tag.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/TagPoint.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/TagSet.pm
    ${CMAKE_CURRENT_SOURCE_DIR}/MNI/Compress.pm
  DESTINATION 
    perl/MNI )


#execute_process(COMMAND 
#  ${PERL_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/autosplit.pl" "${CMAKE_CURRENT_SOURCE_DIR}/MNI.pm" "${CMAKE_CURRENT_BINARY_DIR}/perl/auto" )

# now create autoload files
foreach(pm  ${MNI_PERLLIB_MODULES} )
  get_filename_component(mod_name ${pm} NAME)
  add_custom_target("auto_${mod_name}" ALL
    COMMAND  ${PERL_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/autosplit.pl" "${pm}" "${CMAKE_CURRENT_BINARY_DIR}/perl/auto" )
endforeach(pm)

# install generated auto
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/perl DESTINATION .)
