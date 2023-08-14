PROJECT(ILT)


# install Getopt::Tabular
INSTALL(FILES
    ILT/Executables.pm
    ILT/ImageInfo.pm
    ILT/ImageLayout.pm
    ILT/LayoutInclude.pm
    ILT/LayoutUtils.pm
    ILT/ProgUtils.pm
    ILT/View.pm
    ILT/SceneObject.pm
  DESTINATION
   perl/ILT)

INSTALL(FILES
   ILT/SceneObject/ColourObject.pm
   ILT/SceneObject/GeometricObject.pm
   ILT/SceneObject/IntersectionObject.pm
   ILT/SceneObject/OneSubObject.pm
   ILT/SceneObject/PlaneObject.pm
   ILT/SceneObject/RenderObject.pm
   ILT/SceneObject/TextObject.pm
   ILT/SceneObject/TransformObject.pm
   ILT/SceneObject/UnionObject.pm
   ILT/SceneObject/VolumeObject.pm
 DESTINATION
  perl/ILT/SceneObject)

INSTALL(PROGRAMS
      script/create_verify_image
    PERMISSIONS
      OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ  WORLD_EXECUTE WORLD_READ
    DESTINATION
      bin )

install(DIRECTORY
  share/ILT
  USE_SOURCE_PERMISSIONS
  DESTINATION  share)
