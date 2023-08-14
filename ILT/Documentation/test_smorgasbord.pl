#!/usr/local/bin/perl5 -w

#-------------------------------------------------------
#
#  This script configures an image of M volumes by N slices
#  using the ImageLayout object and related objects.  As such,
#  it is simply a demonstration of the interface to the ImageLayout object
#
#  Author: David MacDonald        March, 1998
#
#--------------------------------------------------------

    use ILT::LayoutInclude;
    use ILT::LayoutUtils;
    use ILT::ProgUtils;

#-------------------------------------------------------
#
#  Define the information to be displayed
#
#  This could be done in any number of ways, so the simplest is shown here
#  for demonstration purposes
#
#--------------------------------------------------------

    $n_rows = 4;
    $n_cols = 4;
    $layout = ILT::ImageLayout->new_grid( $n_rows, $n_cols );

    $layout->horizontal_white_space( 3 );
    $layout->vertical_white_space( 3 );
    $layout->white_space_colour( "dim_gray" );

    for( $row = 0;  $row < $n_rows;  ++$row )
    {
        for( $col = 0;  $col < $n_cols;  ++$col )
        {
            $view = ILT::View->new_canonical( Top_view );
            $surface_object = ILT::GeometricObject->new( "Data/surf1.obj" );
            $image_info = ILT::ImageInfo->new( $surface_object, $view );
            $layout->image_info( $layout->row_col_to_index($row,$col),
                                 $image_info );
        }
    }

#--------------------------------------------------------------------------
#
#   image[0][0] is an x-slice through volume without any colour coding
#
#--------------------------------------------------------------------------

    $plane_object = ILT::PlaneObject->new_canonical( Sagittal_axis, -10 );
    $view = ILT::View->new_canonical( Left_view );
    $volume_object = ILT::VolumeObject->new( "Data/volume1.mnc" );
    $scene_object = ILT::IntersectionObject->new( $plane_object,
                                                  $volume_object );
    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(0,0), $image_info );

#--------------------------------------------------------------------------
#
#   image[0][1] is an x-slice through a volume with colour coding
#
#--------------------------------------------------------------------------

    $plane_object = ILT::PlaneObject->new_canonical( Sagittal_axis, -10 );
    $view = ILT::View->new_canonical( Left_view );
    $volume_object = ILT::VolumeObject->new( "Data/volume1.mnc" );

    $scene_object = 
                  ILT::ColourObject->new_volume_colouring(
                       ILT::IntersectionObject->new( $plane_object,
                                                     $volume_object ),
                       $volume_object, Hot_metal_scale, 50, 150 );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(0,1), $image_info );

#--------------------------------------------------------------------------
#
#   image[0][2] is three colour-coded slices with a three-d view
#
#--------------------------------------------------------------------------

    $plane_object1 = ILT::PlaneObject->new_canonical( Sagittal_axis, 0 );
    $plane_object2 = ILT::PlaneObject->new_canonical( Coronal_axis, 0 );
    $plane_object3 = ILT::PlaneObject->new_canonical( Transverse_axis, 0 );
    $view = ILT::View->new_arbitrary( -1, -2, -1, 0, 0, 1 );
    $volume_object = ILT::VolumeObject->new( "Data/volume1.mnc" );

    $coloured_slice1 = ILT::ColourObject->new_volume_colouring(
                           ILT::IntersectionObject->new( $plane_object1,
                                                         $volume_object ),
                           $volume_object, Hot_metal_scale, 50, 150 );
    $coloured_slice2 = ILT::ColourObject->new_volume_colouring(
                           ILT::IntersectionObject->new( $plane_object2,
                                                         $volume_object ),
                           $volume_object, Hot_metal_scale, 50, 150 );
    $coloured_slice3 = ILT::ColourObject->new_volume_colouring(
                           ILT::IntersectionObject->new( $plane_object3,
                                                         $volume_object ),
                           $volume_object, Hot_metal_scale, 50, 150 );

    $union = ILT::UnionObject->new( $coloured_slice1,
                               $coloured_slice2,
                               $coloured_slice3 );

    $scene_object = ILT::RenderObject->new( $union );
    $scene_object->lighting_state( True );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(0,2), $image_info );

#--------------------------------------------------------------------------
#
#   image[0][3] is an x-slice through a PET volume overlaid on average MRI
#
#--------------------------------------------------------------------------

    $plane_object = ILT::PlaneObject->new_canonical( Sagittal_axis, -10 );
    $view = ILT::View->new_canonical( Left_view );
    $mri_volume = ILT::VolumeObject->new( "Data/volume1.mnc" );
    $pet_volume = ILT::VolumeObject->new(
           "/avgbrain/brain/images/mni_demo_fdg_normal_pet_tal.mnc.gz" );
    $slice_object = ILT::IntersectionObject->new( $plane_object,
                                             $mri_volume ),


    $coloured_mri = ILT::ColourObject->new_volume_colouring( $slice_object,
                           $mri_volume, Gray_scale, 50, 150 );

    $scene_object = ILT::ColourObject->new_volume_colouring( $coloured_mri,
                                      $pet_volume, Spectral_scale,
                                      400, 500 );
    $scene_object->under_colour( "transparent" );
    $scene_object->opacity( .5 );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(0,3), $image_info );

#--------------------------------------------------------------------------
#
#   image[1][0] is a right view of a lit surface
#
#--------------------------------------------------------------------------

    $surface_object = ILT::GeometricObject->new( "Data/surf2.obj" );
    $view = ILT::View->new_canonical( Right_view );

    $scene_object = ILT::RenderObject->new( $surface_object );
    $scene_object->lighting_state( True );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(1,0), $image_info );

#--------------------------------------------------------------------------
#
#   image[1][1] is a cross section of a surface
#
#--------------------------------------------------------------------------

    $plane_object = ILT::PlaneObject->new_canonical( Sagittal_axis, -40 );
    $surface_object = ILT::GeometricObject->new( "Data/surf2.obj" );
    $view = ILT::View->new_canonical( Right_view );

    $cross = ILT::IntersectionObject->new( $plane_object, $surface_object ),
    $scene_object = $cross;

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(1,1), $image_info );

#--------------------------------------------------------------------------
#
#   image[1][2] is a cross section of a surface superimposed on a colour_coded
#               slice
#
#--------------------------------------------------------------------------

    $plane_object = ILT::PlaneObject->new_canonical( Sagittal_axis, -40 );
    $surface_object = ILT::GeometricObject->new( "Data/surf2.obj" );
    $volume_object = ILT::VolumeObject->new( "Data/volume1.mnc" );
    $view = ILT::View->new_canonical( Right_view );

    $cross = ILT::IntersectionObject->new( $plane_object, $surface_object );

    $coloured_slice = ILT::ColourObject->new_volume_colouring(
                           ILT::IntersectionObject->new( $plane_object,
                                                    $volume_object ),
                           $volume_object, Gray_scale, 50, 150 );

    $scene_object = ILT::UnionObject->new( $cross, $coloured_slice );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(1,2), $image_info );

#--------------------------------------------------------------------------
#
#   image[1][3] is a surface colour coded with a volume
#
#--------------------------------------------------------------------------

    $surface_object = ILT::GeometricObject->new( "Data/surf2.obj" );
    $volume_object = ILT::VolumeObject->new(
           "/avgbrain/brain/images/mni_demo_fdg_normal_pet_tal.mnc.gz" );
    $view = ILT::View->new_canonical( Left_view );

    $coloured_surface = ILT::ColourObject->new_volume_colouring(
                             $surface_object,
                             $volume_object, Spectral_scale, 300, 350 );
    $coloured_surface->under_colour( "white" );

    $scene_object = ILT::RenderObject->new( $coloured_surface );
    $scene_object->lighting_state( True );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(1,3), $image_info );

#--------------------------------------------------------------------------
#
#   image[2][0] is a set of fat cross sections of a surface
#
#--------------------------------------------------------------------------

    $plane_object1 = ILT::PlaneObject->new_canonical( Sagittal_axis, -50 );
    $plane_object2 = ILT::PlaneObject->new_canonical( Sagittal_axis, -25 );
    $plane_object3 = ILT::PlaneObject->new_canonical( Sagittal_axis, 0 );
    $plane_object4 = ILT::PlaneObject->new_canonical( Coronal_axis, 0 );
    $plane_object5 = ILT::PlaneObject->new_canonical( Transverse_axis, 0 );
    $surface_object = ILT::GeometricObject->new( "Data/surf2.obj" );
    $view = ILT::View->new_arbitrary( 1, .5, .3, 0, 0, 1 );

    $cross1 = ILT::IntersectionObject->new( $plane_object1, $surface_object ),
    $cross2 = ILT::IntersectionObject->new( $plane_object2, $surface_object ),
    $cross3 = ILT::IntersectionObject->new( $plane_object3, $surface_object ),
    $cross4 = ILT::IntersectionObject->new( $plane_object4, $surface_object ),
    $cross5 = ILT::IntersectionObject->new( $plane_object5, $surface_object ),

    $union = ILT::UnionObject->new( $cross1, $cross2, $cross3, $cross4,
                                    $cross5 );
    $render_object = ILT::RenderObject->new( $union );
    $render_object->lighting_state( True );
    $render_object->line_width( 2 );

    $image_info = ILT::ImageInfo->new( $render_object, $view );
    $layout->image_info( $layout->row_col_to_index(2,0), $image_info );

#--------------------------------------------------------------------------
#
#   image[2][1] is several cross sections of a surface superimposed on a
#   colour_coded slice, where a user defined colour map is used to create
#   a contour-type image
#
#--------------------------------------------------------------------------

    $plane_object1 = ILT::PlaneObject->new_canonical( Transverse_axis, 10 );
    $plane_object2 = ILT::PlaneObject->new_canonical( Transverse_axis, 30 );
    $plane_object3 = ILT::PlaneObject->new_canonical( Transverse_axis, 50 );
    $surface_object = ILT::GeometricObject->new( "Data/surf2.obj" );
    $volume_object = ILT::VolumeObject->new( "Data/volume1.mnc" );
    $view = ILT::View->new_canonical( Top_view );

    $cross1 = ILT::IntersectionObject->new( $plane_object1,
                                            $surface_object );
    $cross2 = ILT::IntersectionObject->new( $plane_object2,
                                            $surface_object );
    $cross3 = ILT::IntersectionObject->new( $plane_object3,
                                            $surface_object );

    $coloured_slice = ILT::ColourObject->new_volume_colouring(
                           ILT::IntersectionObject->new( $plane_object1,
                                                    $volume_object ),
                           $volume_object, Usercc_scale, 0, 180 );
    $coloured_slice->usercc_filename( "Data/contour.map" );

    $scene_object = ILT::UnionObject->new( $cross1, $cross2, $cross3,
                                           $coloured_slice );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(2,1), $image_info );

#--------------------------------------------------------------------------
#
#   image[2][2] is three colour-coded slices with a transparent under colour
#   and a surface in a three d view
#
#--------------------------------------------------------------------------

    $plane_object1 = ILT::PlaneObject->new_canonical( Sagittal_axis, 0 );
    $plane_object2 = ILT::PlaneObject->new_canonical( Coronal_axis, 0 );
    $plane_object3 = ILT::PlaneObject->new_canonical( Transverse_axis, 0 );
    $view = ILT::View->new_arbitrary( -1, -2, -1, 0, 0, 1 );
    $volume_object = ILT::VolumeObject->new( "Data/volume1.mnc" );

    $coloured_slice1 = ILT::ColourObject->new_volume_colouring(
                           ILT::IntersectionObject->new( $plane_object1,
                                                         $volume_object ),
                           $volume_object, Hot_metal_scale, 50, 150 );
    $coloured_slice1->under_colour( "transparent" );

    $coloured_slice2 = ILT::ColourObject->new_volume_colouring(
                           ILT::IntersectionObject->new( $plane_object2,
                                                         $volume_object ),
                           $volume_object, Hot_metal_scale, 50, 150 );
    $coloured_slice2->under_colour( "transparent" );

    $coloured_slice3 = ILT::ColourObject->new_volume_colouring(
                           ILT::IntersectionObject->new( $plane_object3,
                                                         $volume_object ),
                           $volume_object, Hot_metal_scale, 50, 150 );
    $coloured_slice3->under_colour( "transparent" );

    $surface_object = ILT::GeometricObject->new( "Data/surf1.obj" );

    $union = ILT::UnionObject->new( $coloured_slice1,
                                    $coloured_slice2,
                                    $coloured_slice3,
                                    $surface_object );

    $scene_object = ILT::RenderObject->new( $union );
    $scene_object->lighting_state( True );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(2,2), $image_info );

#--------------------------------------------------------------------------
#
#   image[2][3] is a flat shaded view of the ellipsoid and a transformed
#   (rotated by 45 degrees and shifted slightly) smooth shaded copy of it
#
#--------------------------------------------------------------------------

    $view = ILT::View->new_canonical( Top_view );
    $surface_object = ILT::GeometricObject->new( "Data/surf1.obj" );
    $transformed_surface = ILT::TransformObject->new( $surface_object, 
                                                      "Data/rotate.xfm" );

    $render1 = ILT::RenderObject->new( $surface_object );
    $render1->lighting_state( True );
    $render1->surface_shading( Flat_shading );

    $render2 = ILT::RenderObject->new( $transformed_surface );
    $render2->lighting_state( True );
    $render2->surface_shading( Smooth_shading );

    $scene_object = ILT::UnionObject->new( $render1,$render2 );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(2,3), $image_info );

#--------------------------------------------------------------------------
#
#   image[3][0] is an arbitrarily oriented slice through a volume
#
#--------------------------------------------------------------------------

    $plane_object = ILT::PlaneObject->new( 1.0, 2.0, 3.0, 0, 0, 0 );
    $view = ILT::View->new_arbitrary( 1, 2, 3, 0, 1, 0 );
    $volume_object = ILT::VolumeObject->new( "Data/volume1.mnc" );

    $scene_object = 
                  ILT::ColourObject->new_volume_colouring(
                       ILT::IntersectionObject->new( $plane_object,
                                                     $volume_object ),
                       $volume_object, Hot_metal_scale, 50, 150 );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(3,0), $image_info );

#--------------------------------------------------------------------------
#
#   image[3][1] is an arbitrarily oriented slice through a volume,
#   with text annotation added
#
#--------------------------------------------------------------------------

    $plane_object = ILT::PlaneObject->new( 1.0, 2.0, 3.0, 0, 0, 0 );
    $view = ILT::View->new_arbitrary( 1, 2, 3, 0, 1, 0 );
    $volume_object = ILT::VolumeObject->new( "Data/volume1.mnc" );

    $scene_object = 
                  ILT::ColourObject->new_volume_colouring(
                       ILT::IntersectionObject->new( $plane_object,
                                                     $volume_object ),
                       $volume_object, Hot_metal_scale, 50, 150 );

    $text_object = ILT::TextObject->new( "Some Text Annotation", .5, .05 );
    $text_object->font( "-adobe-helvetica-medium-r-normal--14-140-75-75-p-77-iso8859-1" );
    $text_object->colour( "green" );
    $text_object->horizontal_alignment( Align_centre );
    $text_object->vertical_alignment( Align_bottom );

    $scene_object = ILT::UnionObject->new( $scene_object, $text_object );

    $image_info = ILT::ImageInfo->new( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(3,1), $image_info );

#-----------------------------------------------------------------
#   Create a header object
#-----------------------------------------------------------------

    $header = ILT::TextObject->new( "Output of $0", .5, .5 );
    $header->font( "-adobe-helvetica-medium-r-normal--34-240-100-100-p-176-iso8859-1" );
    $header->colour( "yellow" );
    $header->horizontal_alignment( Align_centre );
    $header->vertical_alignment( Align_centre );
    $layout->header( $header );

#-----------------------------------------------------------------
#   Create a footer object
#-----------------------------------------------------------------

    $footer = ILT::TextObject->new( "This is an example footer", .5, .5 );
    $footer->colour( "yellow" );
    $footer->horizontal_alignment( Align_centre );
    $footer->vertical_alignment( Align_centre );
    $layout->footer( $footer );

#--------------------------------------------------------------------------
#
#   now generate the images
#
#--------------------------------------------------------------------------

    $layout->generate_image( "test_smorgasbord.rgb", 1100, 0 );
