#!/usr/local/bin/perl5 -w

    use ILT::LayoutInclude;
    use ILT::LayoutUtils;

    $layout = new_grid ILT::ImageLayout ( 2, 2 );

#---------------------------------  image 0 0

    $plane_object = new_canonical ILT::PlaneObject ( Sagittal_axis, -10 );
    $view = new_canonical ILT::View ( Left_view );
    $volume_object = new ILT::VolumeObject ( "Data/volume1.mnc" );

    $scene_object =
                  new_volume_colouring ILT::ColourObject (
                       new ILT::IntersectionObject ( $plane_object,
                                                     $volume_object ),
                       $volume_object, Hot_metal_scale, 50, 150 );

    $image_info = new ILT::ImageInfo ( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(0,0), $image_info );

#---------------------------------  image 0 1

    $plane_object = new_canonical ILT::PlaneObject ( Sagittal_axis, -40 );
    $view = new_canonical ILT::View ( Left_view );
    $volume_object = new ILT::VolumeObject ( "Data/volume1.mnc" );

    $scene_object =
                  new_volume_colouring ILT::ColourObject (
                       new ILT::IntersectionObject ( $plane_object,
                                                     $volume_object ),
                       $volume_object, Hot_metal_scale, 50, 150 );

    $image_info = new ILT::ImageInfo ( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(0,1), $image_info );

#---------------------------------  image 1 0

    $plane_object = new_canonical ILT::PlaneObject ( Transverse_axis, 10 );
    $view = new_canonical ILT::View ( Top_view );
    $volume_object = new ILT::VolumeObject ( "Data/volume1.mnc" );

    $scene_object =
                  new_volume_colouring ILT::ColourObject (
                       new ILT::IntersectionObject ( $plane_object,
                                                     $volume_object ),
                       $volume_object, Hot_metal_scale, 50, 150 );

    $image_info = new ILT::ImageInfo ( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(1,0), $image_info );

#---------------------------------  image 1 1

    $plane_object = new_canonical ILT::PlaneObject ( Transverse_axis, 40 );
    $view = new_canonical ILT::View ( Top_view );
    $volume_object = new ILT::VolumeObject ( "Data/volume1.mnc" );

    $scene_object =
                  new_volume_colouring ILT::ColourObject (
                       new ILT::IntersectionObject ( $plane_object,
                                                     $volume_object ),
                       $volume_object, Hot_metal_scale, 50, 150 );

    $image_info = new ILT::ImageInfo ( $scene_object, $view );
    $layout->image_info( $layout->row_col_to_index(1,1), $image_info );

#--------------- create the header

    $text_object = ILT::TextObject->new( "A Header", .5, .5 );
    $text_object->colour( "blue" );
    $text_object->horizontal_alignment( Align_centre );
    $text_object->vertical_alignment( Align_centre );
    $text_object->font( "12x24" );
    $layout->header( $text_object );

#--------------- create the footer

    $text_object = ILT::TextObject->new( "This is the footer", .5, .5 );
    $text_object->colour( "red" );
    $text_object->horizontal_alignment( Align_centre );
    $text_object->vertical_alignment( Align_centre );
    $layout->footer( $text_object );

#--------------- render the images

    $layout->generate_image( "ILT_tutorial_example_output/example5.gif", 200, 0 );
