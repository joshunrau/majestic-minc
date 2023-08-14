#!/usr/local/bin/perl5 -w

    use ILT::LayoutInclude;
    use ILT::LayoutUtils;

    $layout = new_grid ILT::ImageLayout ( 1, 1 );

    $plane_object = new_canonical ILT::PlaneObject ( Sagittal_axis, -10 );
    $view = new_canonical ILT::View ( Left_view );
    $volume_object = new ILT::VolumeObject ( "Data/volume1.mnc" );

    $scene_object =
                  new_volume_colouring ILT::ColourObject (
                       new ILT::IntersectionObject ( $plane_object,
                                                     $volume_object ),
                       $volume_object, Hot_metal_scale, 50, 150 );

    $image_info = new ILT::ImageInfo ( $scene_object, $view );
    $layout->image_info( 0, $image_info );

    $layout->generate_image( "ILT_tutorial_example_output/example2.gif", 200, 0 );
