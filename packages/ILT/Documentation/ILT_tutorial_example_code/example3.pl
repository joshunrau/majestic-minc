#!/usr/local/bin/perl5 -w

    use ILT::LayoutInclude;
    use ILT::LayoutUtils;

    $layout = new_grid ILT::ImageLayout ( 1, 1 );

    $plane_object = new_canonical ILT::PlaneObject ( Sagittal_axis, -10 );
    $view = new_canonical ILT::View ( Left_view );
    $volume_object = new ILT::VolumeObject ( "Data/volume1.mnc" );
    $surface_object = new ILT::GeometricObject ( "Data/surf2.obj" );

    $cross = new ILT::IntersectionObject ( $plane_object, $surface_object );

    $slice = new_volume_colouring ILT::ColourObject (
                       new ILT::IntersectionObject ( $plane_object,
                                                     $volume_object ),
                       $volume_object, Gray_scale, 50, 150 );

    $scene_object = new ILT::UnionObject ( $cross, $slice );

    $image_info = new ILT::ImageInfo ( $scene_object, $view );
    $layout->image_info( 0, $image_info );

    $layout->generate_image( "ILT_tutorial_example_output/example3.gif", 200, 0 );
