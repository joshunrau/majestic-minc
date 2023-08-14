#!/usr/local/bin/perl5 -w

    use ILT::LayoutInclude;
    use ILT::LayoutUtils;

    $layout = new_grid ILT::ImageLayout ( 1, 1 );

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

    $image_info = new ILT::ImageInfo ( $scene_object, $view );
    $layout->image_info( 0, $image_info );

    $layout->generate_image( "ILT_tutorial_example_output/example4.gif", 200, 0 );
