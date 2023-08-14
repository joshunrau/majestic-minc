#!/usr/local/bin/perl5 -w

    use ILT::LayoutInclude;
    use ILT::LayoutUtils;

    $layout = new_grid ILT::ImageLayout ( 1, 1 );

    $surface = new ILT::GeometricObject ( "Data/surf2.obj" );
    $view = new_canonical ILT::View ( Left_view );
    $image_info = new ILT::ImageInfo ( $surface, $view );
    $layout->image_info( 0, $image_info );

    $layout->generate_image( "ILT_tutorial_example_output/example1.gif", 200, 0 );
