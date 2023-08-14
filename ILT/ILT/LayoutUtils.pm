#!/usr/local/bin/perl5 -w

    package  ILT::LayoutUtils;
    use  strict;
    use  Carp;
    use  ILT::ProgUtils;
    use  ILT::Executables;
    use UNIVERSAL;


    BEGIN {
        use Exporter   ();
        use vars       qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

        @ISA         = qw(Exporter);
        %EXPORT_TAGS = ( );     # eg: TAG => [ qw!name1 name2! ],

        # your exported package globals go here,
        # as well as any optionally exported functions

        @EXPORT_OK      = ();
        @EXPORT         = qw( &ILT_version
                          &compute_geometry_file_bounding_view
                          &add_vectors
                          &cross_vectors
                          &normalize_vector
                          &scale_vector
                          &Nearest_neighbour_interpolation
                          &Linear_interpolation
                          &Cubic_interpolation
                          &N_interpolation_enums
                          &Flat_shading
                          &Smooth_shading
                          &N_shading_enums
                          &Gray_scale
                          &Hot_metal_scale
                          &Spectral_scale
                          &Red_scale
                          &Green_scale
                          &Blue_scale
                          &Over_colour_scale
                          &Usercc_scale
                          &N_colour_coding_enums
                          &Sagittal_axis
                          &Coronal_axis
                          &Transverse_axis
                          &N_axis_enums
                          &Left_view
                          &Right_view
                          &Top_view
                          &Bottom_view
                          &Back_view
                          &Front_view
                          &N_view_enums
                          &Align_left
                          &Align_centre
                          &Align_right
                          &Align_top
                          &Align_bottom
                          &N_align_enums
                        );
    }

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/LayoutUtils.pm,v 1.8 2006-06-20 21:51:11 claude Exp $';

#---------------

sub ILT_version()
{
    return( "1.2" );
}

#---------------

sub compute_geometry_file_bounding_view( $$$@ )
{
    my( $filename )           = arg_string( shift );
    my( $view_direction_ref ) = arg_array_ref( shift, 3 );
    my( $up_direction_ref )   = arg_array_ref( shift, 3 );
    my( $transform )          = arg_string( shift );
    end_args( @_ );

    my( @view_direction, @up_direction, $args,
        $x_centre, $y_centre, $z_centre, $out,
        $x_min, $x_max, $y_min, $y_max, $z_min, $z_max );

    @view_direction = @$view_direction_ref;
    @up_direction = @$up_direction_ref;

    if( !defined($transform) )
        { $transform = ""; }

    $args = sprintf( "%s %g %g %g %g %g %g %s",
                        $filename,
                        $view_direction[0],
                        $view_direction[1],
                        $view_direction[2],
                        $up_direction[0],
                        $up_direction[1],
                        $up_direction[2],
                        $transform );

    $out = get_output_of_command( "compute_bounding_view", $args );

    $out =~ /Bounding_box:\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/m;

    $x_min = $1;
    $x_max = $2;
    $y_min = $3;
    $y_max = $4;
    $z_min = $5;
    $z_max = $6;

    return( $x_min, $x_max, $y_min, $y_max, $z_min, $z_max );
}

#------ Enumerated types for interpolation methods

sub Nearest_neighbour_interpolation() { return( 0 ); }
sub Linear_interpolation()            { return( 1 ); }
sub Cubic_interpolation()             { return( 2 ); }
sub N_interpolation_enums()           { return( 3 ); }

#------ Enumerated types for interpolation methods

sub Flat_shading()                    { return( 0 ); }
sub Smooth_shading()                  { return( 1 ); }
sub N_shading_enums()                 { return( 2 ); }

#----- Enumerated types for colour coding methods

sub Gray_scale()               { return( 0 ); }
sub Hot_metal_scale()          { return( 1 ); }
sub Spectral_scale()           { return( 2 ); }
sub Red_scale()                { return( 3 ); }
sub Green_scale()              { return( 4 ); }
sub Blue_scale()               { return( 5 ); }
sub Over_colour_scale()        { return( 6 ); }
sub Usercc_scale()             { return( 7 ); }
sub N_colour_coding_enums()    { return( 8 ); }

#----- Enumerated types for views

sub Left_view()       { return( 0 ); }
sub Right_view()      { return( 1 ); }
sub Top_view()        { return( 2 ); }
sub Bottom_view()     { return( 3 ); }
sub Back_view()       { return( 4 ); }
sub Front_view()      { return( 5 ); }
sub N_view_enums()    { return( 6 ); }

#----- Enumerated types for axes

sub Sagittal_axis()       { return( 0 ); }
sub Coronal_axis()        { return( 1 ); }
sub Transverse_axis()     { return( 2 ); }
sub N_axis_enums()        { return( 3 ); }

#----- Enumerated types for axes

sub Align_left()          { return( 0 ); }
sub Align_centre()        { return( 1 ); }
sub Align_right()         { return( 2 ); }
sub Align_bottom()        { return( 0 ); }
sub Align_top()           { return( 2 ); }
sub N_align_enums()       { return( 3 ); }

#---------- vector operations

sub add_vectors( $$ )
{
    my( $v1 ) = arg_array_ref( shift, 3 );
    my( $v2 ) = arg_array_ref( shift, 3 );
    end_args( @_ );

    return( $$v1[0] + $$v2[0], $$v1[1] + $$v2[1], $$v1[2] + $$v2[2] );
}

sub cross_vectors( $$ )
{
    my( $v1 ) = arg_array_ref( shift, 3 );
    my( $v2 ) = arg_array_ref( shift, 3 );
    end_args( @_ );

    return( $$v1[1] * $$v2[2] - $$v1[2] * $$v2[1],
            $$v1[2] * $$v2[0] - $$v1[0] * $$v2[2],
            $$v1[0] * $$v2[1] - $$v1[1] * $$v2[0] );
}

sub scale_vector( $@ )
{
    my( $scale ) = arg_real( shift );
    my( $x )     = arg_real( shift );
    my( $y )     = arg_real( shift );
    my( $z )     = arg_real( shift );
    end_args( @_ );
 
    return( $scale * $x, $scale * $y, $scale * $z );
}

sub normalize_vector( @ )
{
    my( $x )     = arg_real( shift );
    my( $y )     = arg_real( shift );
    my( $z )     = arg_real( shift );
    end_args( @_ );

    my( $length );
 
    $length = $x * $x + $y * $y + $z * $z;
    $length = sqrt( $length );

    return( scale_vector( 1/$length, $x, $y, $z ) );
}

#---------------- return value for 'use'

1;

