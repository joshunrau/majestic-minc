#!/usr/local/bin/perl5 -w

    package  ILT::SceneObject;

    use      strict;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/SceneObject.pm,v 1.6 2011-02-04 16:48:13 alex Exp $';

my( $this_class ) = "ILT::SceneObject";

sub new( $ )
{
    my( $proto ) = arg_any( shift );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = {};

    bless ($self, $class);
    return $self;
}

sub  create_temp_geometry_file( $ )
{
    my( $self )   = arg_object( shift, $this_class );
    end_args( @_ );

#--- nothing to do for base class
}

sub  make_ray_trace_args()
{
    my( $self )   = arg_object( shift, $this_class );
    end_args( @_ );

    fatal_error( "called make_ray_trace_args() for base class\n" );
}

sub  get_plane_intersection( $$$$ )
{
    my( $self )              = arg_object( shift, $this_class );
    my( $plane_origin_ref )  = arg_array_ref( shift, 3 );
    my( $plane_normal_ref )  = arg_array_ref( shift, 3 );
    my( $output_file )       = arg_string( shift );
    end_args( @_ );

    fatal_error( "called get_plane_intersection() for base class\n" );
}

sub  compute_bounding_view( $$$$ )
{
    my( $self )                = arg_object( shift, $this_class );
    my( $view_direction_ref )  =  arg_array_ref( shift, 3 );
    my( $up_direction_ref )    =  arg_array_ref( shift, 3 );
    my( $transform )           =  arg_string( shift );
    end_args( @_ );

    fatal_error( "called compute_bounding_view() for base class\n" );
}

sub  get_default_view( $ )
{
    my( $self )                = arg_object( shift, $this_class );
    end_args( @_ );

    fatal_error( "called get_default_view() for base class\n" );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_text_image_magick_args
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : a string
#@DESCRIPTION: Returns a string which are arguments to pass to the
#              ImageMagick program to do text rendering.  For generic
#              objects, does nothing.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : July 21, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub get_text_image_magick_args( $$$ )
{
    my( $self )             =  arg_object( shift, $this_class );
    my( $viewport_x_size )  =  arg_real( shift, 0, 1e30 );
    my( $viewport_y_size )  =  arg_real( shift, 0, 1e30 );
    end_args( @_ );

    return( "" );
}

1;
