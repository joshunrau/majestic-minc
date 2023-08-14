#!/usr/local/bin/perl5 -w
# ----------------------------------------------------------------------------
#@COPYRIGHT  :
#              Copyright 1993,1994,1995,1996,1997,1998 David MacDonald,
#              McConnell Brain Imaging Centre,
#              Montreal Neurological Institute, McGill University.
#              Permission to use, copy, modify, and distribute this
#              software and its documentation for any purpose and without
#              fee is hereby granted, provided that the above copyright
#              notice appear in all copies.  The author and McGill University
#              make no representations about the suitability of this
#              software for any purpose.  It is provided "as is" without
#              express or implied warranty.
#-----------------------------------------------------------------------------

#----------------------------- MNI Header -----------------------------------
#@NAME       : ILT::ImageInfo
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: An object class to define how to render an image.  This consists
#              of a scene object, a view object, and auxiliary information such
#              as background colour, and later other global rendering
#              information such as supersampling.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::ImageInfo;

    use      strict;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;
    use      ILT::Executables;
    use      ILT::LayoutInclude;

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/ImageInfo.pm,v 1.6 1998-09-18 13:30:02 david Exp $';

#--------------------------------------------------------------------------
# define the class of this package
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::ImageInfo";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : prototype
#              scene_object
#              view
#@OUTPUT     : 
#@RETURNS    : ImageInfo object
#@DESCRIPTION: Creates an image info object which contains a scene object and
#              a view.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $$$ )
{
    my( $proto )        = arg_any( shift );
    my( $scene_object ) = arg_object( shift, "ILT::SceneObject" );
    my( $view )         = arg_object( shift, "ILT::View" );

    my $class = ref($proto) || $proto;
    my $self  = {};

    bless ($self, $class);

    $self->scene_object( $scene_object );
    $self->scene_view( $view );

    return( $self );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : scene_object
#@INPUT      : self
#              scene_object  optional
#@OUTPUT     : 
#@RETURNS    : scene_object
#@DESCRIPTION: Sets or gets the value of the scene object stored in the
#              ImageInfo object.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub scene_object( $@ )
{
    my( $self )         = arg_object( shift, $this_class );
    my( $scene_object ) = opt_arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    #--------------------------------------------------------------------------
    # set the value if the optional parameter is present
    #--------------------------------------------------------------------------

    if( defined($scene_object) )
        { $self->{OBJECT} = $scene_object; }

    return( $self->{OBJECT} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : scene_view
#@INPUT      : self
#              scene_view  optional
#@OUTPUT     : 
#@RETURNS    : scene_view
#@DESCRIPTION: Sets or gets the value of the scene view stored in the
#              ImageInfo object.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub scene_view( $@ )
{
    my( $self ) = arg_object( shift, $this_class );
    my( $view ) = opt_arg_object( shift, "ILT::View" );
    end_args( @_ );

    #--------------------------------------------------------------------------
    # set the value if the optional parameter is present
    #--------------------------------------------------------------------------

    if( defined($view) )
        { $self->{VIEW} = $view; }

    return( $self->{VIEW} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : background_colour
#@INPUT      : self
#              background_colour  optional
#@OUTPUT     : 
#@RETURNS    : background_colour
#@DESCRIPTION: Sets or gets the value of the background colour stored in the
#              ImageInfo object.  The colour is a string containing a colour 
#              name.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub background_colour( $@ )
{
    my( $self )              = arg_object( shift, $this_class );
    my( $background_colour ) = opt_arg_string( shift );
    end_args( @_ );

    #--------------------------------------------------------------------------
    # set the value if the optional parameter is present
    #--------------------------------------------------------------------------

    if( defined($background_colour) )
        { $self->{BACKGROUND_COLOUR} = $background_colour; }

    return( $self->{BACKGROUND_COLOUR} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : create_image
#@INPUT      : self
#              filename
#              x_size
#              y_size
#              view
#@OUTPUT     : 
#@RETURNS    : scene_object
#@DESCRIPTION: Sets or gets the value of the scene object stored in the
#              ImageInfo object.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub create_image( $$$$$ )
{
    my( $self )       = arg_object( shift, $this_class );
    my( $filename )   = arg_string( shift );
    my( $x_size )     = arg_int( shift, 1, 1e30 );
    my( $y_size )     = arg_int( shift, 1, 1e30 );
    end_args( @_ );

    my( $view, $geom_args, $view_args, $args,
        @view_direction, @up_direction, @eye, $window_width, $bg,
        $text_args );

    #--------------------------------------------------------------------------
    # create the ray trace arguments to render the scene geometry
    #--------------------------------------------------------------------------

    $geom_args = $self->{OBJECT}->make_ray_trace_args();

    #--------------------------------------------------------------------------
    # create the ray trace arguments to define the view
    #--------------------------------------------------------------------------

    $view = $self->scene_view();
    @view_direction = $view->view_direction();
    @up_direction = $view->up_direction();
    @eye = $view->eye_position();
    $window_width = $view->window_width( $y_size / $x_size );

    $view_args = sprintf( "-ortho -view %g %g %g %g %g %g " .
                          " -eye %g %g %g " .
                          " -window_width %g" ,
                          $view_direction[0],
                          $view_direction[1],
                          $view_direction[2],
                          $up_direction[0],
                          $up_direction[1],
                          $up_direction[2],
                          $eye[0], $eye[1], $eye[2], $window_width );

    #--------------------------------------------------------------------------
    # if background colour specified, create the argument
    #--------------------------------------------------------------------------

    if( defined($self->background_colour() ) )
        { $bg = sprintf( "-bg '%s'", $self->background_colour() ); }
    else
        { $bg = "-bg black"; }

    #--------------------------------------------------------------------------
    # assemble all ray_trace arguments into a single string
    #--------------------------------------------------------------------------

    $args = "-output $filename -nolight $bg " .
               " -size $x_size $y_size " .
               " $view_args " .
               " $geom_args ";

    #--------------------------------------------------------------------------
    # execute the command to render the image
    #--------------------------------------------------------------------------

    run_executable( "ray_trace", $args );

    #--------------------------------------------------------------------------
    # Now add any text to the image
    #--------------------------------------------------------------------------

    $text_args = $self->{OBJECT}->get_text_image_magick_args( $x_size, $y_size);

    if( $text_args ne "" )
    {
        run_executable( "mogrify", $text_args . " $filename" );
    }
}

1;
