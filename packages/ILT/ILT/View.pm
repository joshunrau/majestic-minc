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
#@NAME       : ILT::View
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: A package which implements a 3D view definition
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::View;

    use      strict;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/View.pm,v 1.4 1998-05-22 14:44:35 david Exp $';

#--------------------------------------------------------------------------
# define the class name
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::View";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : proto
#@OUTPUT     : 
#@RETURNS    : a View instance
#@DESCRIPTION: Creates a view instance
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $ )
{
    my( $proto ) = arg_any( shift );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = {};

    bless ($self, $class);

    #--------------------------------------------------------------------------
    # The View object can automatically assign the view direction, up
    # direction, and the bounding box of the view, if they have not been set
    # by the user of the View.  The following three variables indicate whether
    # the user has set these values or not
    #--------------------------------------------------------------------------

    $self->{VIEW_DIRECTION_ASSIGNED} = 0;
    $self->{UP_DIRECTION_ASSIGNED} = 0;
    $self->{BOUNDING_BOX_ASSIGNED} = 0;

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : new_canonical
#@INPUT      : self
#              view_type  : a view enum : Left_view, Right_view, etc
#@OUTPUT     : 
#@RETURNS    : a View instance
#@DESCRIPTION: An alternate constructor for the view which takes one of the
#              enumerated View types
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new_canonical( $$ )
{
    my( $proto )    = arg_any( shift );
    my( $view_type )     = arg_enum( shift, N_view_enums );
    end_args( @_ );

    my( $self );

    $self = new ILT::View;

    #--------------------------------------------------------------------------
    # define the view and up direction based on the view name
    #--------------------------------------------------------------------------

    if( $view_type == Left_view )
    {
        $self->view_direction( 1, 0, 0 );
        $self->up_direction( 0, 0, 1 );
    }
    elsif( $view_type == Right_view )
    {
        $self->view_direction( -1, 0, 0 );
        $self->up_direction( 0, 0, 1 );
    }
    elsif( $view_type == Front_view )
    {
        $self->view_direction( 0, -1, 0 );
        $self->up_direction( 0, 0, 1 );
    }
    elsif( $view_type == Back_view )
    {
        $self->view_direction( 0, 1, 0 );
        $self->up_direction( 0, 0, 1 );
    }
    elsif( $view_type == Top_view )
    {
        $self->view_direction( 0, 0, -1 );
        $self->up_direction( 0, 1, 0 );
    }
    elsif( $view_type == Bottom_view )
    {
        $self->view_direction( 0, 0, 1 );
        $self->up_direction( 0, 1, 0 );
    }
    else
    {
        fatal_error( "Error in View new\n" );
    }

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : new_arbitrary
#@INPUT      : prototype
#              x_view    : 3 components of the view direction
#              y_view
#              z_view
#              x_up      : 3 components of the up direction
#              y_up
#              z_up
#@OUTPUT     : 
#@RETURNS    : a View instance
#@DESCRIPTION: An alternate constructor which takes two arbitrary vectors.
#              It should check for non-colinearity of the vectors, but this
#              is not implemented yet.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new_arbitrary( $$$$$$$ )
{
    my( $proto )    = arg_any( shift );
    my( $x_view )   = arg_real( shift );
    my( $y_view )   = arg_real( shift );
    my( $z_view )   = arg_real( shift );
    my( $x_up )     = arg_real( shift );
    my( $y_up )     = arg_real( shift );
    my( $z_up )     = arg_real( shift );
    end_args( @_ );

    my( $self );

    $self = new ILT::View;

    $self->view_direction( $x_view, $y_view, $z_view );
    $self->up_direction( $x_up, $y_up, $z_up );

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : copy
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : a View instance
#@DESCRIPTION: Creates a copy of the View.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub copy( $ )
{
    my( $self ) = arg_any( shift );
    end_args( @_ );

    my( $copy );

    $copy = new ILT::View;

    if( $self->{VIEW_DIRECTION_ASSIGNED} )
        { $copy->view_direction( $self->view_direction() ); }

    if( $self->{UP_DIRECTION_ASSIGNED} )
        { $copy->up_direction( $self->up_direction() ); }

    if( $self->{BOUNDING_BOX_ASSIGNED} )
        { $copy->bounding_box( $self->bounding_box() ); }

    return( $copy );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : _set_view_direction
#@INPUT      : self
#              x_view   : 3 component vector
#              y_view
#              z_view
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: Private routine to sets the view direction.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  _set_view_direction( $$$$ )
{
    my( $self )   = arg_object( shift, $this_class );
    my( $x_view ) = arg_real( shift );
    my( $y_view ) = arg_real( shift );
    my( $z_view ) = arg_real( shift );
    end_args( @_ );

    $self->{VIEW_DIRECTION} = [ $x_view, $y_view, $z_view ];
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : view_direction
#@INPUT      : self
#              x_view   : 3 component vector (OPTIONAL)
#              y_view
#              z_view
#@OUTPUT     : 
#@RETURNS    : 3 component array
#@DESCRIPTION: Sets or gets the view direction, depending on whether the
#              optional arguments are present.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub view_direction
{
    my( $self )           = arg_object( shift, $this_class );
    my( @view_direction ) = opt_arg_array_of_reals( \@_, 3 );
    end_args( @_ );

    if( @view_direction )
    {
        $self->_set_view_direction( @view_direction );
        $self->{VIEW_DIRECTION_ASSIGNED} = 1;
    }

    return( @{$self->{VIEW_DIRECTION}} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : _set_view_direction
#@INPUT      : self
#              x_up   : 3 component vector
#              y_up
#              z_up
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: Private routine to sets the up direction.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  _set_up_direction( $$$$ )
{
    my( $self )   = arg_object( shift, $this_class );
    my( $x_up )   = arg_real( shift );
    my( $y_up )   = arg_real( shift );
    my( $z_up )   = arg_real( shift );
    end_args( @_ );

    $self->{UP_DIRECTION} = [ $x_up, $y_up, $z_up ];
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : up_direction
#@INPUT      : self
#              x_up   : 3 component vector (OPTIONAL)
#              y_up
#              z_up
#@OUTPUT     : 
#@RETURNS    : 3 component array
#@DESCRIPTION: Sets or gets the up direction, depending on whether the optional
#              arguments are present.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub up_direction
{
    my( $self )         = arg_object( shift, $this_class );
    my( @up_direction ) = opt_arg_array_of_reals( \@_, 3 );
    end_args( @_ );

    if( @up_direction )
    {
        $self->_set_up_direction( @up_direction );
        $self->{UP_DIRECTION_ASSIGNED} = 1;
    }

    return( @{$self->{UP_DIRECTION}} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : _set_bounding_box
#@INPUT      : self
#              x_min   : 6 component vector defining bounding box
#              x_max
#              y_min
#              y_max
#              z_min
#              z_max
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: Private routine to set the view bounding box (which is oriented
#              with the view_direction and up_direction).
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  _set_bounding_box( $$$$ )
{
    my( $self )  = arg_object( shift, $this_class );
    my( $x_min ) = arg_real( shift );
    my( $x_max ) = arg_real( shift );
    my( $y_min ) = arg_real( shift );
    my( $y_max ) = arg_real( shift );
    my( $z_min ) = arg_real( shift );
    my( $z_max ) = arg_real( shift );
    end_args( @_ );

    $self->{X_MIN} = $x_min;
    $self->{X_MAX} = $x_max;
    $self->{Y_MIN} = $y_min;
    $self->{Y_MAX} = $y_max;
    $self->{Z_MIN} = $z_min;
    $self->{Z_MAX} = $z_max;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : bounding_box
#@INPUT      : self
#              x_min   : 6 component vector defining bounding box (OPTIONAL)
#              x_max
#              y_min
#              y_max
#              z_min
#              z_max
#@OUTPUT     : 
#@RETURNS    : 3 component array
#@DESCRIPTION: Sets or gets the bounding box, depending on whether the optional
#              arguments are present.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub bounding_box
{
    my( $self )   = arg_object( shift, $this_class );
    my( @bbox )   = opt_arg_array_of_reals( \@_, 6 );
    end_args( @_ );

    if( @bbox )
    {
        $self->_set_bounding_box( @bbox );
        $self->{BOUNDING_BOX_ASSIGNED} = 1;
    }

    return( $self->{X_MIN}, $self->{X_MAX},
            $self->{Y_MIN}, $self->{Y_MAX},
            $self->{Z_MIN}, $self->{Z_MAX} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_bounding_box_x_min
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : x_min
#@DESCRIPTION: Returns the bounding box min x coordinate
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_bounding_box_x_min
{
    my( $self )   = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->{X_MIN} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_bounding_box_x_max
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : x_max
#@DESCRIPTION: Returns the bounding box min x coordinate
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_bounding_box_x_max
{
    my( $self )   = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->{X_MAX} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_bounding_box_y_min
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : y_min
#@DESCRIPTION: Returns the bounding box min y coordinate
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_bounding_box_y_min
{
    my( $self )   = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->{Y_MIN} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_bounding_box_y_max
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : y_max
#@DESCRIPTION: Returns the bounding box max y coordinate
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_bounding_box_y_max
{
    my( $self )   = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->{Y_MAX} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_bounding_box_z_min
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : z_min
#@DESCRIPTION: Returns the bounding box min z coordinate
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_bounding_box_z_min
{
    my( $self )   = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->{Z_MIN} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_bounding_box_z_max
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : z_max
#@DESCRIPTION: Returns the bounding box max z coordinate
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_bounding_box_z_max
{
    my( $self )   = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->{Z_MAX} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_bounding_box_width
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : width
#@DESCRIPTION: Returns the width of the bounding box
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_bounding_box_width
{
    my( $self )   = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->get_bounding_box_x_max() - $self->get_bounding_box_x_min() );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_bounding_box_height
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : height
#@DESCRIPTION: Returns the height of the bounding box
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_bounding_box_height
{
    my( $self )   = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->get_bounding_box_y_max() - $self->get_bounding_box_y_min() );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : _check_compute_view_orientation
#@INPUT      : self
#              scene_object
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: If the view directions have been explicitly assigned, does
#              nothing.  Otherwise asks the scene_object to give back a default
#              view.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub _check_compute_view_orientation
{
    my( $self )           = arg_object( shift, $this_class );
    my( $scene_object )   = arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    my( $x_view, $y_view, $z_view, $x_up, $y_up, $z_up );

    if( $self->{VIEW_DIRECTION_ASSIGNED} && $self->{UP_DIRECTION_ASSIGNED} )
        { return; }

    if( !$self->{VIEW_DIRECTION_ASSIGNED} &&
        !$self->{UP_DIRECTION_ASSIGNED} )
    {
        #----------------------------------------------------------------------
        # this scene_object method has not been implemented yet.  Will die here.
        #----------------------------------------------------------------------

        ( $x_view, $y_view, $z_view, $x_up, $y_up, $z_up ) =
                                       $scene_object->get_default_view();

        if( $x_view == 0 && $y_view == 0 && $z_view == 0 )
        {
            fatal_error( "Cannot compute a default view for this scene\n" );
        }

        $self->_set_view_direction( $x_view, $y_view, $z_view );
        $self->_set_up_direction( $x_up, $y_up, $z_up );
    }
    else
    {
        fatal_error( "Case not implemented in _check_compute_view_orientation" );
    }
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : _check_compute_bounding_box
#@INPUT      : self
#              scene_object
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: Checks if the bounding box for the current view has been
#              explicitly assigned.  If not computes a bounding box
#              oriented along the view_direction and the horizontal and
#              vertical axes implied by the up_direction
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub _check_compute_bounding_box
{
    my( $self )           = arg_object( shift, $this_class );
    my( $scene_object )   = arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    my( @view_direction, @up_direction,
        $x_min, $x_max, $y_min, $y_max, $z_min, $z_max );

    if( $self->{BOUNDING_BOX_ASSIGNED} )
        { return; }

    @view_direction = $self->view_direction();
    @up_direction = $self->up_direction();

    #--------------------------------------------------------------------------
    # compute the view for the given object
    #--------------------------------------------------------------------------

    ( $x_min, $x_max, $y_min, $y_max, $z_min, $z_max ) =
           $scene_object->compute_bounding_view( \@view_direction,
                                                 \@up_direction, "" );

    #--------------------------------------------------------------------------
    # assign the bounding box to the view
    #--------------------------------------------------------------------------

    $self->_set_bounding_box( $x_min, $x_max, $y_min, $y_max, $z_min, $z_max );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : compute_view_for_object
#@INPUT      : self
#              scene_object
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: Checks if the view direction, up direction, and bounding box
#              for the current view have been explicitly assigned.  If not,
#              computes reasonable values, based on the supplied scene_object.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  compute_view_for_object
{
    my( $self )           = arg_object( shift, $this_class );
    my( $scene_object )   = arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    $self->_check_compute_view_orientation( $scene_object );

    $self->_check_compute_bounding_box( $scene_object );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : eye_position
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : vector of 3 reals
#@DESCRIPTION: Computes the eye_position of the viewing system based on
#              the bounding box and the view orientation.  This routine
#              really should let the user set and get the eye_position,
#              and if the user does not explicitly set the eye_position,
#              it can be automatically computed, similar to the view 
#              direction.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  eye_position
{
    my( $self )           = arg_object( shift, $this_class );
    end_args( @_ );

    my( @view_direction, @up_direction, @horizontal_direction,
        $x_offset, $y_offset, $z_offset, @x_offset, @y_offset, @z_offset,
        @eye, $step_back );

    #--------------------------------------------------------------------------
    # get the three orthogonal vectors from the view and up direction
    #--------------------------------------------------------------------------

    @view_direction = $self->view_direction();
    @up_direction = $self->up_direction();
    @horizontal_direction = cross_vectors( \@view_direction, \@up_direction );
    @up_direction = cross_vectors( \@horizontal_direction,
                                         \@view_direction );

    #--------------------------------------------------------------------------
    # normalize the three orthogonal vectors
    #--------------------------------------------------------------------------

    @view_direction = normalize_vector( @view_direction );
    @horizontal_direction = normalize_vector( @horizontal_direction );
    @up_direction = normalize_vector( @up_direction );

    #--------------------------------------------------------------------------
    # move the eye position back from the centre of the front face (zmin)
    # of the bounding box 
    #--------------------------------------------------------------------------

    $step_back = max( $self->get_bounding_box_x_max() -
                      $self->get_bounding_box_x_min(),
                      $self->get_bounding_box_y_max() -
                      $self->get_bounding_box_y_min(),
                      $self->get_bounding_box_z_max() -
                      $self->get_bounding_box_z_min() );

    $x_offset = ($self->get_bounding_box_x_min() +
                 $self->get_bounding_box_x_max()) / 2;
    $y_offset = ($self->get_bounding_box_y_min() +
                 $self->get_bounding_box_y_max()) / 2;
    $z_offset = $self->get_bounding_box_z_min() - $step_back;

    @x_offset = scale_vector( $x_offset, @horizontal_direction );
    @y_offset = scale_vector( $y_offset, @up_direction );
    @z_offset = scale_vector( $z_offset, @view_direction );

    #--------------------------------------------------------------------------
    # add the three offsets along the view axes to get the eye point
    #--------------------------------------------------------------------------

    @eye = add_vectors( \@x_offset, \@y_offset );
    @eye = add_vectors( \@eye, \@z_offset );

    return( @eye );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : window_width
#@INPUT      : self
#              image_aspect
#@OUTPUT     : 
#@RETURNS    : width of view window
#@DESCRIPTION: Returns the width of the view window.  If the image aspect
#              desired matches the aspect of the bounding box, then the
#              bounding box width is returned.  Otherwise, a width is
#              returned that places the bounding box within the image_aspect
#              shaped window, with either wasted space in the horizontal
#              or vertical direction.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  window_width( $$ )
{
    my( $self )           = arg_object( shift, $this_class );
    my( $image_aspect )   = arg_real( shift, 0.0, 1.0e30 );
    end_args( @_ );

    my( $world_x_size, $world_y_size, $window_width );

    $world_x_size = $self->get_bounding_box_x_max() -
                    $self->get_bounding_box_x_min();
    $world_y_size = $self->get_bounding_box_y_max() -
                    $self->get_bounding_box_y_min();

    #--------------------------------------------------------------------------
    # if fitting world window into a taller, thinner window
    #--------------------------------------------------------------------------

    if( $image_aspect > $world_y_size / $world_x_size )
    {
        $window_width = $world_x_size;
    }
    else        # fitting into a shorter, fatter window
    {
        $window_width = $world_y_size / $image_aspect;
    }

    return( $window_width );
}

#--------------------------------------------------------------------------

1;
