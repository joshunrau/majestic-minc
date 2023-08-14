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
#@NAME       : ILT::RenderObject
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Object class to represent render parameters.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::RenderObject;

    use      strict;
    use      vars  qw(@ISA);
    use      ILT::LayoutInclude;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;
    use      ILT::SceneObject::OneSubObject;
    @ISA =   ( "ILT::OneSubObject" );

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/SceneObject/RenderObject.pm,v 1.6 1998-05-22 14:44:45 david Exp $';

#--------------------------------------------------------------------------
# define the name of this class
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::RenderObject";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : self
#              sub-object
#@OUTPUT     : 
#@RETURNS    : instance of a render object
#@DESCRIPTION: Constructs a render object, given another object
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $$ )
{
    my( $proto )      = arg_any( shift );
    my( $sub_object ) = arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = $class->SUPER::new( $sub_object );

    bless ($self, $class);

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : line_width
#@INPUT      : self
#              line_width  (OPTIONAL)
#@OUTPUT     : 
#@RETURNS    : line_width
#@DESCRIPTION: Gets and optionally sets (if optional argument is present) the
#              value of this parameter.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub line_width( $@ )
{
    my( $self )       = arg_object( shift, $this_class );
    my( $line_width ) = opt_arg_real( shift, 0 );
    end_args( @_ );

    if( defined($line_width) )
        { $self->{LINE_WIDTH} = $line_width; }

    return( $self->{LINE_WIDTH} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : lighting_state
#@INPUT      : self
#              lighting_state  (OPTIONAL)
#@OUTPUT     : 
#@RETURNS    : lighting_state (True or False)
#@DESCRIPTION: Gets and optionally sets (if optional argument is present) the
#              value of this parameter.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub lighting_state( $@ )
{
    my( $self )           = arg_object( shift, $this_class );
    my( $lighting_state ) = opt_arg_enum( shift, N_boolean_enums );
    end_args( @_ );

    if( defined($lighting_state) )
        { $self->{LIGHTING_STATE} = $lighting_state; }

    return( $self->{LIGHTING_STATE} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : surface_shading
#@INPUT      : self
#              surface_shading  (OPTIONAL)
#@OUTPUT     : 
#@RETURNS    : surface_shading (Flat_shading or Smooth_shading)
#@DESCRIPTION: Gets and optionally sets (if optional argument is present) the
#              value of this parameter.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub surface_shading( $@ )
{
    my( $self )            = arg_object( shift, $this_class );
    my( $surface_shading ) = opt_arg_enum( shift, N_shading_enums );
    end_args( @_ );

    if( defined($surface_shading) )
        { $self->{SURFACE_SHADING} = $surface_shading; }

    return( $self->{SURFACE_SHADING} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : make_ray_trace_args
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : string
#@DESCRIPTION: Constructs a string to represent the ray_trace arguments 
#              corresponding to the render parameters, and returns this
#              string concatenated with the sub-object ray_trace string.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub make_ray_trace_args( $ )
{
    my( $self )            = arg_object( shift, $this_class );
    end_args( @_ );

    my( $pre_args, $post_args, $args, $object_args, $full_args );

    #--------------------------------------------------------------------------
    # initialize the arguments to go before and after the sub-objects's
    #--------------------------------------------------------------------------

    $pre_args = "";
    $post_args = "";

    #--------------------------------------------------------------------------
    # if necessary, add line width arguments
    #--------------------------------------------------------------------------

    if( defined($self->line_width()) )
    {
        $pre_args = $pre_args . sprintf( " -line_width %g",
                                         $self->line_width() );
        $post_args = $post_args . " -line_width -1";
    }

    #--------------------------------------------------------------------------
    # if necessary, add lighting state arguments
    #--------------------------------------------------------------------------

    if( defined($self->line_width()) )
    {
        $pre_args = $pre_args . sprintf( " -line_width %g",
                                         $self->line_width() );
        $post_args = $post_args . " -line_width -1";
    }

    if( defined($self->lighting_state()) )
    {
        if( $self->lighting_state() )
        {
            $pre_args = $pre_args . " -light ";
            $post_args = $post_args . " -nolight ";
        }
        else
        {
            $pre_args = $pre_args . " -nolight ";
            $post_args = $post_args . " -nolight ";
        }
    }

    #--------------------------------------------------------------------------
    # if necessary, add surface shading arguments
    #--------------------------------------------------------------------------

    if( defined($self->surface_shading()) )
    {
        if( $self->surface_shading() == Linear_interpolation )
        {
            $pre_args = $pre_args . " -smooth ";
            $post_args = $post_args . " -smooth ";
        }
        else
        {
            $pre_args = $pre_args . " -flat ";
            $post_args = $post_args . " -smooth ";
        }
    }

    #--------------------------------------------------------------------------
    # get the arguments of the sub-object
    #--------------------------------------------------------------------------

    $object_args  = $self->SUPER::make_ray_trace_args();

    #--------------------------------------------------------------------------
    # assemble the arguments and return them
    #--------------------------------------------------------------------------

    $full_args = "$pre_args $object_args $post_args" ;

    return( $full_args );
}

1;
