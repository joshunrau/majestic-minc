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
#@NAME       : ColourObject
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: A SceneObject derived class which provides colouring of
#              a sub-object, based on a volume or a solid colour.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr.  3, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::ColourObject;

    use      strict;
    use      vars  qw(@ISA);
    use      ILT::LayoutInclude;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;
    use      ILT::SceneObject::OneSubObject;
    @ISA =   ( "ILT::OneSubObject" );

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/SceneObject/ColourObject.pm,v 1.6 1998-05-22 14:44:45 david Exp $';

#--------------------------------------------------------------------------
# define the name of this package
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::ColourObject";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : prototype
#              sub_object
#@OUTPUT     : 
#@RETURNS    : ColourObject instance
#@DESCRIPTION: Creates an instance of a ColourObject.  At present, this
#              constructor is not called directly.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $$ )
{
    my( $proto )       = arg_string( shift );
    my( $sub_object )  = arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = $class->SUPER::new( $sub_object );

    bless ($self, $class);

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : new_volume_colouring
#@INPUT      : object
#              volume_object
#              method
#              low_limit
#              high_limit
#@OUTPUT     : 
#@RETURNS    : a ColourObject
#@DESCRIPTION: Creates a ColourObject which assigns colour based on a volume
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr.  3, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new_volume_colouring( $$$$$$ )
{
    my( $proto )         = arg_any( shift );
    my( $object )        = arg_object( shift, "ILT::SceneObject" );
    my( $volume_object ) = arg_object( shift, "ILT::VolumeObject" );
    my( $method )        = arg_enum( shift, N_colour_coding_enums );
    my( $low_limit )     = arg_real( shift );
    my( $high_limit )    = arg_real( shift );
    end_args( @_ );

    my( $self );

    $self = new ILT::ColourObject( $object );

    $self->volume( $volume_object );
    $self->method( $method );
    $self->low_limit( $low_limit );
    $self->high_limit( $high_limit );
    $self->under_colour( "black" );
    $self->over_colour( "white" );
    $self->opacity( 1 );

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : volume
#@INPUT      : self
#              volume  OPTIONAL
#@OUTPUT     : 
#@RETURNS    : volume_object
#@DESCRIPTION: Sets or gets the value of the volume object associated with
#              the colouring.  If the optional argument is present, sets the
#              value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub volume( $@ )
{
    my( $self )   =  arg_object( shift, $this_class );
    my( $volume ) =  opt_arg_object( shift, "ILT::VolumeObject" );
    end_args( @_ );

    if( defined($volume) )
        { $self->{VOLUME} = $volume; }

    return( $self->{VOLUME} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : method
#@INPUT      : self
#              method  OPTIONAL
#@OUTPUT     : 
#@RETURNS    : method  (e.g., Gray_scale, Hot_metal_scale, etc.)
#@DESCRIPTION: Sets or gets the value of the method associated with
#              the colouring.  If the optional argument is present, sets the
#              value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub method( $@ )
{
    my( $self )   =  arg_object( shift, $this_class );
    my( $method ) =  opt_arg_enum( shift, N_colour_coding_enums );
    end_args( @_ );

    if( defined($method) )
        { $self->{METHOD} = $method; }

    return( $self->{METHOD} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : low_limit
#@INPUT      : self
#              low_limit  OPTIONAL
#@OUTPUT     : 
#@RETURNS    : value of low_limit
#@DESCRIPTION: Sets or gets the value of the colour coding low limit
#              associated with the colouring.  If the optional argument is
#              present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub low_limit( $@ )
{
    my( $self )      =  arg_object( shift, $this_class );
    my( $low_limit ) =  opt_arg_real( shift );
    end_args( @_ );

    if( defined($low_limit) )
        { $self->{LOW_LIMIT} = $low_limit; }

    return( $self->{LOW_LIMIT} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : high_limit
#@INPUT      : self
#              high_limit  OPTIONAL
#@OUTPUT     : 
#@RETURNS    : value of high_limit
#@DESCRIPTION: Sets or gets the value of the colour coding high limit
#              associated with the colouring.  If the optional argument is
#              present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub high_limit( $@ )
{
    my( $self )      =  arg_object( shift, $this_class );
    my( $high_limit ) =  opt_arg_real( shift );
    end_args( @_ );

    if( defined($high_limit) )
        { $self->{HIGH_LIMIT} = $high_limit; }

    return( $self->{HIGH_LIMIT} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : over_colour
#@INPUT      : self
#              over_colour  OPTIONAL
#@OUTPUT     : 
#@RETURNS    : value of over_colour
#@DESCRIPTION: Sets or gets the value of the colour coding over_colour
#              associated with the colouring.  If the optional argument is
#              present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub over_colour( $@ )
{
    my( $self )        =  arg_object( shift, $this_class );
    my( $over_colour ) =  opt_arg_string( shift );
    end_args( @_ );

    if( defined($over_colour) )
        { $self->{OVER_COLOUR} = $over_colour; }

    return( $self->{OVER_COLOUR} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : under_colour
#@INPUT      : self
#              under_colour  OPTIONAL
#@OUTPUT     : 
#@RETURNS    : value of under_colour
#@DESCRIPTION: Sets or gets the value of the colour coding under_colour
#              associated with the colouring.  If the optional argument is
#              present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub under_colour( $@ )
{
    my( $self )        =  arg_object( shift, $this_class );
    my( $under_colour ) =  opt_arg_string( shift );
    end_args( @_ );

    if( defined($under_colour) )
        { $self->{UNDER_COLOUR} = $under_colour; }

    return( $self->{UNDER_COLOUR} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : opacity
#@INPUT      : self
#              opacity  OPTIONAL
#@OUTPUT     : 
#@RETURNS    : value of opacity
#@DESCRIPTION: Sets or gets the value of the colour coding opacity 
#              associated with the colouring.  If the optional argument is
#              present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub opacity( $@ )
{
    my( $self )    =  arg_object( shift, $this_class );
    my( $opacity ) =  opt_arg_real( shift, 0, 1 );
    end_args( @_ );

    if( defined($opacity) )
        { $self->{OPACITY} = $opacity; }

    return( $self->{OPACITY} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : usercc_filename
#@INPUT      : self
#              usercc_filename  OPTIONAL
#@OUTPUT     : 
#@RETURNS    : filename
#@DESCRIPTION: Sets or gets the value of the user defined colour coding
#              filename.  When the method is Usercc_scale, the file is
#              read.  The file is an ascii file with a piecewise linear
#              function mapping values to colour, eg.
#                     10   black
#                     20   green
#                     30   1 0 0
#              Note that the first and last values in the file are mapped
#              to the low and high limits, and therefore, in most cases,
#              you should set the low and high limits to the first and
#              last values in the file.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub usercc_filename( $@ )
{
    my( $self )            =  arg_object( shift, $this_class );
    my( $usercc_filename ) =  opt_arg_string( shift );
    end_args( @_ );

    if( defined($usercc_filename) )
        { $self->{USERCC_FILENAME} = $usercc_filename; }

    return( $self->{USERCC_FILENAME} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : make_ray_trace_args
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : string
#@DESCRIPTION: Creates a string of arguments to be passed to ray trace
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub make_ray_trace_args( $ )
{
    my( $self )          =  arg_object( shift, $this_class );
    end_args( @_ );

    my( $args, $method_name, $continuity, $object_args, $full_args );

    #--------------------------------------------------------------------------
    # Convert the enum method to a ray_trace argument string
    #--------------------------------------------------------------------------

    if( $self->method() == Gray_scale )
        { $method_name = "-gray"; }
    elsif( $self->method() == Hot_metal_scale )
        { $method_name = "-hot"; }
    elsif( $self->method() == Spectral_scale )
        { $method_name = "-spectral"; }
    elsif( $self->method() == Red_scale )
        { $method_name = "-red"; }
    elsif( $self->method() == Green_scale )
        { $method_name = "-green"; }
    elsif( $self->method() == Blue_scale )
        { $method_name = "-blue"; }
    elsif( $self->method() == Over_colour_scale )
        { $method_name = "-colour"; }
    elsif( $self->method() == Usercc_scale )
        { $method_name = sprintf( "-usercc %s", $self->usercc_filename() );  }
    else
        { fatal_error( "Incorrect method type in colour object\n" ); }

    #--------------------------------------------------------------------------
    # Convert the volume interpolation method to a ray_trace continuity integer
    #--------------------------------------------------------------------------

    if( $self->volume()->interpolation() == Nearest_neighbour_interpolation )
        { $continuity = -1; }
    elsif( $self->volume()->interpolation() == Linear_interpolation )
        { $continuity = 0; }
    elsif( $self->volume()->interpolation() == Cubic_interpolation )
        { $continuity = 2; }

    #--------------------------------------------------------------------------
    # Assemble the colour coding arguments to ray_trace
    #--------------------------------------------------------------------------

    $args = sprintf( "-under %s -over %s %s %g %g %s %d %g",
                     $self->under_colour(),
                     $self->over_colour(), $method_name,
                     $self->low_limit(),
                     $self->high_limit(),
                     $self->volume()->filename(),
                     $continuity, $self->opacity() );

    #--------------------------------------------------------------------------
    # Get the ray_trace arguments of the sub object
    #--------------------------------------------------------------------------

    $object_args  = $self->SUPER::make_ray_trace_args();

    #--------------------------------------------------------------------------
    # Finally, assemble the colour coding arguments, the sub object arguments,
    # and a delete volume (so subsequent objects are not colour coded).
    #--------------------------------------------------------------------------

    $full_args = "-reverse_order_colouring $args $object_args -delete_volume 0";

    return( $full_args );
}

#--------------------------------------------------------------------------

1;
