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
#@NAME       : ILT::PlaneObject
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Object class to represent an infinite 3D plane.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::PlaneObject;

    use      strict;
    use      vars  qw(@ISA);
    use      ILT::LayoutInclude;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;
    @ISA =   ( "ILT::SceneObject" );

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/SceneObject/PlaneObject.pm,v 1.5 1998-05-22 14:44:45 david Exp $';

#--------------------------------------------------------------------------
# define the name of the package
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::PlaneObject";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : prototype
#              x_origin
#              y_origin
#              z_origin
#              x_normal
#              y_normal
#              z_normal
#@OUTPUT     : 
#@RETURNS    : instance of PlaneObject
#@DESCRIPTION: Creates a new plane, with the given origin and orientation.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $$$$$$$ )
{
    my( $proto ) = arg_any( shift );
    my( $x_normal ) = arg_real( shift );
    my( $y_normal ) = arg_real( shift );
    my( $z_normal ) = arg_real( shift );
    my( $x_origin ) = arg_real( shift );
    my( $y_origin ) = arg_real( shift );
    my( $z_origin ) = arg_real( shift );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = $class->SUPER::new();

    bless ($self, $class);

    $self->plane_origin( $x_origin, $y_origin, $z_origin );
    $self->plane_normal( $x_normal, $y_normal, $z_normal );

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : new_canonical
#@INPUT      : self
#              axis               : enum of axis orientation
#              plane_origin       : origin of plane
#@OUTPUT     : 
#@RETURNS    : instance of plane
#@DESCRIPTION: Creates a plane object oriented along the x, y, or z axis.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new_canonical( $$$ )
{
    my( $proto )           = arg_any( shift );
    my( $axis )            = arg_enum( shift, N_axis_enums );
    my( $plane_origin )    = arg_real( shift );
    end_args( @_ );
 
    my( $self );

    if( $axis == Sagittal_axis )
    {
        $self = new ILT::PlaneObject( 1, 0, 0, $plane_origin, 0, 0 );
    }
    elsif( $axis == Coronal_axis )
    {
        $self = new ILT::PlaneObject( 0, 1, 0, 0, $plane_origin, 0 );
    }
    elsif( $axis == Transverse_axis )
    {
        $self = new ILT::PlaneObject( 0, 0, 1, 0, 0, $plane_origin );
    }

    return( $self );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : plane_normal
#@INPUT      : self
#              normal       : OPTIONAL array of 3 reals
#@OUTPUT     : 
#@RETURNS    : 3 reals
#@DESCRIPTION: Gets and optionally sets (if normal arg present) the normal
#              of the plane
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub plane_normal( $@ )
{
    my( $self ) = arg_object( shift, $this_class );
    my( @normal ) = opt_arg_array_of_reals( \@_, 3 );
    end_args( @_ );

    if( @normal )
        { $self->{PLANE_NORMAL} = [ @normal ]; }

    return( @{$self->{PLANE_NORMAL}} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : plane_origin
#@INPUT      : self
#              origin       : OPTIONAL array of 3 reals
#@OUTPUT     : 
#@RETURNS    : 3 reals
#@DESCRIPTION: Gets and optionally sets (if origin arg present) the origin
#              of the plane
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub plane_origin( $@ )
{
    my( $self ) = arg_object( shift, $this_class );
    my( @origin ) = opt_arg_array_of_reals( \@_, 3 );
    end_args( @_ );

    if( @origin )
        { $self->{PLANE_ORIGIN} = [ @origin ]; }

    return( @{$self->{PLANE_ORIGIN}} );
}

#--------------------------------------------------------------------------

1;
