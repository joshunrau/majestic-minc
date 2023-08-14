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
#@NAME       : ILT::VolumeObject
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Object class to implement volumes.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::VolumeObject;

    use      strict;
    use      vars  qw(@ISA);
    use      ILT::LayoutInclude;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;
    use      ILT::Executables;
    @ISA =   ( "ILT::SceneObject" );

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/SceneObject/VolumeObject.pm,v 1.5 1998-05-22 14:44:45 david Exp $';

#--------------------------------------------------------------------------
# define the name of this class
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::VolumeObject";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : prototype
#              volume_filename
#              interpolation
#@OUTPUT     : 
#@RETURNS    : instance of volume object
#@DESCRIPTION: Constructs a volume object
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $$$ )
{
    my( $proto )           = arg_any( shift );
    my( $volume_filename ) = arg_string( shift );
    my( $interpolation )   = opt_arg_enum( shift, N_interpolation_enums );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = $class->SUPER::new();

    bless ($self, $class);

    if( !defined($interpolation) )
        { $interpolation = Linear_interpolation; }

    $self->filename( $volume_filename );
    $self->interpolation( $interpolation );

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : filename
#@INPUT      : self
#              filename   (OPTIONAL)
#@OUTPUT     : 
#@RETURNS    : filename
#@DESCRIPTION: Returns and optionally sets (if optional argument is present),
#              the name of the file containing the volume.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub filename( $$ )
{
    my( $self )     = arg_object( shift, $this_class );
    my( $filename ) = opt_arg_string( shift );
    end_args( @_ );

    if( defined($filename) )
        { $self->{FILENAME} = $filename; }

    return( $self->{FILENAME} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : interpolation
#@INPUT      : self
#              interpolation   (OPTIONAL)   : Linear_interpolation, etc.
#@OUTPUT     : 
#@RETURNS    : filename
#@DESCRIPTION: Returns and optionally sets (if optional argument is present),
#              the interpolation method of the volume.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub interpolation( $$ )
{
    my( $self )          = arg_object( shift, $this_class );
    my( $interpolation ) = opt_arg_enum( shift, N_interpolation_enums );
    end_args( @_ );

    if( defined($interpolation) )
        { $self->{INTERPOLATION} = $interpolation; }

    return( $self->{INTERPOLATION} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : make_ray_trace_args
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : string
#@DESCRIPTION: Returns the arguments needed for ray_trace to render this
#              object.  At present there is no rendering method for volumes.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub make_ray_trace_args( )
{
    my( $self )          = arg_object( shift, $this_class );
    end_args( @_ );

    return( "" );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_plane_intersection
#@INPUT      : self
#              plane_origin_ref  : ref to array of 3 reals
#              plane_normal_ref  : ref to array of 3 reals
#              output_file
#@OUTPUT     :
#@RETURNS    :
#@DESCRIPTION: Creates the intersection of this object and the specified
#              plane, storing the result in the specified file.
#@METHOD     :
#@GLOBALS    :
#@CALLS      :
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   :
#----------------------------------------------------------------------------

sub  get_plane_intersection( $$$$ )
{
    my( $self )              = arg_object( shift, $this_class );
    my( $plane_origin_ref )  = arg_array_ref( shift, 3 );
    my( $plane_normal_ref )  = arg_array_ref( shift, 3 );
    my( $output_file )       = arg_string( shift );
    end_args( @_ );

    my( @plane_origin, @plane_normal, $args );

    @plane_origin = @$plane_origin_ref;
    @plane_normal = @$plane_normal_ref;

    #--------------------------------------------------------------------------
    # make the slice geometry (quadrilateral) for the intersection
    #--------------------------------------------------------------------------

    $args = sprintf( "%s %s %g %g %g %g %g %g",
                     $self->filename(), $output_file,
                     $plane_origin[0],
                     $plane_origin[1],
                     $plane_origin[2],
                     $plane_normal[0],
                     $plane_normal[1],
                     $plane_normal[2] );

    run_executable( "make_slice", $args );
}

#--------------------------------------------------------------------------

1;
