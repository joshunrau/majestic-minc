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
#@NAME       : ILT::GeometricObject
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: An object class to define a geometric object defined by
#              a .obj filename
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::GeometricObject;

    use      strict;
    use      vars qw( @ISA );
    use      ILT::LayoutInclude;
    use      ILT::Executables;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;
    @ISA = ( "ILT::SceneObject" );

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/SceneObject/GeometricObject.pm,v 1.7 1998-09-18 13:30:00 david Exp $';

#--------------------------------------------------------------------------
# define the name of this class
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::GeometricObject";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : prototype
#              filename          : of .obj file
#@OUTPUT     : 
#@RETURNS    : instance of GeometricObject
#@DESCRIPTION: Creates a geometric object, by storing a filename.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $$ )
{
    my( $proto )    = arg_any( shift );
    my( $filename ) = arg_string( shift );
    end_args();

    my $class = ref($proto) || $proto;
    my $self  = $class->SUPER::new();

    bless ($self, $class);

    $self->filename( $filename );

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : filename
#@INPUT      : self
#              filename   (OPTIONAL
#@OUTPUT     : 
#@RETURNS    : filename
#@DESCRIPTION: Gets and possibly sets (if optional argument present) the
#              filename of the object.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub filename( $@ )
{
    my( $self )     = arg_object( shift, $this_class );
    my( $filename ) = opt_arg_string( shift );
    end_args();

    if( defined($filename) )
        { $self->{FILENAME} = $filename; }

    return( $self->{FILENAME} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : make_ray_trace_args
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Returns a string containing the necessary ray trace arguments
#              for this object, which is simply the filename.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub make_ray_trace_args( $ )
{
    my( $self )     = arg_object( shift, $this_class );
    end_args();

    my( $filename );

    $filename = $self->{FILENAME};

    return( " $filename " );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : compute_bounding_view
#@INPUT      : self
#              view_direction_ref  : ref to array of 3
#              up_direction_ref    : ref to array of 3
#              transform_ref       : filename of a transform file
#@OUTPUT     : 
#@RETURNS    : array of 6 reals
#@DESCRIPTION: Computes the bounding box of the object, given a view
#              direction.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub compute_bounding_view( $$$$ )
{
    my( $self )                =  arg_object( shift, $this_class );
    my( $view_direction_ref )  =  arg_array_ref( shift, 3 );
    my( $up_direction_ref )    =  arg_array_ref( shift, 3 );
    my( $transform )           =  arg_string( shift );
    end_args( @_ );

    return( compute_geometry_file_bounding_view( $self->{FILENAME},
                                                 $view_direction_ref,
                                                 $up_direction_ref,
                                                 $transform ) );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_plane_intersection
#@INPUT      : self
#              plane_origin_ref  : ref to array of 3
#              plane_normal_ref  : ref to array of 3
#              output_file       : filename where to place the result
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: Creates a geometry file (.obj) containing the intersection
#              of the object with the specified plane.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  get_plane_intersection( $$$$ )
{
    my( $self )              =  arg_object( shift, $this_class );
    my( $plane_origin_ref )  =  arg_array_ref( shift, 3 );
    my( $plane_normal_ref )  =  arg_array_ref( shift, 3 );
    my( $output_file )       =  arg_string( shift );
    end_args( @_ );

    my( @plane_origin, @plane_normal, $command );

    @plane_origin = @$plane_origin_ref;
    @plane_normal = @$plane_normal_ref;

#--- the program plane_polygon_intersect needs to be rewritten to be
#--- plane_object_intersect

    $command = sprintf( "%s %s %g %g %g %g %g %g",
                        $self->filename(), $output_file,
                        $plane_normal[0], $plane_normal[1], $plane_normal[2],
                        $plane_origin[0], $plane_origin[1], $plane_origin[2] );

    run_executable( "plane_polygon_intersect", $command );
}

1;
