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
#@NAME       : ILT::IntersectionObject
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Object class for intersecting two objects.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::IntersectionObject;

    use      strict;
    use      vars  qw(@ISA);
    use      UNIVERSAL;

    use      ILT::LayoutInclude;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;

    @ISA = ( "ILT::SceneObject" );

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/SceneObject/IntersectionObject.pm,v 1.7 2011-02-04 16:48:14 alex Exp $';

#--------------------------------------------------------------------------
# name of this class
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::IntersectionObject";

#
#  The intersection object currently intersect 2 objects, not n objects.
#  Initially, it only supports intersection of a SliceObject with a SceneObject.
#

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : self
#              object1
#              object2
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: 
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $$$ )
{
    my( $proto )         = arg_any( shift );
    my( $first_object )  = arg_object( shift, "ILT::SceneObject" );
    my( $second_object ) = arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = $class->SUPER::new();
    my( $swap );

    bless ($self, $class);

    if( !UNIVERSAL::isa($first_object,"ILT::PlaneObject" ) &&
        !UNIVERSAL::isa($second_object,"ILT::PlaneObject" ) )
    {
        fatal_error( "ILT::IntersectionObject->new expects one sub-object to " .
                     " be a plane\n" );
    }

    $self->sub_object( 0, $first_object ); 
    $self->sub_object( 1, $second_object ); 

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : sub_object
#@INPUT      : self
#              index
#              sub_object  (OPTIONAL)
#@OUTPUT     : 
#@RETURNS    : the sub object
#@DESCRIPTION: Gets and possibly sets the index'th sub-object.  If the
#              optional argument is present, sets the sub-object.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub sub_object( $$@ )
{
    my( $self )        =  arg_object( shift, $this_class );
    my( $index )       =  arg_int( shift, 0, 1 );
    my( $sub_object )  =  opt_arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    if( defined( $sub_object ) )
        { $self->{SUB_OBJECTS}[$index] = $sub_object; }

    return( $self->{SUB_OBJECTS}[$index] );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : create_temp_geometry_file
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: Creates the temporary .obj file corresponding to the intersection
#              of the two sub-objects.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  create_temp_geometry_file( $ )
{
    my( $self )        =  arg_object( shift, $this_class );
    end_args( @_ );

    my( $tmp_file, @plane_origin, @plane_normal, $plane_index );

    if( defined( $self->{TEMP_GEOMETRY_FILE} ) )
        { return; }

    #--------------------------------------------------------------------------
    # create a temporary .obj filename
    #--------------------------------------------------------------------------

    $tmp_file = get_tmp_file( "obj" );

    #--------------------------------------------------------------------------
    # determine which of the 2 sub objects is the PlaneObject
    #--------------------------------------------------------------------------

    if( UNIVERSAL::isa( $self->sub_object(0), "ILT::PlaneObject" ) )
        { $plane_index = 0; }
    elsif( UNIVERSAL::isa( $self->sub_object(1), "ILT::PlaneObject" ) )
        { $plane_index = 1; }
    else
    {
        fatal_error( "ILT::IntersectionObject expects one sub-object to be " .
                     " a plane\n" );
    }

    #--------------------------------------------------------------------------
    # Intersect the two sub-objects, placing result in the file: tmp_file
    #--------------------------------------------------------------------------

    @plane_origin = $self->sub_object($plane_index)->plane_origin();
    @plane_normal = $self->sub_object($plane_index)->plane_normal();

    $self->sub_object(1-$plane_index)->get_plane_intersection(
                   \@plane_origin, \@plane_normal, $tmp_file );

    #--------------------------------------------------------------------------
    # Record that the temporary geometry file has been created
    #--------------------------------------------------------------------------

    $self->{TEMP_GEOMETRY_FILE} = $tmp_file;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : compute_bounding_view
#@INPUT      : self
#              view_direction_ref    : ref to array of 3 reals
#              up_direction_ref      : ref to array of 3 reals
#              transform             : filename of transform file or ""
#@OUTPUT     : 
#@RETURNS    : array of 6 reals
#@DESCRIPTION: Computes the bounding view of the object, given the view
#              orientation and transform
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  compute_bounding_view( $$$$ )
{
    my( $self )                =  arg_object( shift, $this_class );
    my( $view_direction_ref )  =  arg_array_ref( shift, 3 );
    my( $up_direction_ref )    =  arg_array_ref( shift, 3 );
    my( $transform )           =  arg_string( shift );
    end_args( @_ );

    #--------------------------------------------------------------------------
    # check that the geometry file has been created
    #--------------------------------------------------------------------------

    $self->create_temp_geometry_file();

    #--------------------------------------------------------------------------
    # compute the bounding box of the view
    #--------------------------------------------------------------------------

    return( compute_geometry_file_bounding_view( $self->{TEMP_GEOMETRY_FILE},
                              $view_direction_ref, $up_direction_ref,
                              $transform ) );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : make_ray_trace_args
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : string
#@DESCRIPTION: Creates a string to be used as an argument to ray_trace.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  make_ray_trace_args( $ )
{
    my( $self )  =  arg_object( shift, $this_class );
    end_args( @_ );

    my( $filename );

    #--------------------------------------------------------------------------
    # check that the geometry file has been created
    #--------------------------------------------------------------------------

    $self->create_temp_geometry_file();

    #--------------------------------------------------------------------------
    # now return the temporary filename as the argument
    #--------------------------------------------------------------------------

    $filename = $self->{TEMP_GEOMETRY_FILE};

    return( " $filename " );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_plane_intersection
#@INPUT      : self
#              plane_origin_ref   : ref to array of 3 reals
#              plane_normal_ref   : ref to array of 3 reals
#              output_file        : where to place the .obj
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: Computes the intersection of a plane and this object.
#              At present, intersection of a plane and an intersection object
#              is not implemented, but presumably it is equivalent to the
#              intersection of the intersections of the plane with each of
#              the two sub-objects.
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

    fatal_error( "IntersectionObject->get_plane_intersection() " .
                 " not implemented yet" );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_text_image_magick_args
#@INPUT      : self
#              viewport_x_size
#              viewport_y_size
#@OUTPUT     : 
#@RETURNS    : text arguments for rendering text to image
#@DESCRIPTION: 
#@METHOD     :
#@GLOBALS    :
#@CALLS      :  
#@CREATED    : Jun. 23, 1998    David MacDonald
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

#--------------------------------------------------------------------------

1;
