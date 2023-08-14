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
#@NAME       : ILT::UnionObject
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Object class to implement unions of objects.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::UnionObject;

    use      strict;
    use      vars  qw(@ISA);
    use      ILT::LayoutInclude;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;
    use      UNIVERSAL;
    @ISA =   ( "ILT::SceneObject" );

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/SceneObject/UnionObject.pm,v 1.7 2011-02-04 16:48:14 alex Exp $';

#--------------------------------------------------------------------------
# define the name of this class
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::UnionObject";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : prototype
#              sub_object1
#              sub_object2
#                 ...
#@OUTPUT     : 
#@RETURNS    : instance of union
#@DESCRIPTION: Creates a union object from a list of 1 or more objects
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $$@ )
{
    my( $proto )       = arg_any( shift );
    my( @sub_objects ) = arg_array( \@_ );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = $class->SUPER::new();
    my $object;

    bless ($self, $class);

    #--------------------------------------------------------------------------
    # check if all arguments are scene objects
    #--------------------------------------------------------------------------

    foreach $object ( @sub_objects )
    {
        if( !UNIVERSAL::isa( $object, "ILT::SceneObject" ) )
        {
            fatal_error( "Non SceneObject passed into ${this_class}::new\n" );
        }
    }

    $self->{SUB_OBJECTS} = [ @sub_objects ];

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : n_sub_objects
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : number of objects 
#@DESCRIPTION: Returns the number of objects in the union
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub n_sub_objects( $ )
{
    my( $self ) = arg_object( shift, $this_class );
    end_args( @_ );

    return( scalar( @{$self->{SUB_OBJECTS}} ) );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : sub_object
#@INPUT      : self
#              index
#              sub_object  (OPTIONAL)
#@OUTPUT     : 
#@RETURNS    : sub object
#@DESCRIPTION: Gets and optionally sets (if optional argument is present) the
#              sub object associated with the given index.
#              Setting an object index, may increase the number of objects
#              in the union, but there is no way to decrease it.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub sub_object( $$@ )
{
    my( $self )       = arg_object( shift, $this_class );
    my( $index )      = arg_int( shift, 0 );
    my( $sub_object ) = opt_arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    if( defined($sub_object) )
        { $self->{SUB_OBJECTS}[$index] = $sub_object; }

    return( $self->{SUB_OBJECTS}[$index] );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : make_ray_trace_args
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : string
#@DESCRIPTION: Creates the ray_trace args for the union.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub make_ray_trace_args( $ )
{
    my( $self )       = arg_object( shift, $this_class );
    end_args( @_ );

    my( $args, $sub_args, $sub_object );

    #--------------------------------------------------------------------------
    # simply concatenate the arguments of all the sub objects
    #--------------------------------------------------------------------------

    $args = "";

    foreach $sub_object ( @{$self->{SUB_OBJECTS}} )
    {
        $sub_args = $sub_object->make_ray_trace_args();
        $args = $args . " " . $sub_args;
    }

    return( $args );
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

    my( $tmp_file, $sub_object );

    $tmp_file = get_tmp_file( "obj" );

    #--------------------------------------------------------------------------
    # initialize the output file to empty
    #--------------------------------------------------------------------------

    unlink( $output_file );
    run_executable( "touch", $output_file );

    foreach $sub_object ( @{$self->{SUB_OBJECTS}} )
    {
        #----------------------------------------------------------------------
        # create the plane intersection of one object
        #----------------------------------------------------------------------

        $sub_object->get_plane_intersection( $plane_origin_ref,
                                             $plane_normal_ref, $tmp_file );

        #----------------------------------------------------------------------
        # add the plane intersection to the output file
        #----------------------------------------------------------------------

        run_executable( "cat", "$tmp_file >> $output_file" );
    }
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : compute_bounding_view
#@INPUT      : self
#              view_direction_ref  : ref to array of 3 reals
#              up_direction_ref    : ref to array of 3 reals
#              transform
#@OUTPUT     : 
#@RETURNS    : 6 reals
#@DESCRIPTION: Computes the bounding view of the object, given the view
#              direction.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub compute_bounding_view( $$$$ )
{
    my( $self )                = arg_object( shift, $this_class );
    my( $view_direction_ref )  =  arg_array_ref( shift, 3 );
    my( $up_direction_ref )    =  arg_array_ref( shift, 3 );
    my( $transform )           =  arg_string( shift );
    end_args( @_ );

    my( $x_min, $x_max, $y_min, $y_max, $z_min, $z_max, $sub_object,
        $sub_x_min, $sub_x_max, $sub_y_min, $sub_y_max,
        $sub_z_min, $sub_z_max, $first );

    $first = 1;

    #--------------------------------------------------------------------------
    # take the minimum and maximums of the bounding view limits of the
    # sub_objects
    #--------------------------------------------------------------------------

    foreach $sub_object ( @{$self->{SUB_OBJECTS}} )
    {
        ( $sub_x_min, $sub_x_max, $sub_y_min, $sub_y_max,
          $sub_z_min, $sub_z_max ) = 
                    $sub_object->compute_bounding_view( $view_direction_ref,
                                                        $up_direction_ref,
                                                        $transform );

        if( $sub_x_min > $sub_x_max || $sub_y_min > $sub_y_max ||
            $sub_z_min > $sub_z_max )
            { next; }

        if( $first )
        {
            ($x_min, $x_max, $y_min, $y_max, $z_min, $z_max) =
            ($sub_x_min, $sub_x_max, $sub_y_min, $sub_y_max,
             $sub_z_min, $sub_z_max );

            $first = 0;
        }
        else
        {
            if( $sub_x_min < $x_min )
                { $x_min = $sub_x_min; }
            if( $sub_y_min < $y_min )
                { $y_min = $sub_y_min; }
            if( $sub_z_min < $z_min )
                { $z_min = $sub_z_min; }

            if( $sub_x_max > $x_max )
                { $x_max = $sub_x_max; }
            if( $sub_y_max > $y_max )
                { $y_max = $sub_y_max; }
            if( $sub_z_max > $z_max )
                { $z_max = $sub_z_max; }
        }
    }

    return( $x_min, $x_max, $y_min, $y_max, $z_min, $z_max );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : create_temp_geometry_file
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Creates the temporary geometry file for this object
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  create_temp_geometry_file( $ )
{
    my( $self )                = arg_object( shift, $this_class );
    end_args( @_ );

    my( $sub_object );

    #--------------------------------------------------------------------------
    # simply ask each object to make its own temporary geometry file
    #--------------------------------------------------------------------------

    foreach $sub_object ( @{$self->{SUB_OBJECTS}} )
    {
        $sub_object->create_temp_geometry_file();
    }
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

    my( $args, $sub_object, $sub_text );

    $args = "";

    foreach $sub_object ( @{$self->{SUB_OBJECTS}} )
    {
        $sub_text = $sub_object->get_text_image_magick_args(
                                      $viewport_x_size, $viewport_y_size );

        if( defined($sub_text) && $sub_text ne "" )
        {
             if( $args ne "" )
                 { $args = $args . " "; }

             $args = $args . $sub_text;
        }
    }

    return( $args );
}

#--------------------------------------------------------------------------

1;
