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
#@NAME       : ILT::TransformObject
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: An object class to implement transforms.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::TransformObject;

    use      strict;
    use      vars  qw(@ISA);
    use      ILT::LayoutInclude;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;
    use      ILT::Executables;
    use      ILT::SceneObject::OneSubObject;
    @ISA =   ( "ILT::OneSubObject" );

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/SceneObject/TransformObject.pm,v 1.6 2011-02-04 16:48:14 alex Exp $';

#--------------------------------------------------------------------------
# define the name of this class
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::TransformObject";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : prototype
#              sub-object
#              transform       : filename
#@OUTPUT     : 
#@RETURNS    : instance of transform object
#@DESCRIPTION: Constructs a transform object consisting of a sub-object and
#              a filename defining the transform to apply to the sub-object.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $$$ )
{
    my( $proto )              = arg_any( shift );
    my( $sub_object )         = arg_object( shift, "ILT::SceneObject" );
    my( $transform_filename ) = arg_string( shift );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = $class->SUPER::new( $sub_object );

    bless ($self, $class);

    $self->transform( $transform_filename );

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : transform
#@INPUT      : self
#              transform  (OPTIONAL)
#@OUTPUT     : 
#@RETURNS    : transform filename
#@DESCRIPTION: Gets and optionally sets (if optional argument is present) the
#              value of this attribute.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub transform( $@ )
{
    my( $self )       = arg_object( shift, $this_class );
    my( $transform )  = opt_arg_string( shift );
    end_args( @_ );

    if( defined($transform) )
        { $self->{TRANSFORM_FILENAME} = $transform; }

    return( $self->{TRANSFORM_FILENAME} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : make_ray_trace_args
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : string
#@DESCRIPTION: Creates a string to be used as argument to ray_trace
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub make_ray_trace_args( $ )
{
    my( $self )        = arg_object( shift, $this_class );
    end_args( @_ );

    my( $args );

    #--------------------------------------------------------------------------
    # prefix a -transform, and append a reset-to-identity transform to 
    # the arguments of the sub object
    #--------------------------------------------------------------------------

    $args = sprintf( "-transform %s %s -transform identity",
                     $self->transform(),
                     $self->SUPER::make_ray_trace_args() );

    return( $args );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_plane_intersection
#@INPUT      : self
#              plane_origin_ref   : ref to array of 3 reals
#              plane_normal_ref   : ref to array of 3 reals
#              output_file   
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Creates the intersection of the given plane and this object,
#              storing the result in the file specified by output_file.
#              Until the function inverse_transform_plane is implemented
#              this function will crap out.
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

    my( $tmp_file, $transform, $trans_origin, $trans_normal );

    $transform = $self->transform();

    #--------------------------------------------------------------------------
    # inverse transform the plane,
    #--------------------------------------------------------------------------

    fatal_error( "inverse_transform_plane() not yet implemented\n" );

    ( $trans_origin, $trans_normal ) = inverse_transform_plane( $transform,
                                                        $plane_origin_ref,
                                                        $plane_normal_ref );

    $tmp_file = get_tmp_file( "obj" );

    #--------------------------------------------------------------------------
    # intersect the inverse transformed plane with the non-transformed
    # sub-object
    #--------------------------------------------------------------------------

    $self->SUPER::get_plane_intersection( $trans_origin,
                                          $trans_normal, $tmp_file );

    #--------------------------------------------------------------------------
    # then transform the resulting geometry by the transform
    #--------------------------------------------------------------------------

    run_executable( "transform_objects", "$tmp_file $transform $output_file" );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : compute_bounding_view
#@INPUT      : self
#              view_direction_ref   : ref to array of 3 reals
#              up_direction_ref     : ref to array of 3 reals
#              transform   
#@OUTPUT     : 
#@RETURNS    : 6 reals
#@DESCRIPTION: Computes the bounding view of the transform object
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

    my( $x_min, $x_max, $y_min, $y_max, $z_min, $z_max,
        $tmp_transform, $command );

    #--------------------------------------------------------------------------
    # if a transform is passed in, the effective transform is the
    # concatenation of the passed transform, and the transform of this object
    #--------------------------------------------------------------------------

    if( defined( $transform ) && $transform ne "" )
    {
        $tmp_transform = get_tmp_file( "xfm" );
        $command = sprintf( "-clobber %s %s %s",
                            $transform, $self->transform(), $tmp_transform );
        run_executable( "xfmconcat", $command );
    }
    else
    {
        $tmp_transform = $self->transform();
    }

    #--------------------------------------------------------------------------
    # compute the bounding box of the view, given the computed transform
    #--------------------------------------------------------------------------

    ( $x_min, $x_max, $y_min, $y_max, $z_min, $z_max ) = 
           $self->SUPER::compute_bounding_view( $view_direction_ref,
                                                       $up_direction_ref,
                                                       $tmp_transform );

    return( $x_min, $x_max, $y_min, $y_max, $z_min, $z_max );
}

#--------------------------------------------------------------------------

1;
