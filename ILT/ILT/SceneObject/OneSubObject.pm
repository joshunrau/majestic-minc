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
#@NAME       : ILT::OneSubObject
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Object class to be used as a base class for objects that
#              are modifiers of a sub-object, such as RenderObject, and
#              ColourObject
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::OneSubObject;

    use      strict;
    use      vars  qw(@ISA);
    use      ILT::LayoutInclude;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;
    @ISA =   ( "ILT::SceneObject" );

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/SceneObject/OneSubObject.pm,v 1.5 2011-02-04 16:48:14 alex Exp $';

#--------------------------------------------------------------------------
# the name of this class
#--------------------------------------------------------------------------

my( $this_class ) = "ILT::OneSubObject";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : prototype
#              sub-object
#@OUTPUT     : 
#@RETURNS    : instance of this class
#@DESCRIPTION: Constructs an object of this class, which references a
#              subobject 
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new
{
    my( $proto )      = arg_any( shift );
    my( $sub_object ) = arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = $class->SUPER::new();

    bless ($self, $class);

    $self->sub_object( $sub_object );

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : sub_object
#@INPUT      : self
#              sub_object (OPTIONAL)
#@OUTPUT     : 
#@RETURNS    : sub_object
#@DESCRIPTION: Gets and possibly sets the sub_object of this object.  If
#              the optional argument is present, sets the sub_object.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub sub_object( $@ )
{
    my( $self )       = arg_object( shift, $this_class );
    my( $sub_object ) = opt_arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    if( defined($sub_object) )
        { $self->{SUB_OBJECT} = $sub_object; }

    return( $self->{SUB_OBJECT} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : make_ray_trace_args
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : string
#@DESCRIPTION: Returns a string containing the arguments for ray_trace to
#              render this object.
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

    my( $args );

    #--------------------------------------------------------------------------
    # just return the string for the sub object
    #--------------------------------------------------------------------------

    $args = $self->sub_object()->make_ray_trace_args();

    return( $args );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : get_plane_intersection
#@INPUT      : self
#              plane_origin_ref   : ref to 3 reals
#              plane_normal_ref   : ref to 3 reals
#              output_file
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: Intersects the given plane with this object, placing the
#              result in the specified file.
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

    #--------------------------------------------------------------------------
    # simply pass on the request to the sub object
    #--------------------------------------------------------------------------

    $self->sub_object()->get_plane_intersection( $plane_origin_ref,
                                                 $plane_normal_ref,
                                                 $output_file );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : compute_bounding_view
#@INPUT      : self
#              view_direction_ref   : ref to 3 reals
#              up_direction_ref     : ref to 3 reals
#              transform
#@OUTPUT     : 
#@RETURNS    : 6 reals
#@DESCRIPTION: Computes the bounding box of the object based on the view
#              orientation and transform specified.
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

    #--------------------------------------------------------------------------
    # simply pass on the request to the sub object
    #--------------------------------------------------------------------------

    return( $self->sub_object()->compute_bounding_view(
                                            $view_direction_ref,
                                            $up_direction_ref, $transform ) );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : create_temp_geometry_file
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Creates the temporary geometry file for this object.
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

    #--------------------------------------------------------------------------
    # simply pass on the request to the sub object
    #--------------------------------------------------------------------------

    $self->sub_object()->create_temp_geometry_file();
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

    return( $self->sub_object()->get_text_image_magick_args(
                       $viewport_x_size, $viewport_y_size ) );
}

#--------------------------------------------------------------------------

1;
