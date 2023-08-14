#!/usr/local/bin/perl5 -w

#-------------------------------------------------------
#
#  This script configures an image of M volumes by N slices
#  using the ImageLayout object and related objects.  As such,
#  it is simply a demonstration of the interface to the ImageLayout object
#
#  Author: David MacDonald        March, 1998
#
#--------------------------------------------------------
    use strict;

    use ILT::LayoutInclude;
    use ILT::LayoutUtils;
#    use ILT::ProgUtils;

    my( $layout, $n_rows, $n_cols, @volume_filenames, @slice_pos,
        $scene_object, $view, $image_info, $row, $col, @obj_files, @view_dirs,
        @view, @plane_object, @surface_object, @volume_object,
        $colour_object, $surface_cross_section, @colour_codes,
        $bounding_object, $clipped_view, $extra_space_around_bounding_object,
        $object_defining_view, $header );

#-------------------------------------------------------
#
#  Define the information to be displayed
#
#  This could be done in any number of ways, so the simplest is shown here
#  for demonstration purposes
#
#--------------------------------------------------------

#-----------------------------------------------------------------
#   Row information
#-----------------------------------------------------------------

    @volume_filenames = ( "Data/volume1.mnc.gz",
                          "Data/volume1.mnc.gz",
                          "Data/volume2.mnc.gz",
                          "Data/volume3.mnc.gz" );

    @colour_codes = ( [ Spectral_scale, 50, 150 ],
                      [ Hot_metal_scale, 50, 150 ],
                      [ Usercc_scale, 0, 3, "Data/class.map" ],
                      [ Hot_metal_scale, 100000, 500000 ] );

    @obj_files = ( "Data/surf1.obj",
                   "Data/surf2.obj",
                   "Data/surf2.obj",
                   "Data/surf3.obj" );

#-----------------------------------------------------------------
#   Column information
#-----------------------------------------------------------------

    @slice_pos = ( [ Sagittal_axis,   -50 ],
                   [ Sagittal_axis,   -30 ],
                   [ Coronal_axis,    -24 ],
                   [ Coronal_axis,     24 ],
                   [ Transverse_axis,   0 ],
                   [ Transverse_axis,  50 ] );

    @view_dirs = ( Right_view,
                   Right_view,
                   Back_view,
                   Back_view,
                   Top_view,
                   Top_view );

#-----------------------------------------------------------------
#   Create an image layout that is a grid of $n_rows by $n_cols
#-----------------------------------------------------------------

    $n_rows = @volume_filenames;
    $n_cols = @view_dirs;

    $layout = ILT::ImageLayout->new_grid( $n_rows, $n_cols );
    $layout->horizontal_white_space( 3 );
    $layout->vertical_white_space( 3 );
    $layout->white_space_colour( "black" );

#-----------------------------------------------------------------
#   Since clipping is not yet implemented, the following bounding
#   object specificiation is a temporary hack to allow specifying
#   a narrower region of interest than the whole volume
#
#   Actually, it's not really that much of a hack.  Rather than
#   letting the view define itself based on the scene to be rendered
#   it is letting the view defined itself based on some unrelated
#   scene object which represents the region of interest.  This seems
#   like a reasonably valid thing to do.
#-----------------------------------------------------------------

    $bounding_object = ILT::GeometricObject->new( "Data/surf2.obj" );
    $extra_space_around_bounding_object = 5;

#-----------------------------------------------------------------
#   Create a slice and view for each column
#-----------------------------------------------------------------

    for( $col = 0;  $col < $n_cols;  ++$col )
    {
        $plane_object[$col] = ILT::PlaneObject->new_canonical(
                                        $slice_pos[$col][0],
                                        $slice_pos[$col][1] );
        $view[$col] = ILT::View->new_canonical( $view_dirs[$col] );

#   Temporary modification of view to clip it, this will disappear when
#   clipping is implemented.  You may comment out this line if not desired

        $object_defining_view = ILT::IntersectionObject->new(
                                       $plane_object[$col],
                                       $bounding_object );

        set_clipped_view( $view[$col],
                          $object_defining_view,
                          $extra_space_around_bounding_object );

    }

#-----------------------------------------------------------------
#   Create a volume and surface object for each row
#-----------------------------------------------------------------

    for( $row = 0;  $row < $n_rows;  ++$row )
    {
        $volume_object[$row] = ILT::VolumeObject->new(
                                       $volume_filenames[$row] );
        $surface_object[$row] = ILT::GeometricObject->new( $obj_files[$row] );
    }

#-----------------------------------------------------------------
#   Now describe each image in the grid
#-----------------------------------------------------------------

    for( $row = 0;  $row < $n_rows;  ++$row )
    {
        for( $col = 0;  $col < $n_cols;  ++$col )
        {
            #-----------------------------------------------------------------
            #   Create the scene hierarchy for one image as a union of:
            #       a colour-coded intersection of a slice and volume, and
            #       a cross section of a surface (intersection of slice and
            #                                     surface)
            #-----------------------------------------------------------------

            $colour_object =
                  ILT::ColourObject->new_volume_colouring( 
                       ILT::IntersectionObject->new( $plane_object[$col],
                                                $volume_object[$row] ),
                       $volume_object[$row], $colour_codes[$row][0],
                       $colour_codes[$row][1], $colour_codes[$row][2] );

            $colour_object->usercc_filename( $colour_codes[$row][3] );

            $surface_cross_section = ILT::IntersectionObject->new(
                                           $plane_object[$col],
                                           $surface_object[$row] );

            #-----------------------------------------------------------------
            #   Create the image info for the scene,
            #    which is essentially the scene object and the view
            #-----------------------------------------------------------------

            $scene_object =
               ILT::UnionObject->new( $colour_object, $surface_cross_section );

            $image_info = ILT::ImageInfo->new( $scene_object, $view[$col] );
            $image_info->background_colour( "black" );

            #-----------------------------------------------------------------
            #   Attach the image info to the layout manager
            #-----------------------------------------------------------------

            $layout->image_info( $layout->row_col_to_index($row,$col),
                                 $image_info );
        }
    }

    #-----------------------------------------------------------------
    #   Create a header object
    #-----------------------------------------------------------------

    $header = ILT::TextObject->new( "Output of test_layout.pl",
                                    .5, .5 );
    $header->font( "-adobe-helvetica-medium-r-normal--34-240-100-100-p-176-iso8859-1" );
    $header->colour( "yellow" );
    $header->horizontal_alignment( Align_centre );
    $header->vertical_alignment( Align_centre );
    $layout->header( $header );

    #-----------------------------------------------------------------
    #   Render the images and create the output file
    #-----------------------------------------------------------------

    $layout->generate_image( "test_layout.rgb", 1100, 0 );

#-----------------------------------------------------------------
#   Define the view based on some scene object other than that being
#   rendered.
#-----------------------------------------------------------------

sub  set_clipped_view( $$$ )
{
    my( $view )                               = shift;
    my( $object_defining_view )               = shift;
    my( $extra_space_around_bounding_object ) = shift;

    my( $view_copy, @bbox );

#--- make a temporary copy of the view

    $view_copy = $view->copy();

#--- define the default view based on the object

    $view_copy->compute_view_for_object( $object_defining_view );

#--- explicitly set bounding box of the view, so that it doesn't get
#--- automatically computed

    @bbox = $view_copy->bounding_box();

    $bbox[0] -= $extra_space_around_bounding_object;
    $bbox[1] += $extra_space_around_bounding_object;
    $bbox[2] -= $extra_space_around_bounding_object;
    $bbox[3] += $extra_space_around_bounding_object;
    $bbox[4] -= $extra_space_around_bounding_object;
    $bbox[5] += $extra_space_around_bounding_object;

    $view->bounding_box( @bbox );
}
