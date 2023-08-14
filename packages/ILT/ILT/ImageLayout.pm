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
#@NAME       : ImageLayout.pm
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Object class for creating image files which contain several
#              sub-images which are renderings of geometry and data
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

    package  ILT::ImageLayout;

    use      strict;
    use      ILT::LayoutUtils;
    use      ILT::ProgUtils;
    use      ILT::Executables;
    use      ILT::LayoutInclude;

    my( $rcsid ) = '$Header: /private-cvsroot/libraries/ILT/ILT/ImageLayout.pm,v 1.13 2011-02-04 16:48:13 alex Exp $';

#--- Define the name of this class

my( $this_class ) = "ILT::ImageLayout";

#----------------------------- MNI Header -----------------------------------
#@NAME       : new
#@INPUT      : prototype
#@OUTPUT     : 
#@RETURNS    : an instance of ImageLayout class
#@DESCRIPTION: Creates an instance of the ImageLayout class.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new( $ )
{
    my( $proto )  = arg_any( shift );
    end_args( @_ );

    my $class = ref($proto) || $proto;
    my $self  = {};

    bless ($self, $class);

    $self->n_images( 0 );
    $self->{N_ROWS} = 0;
    $self->{N_COLS} = 0;
    $self->{IMAGES} = [];
    $self->{HEADER} = new ILT::TextObject( "", 0, 0 );
    $self->{FOOTER} = new ILT::TextObject( "", 0, 0 );
    $self->horizontal_white_space( 0 );
    $self->vertical_white_space( 0 );
    $self->white_space_colour( "black" );

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : new_grid
#@INPUT      : prototype
#              n_rows
#              n_cols
#@OUTPUT     : 
#@RETURNS    : an instance of an ImageLayout
#@DESCRIPTION: Creates an ImageLayout object with a grid oriented layout.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub new_grid( $$$ )
{
    my( $proto )  = arg_any( shift );
    my( $n_rows ) = arg_int( shift, 1, 1e30 );
    my( $n_cols ) = arg_int( shift, 1, 1e30 );
    end_args( @_ );

    my( $self );

    $self = new ILT::ImageLayout;

    $self->{N_ROWS} = $n_rows;
    $self->{N_COLS} = $n_cols;
    $self->{N_IMAGES} = $n_rows * $n_cols;

    return $self;
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : n_images
#@INPUT      : self
#              n_images   (optional)
#@OUTPUT     : 
#@RETURNS    : current number of images
#@DESCRIPTION: Returns the number of images.  With the optional argument, sets
#              the number of images.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub n_images( $@ )
{
    my( $self )     = arg_object( shift, $this_class );
    my( $n_images ) = opt_arg_int( shift, 1, 1e30 );
    end_args( @_ );

    if( defined( $n_images ) )
        { $self->{N_IMAGES} = $n_images; }

    return( $self->{N_IMAGES} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : n_rows
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : n rows
#@DESCRIPTION: Returns the number of rows.  Note that you cannot set the number
#              of rows with this routine;  it can only be set on creation by
#              new_grid().
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub n_rows( $ )
{
    my( $self ) = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->{N_ROWS} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : n_cols
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : n cols
#@DESCRIPTION: Returns the number of columns.  Note that you cannot set the
#              number of columns with this routine;  it can only be set on
#              creation by new_grid().
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub n_cols( $ )
{
    my( $self ) = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->{N_COLS} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : horizontal_white_space
#@INPUT      : self
#              horizontal_white_space   (optional)
#@OUTPUT     : 
#@RETURNS    : current value
#@DESCRIPTION: Returns the current value of this parameter.  With the optional
#              parameter present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub horizontal_white_space( $@ )
{
    my( $self )     = arg_object( shift, $this_class );
    my( $horizontal_white_space ) = opt_arg_real( shift, 0.0, 1e30 );
    end_args( @_ );

    if( defined( $horizontal_white_space ) )
        { $self->{HORIZONTAL_WHITE_SPACE} = $horizontal_white_space; }

    return( $self->{HORIZONTAL_WHITE_SPACE} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : vertical_white_space
#@INPUT      : self
#              vertical_white_space   (optional)
#@OUTPUT     : 
#@RETURNS    : current value
#@DESCRIPTION: Returns the current value of this parameter.  With the optional
#              parameter present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub vertical_white_space( $@ )
{
    my( $self )     = arg_object( shift, $this_class );
    my( $vertical_white_space ) = opt_arg_real( shift, 0.0, 1e30 );
    end_args( @_ );

    if( defined( $vertical_white_space ) )
        { $self->{VERTICAL_WHITE_SPACE} = $vertical_white_space; }

    return( $self->{VERTICAL_WHITE_SPACE} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : white_space_colour
#@INPUT      : self
#              white_space_colour   (optional)
#@OUTPUT     : 
#@RETURNS    : current value  (a string consisting of the named of the colour)
#@DESCRIPTION: Returns the current value of this parameter.  With the optional
#              parameter present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub white_space_colour( $@ )
{
    my( $self )     = arg_object( shift, $this_class );
    my( $white_space_colour ) = opt_arg_string( shift );
    end_args( @_ );

    if( defined( $white_space_colour ) )
        { $self->{WHITE_SPACE_COLOUR} = $white_space_colour; }

    return( $self->{WHITE_SPACE_COLOUR} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : header
#@INPUT      : self
#              text_object  (optional)
#@OUTPUT     : 
#@RETURNS    : current value (text_object)
#@DESCRIPTION: Returns the current value of the text header.  With the optional
#              parameter present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Jun. 19, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  header( $@ )
{
    my( $self ) = arg_object( shift, $this_class );
    my( $header_text ) = opt_arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    if( defined( $header_text ) )
        { $self->{HEADER} = $header_text; }

    return( $self->{HEADER} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : footer
#@INPUT      : self
#              text_object  (optional)
#@OUTPUT     : 
#@RETURNS    : current value (text_object)
#@DESCRIPTION: Returns the current value of the text footer.  With the optional
#              parameter present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Jun. 19, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  footer( $@ )
{
    my( $self ) = arg_object( shift, $this_class );
    my( $footer_text ) = opt_arg_object( shift, "ILT::SceneObject" );
    end_args( @_ );

    if( defined( $footer_text ) )
        { $self->{FOOTER} = $footer_text; }

    return( $self->{FOOTER} );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : is_in_grid_mode
#@INPUT      : self
#@OUTPUT     : 
#@RETURNS    : whether in grid mode
#@DESCRIPTION: Returns True if the ImageLayout is in a grid orientation.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  is_in_grid_mode( $ )
{
    my( $self )     = arg_object( shift, $this_class );
    end_args( @_ );

    return( $self->n_rows() > 0 && $self->n_cols() > 0 );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : row_col_to_index
#@INPUT      : self
#              row
#              col
#@OUTPUT     : 
#@RETURNS    : image index
#@DESCRIPTION: Computes the image index based on the row and column index.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub  row_col_to_index( $$$ )
{
    my( $self )    = arg_object( shift, $this_class );
    my( $row )     = arg_int( shift, 0, $self->n_rows()-1 );
    my( $col )     = arg_int( shift, 0, $self->n_cols()-1 );
    end_args( @_ );

    if( ! $self->is_in_grid_mode() )
        { fatal_error( "row_col_to_index(): not in grid mode\n" ); }

    return( $row * $self->{N_COLS} + ($self->{N_COLS}-1-$col) );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : image_info
#@INPUT      : self
#              image_index
#              image_info   (optional)
#@OUTPUT     : 
#@RETURNS    : current value
#@DESCRIPTION: Returns the current image_info object associated with the
#              given image index.  If the optional image_info parameter is
#              present, sets the value.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub image_info( $$@ )
{
    my( $self )         = arg_object( shift, $this_class );
    my( $image_index )  = arg_int( shift, 0, $self->n_images()-1 );
    my( $value_to_set ) = opt_arg_object( shift, "ILT::ImageInfo" );
    end_args( @_ );

    if( defined( $value_to_set ) )
        { $self->{IMAGES}[$image_index] = $value_to_set; }

    return( $self->{IMAGES}[$image_index] );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : compute_image_sizes_and_positions
#@INPUT      : self
#              full_x_size     : final image size
#              full_y_size
#@OUTPUT     : x_pos           : refs to arrays of positions for each sub-image
#              y_pos
#              x_sizes         : refs to arrays of sizes for each sub-image
#              y_sizes
#@RETURNS    : void
#@DESCRIPTION: Assigns the size and position of each sub image.  Assumes all
#              image_info's have been assigned.  Currently, only grid mode
#              is implemented
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub   compute_image_sizes_and_positions( $$$$$$$$$ )
{
    my( $self       )    = arg_object( shift, $this_class );
    my( $full_x_size )   = arg_int( shift, 0, 1e30 );
    my( $full_y_size )   = arg_int( shift, 0, 1e30 );
    my( $header_size )   = arg_int( shift, 0, 1e30 );
    my( $footer_size )   = arg_int( shift, 0, 1e30 );
    my( $x_pos )         = arg_array_ref( shift );
    my( $y_pos )         = arg_array_ref( shift );
    my( $x_sizes )       = arg_array_ref( shift );
    my( $y_sizes )       = arg_array_ref( shift );
    end_args( @_ );

    my( $x_size, $y_size, $row, $col, $n_rows, $n_cols, $image_index,
        $hor_white_space, $vert_white_space,
        $width, $height, @max_world_height, @max_world_width,
        $total_height, $total_width, $x_scale, $y_scale, $scale,
        $x_offset, $y_offset, $current_x, $current_y,
        $used_x_size, $used_y_size );

    if( $self->is_in_grid_mode() )
    {
        $hor_white_space = $self->horizontal_white_space();
        $vert_white_space = $self->vertical_white_space();
        $n_rows = $self->n_rows();
        $n_cols = $self->n_cols();
        if( $header_size <= 0.0 ) {
          $header_size = $vert_white_space;
        }
        if( $footer_size <= 0.0 ) {
          $footer_size = $vert_white_space;
        }

        $used_x_size = $full_x_size - ($n_cols-1) * $hor_white_space;
        $used_y_size = $full_y_size - $header_size - $footer_size -
			($n_rows-1) * $vert_white_space;

        #-----------------------------------------------------------------------
        # compute the maximum world height of an image for each row and the
        # maximum world width for each column
        #-----------------------------------------------------------------------

        for( $row = 0;  $row < $n_rows;  ++$row )
            { $max_world_height[$row] = 0; }

        for( $col = 0;  $col < $n_cols;  ++$col )
            { $max_world_width[$col] = 0; }

        for( $row = 0;  $row < $n_rows;  ++$row )
        {
            for( $col = 0;  $col < $n_cols;  ++$col )
            {
                $image_index = $self->row_col_to_index( $row, $col );

                if( $self->{IMAGES}[$image_index]->scene_view()->
                                                get_bounding_box_width() >
                    $max_world_width[$col] )
                {
                    $max_world_width[$col] =
                           $self->{IMAGES}[$image_index]->scene_view()->
                                               get_bounding_box_width();
                }

                if( $self->{IMAGES}[$image_index]->scene_view()->
                                                 get_bounding_box_height() >
                    $max_world_height[$row] )
                {
                    $max_world_height[$row] =
                           $self->{IMAGES}[$image_index]->scene_view()->
                                                 get_bounding_box_height();
                }
            }
        }

        #----------------------------------------------------------------------
        # sum up the total height across rows and width across columns
        #----------------------------------------------------------------------

        $total_height = 0;
        for( $row = 0;  $row < $n_rows;  ++$row )
            { $total_height += $max_world_height[$row]; }

        $total_width = 0;
        for( $col = 0;  $col < $n_cols;  ++$col )
            { $total_width += $max_world_width[$col]; }

        #----------------------------------------------------------------------
        # now the scale factor that converts world coordinates to pixel
        # coordinates can be computed
        #----------------------------------------------------------------------

        if( $used_x_size < 1 )   #--- no x size specified
        {
            $scale = $used_y_size / $total_height;
            $used_x_size = int( $scale * $total_width + 0.5 );
        }
        elsif( $used_y_size < 1 )
        {
            $scale = $used_x_size / $total_width;
            $used_y_size = int( $scale * $total_height + 0.5 );
        }
        else
        {
            $x_scale = $used_x_size / $total_width;
            $y_scale = $used_y_size / $total_height;

            #-----------------------------------------------------------------
            # due to aspect differences, choose the smaller scale
            #-----------------------------------------------------------------

            if( $x_scale < $y_scale )
                { $scale = $x_scale; }
            else
                { $scale = $y_scale; }
        }

        #----------------------------------------------------------------------
        # Depending on the aspect, one of these offsets will be 0 and the
        # other positive
        #----------------------------------------------------------------------

        $x_offset = ($used_x_size - $scale * $total_width) / 2;
        $y_offset = ($used_y_size - $scale * $total_height) / 2;

        #----------------------------------------------------------------------
        # Now we can step through the rows and columns and assign positions
        # and sizes
        #----------------------------------------------------------------------

        $current_y = $footer_size + $y_offset;
        for( $row = 0;  $row < $n_rows;  ++$row )
        {
            #------------------------------------------------------------------
            # get the pixel height of this row
            #------------------------------------------------------------------

            $height = $scale * $max_world_height[$row];

            $current_x = $x_offset;
            for( $col = 0;  $col < $n_cols;  ++$col )
            {
                #--------------------------------------------------------------
                # get the pixel width of this column
                #--------------------------------------------------------------

                $width = $scale * $max_world_width[$col];

                $image_index = $self->row_col_to_index( $row, $col );

                #--------------------------------------------------------------
                # compute the sizes as the scale factor times the world
                # coordinate size of the bounding box
                #--------------------------------------------------------------

                $$x_sizes[$image_index] = int( $scale *
                             $self->{IMAGES}[$image_index]->scene_view()->
                                       get_bounding_box_width() );

                $$y_sizes[$image_index] = int( $scale *
                             $self->{IMAGES}[$image_index]->scene_view()->
                                 get_bounding_box_height() );


                #--------------------------------------------------------------
                # compute the position, by centering the image within 
                # the space alloted for the row or column.
                #--------------------------------------------------------------

                $$x_pos[$image_index] = int( ($current_x +
                              ($width - $$x_sizes[$image_index]) / 2) + 0.5 );

                $$y_pos[$image_index] = int( ($current_y +
                              ($height - $$y_sizes[$image_index]) / 2) + 0.5 );

                #--------------------------------------------------------------
                # advance to the next column
                #--------------------------------------------------------------

                $current_x += $width + $hor_white_space;
            }

            #--------------------------------------------------------------
            # advance to the next row
            #--------------------------------------------------------------

            $current_y += $height + $vert_white_space;
        }

        $full_x_size = $used_x_size + ($n_cols-1) * $hor_white_space;
        $full_y_size = $used_y_size + ($n_rows-1) * $vert_white_space +
                       $header_size + $footer_size;
    }
    else
    {
        fatal_error( "compute_image_sizes_and_positions():  " .
                     " not implemented yet\n" );
    }

    return( $full_x_size, $full_y_size );
}

#----------------------------- MNI Header -----------------------------------
#@NAME       : generate_image
#@INPUT      : self
#              filename
#              full_x_size    :  one of these is allowed to be 0
#              full_y_size    : 
#@OUTPUT     : 
#@RETURNS    : void
#@DESCRIPTION: Positions, renders, and assembles the montage of sub-images
#              into an image of size full_x_size by full_y_size, storing the
#              result in the file specified by filename.  If either
#              full_x_size or full_y_size is 0, then it will
#              be automatically assigned based on the computed aspect ratio
#              and the other full_*_size specified.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      :  
#@CREATED    : Apr. 16, 1998    David MacDonald
#@MODIFIED   : 
#----------------------------------------------------------------------------

sub generate_image
{
    my( $self     )    = arg_object( shift, $this_class );
    my( $filename )    = arg_string( shift );
    my( $full_x_size ) = arg_int( shift, 0, 1e30 );
    my( $full_y_size ) = arg_int( shift, (($full_x_size < 1) ? 1 : 0), 1e30 );
    end_args( @_ );

    my( $image_index, @x_pos, @y_pos, @x_size, @y_size,
        @tmp_image_files, $white_space_colour, $layout_args,
        $tmp_filename, $header, $footer, $tmp_file, $pos );

    #-------------------------------------------------------------------------
    # foreach sub-image, ask the associated view to compute any view
    # parameters not explicitly set already.  This can be the view
    # orientation, the view direction, and the bounding box of the view
    #-------------------------------------------------------------------------

    for( $image_index = 0;  $image_index < $self->{N_IMAGES};  ++$image_index )
    {
        $self->{IMAGES}[$image_index]->scene_view()->
                        compute_view_for_object(
                                $self->{IMAGES}[$image_index]->scene_object );
    }

    #-------------------------------------------------------------------------
    # assign the sizes and positions of the sub-images within the montage
    #-------------------------------------------------------------------------

    $header = $self->header();
    $footer = $self->footer();

    ($full_x_size, $full_y_size ) = $self->compute_image_sizes_and_positions(
                             $full_x_size, $full_y_size,
                             2 * $header->height(), 2 * $footer->height(),
                             \@x_pos, \@y_pos, \@x_size, \@y_size );

    #-------------------------------------------------------------------------
    # generate all the sub-images into a set of temporary files
    #-------------------------------------------------------------------------

    @tmp_image_files = ();

    $layout_args = "";

    for( $image_index = 0;  $image_index < $self->{N_IMAGES};  ++$image_index )
    {
        #----------------------------------------------------------------------
        # generate a temporary rgb filename
        #----------------------------------------------------------------------

        $tmp_image_files[$image_index] = get_tmp_file( "rgb" );

        #----------------------------------------------------------------------
        # ask the i'th image_info to render itself to the temporary file
        #----------------------------------------------------------------------

        $self->{IMAGES}[$image_index]->create_image(
                                      $tmp_image_files[$image_index],
                                      $x_size[$image_index],
                                      $y_size[$image_index] );

        #----------------------------------------------------------------------
        # add the desired position to the string containing layout arguments
        #----------------------------------------------------------------------

        $layout_args = $layout_args . " " . $tmp_image_files[$image_index] .
                       " " . $x_pos[$image_index] .
                       " " . $y_pos[$image_index];
    }

    #-------------------------------------------------------------------------
    # generate the header if necessary
    #-------------------------------------------------------------------------

    $white_space_colour = $self->{WHITE_SPACE_COLOUR};

    if( $header->height() > 0 )
    {
        $tmp_file = get_tmp_file( "rgb" );
        @tmp_image_files = ( @tmp_image_files, $tmp_file );
        $header->create_text_image_file( $tmp_file, $white_space_colour,
                                         $full_x_size, 2 * $header->height() );
        $pos = $full_y_size - 2 * $header->height();
        $layout_args = $layout_args . " $tmp_file 0 $pos ";
    }

    #-------------------------------------------------------------------------
    # generate the footer if necessary
    #-------------------------------------------------------------------------

    if( $footer->height() > 0 )
    {
        $tmp_file = get_tmp_file( "rgb" );
        @tmp_image_files = ( @tmp_image_files, $tmp_file );
        $footer->create_text_image_file( $tmp_file, $white_space_colour,
                                         $full_x_size, 2 * $footer->height() );
        $layout_args = $layout_args . " $tmp_file 0 0 ";
    }

    #-------------------------------------------------------------------------
    # determine if the file is an rgb or must be converted from another type
    #-------------------------------------------------------------------------

    if( substr( $filename, -4 ) eq ".rgb" )
        { $tmp_filename = $filename; }
    else
        { $tmp_filename = get_tmp_file( "rgb" ); }

    #-------------------------------------------------------------------------
    # place the N_IMAGES sub-images into a single file: $filename
    #-------------------------------------------------------------------------

    run_executable( "place_images", "$tmp_filename '$white_space_colour' " .
                    " -size $full_x_size $full_y_size $layout_args" );

    if( $tmp_filename ne $filename )
    {
        run_executable( "convert", "$tmp_filename $filename " );
    }
}

1;
