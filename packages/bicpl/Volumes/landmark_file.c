#include  <mni.h>

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_default_landmark_file_suffix
@INPUT      : 
@OUTPUT     : 
@RETURNS    : "lmk"
@DESCRIPTION: Returns the default suffix for landmark files.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI char  *get_default_landmark_file_suffix()
{
    return( "lmk" );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_landmark_file
@INPUT      : volume
              filename
              colour
              size
              type
@OUTPUT     : 
@RETURNS    : VIO_OK or VIO_ERROR
@DESCRIPTION: Reads a landmark.lmk file and creates an object list consisting
              of markers, with the given colour, size, and type.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI Status   input_landmark_file(
    Volume         volume,
    char           filename[],
    Colour         colour,
    VIO_Real           size,
    Marker_types   type,
    int            *n_objects,
    object_struct  **object_list[] )
{
    Status                  status;
    object_struct           *object;
    marker_struct           marker;
    FILE                    *file;

    status = open_file_with_default_suffix( filename, "lmk", READ_FILE,
                                            ASCII_FORMAT, &file );

    *n_objects = 0;

    if( status == VIO_OK )
    {
        while( io_tag_point( file, READ_FILE, volume, size, &marker ) == VIO_OK )
        {
            marker.colour = colour;
            marker.type = type;
            object = create_object( MARKER );
            *(get_marker_ptr(object)) = marker;

            add_object_to_list( n_objects, object_list, object );
        }

        status = close_file( file );
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_tag_point
@INPUT      : file
              io_direction
              volume
              size
@OUTPUT     : marker
@RETURNS    : VIO_OK or VIO_ERROR
@DESCRIPTION: Inputs one marker in landmark format from the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI Status  io_tag_point(
    FILE            *file,
    IO_types        io_direction,
    Volume          volume,
    VIO_Real            size,
    marker_struct   *marker )
{
    Status   status;
    STRING   line;
    Point    position;
    int      sizes[MAX_DIMENSIONS];
    int      len, offset;
    VIO_Real     voxel[MAX_DIMENSIONS];
    VIO_Real     x, y, z;
    VIO_Real     x_w, y_w, z_w;

    status = VIO_OK;

    if( volume != (Volume) NULL && get_volume_n_dimensions(volume) != 3 )
    {
        print( "Error:  volume must be 3d to use for input landmarks.\n" );
        volume = (Volume) NULL;
    }

    if( io_direction == WRITE_FILE )
    {
        if( volume == (Volume) NULL )
        {
            position = marker->position;
        }
        else
        {
            convert_world_to_voxel( volume,
                                    Point_x(marker->position),
                                    Point_y(marker->position),
                                    Point_z(marker->position),
                                    voxel );

            get_volume_sizes( volume, sizes );

            convert_voxel_to_talairach( voxel[VIO_X], voxel[VIO_Y], voxel[VIO_Z],
                                        sizes[VIO_X], sizes[VIO_Y], sizes[VIO_Z],
                                        &x, &y, &z );

            fill_Point( position, x, y, z );
        }
    }

    if( status == VIO_OK )
        status = io_point( file, io_direction, ASCII_FORMAT, &position );

    if( io_direction == READ_FILE )
    {
        marker->colour = WHITE;
        marker->type = BOX_MARKER;

        if( volume == (Volume) NULL )
        {
            marker->position = position;
        }
        else
        {
            get_volume_sizes( volume, sizes );

            convert_talairach_to_voxel( Point_x(position),
                                        Point_y(position),
                                        Point_z(position),
                                        sizes[VIO_X], sizes[VIO_Y], sizes[VIO_Z],
                                        &voxel[VIO_X], &voxel[VIO_Y], &voxel[VIO_Z] );

            convert_voxel_to_world( volume, voxel, &x_w, &y_w, &z_w );
            fill_Point( marker->position, x_w, y_w, z_w );
        }
    }

#define USE_X_POSITION_FOR_WEIGHT
#ifdef  USE_X_POSITION_FOR_WEIGHT
    if( status == VIO_OK )
    {
        if( io_direction == WRITE_FILE )
            status = io_float( file, io_direction, ASCII_FORMAT,
                               &Point_x(position));
        else
        {
            status = io_real( file, io_direction, ASCII_FORMAT, &marker->size );
            marker->size = size;
        }
    }
#else
    if( status == VIO_OK )
        status = io_real( file, io_direction, ASCII_FORMAT, &marker->size );
#endif

    if( status == VIO_OK )
        status = io_int( file, io_direction, ASCII_FORMAT,
                         &marker->structure_id );

    if( status == VIO_OK )
        status = io_int( file, io_direction, ASCII_FORMAT,
                         &marker->patient_id );

    if( io_direction == WRITE_FILE )
    {
        if( status == VIO_OK && strlen(marker->label) > 0 )
            status = io_quoted_string( file, io_direction, ASCII_FORMAT,
                                       marker->label, MAX_STRING_LENGTH );
    }
    else
    {
        if( status == VIO_OK )
            status = input_line( file, line, MAX_STRING_LENGTH );

        if( status == VIO_OK )
        {
            strip_blanks( line, line );

            if( line[0] == '"' )
                offset = 1;
            else
                offset = 0;

            (void) strcpy( marker->label, &line[offset] );

            len = strlen( marker->label );

            if( len > 0 && marker->label[len-1] == '"' )
                 marker->label[len-1] = (char) 0;
        }
    }

    if( status == VIO_OK )
        status = io_newline( file, io_direction, ASCII_FORMAT );

    return( status );
}
