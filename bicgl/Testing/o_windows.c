 
#include  <internal_volume_io.h>
#include  <GS_graphics.h>

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_create_window
@INPUT      : title
              x_pos
              y_pos
              width
              height
@OUTPUT     : window
@RETURNS    : 
@DESCRIPTION: Creates a window of the specified size and shape.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

  Status  GS_create_window(
    GSwindow       window,
    STRING         title,
    int            x_pos,
    int            y_pos,
    int            width,
    int            height,
    VIO_BOOL        colour_map_desired,
    VIO_BOOL        double_buffer_desired,
    VIO_BOOL        depth_buffer_desired,
    int            n_overlay_planes_desired,
    VIO_BOOL        *actual_colour_map_flag,
    VIO_BOOL        *actual_double_buffer_flag,
    VIO_BOOL        *actual_depth_buffer_flag,
    int            *actual_n_overlay_planes )
{
    Status       status;

    status = WS_create_window( title, x_pos, y_pos, width, height,
                               colour_map_desired, double_buffer_desired,
                               depth_buffer_desired,
                               0, actual_colour_map_flag,
                               actual_double_buffer_flag,
                               actual_depth_buffer_flag,
                               actual_n_overlay_planes, &window->WS_window);

    return( status );
}
