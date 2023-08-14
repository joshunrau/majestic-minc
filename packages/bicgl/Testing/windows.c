#include  <internal_volume_io.h>
#include  <graphics.h>

  Status  G_create_window(
    STRING         title,
    int            x_pos,
    int            y_pos,
    int            width,
    int            height,
    VIO_BOOL        colour_map_desired,
    VIO_BOOL        double_buffer_desired,
    VIO_BOOL        depth_buffer_desired,
    int            n_overlay_planes_desired,
    Gwindow        *window )
{
    Status   status;
    VIO_BOOL  actual_colour_map_flag;
    VIO_BOOL  actual_double_buffer_flag;
    VIO_BOOL  actual_depth_buffer_flag;
    int      actual_n_overlay_planes;

    ALLOC( *window, 1 );

    ALLOC( (*window)->GS_window, 1 );

    status = GS_create_window( (*window)->GS_window,
                               title, x_pos, y_pos, width, height,
                               colour_map_desired,
                               double_buffer_desired,
                               depth_buffer_desired,
                               0,
                               &actual_colour_map_flag,
                               &actual_double_buffer_flag,
                               &actual_depth_buffer_flag,
                               &actual_n_overlay_planes );

    return( status );
}
