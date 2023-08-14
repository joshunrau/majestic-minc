#include  <internal_volume_io.h>
#include  <x_windows.h>

/* ARGSUSED */

  Status  X_create_window(
    STRING           title,
    int              initial_x_pos,
    int              initial_y_pos,
    int              initial_x_size,
    int              initial_y_size,
    VIO_BOOL          colour_map_mode,  
    VIO_BOOL          double_buffer_flag,  
    VIO_BOOL          depth_buffer_flag,  
    int              n_overlay_planes,
    X_window_struct  *window )
{
    int              depth, class;
    Status           status;
    XVisualInfo      visual;
    Colormap         cmap;
    XGCValues        gc_values;

    if( colour_map_mode )
        class = PseudoColor;
    else
        class = TrueColor;

    depth = 8;

    if( !XMatchVisualInfo( X_get_display(),  X_get_screen(),
                           depth, class, &visual ) )
    {
        print_error( "Unable to get needed visual\n" );
        return( ERROR );
    }

    cmap = DefaultColormap( X_get_display(), X_get_screen() );

    status = X_create_window_with_visual( title,
                                          initial_x_pos, initial_y_pos,
                                          initial_x_size, initial_y_size,
                                          colour_map_mode, &visual, cmap,
                                          window );

    if( status == OK )
    {
        gc_values.foreground = 0;
        gc_values.background = 0;
        window->graphics_context = XCreateGC( X_get_display(),
                                              window->window_id,
                                              (GCForeground | GCBackground),
                                              &gc_values );
    }

    return( status );
}

  void  X_draw_pixels(
    X_window_struct   *window,
    pixels_struct     *pixels )
{
    int       format, bitmap_pad, bytes_per_line;
    char      *data;
    XImage    *image;

    if( pixels->x_size <= 0 || pixels->y_size <= 0 )
        return;

    if( window->colour_map_mode != (pixels->pixel_type != RGB_PIXEL) )
    {
        print_error(
          "Invalid pixel type in X_draw_pixels for rgb/colourmap mode.\n");
        return;
    }

    switch( pixels->pixel_type )
    {
    case COLOUR_INDEX_8BIT_PIXEL:
        bitmap_pad = 8;
        bytes_per_line = pixels->x_size;
        data = (void *) &PIXEL_COLOUR_INDEX_8(*pixels,0,0);
        break;
    }

    format = XYPixmap;
    image = XCreateImage( X_get_display(), window->visual.visual,
                          window->visual.depth, format, 0,
                          data, pixels->x_size, pixels->y_size, bitmap_pad,
                          bytes_per_line );

    XPutImage( X_get_display(), X_get_screen(), window->graphics_context,
               image, 0, 0,
               pixels->x_position, pixels->x_position,
               pixels->x_size, pixels->y_size );

    XDestroyImage( image );
}

  void  X_update_window(
    X_window_struct   *window )
{
    XFlush( X_get_display() );
}
