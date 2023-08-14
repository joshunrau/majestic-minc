#include  <internal_volume_io.h>
#include  <x_windows.h>

#define  SIZE   100

int  main(
    int    argc,
    char   *argv[] )
{
    int                   i, j, n_colours, min_colour, max_colour, intensity, n;
    VIO_BOOL               done;
    Status                status;
    X_window_struct       window;
    pixels_struct         pixels;
    Event_types           event_type;
    Window_id             event_window_id;
    event_info_struct     info;

    status = X_create_window( "Test X Window", 100, 100, 300, 300,
                              TRUE, FALSE, FALSE, 0, &window );

    if( status != OK )
        return( 1 );

    done = FALSE;
    while( !done )
    {
        n = 0;
        while( !X_get_event( &event_type, &event_window_id, &info ) )
        {
            ++n;
        }

        if( event_type != NO_EVENT && n != 0 )
            print( "[%d]\t", n );

        if( window.window_id != event_window_id )
            print( "Window id: %d\t", event_window_id );

        switch( event_type )
        {
        case NO_EVENT:   break;

        case KEY_DOWN_EVENT:
            print( "<%c> %d down\n", info.key_pressed, info.key_pressed );
            if( info.key_pressed == '' )
                done = TRUE;
            break;

        case KEY_UP_EVENT:
            print( "<%c> %d up\n", info.key_pressed, info.key_pressed );
            break;

        case MOUSE_MOVEMENT_EVENT:
            /* print( "Mouse %d %d\n", info.x_mouse, info.y_mouse ); */
            break;

        case LEFT_MOUSE_DOWN_EVENT:
            print( "Left Mouse Down %d %d\n", info.x_mouse, info.y_mouse );
            break;

        case MIDDLE_MOUSE_DOWN_EVENT:
            print( "Middle Mouse Down %d %d\n", info.x_mouse, info.y_mouse );
            break;

        case RIGHT_MOUSE_DOWN_EVENT:
            print( "Right Mouse Down %d %d\n", info.x_mouse, info.y_mouse );
            break;

        case LEFT_MOUSE_UP_EVENT:
            print( "Left Mouse Up %d %d\n", info.x_mouse, info.y_mouse );
            break;

        case MIDDLE_MOUSE_UP_EVENT:
            print( "Middle Mouse Up %d %d\n", info.x_mouse, info.y_mouse );
            break;

        case RIGHT_MOUSE_UP_EVENT:
            print( "Right Mouse Up %d %d\n", info.x_mouse, info.y_mouse );
            break;

        case REDRAW_OVERLAY_EVENT:
            print( "Redraw Overlay.\n" );
            break;

        case WINDOW_REDRAW_EVENT:
            print( "Redraw Window.\n" );
            break;

        case WINDOW_RESIZE_EVENT:
            print( "Resize Window %d %d \t %d %d.\n",
                   info.x_position, info.y_position,
                   info.x_size, info.y_size );
            break;

        case WINDOW_ICONIZED_EVENT:
            print( "Iconize Window.\n" );
            break;

        case WINDOW_DEICONIZED_EVENT:
            print( "Deiconize Window.\n" );
            break;

        case WINDOW_LEAVE_EVENT:
            print( "Leave Window %d.\n", event_window_id );
            break;

        case WINDOW_ENTER_EVENT:
            print( "Enter Window %d.\n", event_window_id );
            break;

        case WINDOW_QUIT_EVENT:
            print( "Quit Window.\n" );
            break;
        }
    }

#ifdef  DRAW
    n_colours = X_get_n_colours( &window );

    min_colour = MAX( 0, n_colours - 64 );
    max_colour = n_colours-1;

    for_less( i, min_colour, max_colour )
    {
        intensity = ((VIO_Real) i - (VIO_Real) min_colour) /
                    ((VIO_Real) max_colour - (VIO_Real) min_colour) * 255.0;

        X_set_colour_map_entry( &window, i, intensity, intensity, intensity );
    }

    initialize_pixels( &pixels, 10, 10, SIZE, SIZE, 1.0, 1.0,
                       COLOUR_INDEX_8BIT_PIXEL );

    for_less( i, 0, SIZE )
    {
        for_less( j, 0, SIZE )
        {
            intensity = min_colour + (max_colour - min_colour) *
                        (VIO_Real) i * (VIO_Real) j / (VIO_Real) SIZE / (VIO_Real) SIZE;
        }
    }

    X_draw_pixels( &window, &pixels );

    X_update_window( &window );
#endif

    X_delete_window( &window );

    return( 0 );
}
