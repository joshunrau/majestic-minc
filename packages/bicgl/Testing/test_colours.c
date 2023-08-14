
#include  <graphics.h>

#define  LIGHT_INDEX         0

#define  PIXELS_X_ZOOM       2.0
#define  PIXELS_Y_ZOOM       2.0

main()
{
    Status            status;
    window_struct     *window;
    int               i, j;

    status = G_create_window( "Test Window", -1, -1, -1, -1, &window );

    G_set_double_buffer_state( window, TRUE );
    G_set_colour_map_state( window, TRUE );
    G_set_background_colour( window, 256 );

    G_set_colour_map_entry( window, 256, DARK_SLATE_GREY );

    G_clear_window( window );
    G_update_window( window );

    for_less( j, 0, 100000 )
    for_less( i, 288, 600 )
        G_set_colour_map_entry( window, i, WHITE );

    (void) getchar();

    status = G_delete_window( window );

    return( status != OK );
}
