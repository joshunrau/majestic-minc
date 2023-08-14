
#include  <graphics.h>

main()
{
    Status            status;
    window_struct     *window;
    lines_struct      lines;

    status = G_create_window( "Test Window", -1, -1, -1, -1, &window );

    G_set_background_colour( window, DARK_SLATE_GREY );
    G_set_automatic_clear_state( window, FALSE );

    /* ------------ define line to be drawn  ------------- */

    initialize_lines( &lines, 1 );

    lines.n_points = 4;
    ALLOC( lines.points, 4 );
    fill_Point( lines.points[0],  10.0,  10.0, 0.0 );
    fill_Point( lines.points[1], 100.0,  10.0, 0.0 );
    fill_Point( lines.points[2], 100.0, 100.0, 0.0 );
    fill_Point( lines.points[3],  10.0, 100.0, 0.0 );

    lines.n_items = 1;
    ALLOC( lines.end_indices, lines.n_items );
    lines.end_indices[0] = 5;

    ALLOC( lines.indices, lines.end_indices[lines.n_items-1] );
    lines.indices[0] = 0;
    lines.indices[1] = 1;
    lines.indices[2] = 2;
    lines.indices[3] = 3;
    lines.indices[4] = 0;

    G_set_bitplanes( window, NORMAL_PLANES );
    G_clear_window( window );
    G_update_window( window );

    do
    {
        G_set_bitplanes( window, OVERLAY_PLANES );
        G_set_overlay_colour_map( window, 1, GREEN );
        G_set_overlay_colour_map( window, 2, GREEN );
        G_set_overlay_colour_map( window, 3, GREEN );
        G_clear_window( window );
        G_set_view_type( window, PIXEL_VIEW );
        G_draw_lines( window, &lines );
        G_update_window( window );
    }
    while( getchar() == '\n' );

    /* delete drawing objects and window (text does not need to be deleted */

    delete_lines( &lines );

    status = G_delete_window( window );

    return( status != OK );
}
