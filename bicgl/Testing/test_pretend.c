#include  <graphics.h>

main()
{
    int               line_colour_index = 256;
    int               polygon_colour_index = 257;
    Status            status;
    window_struct     *window;
    text_struct       text;
    lines_struct      lines;
    polygons_struct   polygons;
    VIO_Point             point;
    int               i, j, x_position, y_position;

    status = G_create_window( "Test Window1", -1, -1, -1, -1, &window );

    G_set_pretend_colour_map_state( window, TRUE );
    G_set_pretend_colour_map( window, line_colour_index, RED );
    G_set_pretend_colour_map( window, polygon_colour_index,
                              GREEN );

    /* ------------ define line to be drawn  ------------- */

    initialize_lines( &lines, line_colour_index );

    lines.n_points = 4;
    ALLOC( lines.points, 4 );
    fill_Point( lines.points[0], 0.0, 0.0, 0.0 );
    fill_Point( lines.points[1], 100.0, 100.0, 0.0 );
    fill_Point( lines.points[2], 30.0, 100.0, 0.0 );
    fill_Point( lines.points[3], 100.0, 30.0, 0.0 );

    lines.n_items = 3;
    ALLOC( lines.end_indices, lines.n_items );
    lines.end_indices[0] = 2;
    lines.end_indices[1] = 4;
    lines.end_indices[2] = 6;

    ALLOC( lines.indices, lines.end_indices[lines.n_items-1] );
    lines.indices[0] = 0;
    lines.indices[1] = 1;

    lines.indices[2] = 0;
    lines.indices[3] = 2;

    lines.indices[4] = 0;
    lines.indices[5] = 3;

    /* ------------ define polygons to be drawn  ------------- */

    initialize_polygons( &polygons, polygon_colour_index, (Surfprop *) NULL );

    fill_Point( point, 10.0, 10.0, 0.0 );
    add_point_to_polygon( &polygons, &point, (VIO_Vector *) 0 );
    fill_Point( point, 90.0, 10.0, 0.0 );
    add_point_to_polygon( &polygons, &point, (VIO_Vector *) 0 );
    fill_Point( point, 90.0, 90.0, 0.0 );
    add_point_to_polygon( &polygons, &point, (VIO_Vector *) 0 );
    fill_Point( point, 10.0, 90.0, 0.0 );
    add_point_to_polygon( &polygons, &point, (VIO_Vector *) 0 );

    G_set_bitplanes( window, NORMAL_PLANES );

    while( TRUE )
    {
        /* draw window */

        if( G_get_colour_map_state(window) )
        {
            G_set_colour_map_state( window, FALSE );
        }
        else
        {
            G_set_colour_map_state( window, TRUE );
            G_set_colour_map_entry( window, line_colour_index, RED );
            G_set_colour_map_entry( window, polygon_colour_index, GREEN );
        }

        G_set_viewport( window, 100, 300, 100, 300 );

        G_set_view_type( window, PIXEL_VIEW );

        G_draw_polygons( window, &polygons );
        G_draw_lines( window, &lines );

        G_update_window( window );

        print( "Hit return: " );
        (void) getchar();
    }

    /* delete drawing objects and window (text does not need to be deleted */

    delete_lines( &lines );

    delete_polygons( &polygons );

    status = G_delete_window( window );

    return( status != OK );
}
