#include  <graphics.h>

main()
{
    Status            status;
    window_struct     *window1, *window2;
    text_struct       text;
    lines_struct      lines;
    polygons_struct   polygons;
    VIO_Point             point;
    int               i, j, x_position, y_position;

    status = G_create_window( "Test Window2", -1, -1, -1, -1, &window2 );

    status = G_create_window( "Test Window1", -1, -1, -1, -1, &window1 );

    /* ------------ define line to be drawn  ------------- */

    initialize_lines( &lines, make_Colour(255,255,0) );

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

    initialize_polygons( &polygons, make_Colour(0,255,255), (Surfprop *) NULL );

    fill_Point( point, 10.0, 10.0, 0.0 );
    add_point_to_polygon( &polygons, &point, (VIO_Vector *) 0 );
    fill_Point( point, 90.0, 10.0, 0.0 );
    add_point_to_polygon( &polygons, &point, (VIO_Vector *) 0 );
    fill_Point( point, 90.0, 90.0, 0.0 );
    add_point_to_polygon( &polygons, &point, (VIO_Vector *) 0 );
    fill_Point( point, 10.0, 90.0, 0.0 );
    add_point_to_polygon( &polygons, &point, (VIO_Vector *) 0 );

    /* draw window 1 */

    G_set_viewport( window1, 100, 300, 100, 300 );

    G_set_view_type( window1, PIXEL_VIEW );

    G_draw_polygons( window1, &polygons );
    G_draw_lines( window1, &lines );

    G_update_window( window1 );

    /* draw window 2 */

    G_set_viewport( window2, 200, 400, 000, 200 );

    G_set_view_type( window2, PIXEL_VIEW );

    G_draw_polygons( window2, &polygons );
    G_draw_lines( window2, &lines );

    G_update_window( window2 );

    print( "Hit return: " );
    (void) getchar();

    /* delete drawing objects and window (text does not need to be deleted */

    delete_lines( &lines );

    delete_polygons( &polygons );

    status = G_delete_window( window1 );
    status = G_delete_window( window2 );

    return( status != OK );
}
