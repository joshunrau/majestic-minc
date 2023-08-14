/*------------------------------------------------------

    example_graphics.c

This program demonstrates the use of the BIC Graphics
Library by creating a window and displaying some objects
through a simple interaction mechanism.

------------------------------------------------------ */

#include  <graphics.h>

#define   LIGHT_INDEX   0

int main(
    int    argc,
    char   *argv[] )
{
    Status            status;
    window_struct     *window, *event_window;
    text_struct       text;
    lines_struct      lines;
    polygons_struct   polygons;
    pixels_struct     pixels;
    static Surfprop   spr = { 0.2, 0.5, 0.5, 20.0, 1.0 };
    VIO_Point             point;
    VIO_Vector            normal, light_direction;
    Event_types       event_type;
    VIO_BOOL           update_required, done;
    int               key_pressed;
    int               x_position, y_position, x_size, y_size;
    int               i, j, pixels_x_size, pixels_y_size;
    static VIO_Point      origin = { 0.0, 0.0, 2.0 };
    static VIO_Vector     up_direction = { 0.0, 1.0, 0.0 };
    static VIO_Vector     line_of_sight = { 0.0, 0.0, -1.0 };

    /* ------- create a graphics window ------ */

    status = G_create_window( "Test Window", -1, -1, -1, -1,
                              &window );

    /*--- define the 3D view */

    G_set_3D_view( window, &origin, &line_of_sight,
                   &up_direction, 0.01, 4.0, TRUE, 2.0,
                   FALSE, 0.0, 2.0, 2.0 );

    /* -------- define text to be drawn
                (text.string will be filled in later) ----- */

    fill_Point( point, 10.0, 10.0, 0.0 );
    initialize_text( &text, &point, make_Colour(255,0,255),
                     SIZED_FONT, 14.0 );
    text.string[0] = (char) 0;

    /* ------------ define line to be drawn  ------------- */

    initialize_lines( &lines, make_Colour(255,255,0) );

    lines.n_points = 4;
    ALLOC( lines.points, 4 );
    fill_Point( lines.points[0], 0.0, 0.0, 0.0 );
    fill_Point( lines.points[1], 1.0, 1.0, 0.0 );
    fill_Point( lines.points[2], -0.3, 1.0, -1.0 );
    fill_Point( lines.points[3], 0.5, 0.4, 0.0 );

    lines.n_items = 3;
    ALLOC( lines.end_indices, lines.n_items );
    lines.end_indices[0] = 2;
    lines.end_indices[1] = 4;
    lines.end_indices[2] = 5;

    ALLOC( lines.indices, lines.end_indices[lines.n_items-1] );
    lines.indices[0] = 0;
    lines.indices[1] = 1;

    lines.indices[2] = 0;
    lines.indices[3] = 2;

    lines.indices[4] = 3;

    /* ------------ define pixels to be drawn  ------------- */

    pixels_x_size = 256;
    pixels_y_size = 256;

    x_position = 10;
    y_position = 10;
    initialize_pixels( &pixels, x_position, y_position,
                       pixels_x_size, pixels_y_size,
                       1.0, 1.0, RGB_PIXEL );

    for_less( i, 0, pixels_x_size )
    {
        for_less( j, 0, pixels_y_size )
        {
            PIXEL_RGB_COLOUR(pixels,i,j) =
                   make_Colour( i % 256, j % 256, 0 );
        }
    }

    /* ------------ define polygons to be drawn  ------------- */

    initialize_polygons( &polygons, make_Colour(0,255,255),
                         &spr );

    start_new_polygon( &polygons );

    fill_Point( point, -0.3, -0.3, 0.0 );
    fill_Vector( normal, 0.0, 0.0, 1.0 );
    add_point_to_polygon( &polygons, &point, &normal );
    fill_Point( point, 0.3, -0.3, 0.0 );
    add_point_to_polygon( &polygons, &point, &normal );
    fill_Point( point, 0.3, 0.3, 0.0 );
    add_point_to_polygon( &polygons, &point, &normal );
    fill_Point( point, -0.3, 0.3, 0.0 );
    add_point_to_polygon( &polygons, &point, &normal );

    /* ------------ define a directional light from over
                    the left shoulder ----------------- */

    fill_Vector( light_direction, 1.0, 1.0, -1.0 );

    G_define_light( window, LIGHT_INDEX, DIRECTIONAL_LIGHT,
                    make_Colour(255,255,255),
                    &light_direction, (VIO_Point *) 0, 0.0, 0.0 );
    G_set_light_state( window, LIGHT_INDEX, TRUE );

    G_set_lighting_state( window, TRUE );

    /* --------------------------------------- */
    /* ------------ do main loop ------------- */
    /* --------------------------------------- */

    update_required = TRUE;

    print( "Hit left mouse button to exit\n" );

    done = FALSE;

    do
    {
        /* ----- get events until none found ------- */

        do
        {
            event_type = G_get_event( &event_window,
                                      &key_pressed );

            if( event_window == window )
            {
                switch( event_type )
                {
                case KEY_DOWN_EVENT:
                    print( "Key pressed: \"%c\"\n",
                           key_pressed );
                    break;

                case KEY_UP_EVENT:
                    print( "Key released: \"%c\"\n",
                           key_pressed );
                    break;

                case LEFT_MOUSE_DOWN_EVENT:
                    print( "Left mouse DOWN\n" );
                    done = TRUE;
                    break;

                case WINDOW_REDRAW_EVENT:
                    print( "Window will be redrawn.\n" );
                    update_required = TRUE;
                    break;

                case WINDOW_RESIZE_EVENT:
                    G_get_window_position( window,
                                           &x_position,
                                           &y_position );
                    G_get_window_size( window,
                                       &x_size, &y_size );
                    print( "Window resized, " );
                    print( " new position: %d %d\n",
                           x_position, y_position );
                    print( " new size: %d %d\n",
                           x_size, y_size );
                    update_required = TRUE;
                    break;
                }
            }
        }     /* break to do update when no events */
        while( event_type != NO_EVENT );

        /* ------- if one or more events caused an update,
                   redraw the screen */

        if( update_required )
        {
            /* ----- clear the window */

            G_clear_window( window );

            /* ----- draw the objects */

            G_draw_pixels( window, &pixels );

            G_set_view_type( window, MODEL_VIEW );
            G_draw_lines( window, &lines );
            G_draw_polygons( window, &polygons );

            G_set_view_type( window, PIXEL_VIEW );
            G_draw_text( window, &text );

            /* ----- finished drawing scene, now update */

            G_update_window( window );

            update_required = FALSE;
        }
    }
    while( !done );

    /* ----- delete drawing objects and window 
                 (text does not need to be deleted */

    delete_lines( &lines );

    delete_polygons( &polygons );

    delete_pixels( &pixels );

    status = G_delete_window( window );

    G_terminate();

    return( status != VIO_OK );
}
