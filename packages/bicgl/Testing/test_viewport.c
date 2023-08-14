
#include  <volume_io.h>
#include  <graphics.h>

#define  X_SIZE  300
#define  Y_SIZE  300

#define  X_MIN_VIEWPORT    10
#define  X_MAX_VIEWPORT   290
#define  Y_MIN_VIEWPORT    10
#define  Y_MAX_VIEWPORT   290

static  void  draw_rectangle(
    Gwindow   window,
    VIO_Colour    colour,
    int       x_min,
    int       x_max,
    int       y_min,
    int       y_max );

static  void  draw_region(
    Gwindow   window,
    VIO_Colour    colour,
    int       x_min,
    int       x_max,
    int       y_min,
    int       y_max );

static  void  draw_pixels(
    Gwindow   window,
    VIO_Colour    colour,
    int       x1,
    int       y1,
    int       x2,
    int       y2 );

int main(
    int    argc,
    char   *argv[] )
{
    Gwindow           window;

    if( G_create_window( "Test Window",
                         100, 600, X_SIZE, Y_SIZE,
                         FALSE, FALSE, FALSE, 0, &window ) != OK )
        return( 1 );
    G_set_zbuffer_state( window, FALSE );
    G_set_lighting_state( window, FALSE );
    G_set_view_type( window, PIXEL_VIEW );
    G_set_automatic_clear_state( window, FALSE );

    G_clear_window( window );

    draw_rectangle( window, YELLOW, 1, X_SIZE-2, 1, Y_SIZE-2 );
    draw_rectangle( window, ORANGE, 0, X_SIZE-1, 0, Y_SIZE-1 );
    draw_rectangle( window, BLUE, -1, X_SIZE, -1, Y_SIZE );

    draw_rectangle( window, RED, X_MIN_VIEWPORT-1, X_MAX_VIEWPORT+1,
                            Y_MIN_VIEWPORT-1, Y_MAX_VIEWPORT+1 );
    draw_rectangle( window, GREEN, X_MIN_VIEWPORT, X_MAX_VIEWPORT,
                            Y_MIN_VIEWPORT, Y_MAX_VIEWPORT );
    draw_rectangle( window, CYAN, X_MIN_VIEWPORT+1, X_MAX_VIEWPORT-1,
                            Y_MIN_VIEWPORT+1, Y_MAX_VIEWPORT-1 );
    draw_rectangle( window, BLUE, X_MIN_VIEWPORT-2, X_MAX_VIEWPORT+2,
                            Y_MIN_VIEWPORT-2, Y_MAX_VIEWPORT+2 );

    draw_pixels( window, RED, ROUND( X_SIZE / 2), 0,
                              ROUND( X_SIZE / 2) + 2, Y_SIZE - 1 );

    G_update_window( window );

    print( "Hit return: " );
    (void) getchar();

    G_set_viewport( window, X_MIN_VIEWPORT, X_MAX_VIEWPORT,
                            Y_MIN_VIEWPORT, Y_MAX_VIEWPORT );

    G_clear_viewport( window, make_Colour_0_1( .2, .2, .2 ) );

    draw_rectangle( window, BLUE, -10, X_MAX_VIEWPORT - X_MIN_VIEWPORT + 10,
                                  -10, (Y_MAX_VIEWPORT + Y_MIN_VIEWPORT)/2 );
    draw_rectangle( window, BLUE, -10, (X_MAX_VIEWPORT + X_MIN_VIEWPORT)/2,
                                  -10, Y_MAX_VIEWPORT - Y_MIN_VIEWPORT + 10 );

    G_update_window( window );

    print( "Viewport cleared:  Hit return: " );
    (void) getchar();

    draw_rectangle( window, WHITE, -1, X_MAX_VIEWPORT - X_MIN_VIEWPORT + 1,
                                   -1, Y_MAX_VIEWPORT - Y_MIN_VIEWPORT + 1 );
    draw_rectangle( window, MAGENTA, 0, X_MAX_VIEWPORT - X_MIN_VIEWPORT,
                                     0, Y_MAX_VIEWPORT - Y_MIN_VIEWPORT );
    draw_rectangle( window, YELLOW, 1, X_MAX_VIEWPORT - X_MIN_VIEWPORT-1,
                                    1, Y_MAX_VIEWPORT - Y_MIN_VIEWPORT-1 );
    draw_pixels( window, RED, ROUND( (X_MAX_VIEWPORT+X_MIN_VIEWPORT)/ 2), 0,
                              ROUND( (X_MAX_VIEWPORT+X_MIN_VIEWPORT)/ 2)+2,
                              Y_MAX_VIEWPORT - Y_MIN_VIEWPORT );
    G_update_window( window );

    print( "Viewport redrawn:  Hit return: " );
    (void) getchar();

    (void) G_delete_window( window );

    return( 0 );
}

static  void  draw_rectangle(
    Gwindow   window,
    VIO_Colour    colour,
    int       x_min,
    int       x_max,
    int       y_min,
    int       y_max )
{
    lines_struct  lines;

    initialize_lines( &lines, colour );

    lines.n_points = 4;
    ALLOC( lines.points, 4 );
    fill_Point( lines.points[0], (VIO_Real) x_min, (VIO_Real) y_min, 0.0 );
    fill_Point( lines.points[1], (VIO_Real) x_max, (VIO_Real) y_min, 0.0 );
    fill_Point( lines.points[2], (VIO_Real) x_max, (VIO_Real) y_max, 0.0 );
    fill_Point( lines.points[3], (VIO_Real) x_min, (VIO_Real) y_max, 0.0 );

    lines.n_items = 1;
    ALLOC( lines.end_indices, lines.n_items );
    lines.end_indices[0] = 5;

    ALLOC( lines.indices, lines.end_indices[lines.n_items-1] );
    lines.indices[0] = 0;
    lines.indices[1] = 1;
    lines.indices[2] = 2;
    lines.indices[3] = 3;
    lines.indices[4] = 0;

    G_draw_lines( window, &lines );

    delete_lines( &lines );
}

static  void  draw_region(
    Gwindow   window,
    VIO_Colour    colour,
    int       x_min,
    int       x_max,
    int       y_min,
    int       y_max )
{
    polygons_struct  polygons;

    initialize_polygons( &polygons, colour, NULL );

    polygons.n_points = 4;
    ALLOC( polygons.points, 4 );
    fill_Point( polygons.points[0], (VIO_Real) x_min, (VIO_Real) y_min, 0.0 );
    fill_Point( polygons.points[1], (VIO_Real) x_max, (VIO_Real) y_min, 0.0 );
    fill_Point( polygons.points[2], (VIO_Real) x_max, (VIO_Real) y_max, 0.0 );
    fill_Point( polygons.points[3], (VIO_Real) x_min, (VIO_Real) y_max, 0.0 );

    polygons.n_items = 1;
    ALLOC( polygons.end_indices, polygons.n_items );
    polygons.end_indices[0] = 4;

    ALLOC( polygons.indices, polygons.end_indices[polygons.n_items-1] );
    polygons.indices[0] = 0;
    polygons.indices[1] = 1;
    polygons.indices[2] = 2;
    polygons.indices[3] = 3;

    G_draw_polygons( window, &polygons );

    delete_polygons( &polygons );
}

static  void  draw_pixel(
    Gwindow   window,
    VIO_Colour    colour,
    VIO_Real      x,
    VIO_Real      y )
{
    lines_struct  lines;

    initialize_lines( &lines, colour );

    lines.n_points = 1;
    ALLOC( lines.points, 1 );
    fill_Point( lines.points[0], x, y, 0.0 );

    lines.n_items = 1;
    ALLOC( lines.end_indices, lines.n_items );
    lines.end_indices[0] = 1;

    ALLOC( lines.indices, lines.end_indices[lines.n_items-1] );
    lines.indices[0] = 0;

    G_draw_lines( window, &lines );

    delete_lines( &lines );
}

static  void  draw_pixels(
    Gwindow   window,
    VIO_Colour    colour,
    int       x1,
    int       y1,
    int       x2,
    int       y2 )
{
    int   i;
    VIO_Real  x;

    for_inclusive( i, y1, y2 )
    {
        x = x1 + (VIO_Real) (x2 - x1) *
            ((VIO_Real) i - (VIO_Real) y1) / ((VIO_Real) y2 - (VIO_Real) y1);
        draw_pixel( window, colour, x, (VIO_Real) i );
    }
}
