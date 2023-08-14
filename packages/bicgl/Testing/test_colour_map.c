
#include  <graphics.h>

#define  HEIGHT   10
#define  MULTIPLES   100

int main(
    int    argc,
    char   *argv[] )
{
    Status            status;
    window_struct     *window;
    pixels_struct     pixels;
    int               x_position, y_position, n_colours;
    int               x_pixel, y_pixel;
    int               i, j, pixels_x_size, pixels_y_size;
    int               n_c;
    VIO_Colour            colours_list[7];

    status = G_create_window( "Test Window",
                              100, 600, 300, 300,
                              TRUE, FALSE, FALSE, 0, &window );

    /* ------------ define pixels to be drawn  ------------- */

    pixels_x_size = 256;
    pixels_y_size = 256;

    x_position = 10;
    y_position = 10;

    initialize_pixels( &pixels, x_position, y_position,
                       pixels_x_size, pixels_y_size,
                       1.0, 1.0, COLOUR_INDEX_16BIT_PIXEL );

    n_colours = G_get_n_colour_map_entries( window );
    n_c = VIO_SIZEOF_STATIC_ARRAY( colours_list );
    colours_list[0] = RED;
    colours_list[1] = GREEN;
    colours_list[2] = BLUE;
    colours_list[3] = CYAN;
    colours_list[4] = MAGENTA;
    colours_list[5] = YELLOW;
    colours_list[6] = WHITE;

    for_less( i, 1, n_colours )
    {
        G_set_colour_map_entry( window, i, colours_list[i%n_c] );
    }

    for_less( i, 0, pixels_x_size )
    {
        for_less( j, 0, pixels_y_size )
        {
            PIXEL_COLOUR_INDEX_16(pixels,i,j) = 
                       IJ(j/HEIGHT,i/HEIGHT,pixels_x_size) % n_colours;
        }
    }

    G_clear_window( window );

    G_set_zbuffer_state( window, FALSE );
    G_set_lighting_state( window, FALSE );
    G_set_view_type( window, PIXEL_VIEW );

    G_draw_pixels( window, &pixels );

    G_update_window( window );

    (void) getchar();

    delete_pixels( &pixels );

    status = G_delete_window( window );

    return( status != OK );
}
