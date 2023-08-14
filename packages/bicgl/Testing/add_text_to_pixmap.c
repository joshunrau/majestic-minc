
#include  <internal_volume_io.h>
#include  <graphics.h>
#include  <images.h>

static  void  usage(
    char   executable[] )
{
    char   message[] = "\n\
Usage: %s  input.rgb  output.rgb   [x y size colour text] ...\n\n";

    print_error( message, executable );
}


int main(
    int    argc,
    char   *argv[] )
{
    window_struct     *window;
    text_struct       text;
    Font_types        font_type;
    pixels_struct     pixels;
    STRING            input_filename, output_filename, colour_name, string;
    VIO_Colour            colour;
    int               x, y, font_size, x_size, y_size;
    VIO_Point             point;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    if( input_rgb_file( input_filename, &pixels ) != VIO_OK )
        return( 1 );

    x_size = pixels.x_size;
    y_size = pixels.y_size;

    if( G_create_window( "Add Text To Pixmap",
                         100, 100, x_size, y_size,
                         FALSE, FALSE, FALSE, 0, &window ) != VIO_OK )
        return( 1 );


    G_clear_window( window );

    G_set_zbuffer_state( window, FALSE );
    G_set_lighting_state( window, FALSE );
    G_set_view_type( window, PIXEL_VIEW );
    G_draw_pixels( window, &pixels );

    while( get_int_argument( 0, &x ) &&
           get_int_argument( 0, &y ) &&
           get_int_argument( 0, &font_size ) &&
           get_string_argument( NULL, &colour_name ) &&
           get_string_argument( NULL, &string ) )
    {
        if( !lookup_colour( colour_name, &colour ) )
            print_error( "Invalid colour: %s\n", colour_name );

        fill_Point( point, (VIO_Real) x, (VIO_Real) y, 0.0 );
        if( font_size <= 0 )
            font_type = FIXED_FONT;
        else
            font_type = SIZED_FONT;

        initialize_text( &text, &point, colour, font_type, (VIO_Real) font_size );
        text.string = create_string( string );

        G_draw_text( window, &text );
        delete_text( &text );
    }

    G_update_window( window );

    G_read_pixels( window, 0, x_size-1, 0, y_size-1,
                   &PIXEL_RGB_COLOUR(pixels,0,0) );

    if( output_rgb_file( output_filename, &pixels ) != VIO_OK )
        return( 1 );

    delete_pixels( &pixels );

    (void) G_delete_window( window );

    G_terminate();

    return( 0 );
}
