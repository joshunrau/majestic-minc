#include  <graphics.h>

int  main(
    int   argc,
    char  *argv[] )
{
    Status         status;
    int            n_alloced, x_size, y_size, i, sizes[MAX_DIMENSIONS];
    VIO_Real           intensity, separation[MAX_DIMENSIONS];
    VIO_Real           min_value, max_value;
    VIO_Real           min_real, max_real;
    VIO_Real           fraction_done;
    pixels_struct  pixels;
    Volume         volume;
    window_struct  *window;
    Minc_file      file;
    VIO_Real           scale;
    VIO_Colour         *rgb_map;
    char           *filename;
    static STRING  dim_names[N_DIMENSIONS] = { MIxspace, MIyspace, MIzspace };

    initialize_argument_processing( argc, argv );
    (void) get_string_argument( "/nil/david/mr_data/indiv_mr.mnc", &filename );

    volume = create_volume( 3, dim_names, NC_UNSPECIFIED, FALSE,
                            0.0, 0.0 );
    file = initialize_minc_input( filename, volume );

    if( file == (Minc_file) NULL )
    {
        print( "Error reading file.\n" );
        return( 1 );
    }

    while( input_more_minc_file( file, &fraction_done ) )
    {
    }

    (void) close_minc_input( file );

    get_volume_sizes( volume, sizes );
    get_volume_separations( volume, separation );

    get_volume_range( volume, &min_real, &max_real );
    get_volume_voxel_range( volume, &min_value, &max_value );

    print( "Volume %s: %d by %d by %d\n",
            filename, sizes[X], sizes[Y], sizes[Z] );
    print( "Separation: %g %g %g\n",
            separation[X], separation[Y], separation[Z] );
    print( "[%g,%g] == [%g,%g]\n", min_value, max_value, min_real, max_real );

    status = G_create_window( "Volume Browser", -1, -1, -1, -1, &window );

    G_get_window_size( window, &x_size, &y_size );

    n_alloced = 0;

    ALLOC( rgb_map, (int) max_value + 1 );

    for_less( i, 0, (int) max_value + 1 )
    {
        intensity = (VIO_Real) i / max_value;
        rgb_map[i] = make_Colour_0_1( intensity, intensity, intensity );
    }

    scale = (VIO_Real) x_size / (VIO_Real) sizes[X];
    if( (VIO_Real) y_size / (VIO_Real) sizes[Y] < scale )
        scale = (VIO_Real) y_size / (VIO_Real) sizes[Y];

    create_volume_slice( volume, (VIO_Real) sizes[Z] / 2.0,
                         0.0, 0.0, scale, scale,
                         (Volume) NULL, 0.0, 0.0, 0.0, 0.0, 0.0,
                         X, Y, x_size, y_size, RGB_PIXEL, FALSE,
                         (unsigned short **) NULL,
                         &rgb_map, &n_alloced, &pixels );

    G_set_view_type( window, PIXEL_VIEW );
    G_draw_pixels( window, &pixels );
    G_update_window( window );

    print( "Hit return: " );
    (void) getchar();

    return( status != OK );
}
