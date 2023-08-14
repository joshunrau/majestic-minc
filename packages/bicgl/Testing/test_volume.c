#include  <graphics.h>

#define  DEGREES_CONTINUITY  -1

int  main(
    int   argc,
    char  *argv[] )
{
    Status         status;
    int            n_alloced, x_size, y_size, i, sizes[MAX_DIMENSIONS];
    int            iter, n_iters;
    VIO_Real           alpha;
    VIO_Real           intensity, separations[MAX_DIMENSIONS];
    VIO_Real           min_value, max_value;
    VIO_Real           x_axis[N_DIMENSIONS];
    VIO_Real           y_axis[N_DIMENSIONS];
    VIO_Real           origin[N_DIMENSIONS];
    int            used_x_viewport_size, used_y_viewport_size;
    pixels_struct  pixels1, pixels2;
    Volume         volume;
    window_struct  *window;
    VIO_Real           x_scale, y_scale, x_translation, y_translation;
    VIO_Real           slice_fit_oversize = 0.1;
    VIO_Colour         *rgb_map;
    char           *filename;
    static char    *dim_names[] = { MIxspace, MIyspace, MIzspace };

    initialize_argument_processing( argc, argv );
    (void) get_string_argument( "/nil/david/big_data/sphere.fre", &filename );
    (void) get_real_argument( 1.0, &alpha );

    status = input_volume( filename, 3, dim_names, NC_UNSPECIFIED, FALSE,
                           0.0, 0.0, TRUE, &volume,
                           (minc_input_options *) NULL );

    get_volume_voxel_range( volume, &min_value, &max_value );

    print( "%g %g   %g %g\n", min_value, max_value,
           CONVERT_VOXEL_TO_VALUE( volume, min_value ),
           CONVERT_VOXEL_TO_VALUE( volume, max_value ) );

    if( status != OK )
        return( 1 );

    get_volume_sizes( volume, sizes );
    get_volume_separations( volume, separations );
    print( "Volume %s: %d by %d by %d\n",
            filename, sizes[X], sizes[Y], sizes[Z] );
    print( "Thickness: %g %g %g\n",
            separations[X], separations[Y], separations[Z] );

    status = G_create_window( "Volume Browser", -1, -1, -1, -1,
                              FALSE, TRUE, FALSE, 0, &window );
    G_set_double_buffer_state( window, FALSE );
    G_set_transparency_state( window, TRUE );

    G_get_window_size( window, &x_size, &y_size );

    n_alloced = 0;

    get_volume_voxel_range( volume, &min_value, &max_value );
    ALLOC( rgb_map, (int)max_value+1 );

    for_less( i, 0, (int) max_value+1 )
    {
        intensity = (VIO_Real) i / max_value;
        rgb_map[i] = make_rgba_Colour_0_1( intensity, intensity, intensity, 1.0);
    }

    origin[X] = 0.0;
    origin[Y] = 0.0;
    origin[Z] = (VIO_Real) (sizes[Z] - 1) / 2.0;
    x_axis[X] = 1.0;
    x_axis[Y] = 0.0;
    x_axis[Z] = 0.0;
    y_axis[X] = 0.0;
    y_axis[Y] = 1.0;
    y_axis[Z] = 0.0;
    fit_volume_slice_to_viewport( volume, origin, x_axis, y_axis,
                                  x_size, y_size, slice_fit_oversize,
                                  &x_translation, &y_translation,
                                  &x_scale, &y_scale,
                                  &used_x_viewport_size, &used_y_viewport_size);

    n_iters = 100;
    n_iters = 1;
    start_timing();
    for_less( iter, 0, n_iters )
    create_volume_slice( volume, NEAREST_NEIGHBOUR, 0.0,
                         origin, x_axis, y_axis,
                         x_translation, y_translation,
                         x_scale, y_scale,
                         (Volume) NULL, NEAREST_NEIGHBOUR, 0.0,
                         (VIO_Real *) NULL, (VIO_Real *) NULL, (VIO_Real *) NULL,
                         0.0, 0.0, 0.0, 0.0,
                         x_size, y_size, RGB_PIXEL, DEGREES_CONTINUITY,
                         (unsigned short **) NULL,
                         &rgb_map, BLACK, NULL, &n_alloced, &pixels1 );
    end_timing( "Creating Pixmap", n_iters );

    if( alpha != 1.0 )
    {
        origin[X] = 0.0;
        origin[Y] = 0.0;
        origin[Z] = (VIO_Real) (sizes[Z] - 1) / 3.0;

        for_less( i, 0, (int) max_value+1 )
        {
            intensity = (VIO_Real) i / max_value;
            if( i < max_value / 2 )
                rgb_map[i] = make_rgba_Colour_0_1( 0.5, 0.5, 0.5, 0.0 );
            else
                rgb_map[i] = make_rgba_Colour_0_1( 0.0, 0.5, 0.0, 1.0 - alpha );
        }

        create_volume_slice( volume, NEAREST_NEIGHBOUR, 0.0,
                             origin, x_axis, y_axis,
                             x_translation, y_translation,
                             x_scale, y_scale,
                             (Volume) NULL, NEAREST_NEIGHBOUR, 0.0,
                             (VIO_Real *) NULL, (VIO_Real *) NULL, (VIO_Real *) NULL,
                             0.0, 0.0, 0.0, 0.0,
                             x_size, y_size, RGB_PIXEL, DEGREES_CONTINUITY,
                             (unsigned short **) NULL,
                             &rgb_map, BLACK, NULL, &n_alloced, &pixels2 );
    }

    G_set_view_type( window, PIXEL_VIEW );

    n_iters = 1000;
    n_iters = 1;
    start_timing();
    for_less( iter, 0, n_iters )
        G_draw_pixels( window, &pixels1 );
    end_timing( "Drawing Pixmap", n_iters );

    if( alpha != 1.0 )
        G_draw_pixels( window, &pixels2 );

    G_update_window( window );

    print( "Hit return: " );
    (void) getchar();

    return( status != OK );
}
