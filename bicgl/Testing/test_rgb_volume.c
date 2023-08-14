#include  <graphics.h>

int  main(
    int   argc,
    char  *argv[] )
{
    Status         status;
    int            n_alloced, x_size, y_size, i, sizes[MAX_DIMENSIONS];
    VIO_Real           separations[MAX_DIMENSIONS];
    VIO_Real           min_value, max_value;
    VIO_Real           x_axis[N_DIMENSIONS];
    VIO_Real           y_axis[N_DIMENSIONS];
    VIO_Real           origin[N_DIMENSIONS];
    int            used_x_viewport_size, used_y_viewport_size;
    int            v[MAX_DIMENSIONS];
    pixels_struct  pixels1;
    Volume         volume;
    window_struct  *window;
    VIO_Real           x_scale, y_scale, x_translation, y_translation;
    VIO_Real           slice_fit_oversize = 0.1;
    VIO_Colour         colour;
    char           *filename;
    static char    *dim_names[] = { MIxspace, MIyspace, MIzspace };
    minc_input_options  options;

    initialize_argument_processing( argc, argv );
    (void) get_string_argument( "/nil/david/big_data/sphere.fre", &filename );

    set_default_minc_input_options( &options );
    set_minc_input_vector_to_colour_flag( &options, TRUE );

    status = input_volume( filename, 3, dim_names, NC_UNSPECIFIED, FALSE,
                           0.0, 0.0, TRUE, &volume, &options );

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

    status = G_create_window( "Volume Browser", -1, -1, -1, -1, &window );
    G_set_double_buffer_state( window, FALSE );

    G_get_window_size( window, &x_size, &y_size );

    n_alloced = 0;

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

    create_volume_slice( volume, BOX_FILTER, 0.0,
                         origin, x_axis, y_axis,
                         x_translation, y_translation,
                         x_scale, y_scale,
                         (Volume) NULL, BOX_FILTER, 0.0,
                         (VIO_Real *) NULL, (VIO_Real *) NULL, (VIO_Real *) NULL,
                         0.0, 0.0, 0.0, 0.0,
                         x_size, y_size, RGB_PIXEL, FALSE,
                         (unsigned short **) NULL,
                         NULL, BLACK, NULL, &n_alloced, &pixels1 );

    G_set_view_type( window, PIXEL_VIEW );
    G_draw_pixels( window, &pixels1 );

    G_update_window( window );

    BEGIN_ALL_VOXELS( volume, v[0], v[1], v[2], v[3], v[4] )

        GET_VOXEL( colour, volume, v[0], v[1], v[2], v[3], v[4] );

        print( "%d %d %d\n", get_Colour_r(colour),
                             get_Colour_g(colour),
                             get_Colour_b(colour) );
    END_ALL_VOXELS


    print( "Hit return: " );
    (void) getchar();

    return( status != OK );
}
