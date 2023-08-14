#include  <bicpl.h>

int  main(
    int   argc,
    char  *argv[] )
{
    int            i, n_points;
    VIO_Real           x_min, x_max, y_min, y_max;
    VIO_Real           x_axis[N_DIMENSIONS];
    VIO_Real           y_axis[N_DIMENSIONS];
    VIO_Real           origin[N_DIMENSIONS];
    VIO_Real           x_trans, y_trans, x_scale, y_scale;
    int            x_width, y_width;
    Volume         volume;
    char           *filename;
    static char    *dim_names[] = { MIxspace, MIyspace, MIzspace };
    VIO_Real           clipped_voxels[2*MAX_DIMENSIONS][MAX_DIMENSIONS];

    initialize_argument_processing( argc, argv );
    (void) get_string_argument( "$AVG128", &filename );

    if( input_volume( filename, 3, dim_names, NC_UNSPECIFIED, FALSE,
                           0.0, 0.0, TRUE, &volume,
                           (minc_input_options *) NULL ) != VIO_OK )
        return( 1 );

    origin[VIO_X] = 22.766879336003331;
    origin[VIO_Y] = 24.504077531014342;
    origin[VIO_Z] = -24.040681453003668;

    x_axis[VIO_X] = 0.85089358491769951;
    x_axis[VIO_Y] = -0.34908997859790536;
    x_axis[VIO_Z] = 0.3424883455941593;

    y_axis[VIO_X] = 0.0;
    y_axis[VIO_Y] = 0.65193891922728087;
    y_axis[VIO_Z] = 0.66450536576761754;

    origin[VIO_X] = 63.5;
    origin[VIO_Y] = 63.5;
    origin[VIO_Z] = 39.5;

    x_axis[VIO_X] = 1.0;
    x_axis[VIO_Y] = 0.0;
    x_axis[VIO_Z] = 0.0;

    y_axis[VIO_X] = 0.0;
    y_axis[VIO_Y] = 1.0;
    y_axis[VIO_Z] = 0.0;

    x_scale = 1.2924091412463503;
    y_scale = x_scale;
    x_trans = 172.0;
    y_trans = 156.0;

    x_width = 344;
    y_width = 312;

    for_less( i, 0, 100 )
    {
/*
        n_points = get_volume_cross_section( volume, origin, x_axis, y_axis,
                                         clipped_voxels );
*/
        get_volume_mapping_range( volume, origin, x_axis, y_axis,
                                  x_trans, y_trans, x_scale, y_scale,
                                  &x_min, &x_max,
                                  &y_min, &y_max );
        print( "%g %g %g %g\n", x_min, x_max, y_min, y_max );
    }

    return( 0 );
}
