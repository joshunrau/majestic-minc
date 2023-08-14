#include <bicpl.h>

#define  BINTREE_FACTOR  0.3

static  void  usage(
    VIO_STR   executable )
{
    VIO_STR  usage_str = "\n\
Usage: %s  input.obj  values.txt [output.txt]\n\
\n\
     .\n\n";

    fprintf(stderr, usage_str, executable );
}

int  main(
    int    argc,
    char   *argv[] )
{
    VIO_STR           input_filename, output_filename, values_filename;
    FILE             *file;
    int              n_objects, *n_neighbours, **neighbours, p, n;
    int              n_non, n_read;
    VIO_Real             avg;
    VIO_File_formats     format;
    object_struct    **object_list;
    polygons_struct  *polygons;
    VIO_Real             *values;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &values_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    (void) get_string_argument( values_filename, &output_filename );

    if( input_graphics_file( input_filename, &format, &n_objects,
                             &object_list ) != VIO_OK )
        return( 1 );

    if( n_objects != 1 || get_object_type(object_list[0]) != POLYGONS )
    {
        fprintf(stderr, "Must be polygons\n" );
        return( 1 );
    }

    polygons = get_polygons_ptr( object_list[0] );

    if( open_file( values_filename, READ_FILE, ASCII_FORMAT, &file ) != VIO_OK )
        return( 1 );

    ALLOC( values, polygons->n_points );

    n_read = 0;
    while( input_real( file, &values[n_read] ) == VIO_OK )
        ++n_read;

    (void) close_file( file );

    for_less( p, n_read, polygons->n_points )
        values[p] = 1.0e30;

    create_polygon_point_neighbours( polygons, FALSE, &n_neighbours,
                                     &neighbours, NULL, NULL );

    for_less( p, n_read, polygons->n_points )
    {
        n_non = 0;
        avg = 0.0;
        for_less( n, 0, n_neighbours[p] )
        {
            if( neighbours[p][n] < n_read )
            {
                ++n_non;
                avg += values[neighbours[p][n]];
            }
        }

        if( n_non != 2 )
            fprintf(stderr, "n_non: %d\n", n_non );
        values[p] = avg / 2.0;
    }

    delete_polygon_point_neighbours( polygons, n_neighbours, neighbours,
                                     NULL, NULL );

    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != VIO_OK )
        return( 1 );

    for_less( p, 0, polygons->n_points )
    {
        if( output_real( file, values[p] ) != VIO_OK ||
            output_newline( file ) != VIO_OK )
            return( 1 );
    }

    (void) close_file( file );

    return( 0 );
}
