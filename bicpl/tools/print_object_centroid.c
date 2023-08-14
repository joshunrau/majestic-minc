#include  <bicpl.h>

static  void  usage(
    VIO_STR   executable )
{
    VIO_STR  usage_str = "\n\
Usage: %s  input.obj\n\
\n\
     Prints the centroid of the object.\n\n";

    fprintf(stderr, usage_str, executable );
}

int  main(
    int    argc,
    char   *argv[] )
{
    VIO_STR              input_filename;
    int                 i, n_objects, n_points;
    VIO_File_formats        format;
    object_struct       **object_list;
    VIO_Point               *points, centroid;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( "", &input_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    if( input_graphics_file( input_filename, &format, &n_objects,
                             &object_list ) != VIO_OK )
        return( 1 );

    for_less( i, 0, n_objects )
    {
        n_points = get_object_points( object_list[i], &points );
        get_points_centroid( n_points, points, &centroid );
        print( "%g %g %g\n", RPoint_x( centroid ),
                             RPoint_y( centroid ),
                             RPoint_z( centroid ) );
    }

    delete_object_list( n_objects, object_list );

    return( 0 );
}
