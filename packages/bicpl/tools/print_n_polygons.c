#include <bicpl.h>

#define  BINTREE_FACTOR  0.3

static  void  usage(
    VIO_STR   executable )
{
    VIO_STR  usage_str = "\n\
Usage: %s  input.obj\n\
\n\
     Writes out the number of items in the object.\n\n";

    fprintf(stderr, usage_str, executable );
}

int  main(
    int    argc,
    char   *argv[] )
{
    VIO_STR           input_filename;
    int              n_objects;
    VIO_File_formats     format;
    object_struct    **object_list;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    if( input_graphics_file( input_filename, &format, &n_objects,
                             &object_list ) != VIO_OK )
        return( 1 );

    if( n_objects != 1 || get_object_type(object_list[0]) != POLYGONS )
    {
        fprintf(stderr, "File must contain exactly 1 polygons struct.\n" );
        return( 1 );
    }

    print( "%d\n", get_polygons_ptr( object_list[0] )->n_items );

    return( 0 );
}
