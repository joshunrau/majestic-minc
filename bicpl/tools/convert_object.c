/* Convert a surface .obj file from ASCII to BINARY, or vice-versa.
   Simply toggle the format of the input to create the output. 

   Claude Lepage    May 12, 2006 - Initial implementation  
 */

#include  <bicpl.h>

int  main( int    argc, char   *argv[] ) {

    VIO_Status           status;
    VIO_STR           input_filename, output_filename;
    int              n_objects;
    VIO_File_formats     format;
    object_struct    **object_list;

    status = VIO_OK;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ) ) {
        fprintf(stderr, "Toggle format of a .obj file from ASCII to BINARY, or vice-versa\n" );
        fprintf(stderr, "Usage: %s input.obj output.obj\n", argv[0] );
        return( 1 );
    }

    if( input_graphics_file( input_filename, &format, &n_objects, &object_list ) != VIO_OK ) {
        print( "Couldn't read input object %s.\n", input_filename );
        return( 1 );
    }

    if( format == ASCII_FORMAT ) {
      printf( "Converting object from ASCII to BINARY...\n" );
      format = BINARY_FORMAT;
    } else {
      printf( "Converting object from BINARY to ASCII...\n" );
      format = ASCII_FORMAT;
    }

    status = output_graphics_file( output_filename, format, n_objects, object_list );

    return( status != VIO_OK );
}

