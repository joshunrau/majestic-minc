/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include "bicpl_internal.h"

static  VIO_Status  output_texture_values_ascii(
    VIO_STR   filename,
    int      n_values,
    VIO_Real     values[] )
{
    int      v;
    VIO_Status   status;
    FILE     *file;

    status = open_file( filename, WRITE_FILE, ASCII_FORMAT, &file );

    if( status != VIO_OK )
        return( status );

    for_less( v, 0, n_values )
    {
        if( output_real( file, values[v] ) != VIO_OK ||
            output_newline( file ) != VIO_OK )
        {
            print_error( "Error outputting %d'th value out of %d to file %s\n",
                         v+1, n_values, filename );
            return( VIO_ERROR );
        }
    }

    (void) close_file( file );

    return( VIO_OK );
}

static  VIO_Status  input_texture_values_ascii(
    VIO_STR   filename,
    int      *n_values,
    VIO_Real     *values[] )
{
    VIO_Status   status;
    FILE     *file;
    VIO_Real     value;

    status = open_file( filename, READ_FILE, ASCII_FORMAT, &file );

    if( status != VIO_OK )
        return( status );

    *n_values = 0;
    *values = NULL;

    while( input_real( file, &value ) == VIO_OK )
    {
        ADD_ELEMENT_TO_ARRAY( *values, *n_values, value, DEFAULT_CHUNK_SIZE );
    }

    (void) close_file( file );

    return( VIO_OK );
}

static  VIO_Status  output_texture_values_binary(
    VIO_STR   filename,
    int      n_values,
    VIO_Real     values[] )
{
    int      v, sizes[2];
    VIO_Status   status;
    VIO_Volume   volume;
    VIO_STR   dim_names[] = { MIxspace, MIyspace };

    volume = create_volume( 2, dim_names, NC_FLOAT, FALSE, 0.0, 0.0 );
    sizes[0] = 1;
    sizes[1] = n_values;
    set_volume_sizes( volume, sizes );
    alloc_volume_data( volume );

    for_less( v, 0, n_values )
        set_volume_real_value( volume, 0, v, 0, 0, 0, values[v] );

    status = output_volume( filename, NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                            volume, "Texture values.\n", NULL );

    delete_volume( volume );

    return( status );
}

static  VIO_Status  input_texture_values_binary(
    VIO_STR   filename,
    int      *n_values,
    VIO_Real     *values[] )
{
    int      v, sizes[2];
    VIO_Status   status;
    VIO_Volume   volume;
    VIO_STR   dim_names[] = { MIxspace, MIyspace };

    status = input_volume( filename, 2, dim_names,
                           NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                           TRUE, &volume, NULL );

    if( status != VIO_OK )
        return( status );

    get_volume_sizes( volume, sizes );

    *n_values = sizes[1];

    ALLOC( *values, *n_values );

    for_less( v, 0, *n_values )
        (*values)[v] = get_volume_real_value( volume, 0, v, 0, 0, 0 );

    delete_volume( volume );

    return( VIO_OK );
}

/*!
 * \brief Write a set of real values to file.
 *
 * The file will be created, if necessary.
 * BINARY_FORMAT files will not be portable.
 */
BICAPI  VIO_Status  output_texture_values(
    VIO_STR         filename,
    VIO_File_formats   format,
    int            n_values,
    VIO_Real           values[] )
{
    VIO_Status   status;

    if( format == ASCII_FORMAT )
        status = output_texture_values_ascii( filename, n_values, values );
    else
        status = output_texture_values_binary( filename, n_values, values );

    return( status );
}

/*!
 * \brief Read a set of real values from file.
 *
 * Reading a BINARY_FORMAT file that was created
 * on another system will produce undefined behaviour.
 */
BICAPI  VIO_Status  input_texture_values(
    VIO_STR         filename,
    int            *n_values,
    VIO_Real           *values[] )
{
    VIO_Status         status;

    if( filename_extension_matches( filename, MNC_ENDING ) )
        status = input_texture_values_binary( filename, n_values, values );
    else
        status = input_texture_values_ascii( filename, n_values, values );

    return( status );
}
