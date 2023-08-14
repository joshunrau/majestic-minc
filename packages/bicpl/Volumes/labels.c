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

#include  "bicpl_internal.h"

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_label_volume_real_range
@INPUT      : volume
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the real range to be equal to the voxel range.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI void  set_label_volume_real_range(
    VIO_Volume  volume )
{
    if( get_volume_data_type(volume) != VIO_FLOAT &&
        get_volume_data_type(volume) != VIO_DOUBLE )
    {
        set_volume_real_range( volume,
                               get_volume_voxel_min(volume),
                               get_volume_voxel_max(volume) );
    }
    else
        volume->real_range_set = FALSE;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_label_volume
@INPUT      : volume
              type
@OUTPUT     : 
@RETURNS    : a volume
@DESCRIPTION: Creates a label volume with the same tessellation as the volume.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI VIO_Volume  create_label_volume(
    VIO_Volume  volume,
    nc_type type )
{
    VIO_Volume   label_volume;

    if( type == NC_UNSPECIFIED )
        type = NC_BYTE;

    label_volume = copy_volume_definition_no_alloc( volume, type, FALSE,
                                                    0.0, -1.0 );

    set_label_volume_real_range( label_volume );

    return( label_volume );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_alloc_label_data
@INPUT      : volume
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Checks if the label data has been allocated.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static void  check_alloc_label_data(
    VIO_Volume  volume )
{
    if( !volume_is_alloced( volume ) && !volume_is_cached(volume) )
    {
        alloc_volume_data( volume );
        set_all_volume_label_data( volume, 0 );
    }
}

BICAPI VIO_BOOL  is_label_volume_initialized(
    VIO_Volume  volume )
{
    return( volume != NULL && 
            ((volume->is_cached_volume &&
             cached_volume_has_been_modified( &volume->cache )) ||
             (!volume->is_cached_volume &&
              multidim_array_is_alloced( &volume->array ))) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_all_volume_label_data
@INPUT      : volume
              value
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the label value of all voxels to the value specified.
@CREATED    : Mar   1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI void  set_all_volume_label_data(
    VIO_Volume    volume,
    int       value )
{
    VIO_Data_types      type;
    void            *ptr;
    VIO_Real            real_value;
    int             v0, v1, v2, v3, v4;
    unsigned int    n_voxels;

    check_alloc_label_data( volume );

    type = get_volume_data_type( volume );
    if( !volume->is_cached_volume && value == 0 &&
        type != VIO_FLOAT && type != VIO_DOUBLE )
    {
        GET_VOXEL_PTR( ptr, volume, 0, 0, 0, 0, 0 );
        n_voxels = get_volume_total_n_voxels( volume );
        (void)memset( ptr, 0, (size_t) n_voxels * (size_t) get_type_size(type));
    }
    else
    {
        real_value = (VIO_Real) value;
        BEGIN_ALL_VOXELS( volume, v0, v1, v2, v3, v4 )

            set_volume_real_value( volume, v0, v1, v2, v3, v4, real_value );

        END_ALL_VOXELS
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_volume_label_data
@INPUT      : volume
              x
              y
              z
              value
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the label data of the given voxel to the value.
@CREATED    : Mar   1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI void  set_volume_label_data_5d(
    VIO_Volume          volume,
    int             v0,
    int             v1,
    int             v2,
    int             v3,
    int             v4,
    int             value )
{
    check_alloc_label_data( volume );

    set_volume_real_value( volume, v0, v1, v2, v3, v4, (VIO_Real) value );
}

BICAPI void  set_volume_label_data(
    VIO_Volume          volume,
    int             voxel[],
    int             value )
{
    set_volume_label_data_5d( volume,
                              voxel[0], voxel[1], voxel[2], voxel[3], voxel[4],
                              value );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_volume_label_data
@INPUT      : volume
              x
              y
              z
@OUTPUT     : 
@RETURNS    : value
@DESCRIPTION: Returns the label data of the given voxel.
@CREATED    : Mar   1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI int  get_volume_label_data(
    VIO_Volume          volume,
    int             voxel[] )
{
    return( get_volume_label_data_5d( volume,
                       voxel[0], voxel[1], voxel[2], voxel[3], voxel[4] ) );
}

BICAPI int  get_volume_label_data_5d(
    VIO_Volume          volume,
    int             v0,
    int             v1,
    int             v2,
    int             v3,
    int             v4 )
{
    int    label;

    if( volume == (VIO_Volume) NULL || !volume_is_alloced( volume ) )
        return( 0 );
    else
    {
        label = (int) get_volume_real_value( volume, v0, v1, v2, v3, v4 );
        return( label );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_3D_volume_label_data
@INPUT      : volume
              x
              y
              z
@OUTPUT     : 
@RETURNS    : label
@DESCRIPTION: Gets the label value of a 3D label volume.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI int  get_3D_volume_label_data(
    VIO_Volume          volume,
    int             x,
    int             y,
    int             z )
{
    int    label;

    if( volume == (VIO_Volume) NULL || !volume_is_alloced( volume ) )
        return( 0 );
    else
    {
        label = (int) get_volume_real_value( volume, x, y, z, 0, 0 );
        return( label );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_voxel_label_bit
@INPUT      : volume
              voxel
              bit
@OUTPUT     : 
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Returns the label bit.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI VIO_BOOL  get_voxel_label_bit(
    VIO_Volume          volume,
    int             voxel[],
    int             bit )
{
    return( (get_volume_label_data( volume, voxel ) & bit) == 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_voxel_label_bit
@INPUT      : volume
              voxel
              bit
              value
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the voxel label bit.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI void  set_voxel_label_bit(
    VIO_Volume          volume,
    int             voxel[],
    int             bit,
    VIO_BOOL         value )
{
    int     i, n_dims, v[VIO_MAX_DIMENSIONS], label, anded, new_label;

    check_alloc_label_data( volume );

    n_dims = get_volume_n_dimensions(volume);

    for_less( i, 0, n_dims )
        v[i] = voxel[i];

    label = (int) get_volume_real_value( volume, v[0], v[1], v[2], v[3], v[4]);

    anded = (label & bit);

    if( value )
    {
        if( anded != bit )
        {
            new_label = label | bit;
            set_volume_real_value( volume, v[0], v[1], v[2], v[3], v[4],
                                    (VIO_Real) new_label );
        }
    }
    else if( anded != 0 )
    {
        new_label = label & (~bit);
        set_volume_real_value( volume, v[0], v[1], v[2], v[3], v[4],
                                (VIO_Real) new_label );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_all_volume_label_data_bit
@INPUT      : volume
              bit
              value - TRUE or FALSE
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets just the given bit of all the voxels' label data to the
              given value.
@CREATED    : Mar   1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI void  set_all_volume_label_data_bit(
    VIO_Volume         volume,
    int            bit,
    VIO_BOOL        value )
{
    int             v[VIO_MAX_DIMENSIONS];

    check_alloc_label_data( volume );

    BEGIN_ALL_VOXELS( volume, v[0], v[1], v[2], v[3], v[4] )

        set_voxel_label_bit( volume, v, bit, value );

    END_ALL_VOXELS
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_volume_voxel_activity
@INPUT      : volume
              voxel
              activity_if_mixed
@OUTPUT     : 
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Returns the voxel activity, by looking at the 2^d corners of the
              voxel.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI VIO_BOOL  get_volume_voxel_activity(
    VIO_Volume     volume,
    VIO_Real       voxel[],
    VIO_BOOL    activity_if_mixed )
{
    VIO_BOOL  active_found, inactive_found;
    int      c, int_index[VIO_MAX_DIMENSIONS], ind[VIO_MAX_DIMENSIONS];
    int      n[VIO_MAX_DIMENSIONS], sizes[VIO_MAX_DIMENSIONS];

    if( volume == (VIO_Volume) NULL || !volume_is_alloced( volume ) )
        return( TRUE );

    get_volume_sizes( volume, sizes );

    for_less( c, 0, get_volume_n_dimensions(volume) )
        if( voxel[c] < 0.0 || voxel[c] > (VIO_Real) sizes[c]-1.0 )
            return( FALSE );

    for_less( c, 0, get_volume_n_dimensions(volume) )
    {
        int_index[c] = (int) voxel[c];
        if( int_index[c] == sizes[c] - 1 )
            int_index[c] = sizes[c] - 2;
        n[c] = 2;
    }

    for_less( c, get_volume_n_dimensions(volume), VIO_MAX_DIMENSIONS )
    {
        n[c] = 1;
        int_index[c] = 0;
    }

    active_found = FALSE;
    inactive_found = FALSE;

    for_less( ind[VIO_X], int_index[VIO_X], int_index[VIO_X] + n[VIO_X] )
    for_less( ind[VIO_Y], int_index[VIO_Y], int_index[VIO_Y] + n[VIO_Y] )
    for_less( ind[VIO_Z], int_index[VIO_Z], int_index[VIO_Z] + n[VIO_Z] )
    for_less( ind[3], int_index[3], int_index[3] + n[3] )
    for_less( ind[4], int_index[4], int_index[4] + n[4] )
    {
        if( get_volume_label_data( volume, ind ) == 0 )
        {
            if( inactive_found )
                return( activity_if_mixed );
            active_found = TRUE;
        }
        else
        {
            if( active_found )
                return( activity_if_mixed );
            inactive_found = TRUE;
        }
    }

    if( active_found && !inactive_found )
        return( TRUE );
    else if( !active_found && inactive_found )
        return( FALSE );
    else
        return( activity_if_mixed );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_input_volume_label_limits
@INPUT      : volume1
              volume2
              slice
@OUTPUT     : limits
@RETURNS    : 
@DESCRIPTION: Computes the range of overlap of volume2 in volume1.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Aug. 1, 1995    David MacDonald
@MODIFIED   : May. 6, 1997    D. MacDonald     -  one slice of 3d volume
---------------------------------------------------------------------------- */

static void  get_input_volume_label_limits(
    VIO_Volume   volume1,
    VIO_Volume   volume2,
    int      slice,
    int      limits[2][VIO_N_DIMENSIONS] )
{
    int       sizes1[VIO_MAX_DIMENSIONS], sizes2[VIO_MAX_DIMENSIONS];
    int       d, t[VIO_N_DIMENSIONS];
    VIO_Real  voxel1[VIO_MAX_DIMENSIONS], voxel2[VIO_MAX_DIMENSIONS];
    int       pos;
    VIO_Real      xw, yw, zw;
    VIO_BOOL   first;

    get_volume_sizes( volume1, sizes1 );
    get_volume_sizes( volume2, sizes2 );

    first = TRUE;

    for_less( t[0], 0, 2 )
    for_less( t[1], 0, 2 )
    for_less( t[2], 0, 2 )
    {
        voxel2[0] = (VIO_Real) slice - 0.5 + (VIO_Real) t[0];
        voxel2[1] = -0.5 + (VIO_Real) t[1] * (VIO_Real) sizes2[1];
        voxel2[2] = -0.5 + (VIO_Real) t[2] * (VIO_Real) sizes2[2];

        convert_voxel_to_world( volume2, voxel2, &xw, &yw, &zw );
        convert_world_to_voxel( volume1, xw, yw, zw, voxel1 );

        for_less( d, 0, VIO_N_DIMENSIONS )
        {
            pos = VIO_FLOOR( voxel1[d] + 0.5 );

            if( first )
            {
                limits[0][d] = pos;
                limits[1][d] = pos;
            }
            else
            {
                if( pos < limits[0][d] )
                    limits[0][d] = pos;
                else if( pos > limits[1][d] )
                    limits[1][d] = pos;
            }
        }

        first = FALSE;
    }

    for_less( d, 0, VIO_N_DIMENSIONS )
    {
        if( limits[0][d] < 0 )
            limits[0][d] = 0;
        if( limits[1][d] >= sizes1[d] )
            limits[1][d] = sizes1[d] - 1;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_label_volume
@INPUT      : Filename
@OUTPUT     : label_volume
@RETURNS    : VIO_ERROR or VIO_OK
@DESCRIPTION: Loads the label volume.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : Aug. 1, 1995    D. MacDonald  -  loaded volume no longer needs
                                               to be same grid.
@MODIFIED   : May. 6, 1997    D. MacDonald  -  now loads one slice at a time,
                                               and does not overwrite labels
                                               with 0's.
@MODIFIED   : Jun. 6, 1997    D. MacDonald  -  special case, fast, for matching
                                               coordinate systems
---------------------------------------------------------------------------- */

BICAPI VIO_Status  load_label_volume(
    VIO_STR   filename,
    VIO_Volume   label_volume )
{
    int                   slice, n_slices;
    int                   label_voxel[VIO_MAX_DIMENSIONS];
    int                   int_voxel[VIO_N_DIMENSIONS], int_file_value;
    int                   limits[2][VIO_N_DIMENSIONS];
    int                   file_sizes[VIO_N_DIMENSIONS];
    VIO_Real              file_value, xw, yw, zw;
    int                   int_y_voxel[VIO_N_DIMENSIONS];
    VIO_Real                  y_voxel[VIO_N_DIMENSIONS];
    VIO_Real                  voxel[VIO_N_DIMENSIONS], z_dx, z_dy, z_dz;
    VIO_Real                  y_dx, y_dy, y_dz;
    int                   int_z_dx, int_z_dy, int_z_dz;
    int                   int_y_dx, int_y_dy, int_y_dz;
    VIO_Real                  start_voxel[VIO_N_DIMENSIONS];
    VIO_Real                  end_voxel[VIO_N_DIMENSIONS];
    VIO_Volume            file_volume;
    VIO_BOOL               is_linear, is_integer_step;
    VIO_progress_struct       progress;
    static VIO_STR         file_order_dim_names[VIO_MAX_DIMENSIONS] = {
      "", "", "", "", ""
    };
    int                   i;
    minc_input_options    options;

    for_less(i, 0, VIO_MAX_DIMENSIONS)
        label_voxel[i] = 0;       /* Initialize. */

    check_alloc_label_data( label_volume );

    set_default_minc_input_options( &options );
    set_minc_input_vector_to_scalar_flag( &options, FALSE );

    if (input_volume( filename, VIO_N_DIMENSIONS, file_order_dim_names,
                      NC_UNSPECIFIED, FALSE, 0.0, 0.0, TRUE,
                      &file_volume, &options ) != VIO_OK)
        return (VIO_ERROR);

    get_volume_sizes( file_volume, file_sizes );

    n_slices = file_sizes[0];
    /*--- check if the voxel1-to-voxel2 transform is linear */

    is_linear = get_transform_type(
                get_voxel_to_world_transform(label_volume) ) == LINEAR &&
                get_transform_type(
                get_voxel_to_world_transform(file_volume) ) == LINEAR;
    is_integer_step = FALSE;

    if( is_linear )
    {
        convert_3D_voxel_to_world( label_volume, 0.0, 0.0, 0.0,
                                   &xw, &yw, &zw );
        convert_world_to_voxel( file_volume, xw, yw, zw, start_voxel );

        convert_3D_voxel_to_world( label_volume, 0.0, 0.0, 1.0,
                                   &xw, &yw, &zw );
        convert_world_to_voxel( file_volume, xw, yw, zw, end_voxel );

        z_dx = end_voxel[0] - start_voxel[0];
        z_dy = end_voxel[1] - start_voxel[1];
        z_dz = end_voxel[2] - start_voxel[2];

        convert_3D_voxel_to_world( label_volume, 0.0, 1.0, 0.0,
                                   &xw, &yw, &zw );
        convert_world_to_voxel( file_volume, xw, yw, zw, end_voxel );

        y_dx = end_voxel[0] - start_voxel[0];
        y_dy = end_voxel[1] - start_voxel[1];
        y_dz = end_voxel[2] - start_voxel[2];

        is_integer_step = FALSE;
        if( VIO_IS_INT(z_dx) && VIO_IS_INT(z_dy) && VIO_IS_INT(z_dz) &&
            VIO_IS_INT(y_dx) && VIO_IS_INT(y_dy) && VIO_IS_INT(y_dz) )
        {
            is_integer_step = TRUE;
            int_z_dx = (int) z_dx;
            int_z_dy = (int) z_dy;
            int_z_dz = (int) z_dz;
            int_y_dx = (int) y_dx;
            int_y_dy = (int) y_dy;
            int_y_dz = (int) y_dz;
        }
    }

    /*--- input label slices */

    initialize_progress_report( &progress, FALSE, n_slices,
                                "Reading Labels" );

    for_less( slice, 0, n_slices )
    {
        get_input_volume_label_limits( label_volume, file_volume,
                                       slice, limits );
        for_inclusive( label_voxel[VIO_X], limits[0][VIO_X], limits[1][VIO_X] )
        {
            for_inclusive( label_voxel[VIO_Y], limits[0][VIO_Y], limits[1][VIO_Y] )
            {
                if( is_linear )
                {
                    if( label_voxel[VIO_Y] == limits[0][VIO_Y] )
                    {
                        convert_3D_voxel_to_world( label_volume,
                                                   (VIO_Real) label_voxel[VIO_X],
                                                   (VIO_Real) label_voxel[VIO_Y],
                                                   (VIO_Real) limits[0][VIO_Z],
                                                   &xw, &yw, &zw );
                        convert_world_to_voxel( file_volume, xw, yw, zw,
                                                y_voxel );
                        y_voxel[VIO_X] += 0.5;
                        y_voxel[VIO_Y] += 0.5;
                        y_voxel[VIO_Z] += 0.5;

                        if( is_integer_step )
                        {
                            int_y_voxel[VIO_X] = VIO_FLOOR( y_voxel[VIO_X] );
                            int_y_voxel[VIO_Y] = VIO_FLOOR( y_voxel[VIO_Y] );
                            int_y_voxel[VIO_Z] = VIO_FLOOR( y_voxel[VIO_Z] );
                        }
                    }
                    else
                    {
                        if( is_integer_step )
                        {
                            int_y_voxel[VIO_X] += int_y_dx;
                            int_y_voxel[VIO_Y] += int_y_dy;
                            int_y_voxel[VIO_Z] += int_y_dz;
                        }
                        else
                        {
                            y_voxel[VIO_X] += y_dx;
                            y_voxel[VIO_Y] += y_dy;
                            y_voxel[VIO_Z] += y_dz;
                        }
                    }

                    if( is_integer_step )
                    {
                        int_voxel[VIO_X] = int_y_voxel[VIO_X];
                        int_voxel[VIO_Y] = int_y_voxel[VIO_Y];
                        int_voxel[VIO_Z] = int_y_voxel[VIO_Z];
                    }
                    else
                    {
                        voxel[VIO_X] = y_voxel[VIO_X];
                        voxel[VIO_Y] = y_voxel[VIO_Y];
                        voxel[VIO_Z] = y_voxel[VIO_Z];
                    }
                }

                for_inclusive( label_voxel[VIO_Z], limits[0][VIO_Z], limits[1][VIO_Z] )
                {
                    if( !is_linear )
                    {
                        convert_3D_voxel_to_world( label_volume,
                                                   (VIO_Real) label_voxel[VIO_X],
                                                   (VIO_Real) label_voxel[VIO_Y],
                                                   (VIO_Real) label_voxel[VIO_Z],
                                                   &xw, &yw, &zw );
                        convert_world_to_voxel( file_volume, xw, yw, zw,
                                                voxel );
                        voxel[VIO_X] += 0.5;
                        voxel[VIO_Y] += 0.5;
                        voxel[VIO_Z] += 0.5;
                    }
                    else if( label_voxel[VIO_Z] != limits[0][VIO_Z] )
                    {
                        if( is_integer_step )
                        {
                            int_voxel[VIO_X] += int_z_dx;
                            int_voxel[VIO_Y] += int_z_dy;
                            int_voxel[VIO_Z] += int_z_dz;
                        }
                        else
                        {
                            voxel[VIO_X] += z_dx;
                            voxel[VIO_Y] += z_dy;
                            voxel[VIO_Z] += z_dz;
                        }
                    }

                    if( !is_integer_step )
                    {
                        int_voxel[VIO_X] = VIO_FLOOR( voxel[VIO_X] );
                        int_voxel[VIO_Y] = VIO_FLOOR( voxel[VIO_Y] );
                        int_voxel[VIO_Z] = VIO_FLOOR( voxel[VIO_Z] );
                    }

                    if( int_voxel[VIO_X] == slice &&
                        int_voxel[VIO_Y] >= 0 && int_voxel[VIO_Y] < file_sizes[VIO_Y] &&
                        int_voxel[VIO_Z] >= 0 && int_voxel[VIO_Z] < file_sizes[VIO_Z] )
                    {
                        file_value = get_volume_real_value( file_volume,
                                        int_voxel[VIO_X], int_voxel[VIO_Y], int_voxel[VIO_Z], 0, 0 );

                        if( file_value > 0.0 )
                        {
                            int_file_value = VIO_ROUND( file_value );
                            set_volume_label_data( label_volume, label_voxel,
                                                   int_file_value );
                        }
                    }
                }
            }
        }
        update_progress_report( &progress, slice+1 );
    }

    terminate_progress_report( &progress );

    delete_volume( file_volume );

    return( VIO_OK );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : save_label_volume
@INPUT      : filename
              original_filename
              label_volume
@OUTPUT     : 
@RETURNS    : VIO_OK or VIO_ERROR
@DESCRIPTION: Saves the label volume.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : Aug. 1, 1995    D. MacDonald    - crops the volume on output
---------------------------------------------------------------------------- */

BICAPI VIO_Status  save_label_volume(
    VIO_STR   filename,
    VIO_STR   original_filename,
    VIO_Volume   label_volume,
    VIO_Real     crop_threshold )
{
    VIO_Status   status;
    VIO_BOOL  cropping;
    VIO_Volume   cropped_volume;
    int      n_voxels, n_voxels_cropped;
    int      c, limits[2][VIO_MAX_DIMENSIONS], sizes[VIO_MAX_DIMENSIONS];

    check_alloc_label_data( label_volume );

    if( crop_threshold > 0.0 )
    {
        if( find_volume_crop_bounds( label_volume, -1.0, 0.5, limits ) )
        {
            get_volume_sizes( label_volume, sizes );
            n_voxels = 1;
            n_voxels_cropped = 1;
            for_less( c, 0, VIO_N_DIMENSIONS )
            {
                n_voxels *= sizes[c];
                n_voxels_cropped *= limits[1][c] - limits[0][c] + 1;
            }

            cropping = ((VIO_Real) n_voxels_cropped / (VIO_Real) n_voxels <
                        crop_threshold);
        }
        else
        {
            for_less( c, 0, VIO_N_DIMENSIONS )
            {
                limits[0][c] = 0;
                limits[1][c] = 0;
            }
            cropping = TRUE;
        }
    }
    else
        cropping = FALSE;

    if( cropping )
        cropped_volume = create_cropped_volume( label_volume, limits );
    else
        cropped_volume = label_volume;

    if( original_filename != NULL )
    {
        status = output_modified_volume( filename,
                                NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                                cropped_volume, original_filename,
                                "Label volume\n", NULL );
    }
    else
    {
        status = output_volume( filename, NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                                cropped_volume, "Label volume\n", NULL );
    }

    if( cropping )
        delete_volume( cropped_volume );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_tags_as_labels
@INPUT      : file
              volume
@OUTPUT     : label_volume
@RETURNS    : VIO_OK or VIO_ERROR
@DESCRIPTION: Inputs a tag file into a label volume.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993             David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI VIO_Status  input_tags_as_labels(
    FILE    *file,
    VIO_Volume  volume,
    VIO_Volume  label_volume )
{
    VIO_Status          status;
    int             c, label, ind[VIO_MAX_DIMENSIONS];
    VIO_Real            voxel[VIO_MAX_DIMENSIONS];
    int             n_volumes;
    VIO_Real            tag1[VIO_N_DIMENSIONS];
    int             structure_id;
    VIO_Real            min_label, max_label;

    check_alloc_label_data( label_volume );

    get_volume_real_range( label_volume, &min_label, &max_label );

    status = initialize_tag_file_input( file, &n_volumes );

    while( status == VIO_OK &&
           input_one_tag( file, n_volumes,
                          tag1, NULL, NULL, &structure_id, NULL, NULL,
                          &status ) )
    {
        convert_world_to_voxel( volume, tag1[VIO_X], tag1[VIO_Y], tag1[VIO_Z], voxel );

        for_less( c, 0, get_volume_n_dimensions(volume) )
        {
            ind[c] = VIO_ROUND( voxel[c] );
        }

        label = structure_id;
        if( (VIO_Real) label >= min_label && (VIO_Real) label <= max_label &&
            int_voxel_is_within_volume( volume, ind ) )
        {
            set_volume_label_data( label_volume, ind, label);
        }
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_label_volume_from_file
@INPUT      : filename
              volume
@OUTPUT     : label_volume
@RETURNS    : VIO_OK or VIO_ERROR
@DESCRIPTION: Creates a label volume for the given volume from a tag file or
              a minc file.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : Dec.  8, 1995    David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

BICAPI VIO_Status  create_label_volume_from_file(
    VIO_STR   filename,
    VIO_Volume   volume,
    VIO_Volume   *label_volume )
{
    VIO_Status  status;
    VIO_STR  *dim_names;
    VIO_Volume  file_volume;
    FILE    *file;
    VIO_BOOL same_grid;

    status = VIO_OK;

    if( filename_extension_matches( filename, "mnc" ) )
    {
        dim_names = get_volume_dimension_names( volume );
        status = input_volume_header_only( filename,
                                           get_volume_n_dimensions(volume),
                                           dim_names, &file_volume, NULL );

        if( status != VIO_OK )
        {
            delete_dimension_names( volume, dim_names );
            return( status );
        }

        same_grid = volumes_are_same_grid( volume, file_volume );
        delete_volume( file_volume );

        if( same_grid )
        {
            status = input_volume( filename, get_volume_n_dimensions(volume),
                                   dim_names, NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                                   TRUE, label_volume, NULL );
        }
        else
        {
            *label_volume = create_label_volume( volume, NC_UNSPECIFIED );
            status = load_label_volume( filename, *label_volume );
        }

        delete_dimension_names( volume, dim_names );
    }
    else
    {
        *label_volume = create_label_volume( volume, NC_UNSPECIFIED );

        if( open_file( filename, READ_FILE, ASCII_FORMAT, &file )!=VIO_OK)
            return( VIO_ERROR );

        if( input_tags_as_labels( file, volume, *label_volume ) != VIO_OK )
            return( VIO_ERROR );

        (void) close_file( file );
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_labels_as_tags
@INPUT      : file
              volume
              label_volume
              desired_label
              size
              patient_id
@OUTPUT     : 
@RETURNS    : VIO_OK or VIO_ERROR
@DESCRIPTION: Outputs a set of labels as tags.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : Oct. 19, 1995   D. MacDonald - now writes tags 1 at a time to
                                             be more memory efficient.
---------------------------------------------------------------------------- */

BICAPI VIO_Status  output_labels_as_tags(
    FILE    *file,
    VIO_Volume  volume,
    VIO_Volume  label_volume,
    int     desired_label,
    VIO_Real    size,
    int     patient_id )
{
    int             ind[VIO_MAX_DIMENSIONS];
    int             label, sizes[VIO_MAX_DIMENSIONS];
    VIO_Real        real_ind[VIO_N_DIMENSIONS];
    VIO_Real        tags[VIO_N_DIMENSIONS];
    int             n_tags;

    if( get_volume_n_dimensions(volume) != 3 )
    {
        print_error( "output_labels_as_tags:  volume must be 3D\n" );
        return( VIO_ERROR );
    }

    check_alloc_label_data( label_volume );
    get_volume_sizes( label_volume, sizes );

    n_tags = 0;

    for_less( ind[VIO_X], 0, sizes[VIO_X] )
    {
        real_ind[VIO_X] = (VIO_Real) ind[VIO_X];
        for_less( ind[VIO_Y], 0, sizes[VIO_Y] )
        {
            real_ind[VIO_Y] = (VIO_Real) ind[VIO_Y];
            for_less( ind[VIO_Z], 0, sizes[VIO_Z] )
            {
                real_ind[VIO_Z] = (VIO_Real) ind[VIO_Z];
                label = get_volume_label_data( label_volume, ind );

                if( label == desired_label || (desired_label < 0 && label > 0) )
                {
                    convert_voxel_to_world( volume, real_ind,
                                            &tags[VIO_X], &tags[VIO_Y], &tags[VIO_Z] );

                    if( n_tags == 0 &&
                        initialize_tag_file_output( file, NULL, 1 ) != VIO_OK )
                        return( VIO_ERROR );

                    if( output_one_tag( file, 1, tags, NULL,
                               &size, &label, &patient_id, NULL ) != VIO_OK )
                        return( VIO_ERROR );

                    ++n_tags;
                }
            }
        }
    }

    if( n_tags > 0 )
        terminate_tag_file_output( file );

    return( VIO_OK );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_landmarks_as_labels
@INPUT      : file
              volume
@OUTPUT     : label_volume
@RETURNS    : VIO_OK or VIO_ERROR
@DESCRIPTION: Loads a set of landmarks into the label_volume.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : Oct. 19, 1995   D. MacDonald - now reads tags 1 at a time to
                                             be more memory efficient.
---------------------------------------------------------------------------- */

BICAPI VIO_Status  input_landmarks_as_labels(
    FILE    *file,
    VIO_Volume  volume,
    VIO_Volume  label_volume )
{
    int             c, label, ind[VIO_MAX_DIMENSIONS];
    VIO_Real            voxel[VIO_MAX_DIMENSIONS];
    marker_struct   marker;
    VIO_Real            min_label, max_label;

    check_alloc_label_data( label_volume );

    get_volume_real_range( label_volume, &min_label, &max_label );

    while( io_tag_point( file, READ_FILE, volume, 1.0, &marker ) == VIO_OK )
    {
        convert_world_to_voxel( volume,
                                (VIO_Real) Point_x(marker.position),
                                (VIO_Real) Point_y(marker.position),
                                (VIO_Real) Point_z(marker.position), voxel );

        for_less( c, 0, get_volume_n_dimensions(volume) )
            ind[c] = VIO_ROUND( voxel[c] );

        label = marker.structure_id;
        if( (VIO_Real) label >= min_label && (VIO_Real) label <= max_label &&
            int_voxel_is_within_volume( volume, ind ) )
        {
            set_volume_label_data( label_volume, ind, label );
        }
    }

    return( VIO_OK );
}
