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

#define  DEBUG
#undef  DEBUG

#ifdef DEBUG
static VIO_Real  get_correct_amount(
    VIO_Volume   volume,
    int      x,
    int      y,
    int      z,
    VIO_Real     x_width,
    VIO_Real     y_width,
    VIO_Real     z_width );
#endif

#define  INITIALIZE_SAMPLE( half_width, receding, advancing, left_weight, right_weight ) \
    { \
        VIO_Real  min_pos, max_pos; \
 \
        min_pos = -(half_width); \
        max_pos = (half_width); \
 \
        (receding) = VIO_ROUND( min_pos ); \
        (advancing) = VIO_ROUND( max_pos ); \
 \
        (left_weight) = (VIO_Real) (receding) + 0.5 - min_pos; \
        if( (left_weight) >= 1.0 ) \
            --(advancing); \
        (right_weight) = 1.0 - (left_weight); \
    }

#define  GET_FIRST_SAMPLE( advancing, left_weight, size, data, sample ) \
    { \
        int   _I; \
 \
        (sample) = 0.0; \
        for_inclusive( _I, 0, MIN( (size)-1, (advancing)-1) ) \
            (sample) += (data); \
 \
        if( (advancing) < (size) ) \
            (sample) += (data) * (left_weight); \
    }

#define  GET_NEXT_SAMPLE( receding, advancing, left_weight, right_weight, size, data, sample ) \
    { \
        int   _I; \
 \
        if( (advancing) < (size) ) \
        { \
            _I = (advancing); \
            (sample) += (data) * (right_weight); \
        } \
        if( (advancing) < (size)-1 ) \
        { \
            _I = (advancing) + 1; \
            (sample) += (data) * (left_weight); \
        } \
\
        if( (receding) >= 0 ) \
        { \
            _I = (receding); \
            (sample) -= (data) * (left_weight); \
        } \
        if( (receding) >= -1 ) \
        { \
            _I = (receding) + 1; \
            (sample) -= (data) * (right_weight); \
        } \
    }

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_box_filtered_volume
@INPUT      : volume           - volume to filter
              nc_data_type     - if NC_UNSPECIFIED, this and next 3 args ignored
              sign_flag
              real_min_value
              real_max_value
              x_width          - full width of box filter
              y_width
              z_width
@OUTPUT     : 
@RETURNS    : filtered volume
@DESCRIPTION: Filters a volume with a box filter, creating a new volume with
              the same number of samples per dimension.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI VIO_Volume  create_box_filtered_volume(
    VIO_Volume   volume,
    nc_type  nc_data_type,
    VIO_BOOL  sign_flag,
    VIO_Real     real_min_value,
    VIO_Real     real_max_value,
    VIO_Real     x_width,
    VIO_Real     y_width,
    VIO_Real     z_width )
{
    VIO_Volume             resampled_volume;
    int                x, y, z;
    int                sizes[VIO_MAX_DIMENSIONS], nx_required;
#ifdef DEBUG
    VIO_Real               correct_voxel;
#endif
    VIO_Real               total_volume, value, sum, sample;
    VIO_Real               *volume_row, ***volume_cache;
    int                n_alloced, max_alloced;
    int                x_init_recede_index, y_init_recede_index;
    int                z_init_recede_index;
    int                x_init_advance_index, y_init_advance_index;
    int                z_init_advance_index;
    int                x_recede_index, y_recede_index, z_recede_index;
    int                x_advance_index, y_advance_index, z_advance_index;
    VIO_Real               x_left_weight, x_right_weight;
    VIO_Real               y_left_weight, y_right_weight;
    VIO_Real               z_left_weight, z_right_weight;
    float              **slice, *row;
    VIO_progress_struct    progress;

    if( get_volume_n_dimensions(volume) != 3 )
    {
        handle_internal_error(
           "create_box_filtered_volume: volume must be 3D.\n" );
    }

    x_width = VIO_FABS( x_width );
    y_width = VIO_FABS( y_width );
    z_width = VIO_FABS( z_width );

    if( x_width < 1.0 )
        x_width = 1.0;
    if( y_width < 1.0 )
        y_width = 1.0;
    if( z_width < 1.0 )
        z_width = 1.0;

    get_volume_sizes( volume, sizes );

    resampled_volume = copy_volume_definition( volume, nc_data_type, sign_flag,
                                               real_min_value, real_max_value );

    total_volume = x_width * y_width * z_width;

    x_width /= 2.0;
    y_width /= 2.0;
    z_width /= 2.0;

    INITIALIZE_SAMPLE( x_width, x_init_recede_index, x_init_advance_index,
                       x_left_weight, x_right_weight )
    INITIALIZE_SAMPLE( y_width, y_init_recede_index, y_init_advance_index,
                       y_left_weight, y_right_weight )
    INITIALIZE_SAMPLE( z_width, z_init_recede_index, z_init_advance_index,
                       z_left_weight, z_right_weight )

    initialize_progress_report( &progress, FALSE, sizes[VIO_X] * sizes[VIO_Y],
                                "Box Filtering" );

    VIO_ALLOC2D( slice, sizes[VIO_Y], sizes[VIO_Z] );
    ALLOC( row, sizes[VIO_Z] );

    nx_required = MAX( sizes[VIO_X], x_init_advance_index + 1 );
    ALLOC( volume_row, nx_required );

    for_less( y, 0, sizes[VIO_Y] )
    {
        for_less( z, 0, sizes[VIO_Z] )
        {
            get_volume_value_hyperslab_3d( volume, 0, y, z, nx_required, 1, 1,
                                           volume_row );
            GET_FIRST_SAMPLE( x_init_advance_index,
                              x_left_weight, sizes[VIO_X],
                              volume_row[_I], sample );
            slice[y][z] = (float) sample;
        }
    }

    FREE( volume_row );

    ALLOC( volume_cache, sizes[VIO_X] );
    for_less( x, 0, sizes[VIO_X] )
        volume_cache[x] = NULL;
    n_alloced = 0;
    VIO_ALLOC2D( volume_cache[0], sizes[VIO_Y], sizes[VIO_Z] );
    get_volume_value_hyperslab_3d( volume, 0, 0, 0, 1, sizes[VIO_Y], sizes[VIO_Z],
                                   &volume_cache[0][0][0] );
    ++n_alloced;
    if( x_init_recede_index >= 0 )
    {
        VIO_ALLOC2D( volume_cache[1], sizes[VIO_Y], sizes[VIO_Z] );
        get_volume_value_hyperslab_3d( volume, 1, 0, 0, 1, sizes[VIO_Y], sizes[VIO_Z],
                                       &volume_cache[1][0][0] );
        ++n_alloced;
    }

    if( x_init_advance_index < sizes[VIO_X] )
    {
        VIO_ALLOC2D( volume_cache[x_init_advance_index], sizes[VIO_Y], sizes[VIO_Z] );
        get_volume_value_hyperslab_3d( volume, x_init_advance_index, 0, 0,
                                       1, sizes[VIO_Y], sizes[VIO_Z],
                                    &volume_cache[x_init_advance_index][0][0] );
        ++n_alloced;
    }

    if( x_init_advance_index+1 < sizes[VIO_X] )
    {
        VIO_ALLOC2D( volume_cache[x_init_advance_index+1], sizes[VIO_Y], sizes[VIO_Z] );
        get_volume_value_hyperslab_3d( volume, x_init_advance_index+1, 0, 0,
                                 1, sizes[VIO_Y], sizes[VIO_Z],
                                 &volume_cache[x_init_advance_index+1][0][0] );
        ++n_alloced;
    }

    if( x_init_recede_index == x_init_advance_index )
        max_alloced = 2;
    else if( x_init_recede_index+1 == x_init_advance_index )
        max_alloced = 3;
    else
        max_alloced = 4;

    x_recede_index = x_init_recede_index;
    x_advance_index = x_init_advance_index;

    for_less( x, 0, sizes[VIO_X] )
    {
        for_less( z, 0, sizes[VIO_Z] )
        {
            GET_FIRST_SAMPLE( y_init_advance_index, y_left_weight, sizes[VIO_Y],
                              (VIO_Real) slice[_I][z], sample )
            row[z] = (float) sample;
        }

        y_recede_index = y_init_recede_index;
        y_advance_index = y_init_advance_index;

        for_less( y, 0, sizes[VIO_Y] )
        {
            GET_FIRST_SAMPLE( z_init_advance_index, z_left_weight,
                              sizes[VIO_Z], (VIO_Real) row[_I], sum )
            z_recede_index = z_init_recede_index;
            z_advance_index = z_init_advance_index;

            for_less( z, 0, sizes[VIO_Z] )
            {
#ifdef DEBUG
                correct_voxel = get_correct_amount( volume, x, y, z,
                                                    x_width, y_width, z_width );

                if( !numerically_close( sum, correct_voxel, 0.001 ) )
                    handle_internal_error( "Dang" );
#endif

                value = sum / total_volume;
                set_volume_real_value( resampled_volume, x, y, z, 0, 0, value );

                if( z == sizes[VIO_Z]-1 )
                    continue;

                GET_NEXT_SAMPLE( z_recede_index, z_advance_index,
                                 z_left_weight, z_right_weight,
                                 sizes[VIO_Z], (VIO_Real) row[_I], sum )

                ++z_recede_index;
                ++z_advance_index;
            }

            if( y == sizes[VIO_Y]-1 )
                continue;

            for_less( z, 0, sizes[VIO_Z] )
            {
                sample = (VIO_Real) row[z];
                GET_NEXT_SAMPLE( y_recede_index, y_advance_index,
                                 y_left_weight, y_right_weight,
                                 sizes[VIO_Y], (VIO_Real) slice[_I][z], sample )
                row[z] = (float) sample;
            }

            ++y_recede_index;
            ++y_advance_index;

            update_progress_report( &progress, x * sizes[VIO_Y] + y + 1 );
        }

        if( x == sizes[VIO_X] - 1 )
            continue;

        for_less( y, 0, sizes[VIO_Y] )
        {
            for_less( z, 0, sizes[VIO_Z] )
            {
                sample = (VIO_Real) slice[y][z];
                GET_NEXT_SAMPLE( x_recede_index, x_advance_index,
                                 x_left_weight, x_right_weight,
                                 sizes[VIO_X], volume_cache[_I][y][z],
                                 sample )
                slice[y][z] = (float) sample;
            }
        }

        ++x_recede_index;
        ++x_advance_index;

        if( x_recede_index+1 >= 0 && x_recede_index+1 < sizes[VIO_X] &&
            volume_cache[x_recede_index+1] == NULL )
        {
            if( n_alloced < max_alloced )
            {
                VIO_ALLOC2D( volume_cache[x_recede_index+1], sizes[VIO_Y], sizes[VIO_Z] );
                ++n_alloced;
            }
            else
            {
                volume_cache[x_recede_index+1] = volume_cache[x_recede_index-1];
                volume_cache[x_recede_index-1] = NULL;
            }

            get_volume_value_hyperslab_3d( volume, x_recede_index+1, 0, 0,
                                        1, sizes[VIO_Y], sizes[VIO_Z],
                                        &volume_cache[x_recede_index+1][0][0] );
        }

        if( x_advance_index+1 >= 0 && x_advance_index+1 < sizes[VIO_X] &&
            volume_cache[x_advance_index+1] == NULL )
        {
            if( x_advance_index - 1 == x_recede_index ||
                x_advance_index - 1 == x_recede_index+1 )
            {
                VIO_ALLOC2D( volume_cache[x_advance_index+1], sizes[VIO_Y], sizes[VIO_Z] );
                ++n_alloced;
            }
            else
            {
                volume_cache[x_advance_index+1] =
                                            volume_cache[x_advance_index-1];
                volume_cache[x_advance_index-1] = NULL;
            }

            get_volume_value_hyperslab_3d( volume, x_advance_index+1, 0, 0,
                                        1, sizes[VIO_Y], sizes[VIO_Z],
                                      &volume_cache[x_advance_index+1][0][0] );
        }
    }

    terminate_progress_report( &progress );

    for_less( x, 0, sizes[VIO_X] )
    {
        if( volume_cache[x] != NULL )
        {
            VIO_FREE2D( volume_cache[x] );
            --n_alloced;
        }
    }

    FREE( volume_cache );

    VIO_FREE2D( slice );
    FREE( row );

    return( resampled_volume );
}

#ifdef DEBUG
static void  get_voxel_range(
    int      size,
    VIO_Real     min_pos,
    VIO_Real     max_pos,
    int      *min_voxel,
    int      *max_voxel,
    VIO_Real     *start_weight,
    VIO_Real     *end_weight )
{
    if( min_pos < -0.5 )
        min_pos = -0.5;
    if( max_pos > (VIO_Real) size - 0.5 )
        max_pos = (VIO_Real) size - 0.5;

    *min_voxel = VIO_ROUND( min_pos );
    *max_voxel = VIO_ROUND( max_pos );

    if( (VIO_Real) (*max_voxel) - 0.5 == max_pos )
        --(*max_voxel);

    if( *min_voxel == *max_voxel )
        *start_weight = max_pos - min_pos;
    else
    {
        *start_weight = ((VIO_Real) (*min_voxel) + 0.5) - min_pos;
        *end_weight = max_pos - ((VIO_Real) (*max_voxel) - 0.5);
    }
}
    
static VIO_Real  get_amount_in_box(
    VIO_Volume    volume,
    int       x_min_voxel,
    int       x_max_voxel,
    int       y_min_voxel,
    int       y_max_voxel,
    int       z_min_voxel,
    int       z_max_voxel,
    VIO_Real      x_weight_start,
    VIO_Real      x_weight_end,
    VIO_Real      y_weight_start,
    VIO_Real      y_weight_end,
    VIO_Real      z_weight_start,
    VIO_Real      z_weight_end )
{
    int     x, y, z;
    VIO_Real    sum, z_sum, x_sum, value;

    sum = 0.0;

    for_inclusive( z, z_min_voxel, z_max_voxel )
    {
        z_sum = 0.0;

        for_inclusive( x, x_min_voxel, x_max_voxel )
        {
            x_sum = 0.0;

            for_inclusive( y, y_min_voxel, y_max_voxel )
            {
                value = get_volume_real_value( volume, x, y, z, 0, 0 );

                if( y == y_min_voxel )
                    value *= y_weight_start;
                else if( y == y_max_voxel )
                    value *= y_weight_end;

                x_sum += value;
            }

            if( x == x_min_voxel )
                x_sum *= x_weight_start;
            else if( x == x_max_voxel )
                x_sum *= x_weight_end;

            z_sum += x_sum;
        }

        if( z == z_min_voxel )
            z_sum *= z_weight_start;
        else if( z == z_max_voxel )
            z_sum *= z_weight_end;

        sum += z_sum;
    }

    return( sum );
}

static VIO_Real  get_correct_amount(
    VIO_Volume   volume,
    int      x,
    int      y,
    int      z,
    VIO_Real     x_width,
    VIO_Real     y_width,
    VIO_Real     z_width )
{
    int     sizes[VIO_MAX_DIMENSIONS];
    VIO_Real    x_start_weight, x_end_weight;
    VIO_Real    y_start_weight, y_end_weight;
    VIO_Real    z_start_weight, z_end_weight;
    VIO_Real    sum;
    int     x_min_voxel, x_max_voxel;
    int     y_min_voxel, y_max_voxel;
    int     z_min_voxel, z_max_voxel;

    get_volume_sizes( volume, sizes );

    get_voxel_range( sizes[VIO_X], (VIO_Real) x - x_width, (VIO_Real) x + x_width,
                    &x_min_voxel, &x_max_voxel, &x_start_weight, &x_end_weight);
    get_voxel_range( sizes[VIO_Y], (VIO_Real) y - y_width, (VIO_Real) y + y_width,
                    &y_min_voxel, &y_max_voxel, &y_start_weight, &y_end_weight);
    get_voxel_range( sizes[VIO_Z], (VIO_Real) z - z_width, (VIO_Real) z + z_width,
                    &z_min_voxel, &z_max_voxel, &z_start_weight, &z_end_weight);

    sum = get_amount_in_box( volume,
                             x_min_voxel, x_max_voxel,
                             y_min_voxel, y_max_voxel,
                             z_min_voxel, z_max_voxel,
                             x_start_weight, x_end_weight,
                             y_start_weight, y_end_weight,
                             z_start_weight, z_end_weight );

    return( sum );
}
#endif
