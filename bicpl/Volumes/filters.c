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
#include  <bicpl/vols.h>

#define  N_STD_DEVIATIONS  3.0
#define  N_SAMPLES         9

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_slice_weights_for_filter
@INPUT      : volume
              voxel_position
              voxel_direction
              filter_type
              full_width_half_max
@OUTPUT     : positions
              weights
@RETURNS    : 
@DESCRIPTION: Computes the slice weights for the given filter.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI int  get_slice_weights_for_filter(
    VIO_Volume         volume,
    VIO_Real           voxel_position[],
    VIO_Real           voxel_direction[],   /* if filter_type != NEAREST */
    VIO_Filter_types   filter_type,
    VIO_Real           full_width_half_max,
    VIO_Real           ***positions,
    VIO_Real           *weights[] )
{
    int   i, c, n_slices, n_dims, axis;
    int   slice, first_slice, last_slice, s;
    VIO_Real  sum_weights, frac, start, end, weight;
    VIO_Real  half_width, sigma, start_interval, end_interval, x;
    VIO_Real  *origins;

    n_dims = get_volume_n_dimensions( volume );

    n_slices = 0;

    switch( filter_type )
    {
    case NEAREST_NEIGHBOUR:
        ALLOC( origins, 1 );
        ALLOC( *weights, 1 );
        origins[0] = 0.0;
        (*weights)[0] = 1.0;
        n_slices = 1;
        break;

    case LINEAR_INTERPOLATION:
        ALLOC( origins, 2 );
        ALLOC( *weights, 2 );

        axis = -1;
        for_less( c, 0, n_dims )
        {
            if( voxel_direction[c] != 0.0 )
            {
                if( axis == -1 )
                    axis = c;
                else
                {
                    /* Cannot do linear interpolation on non-ortho axis */
                    axis = -1;
                    break;
                }
            }
        }
        frac = (axis >= 0) ? VIO_FRACTION( voxel_position[axis] ) : 0.0;
        origins[0] = 0.0;
        (*weights)[0] = frac;
        n_slices = 1;
        if( frac > 0.0 )
        {
            origins[1] = 1.0;
            (*weights)[1] = 1.0 - frac;
            n_slices = 2;
        }
        break;

    case BOX_FILTER:
    case TRIANGLE_FILTER:
    case GAUSSIAN_FILTER:
        switch( filter_type )
        {
        case BOX_FILTER:       half_width = full_width_half_max / 2.0;  break;
        case TRIANGLE_FILTER:  half_width = full_width_half_max;  break;
        case GAUSSIAN_FILTER:
            sigma = full_width_half_max / 2.0 / sqrt(log(2.0));
            half_width = N_STD_DEVIATIONS * sigma;
            break;
        default:
            break;                /* avoid warning */
        }
        first_slice = (int) (-half_width - 0.5);
        last_slice = (int) (half_width + 0.5);

        n_slices = last_slice - first_slice + 1;
        ALLOC( origins, n_slices );
        ALLOC( *weights, n_slices );

        for_inclusive( slice, first_slice, last_slice )
        {
            origins[slice-first_slice] = (VIO_Real) slice;

            if( slice == first_slice )
                start_interval = -half_width;
            else
                start_interval = (VIO_Real) slice - 0.5;

            if( slice == last_slice )
                end_interval = half_width;
            else
                end_interval = (VIO_Real) slice + 0.5;

            switch( filter_type )
            {
            case BOX_FILTER:
                weight = end_interval - start_interval;
                break;
            case TRIANGLE_FILTER:
                weight = 0.0;
                if( start_interval < 0.0 )
                {
                    end = MIN( end_interval, 0.0);
                    weight = (end - start_interval) *
                             ((start_interval + end) / 2.0 + half_width) /
                             half_width;
                }
                if( end_interval > 0.0 )
                {
                    start = MAX( start_interval, 0.0);
                    weight += (end_interval - start) *
                             (half_width - (start + end_interval) / 2.0) /
                             half_width;
                }
                break;
            case GAUSSIAN_FILTER:
                weight = 0.0;
                if( end_interval > start_interval )
                {
                    for_less( s, 0, N_SAMPLES )
                    {
                        x = start_interval + (end_interval - start_interval) *
                            ((VIO_Real) s + 0.5) / (VIO_Real) N_SAMPLES;
                        weight += (end_interval - start_interval) *
                                  exp ( - x * x / sigma / sigma );
                    }
                }
                break;
            default:
                break;          /* avoid warning */
            }

            (*weights)[slice-first_slice] = weight;
        }

        break;
    }

    if( n_slices < 1 )
    {
        handle_internal_error( "get_slice_weights_for_filter" );
        return( 0 );
    }

    VIO_ALLOC2D( *positions, n_slices, n_dims );
    sum_weights = 0.0;
    for_less( i, 0, n_slices )
        sum_weights += (*weights)[i];

    for_less( i, 0, n_slices )
    {
        if( sum_weights == 0.0 )
            (*weights)[i] = 1.0 / (VIO_Real) n_slices;
        else
            (*weights)[i] /= sum_weights;

        for_less( c, 0, n_dims )
        {
            if( filter_type == NEAREST_NEIGHBOUR )
            {
                 (*positions)[i][c] = voxel_position[c];
            }
            else
            {
                 (*positions)[i][c] = voxel_position[c] +
                                      origins[i] * voxel_direction[c];
            }
        }
    }

    FREE( origins );

    return( n_slices );
}
