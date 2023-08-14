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
@NAME       : interpolate_volume_to_slice
@INPUT      : volume1
              n_dims1
              origin1
              x_axis1
              y_axis1
              volume2
              n_dims2
              origin2
              x_axis2
              y_axis2
              degrees_continuity
              cmode_colour_map
              rgb_colour_map
              empty_colour
              colour_coding,
@OUTPUT     : pixels
@RETURNS    : 
@DESCRIPTION: Interpolates the volume to the given slice, using either
              nearest neighbour, linear or cubic interpolation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

BICAPI void  interpolate_volume_to_slice(
    VIO_Volume      volume1,
    int             n_dims1,
    VIO_Real        origin1[],
    VIO_Real        x_axis1[],
    VIO_Real        y_axis1[],
    VIO_Volume      volume2,
    int             n_dims2,
    VIO_Real        origin2[],
    VIO_Real        x_axis2[],
    VIO_Real        y_axis2[],
    int             x_pixel_start,
    int             x_pixel_end,
    int             y_pixel_start,
    int             y_pixel_end,
    int             degrees_continuity,
    unsigned short  **cmode_colour_map,
    VIO_Colour      **rgb_colour_map,
    VIO_Colour      empty_colour,
    colour_coding_struct *colour_coding,
    pixels_struct   *pixels )
{
    int              dim, x, y;
    int              int_voxel_value1, int_voxel_value2;
    VIO_Real         outside_value1, outside_value2;
    VIO_Real         start_voxel1[VIO_MAX_DIMENSIONS], voxel1[VIO_MAX_DIMENSIONS];
    VIO_Real         start_voxel2[VIO_MAX_DIMENSIONS], voxel2[VIO_MAX_DIMENSIONS];
    VIO_Real         value1, voxel_value1, value2, voxel_value2;
    VIO_Real         min_voxel1, max_voxel1, min_voxel2, max_voxel2;
    unsigned short   *cmode_ptr;
    VIO_Colour       *rgb_ptr;
    VIO_BOOL         inside1, inside2;
    Pixel_types      pixel_type;
#ifdef  REPORT_PROGRESS
    VIO_progress_struct  progress;
#endif

    if( is_an_rgb_volume( volume1 ) ||
        (volume2 != NULL && is_an_rgb_volume( volume2 )) ) {
        colour_coding = NULL;   /* Don't use the colour coding. */
        if (degrees_continuity == 2)
          degrees_continuity = 0;
    }

    get_volume_voxel_range( volume1, &min_voxel1, &max_voxel1 );

    pixel_type = pixels->pixel_type;

    for_less( dim, 0, n_dims1 )
    {
        start_voxel1[dim] = origin1[dim] + (VIO_Real) x_pixel_start * x_axis1[dim] +
                                           (VIO_Real) y_pixel_start * y_axis1[dim];
    }

    outside_value1 = 0.0;

    if( volume2 != NULL )
    {
        get_volume_voxel_range( volume2, &min_voxel2, &max_voxel2 );
        for_less( dim, 0, n_dims2 )
        {
            start_voxel2[dim] = origin2[dim] +
                                (VIO_Real) x_pixel_start * x_axis2[dim] +
                                (VIO_Real) y_pixel_start * y_axis2[dim];
        }

        outside_value2 = 0.0;
    }

#ifdef  REPORT_PROGRESS
    initialize_progress_report( &progress, FALSE,
                                y_pixel_end - y_pixel_start + 1,
                                "Creating Slice" );
#endif

    for_inclusive( y, y_pixel_start, y_pixel_end )
    {
        for_less( dim, 0, n_dims1 )
            voxel1[dim] = start_voxel1[dim];

        if( volume2 != NULL )
        {
            for_less( dim, 0, n_dims2 )
                voxel2[dim] = start_voxel2[dim];
        }

        if( pixel_type == RGB_PIXEL )
            rgb_ptr = &PIXEL_RGB_COLOUR(*pixels,x_pixel_start,y);
        else
            cmode_ptr = &PIXEL_COLOUR_INDEX_16(*pixels,x_pixel_start,y);

        for_inclusive( x, x_pixel_start, x_pixel_end )
        {
            inside1 = voxel_is_within_volume( volume1, voxel1 );

            if( inside1 )
            {
                (void) evaluate_volume( volume1, voxel1, NULL,
                                        degrees_continuity,
                                        FALSE, outside_value1,
                                        &value1, NULL, NULL );

                voxel_value1 = convert_value_to_voxel( volume1, value1 );

                if( voxel_value1 < min_voxel1 )
                    voxel_value1 = min_voxel1;
                else if( voxel_value1 > max_voxel1 )
                    voxel_value1 = max_voxel1;
            }

            for_less( dim, 0, n_dims1 )
                voxel1[dim] += x_axis1[dim];

            if( volume2 != NULL )
            {
                inside2 = voxel_is_within_volume( volume2, voxel2 );

                if( inside2 )
                {
                    (void) evaluate_volume( volume2, voxel2, NULL,
                                            degrees_continuity,
                                            FALSE, outside_value2,
                                            &value2, NULL, NULL );

                    voxel_value2 = convert_value_to_voxel( volume2, value2 );

                    if( voxel_value2 < min_voxel2 )
                        voxel_value2 = min_voxel2;
                    else if( voxel_value2 > max_voxel2 )
                        voxel_value2 = max_voxel2;
                }

                for_less( dim, 0, n_dims2 )
                    voxel2[dim] += x_axis2[dim];

                if( pixel_type == RGB_PIXEL )
                {
                    if( inside1 && inside2 )
                    {
                        int_voxel_value1 = VIO_ROUND( voxel_value1 );
                        int_voxel_value2 = VIO_ROUND( voxel_value2 );
                        *rgb_ptr = rgb_colour_map[int_voxel_value1]
                                                 [int_voxel_value2];
                    }
                    else
                        *rgb_ptr = empty_colour;

                    ++rgb_ptr;
                }
                else
                {
                    if( inside1 && inside2 )
                    {
                        int_voxel_value1 = VIO_ROUND( voxel_value1 );
                        int_voxel_value2 = VIO_ROUND( voxel_value2 );
                        *cmode_ptr = cmode_colour_map[int_voxel_value1]
                                                     [int_voxel_value2];
                    }
                    else
                        *cmode_ptr = (unsigned short) empty_colour;

                    ++cmode_ptr;
                }
            }
            else
            {
                if( pixel_type == RGB_PIXEL )
                {
                    if( inside1 )
                    {
                        if( rgb_colour_map == NULL )
                        {
                            if (colour_coding != NULL)
                            {
                                *rgb_ptr = get_colour_code(colour_coding,
                                                           value1);
                            }
                            else
                                *rgb_ptr = (VIO_Colour) voxel_value1;
                        }
                        else
                        {
                            int_voxel_value1 = VIO_ROUND( voxel_value1 );
                            *rgb_ptr = rgb_colour_map[0][int_voxel_value1];
                        }
                    }
                    else
                        *rgb_ptr = empty_colour;

                    ++rgb_ptr;
                }
                else
                {
                    if( inside1 )
                    {
                        int_voxel_value1 = VIO_ROUND( voxel_value1 );
                        if( cmode_colour_map == NULL )
                            *cmode_ptr = (unsigned short) int_voxel_value1;
                        else
                        {
                            *cmode_ptr = cmode_colour_map[0][int_voxel_value1];
                        }
                    }
                    else
                        *cmode_ptr = (unsigned short) empty_colour;

                    ++cmode_ptr;
                }
            }
        }

        for_less( dim, 0, n_dims1 )
            start_voxel1[dim] += y_axis1[dim];

        if( volume2 != NULL )
        {
            for_less( dim, 0, n_dims2 )
                start_voxel2[dim] += y_axis2[dim];
        }

#ifdef  REPORT_PROGRESS
        update_progress_report( &progress, y-y_pixel_start+1 );
#endif
    }

#ifdef  REPORT_PROGRESS
    terminate_progress_report( &progress );
#endif
}
