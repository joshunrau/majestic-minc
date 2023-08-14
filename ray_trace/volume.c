#if HAVE_CONFIG_H
#include "config.h"
#endif

#include  <ray_trace.h>
#include  <bicpl/deform.h>

#define  BOUNDARY_TOLERANCE   1.0e-4

/* ARGSUSED */

static  VIO_BOOL  search_volume(
    VIO_Volume             volume,
    int                continuity,
    bitlist_3d_struct  *done_bits,
    bitlist_3d_struct  *surface_bits,
    VIO_Real               threshold,
    VIO_Point              *origin,
    VIO_Vector             *direction,
    VIO_Real               t_min,
    VIO_Real               t_max,
    VIO_Real               *t )
{
    boundary_definition_struct   boundary_def;

    boundary_def.min_isovalue = threshold;
    boundary_def.max_isovalue = threshold;
    boundary_def.gradient_threshold = -1.0;
    boundary_def.min_dot_product = -2.0;
    boundary_def.max_dot_product = 2.0;
    boundary_def.normal_direction = ANY_DIRECTION;
    boundary_def.tolerance = BOUNDARY_TOLERANCE;

    return( find_boundary_in_direction( volume, NULL, NULL,
                                        done_bits, surface_bits, 0.0,
                                        origin, direction, direction,
                                        1.0e30, 0.0, continuity,
                                        &boundary_def, t ) );
}

  VIO_BOOL  ray_intersects_a_volume(
    VIO_Point              *origin,
    VIO_Vector             *direction,
    VIO_Volume             volume,
    int                continuity,
    bitlist_3d_struct  *done_bits,
    bitlist_3d_struct  *surface_bits,
    VIO_Real               threshold,
    VIO_Point              *point,
    VIO_Vector             *normal,
    VIO_Real               *dist )
{
    int      sizes[VIO_MAX_DIMENSIONS];
    VIO_Real     t, t_min, t_max, xw, yw, zw;
    VIO_Real     voxel[VIO_MAX_DIMENSIONS];
    VIO_Real     voxel_origin[VIO_MAX_DIMENSIONS];
    VIO_Real     voxel_direction[VIO_MAX_DIMENSIONS];
    VIO_Real     value, dx, dy, dz;
    VIO_Point    p_origin, int_point;
    VIO_Vector   p_direction;
    VIO_BOOL  found;

    convert_world_to_voxel( volume,
                            (VIO_Real) Point_x(*origin),
                            (VIO_Real) Point_y(*origin),
                            (VIO_Real) Point_z(*origin), voxel_origin );

    convert_world_vector_to_voxel( volume,
                                   (VIO_Real) Vector_x(*direction),
                                   (VIO_Real) Vector_y(*direction),
                                   (VIO_Real) Vector_z(*direction),
                                   voxel_direction );

    fill_Point( p_origin, voxel_origin[VIO_X], voxel_origin[VIO_Y], voxel_origin[VIO_Z] );
    fill_Vector( p_direction, voxel_direction[VIO_X], voxel_direction[VIO_Y],
                 voxel_direction[VIO_Z] );

    get_volume_sizes( volume, sizes );

    if( !clip_line_to_box( &p_origin, &p_direction,
                           0.0, (VIO_Real) sizes[VIO_X]-1.0,
                           0.0, (VIO_Real) sizes[VIO_Y]-1.0,
                           0.0, (VIO_Real) sizes[VIO_Z]-1.0,
                           &t_min, &t_max ) )
    {
        t_min = 0.0;
        t_max = 0.0;
    }

    if( t_min < 0.0 )
        t_min = 0.0;
    if( t_max < 0.0 )
        t_max = 0.0;

    if( t_min >= t_max )
        return( FALSE );

    found = search_volume( volume, continuity, done_bits, surface_bits,
                           threshold, origin, direction, t_min, t_max, &t );

    if( found )
    {
        GET_POINT_ON_RAY( int_point, p_origin, p_direction, t );

        voxel[VIO_X] = (VIO_Real) Point_x( int_point );
        voxel[VIO_Y] = (VIO_Real) Point_y( int_point );
        voxel[VIO_Z] = (VIO_Real) Point_z( int_point );

        convert_voxel_to_world( volume, voxel, &xw, &yw, &zw );

        if( point != (VIO_Point *) NULL )
        {
            fill_Point( *point, xw, yw, zw );

            evaluate_volume_in_world( volume, xw, yw, zw, 2, FALSE,
                                      get_volume_real_min(volume),
                                      &value, &dx, &dy, &dz,
                                      NULL, NULL, NULL,
                                      NULL, NULL, NULL );

            fill_Vector( *normal, dx, dy, dz );

            if( DOT_VECTORS( *normal, *direction ) > 0.0 )
                SCALE_VECTOR( *normal, *normal, -1.0 );

            NORMALIZE_VECTOR( *normal, *normal );
        }

        *dist = t;
    }

    return( found );
}

  void  get_range_of_volume(
    VIO_Volume    volume,
    VIO_Point     *min_point,
    VIO_Point     *max_point )
{
    int      i, j, k, sizes[VIO_MAX_DIMENSIONS];
    VIO_Real     voxel[VIO_MAX_DIMENSIONS], xw, yw, zw;
    VIO_Point    pos;

    get_volume_sizes( volume, sizes );

    for_less( i, 0, 2 )
    {
        if( i == 0 )
            voxel[VIO_X] = 0.0;
        else
            voxel[VIO_X] = (VIO_Real) sizes[VIO_X] - 1.0;
        for_less( j, 0, 2 )
        {
            if( j == 0 )
                voxel[VIO_Y] = 0.0;
            else
                voxel[VIO_Y] = (VIO_Real) sizes[VIO_Y] - 1.0;
            for_less( k, 0, 2 )
            {
                if( k == 0 )
                    voxel[VIO_Z] = 0.0;
                else
                    voxel[VIO_Z] = (VIO_Real) sizes[VIO_Z] - 1.0;

                convert_voxel_to_world( volume, voxel, &xw, &yw, &zw );
                fill_Point( pos, xw, yw, zw );

                if( i == 0 && j == 0 && k == 0 )
                {
                    *min_point = pos;
                    *max_point = pos;
                }
                else
                {
                    apply_point_to_min_and_max( &pos, min_point, max_point );
                    apply_point_to_min_and_max( &pos, min_point, max_point );
                }
            }
        }
    }
}
