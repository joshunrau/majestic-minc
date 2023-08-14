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

static  void  recursive_intersect_ray(
    VIO_Point                 *origin,
    VIO_Vector                *direction,
    VIO_Real                  t_min,
    VIO_Real                  t_max,
    bintree_node_struct   *node,
    object_struct         *object,
    int                   *obj_index,
    VIO_Real                  *closest_dist,
    int                   *n_intersections,
    VIO_Real                  *distances[] );

static  int  n_nodes_searched = 0;
static  int  n_objects_searched = 0;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : print_bintree_stats
@INPUT      : n_objects
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Prints information on the bintree search structure.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jun 21, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  print_bintree_stats(
    int   n_objects )
{
    print( "Nodes %g  ", (VIO_Real) n_nodes_searched / (VIO_Real) n_objects );
    print( "Objects %g\n", (VIO_Real) n_objects_searched / (VIO_Real) n_objects );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : intersect_ray_with_bintree
@INPUT      : origin
              direction
              bintree
@OUTPUT     : object
              obj_index
              dist
              distances
@RETURNS    : number of intersections
@DESCRIPTION: Tests if the ray intersects the objects in the bintree.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  intersect_ray_with_bintree(
    VIO_Point               *origin,
    VIO_Vector              *direction,
    bintree_struct_ptr  bintree,
    object_struct       *object,
    int                 *obj_index,
    VIO_Real                *dist,
    VIO_Real                *distances[] )
{
    int       n_intersections;
    VIO_Real      t_min, t_max;

    n_intersections = 0;
    if( obj_index != (int *) NULL )
        *obj_index = -1;

    if( ray_intersects_range( &bintree->range, origin, direction,
                              &t_min, &t_max ) )
    {
        recursive_intersect_ray( origin, direction, t_min, t_max,
                                 bintree->root, object, obj_index, dist,
                                 &n_intersections, distances );
    }

    return( n_intersections );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : recursive_intersect_ray
@INPUT      : origin
              direction
              t_min
              t_max
              node
@OUTPUT     : object
              obj_index
              closest_dist
              n_intersections
              distances
@RETURNS    : 
@DESCRIPTION: Traverses the bintree testing for ray intersection.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jun 21, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  recursive_intersect_ray(
    VIO_Point                 *origin,
    VIO_Vector                *direction,
    VIO_Real                  t_min,
    VIO_Real                  t_max,
    bintree_node_struct   *node,
    object_struct         *object,
    int                   *obj_index,
    VIO_Real                  *closest_dist,
    int                   *n_intersections,
    VIO_Real                  *distances[] )
{
    VIO_BOOL               test_child, searching_left;
    int                   i, n_objects, *object_list, axis_index;
    bintree_node_struct   *left_child, *right_child;
    VIO_Real                  delta, left_limit, right_limit;
    VIO_Real                  t, t_min_child, t_max_child;

    if( distances == NULL && obj_index != NULL &&
        *obj_index >= 0 && *closest_dist < t_min )
        return;

    ++n_nodes_searched;

    if( bintree_node_is_leaf( node ) )
    {
        n_objects = get_bintree_leaf_objects( node, &object_list );

        for_less( i, 0, n_objects )
        {
            ++n_objects_searched;

             intersect_ray_object( origin, direction,
                                   object, object_list[i], obj_index,
                                   closest_dist, n_intersections,
                                   distances );
        }
    }
    else
    {
        axis_index = get_node_split_axis( node );
        delta = (VIO_Real) Vector_coord( *direction, axis_index );

        if( delta > 0.0 )
            searching_left = TRUE;
        else
            searching_left = FALSE;

        for_less( i, 0, 2 )
        {
            t_min_child = t_min;
            t_max_child = t_max;

            if( searching_left && get_bintree_left_child( node, &left_child ) )
            {
                left_limit = get_node_split_position( left_child );

                if( delta == 0.0 )
                {
                    test_child = ((VIO_Real) Point_coord(*origin,axis_index) <=
                                  left_limit);
                }
                else
                {
                    test_child = FALSE;

                    t = (left_limit - (VIO_Real) Point_coord(*origin,axis_index)) /
                        delta;

                    if( delta < 0.0 && t <= t_max_child )
                    {
                        test_child = TRUE;

                        if( t > t_min_child )
                            t_min_child = t;
                    }
                    else if( delta > 0.0 && t >= t_min_child )
                    {
                        test_child = TRUE;

                        if( t < t_max_child )
                            t_max_child = t;
                    }
                }

                if( test_child )
                {
                    recursive_intersect_ray( origin, direction,
                                             t_min_child, t_max_child,
                                             left_child, object,
                                             obj_index, closest_dist,
                                             n_intersections, distances );

                    if( distances == NULL && obj_index != NULL &&
                        *obj_index >= 0 && *closest_dist < t_min )
                        return;
                }
            }
            else if( !searching_left &&
                     get_bintree_right_child( node, &right_child ) )
            {
                right_limit = get_node_split_position( right_child );

                if( delta == 0.0 )
                {
                    test_child = ((VIO_Real) Point_coord(*origin,axis_index) >=
                                  right_limit);
                }
                else
                {
                    test_child = FALSE;

                    t = (right_limit - (VIO_Real) Point_coord(*origin,axis_index)) /
                        delta;

                    if( delta < 0.0 && t >= t_min_child )
                    {
                        test_child = TRUE;

                        if( t < t_max_child )
                            t_max_child = t;
                    }
                    else if( delta > 0.0 && t <= t_max_child )
                    {
                        test_child = TRUE;

                        if( t > t_min_child )
                            t_min_child = t;
                    }
                }

                if( test_child )
                {
                    recursive_intersect_ray( origin, direction,
                                             t_min_child, t_max_child,
                                             right_child, object,
                                             obj_index, closest_dist,
                                             n_intersections, distances );

                    if( distances == NULL && obj_index != NULL &&
                        *obj_index >= 0 && *closest_dist < t_min )
                        return;
                }
            }

            searching_left = !searching_left;
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : ray_intersects_range
@INPUT      : range   - a box
              origin
              direction
@OUTPUT     : t_min
              t_max
@RETURNS    : TRUE if ray intersects box
@DESCRIPTION: Tests if the ray intersects the box and passes back the
              two intersection distances.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jun 21, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  VIO_BOOL  ray_intersects_range(
    range_struct  *range,
    VIO_Point         *origin,
    VIO_Vector        *direction,
    VIO_Real          *t_min,
    VIO_Real          *t_max )
{
    VIO_BOOL  intersects;

    intersects = clip_line_to_box( origin, direction,
                                   (VIO_Real) range->limits[VIO_X][0],
                                   (VIO_Real) range->limits[VIO_X][1],
                                   (VIO_Real) range->limits[VIO_Y][0],
                                   (VIO_Real) range->limits[VIO_Y][1],
                                   (VIO_Real) range->limits[VIO_Z][0],
                                   (VIO_Real) range->limits[VIO_Z][1], t_min, t_max );


    if( intersects )
    {
        if( *t_min < 0.0 )
            *t_min = 0.0;

        intersects = (*t_min <= *t_max);
    }

    return( intersects );
}
