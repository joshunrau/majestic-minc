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

/* ----------------------------- MNI Header -----------------------------------
@NAME       : null_Point
@INPUT      : p
@OUTPUT     : 
@RETURNS    : TRUE if point is ( 0, 0, 0 ).
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  VIO_BOOL  null_Point(
    VIO_Point   *p )
{
    return( Point_x(*p) == 0.0f &&
            Point_y(*p) == 0.0f &&
            Point_z(*p) == 0.0f );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : null_Vector
@INPUT      : 
@OUTPUT     : 
@RETURNS    : TRUE if vector is ( 0, 0, 0 ).
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  VIO_BOOL  null_Vector(
    VIO_Vector   *v )
{
    return( Vector_x(*v) == 0.0f &&
            Vector_y(*v) == 0.0f &&
            Vector_z(*v) == 0.0f );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : distance_between_points
@INPUT      : p1
              p2
@OUTPUT     : 
@RETURNS    : Euclidean distance between points
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  VIO_Real  distance_between_points(
    VIO_Point  *p1,
    VIO_Point  *p2 )
{
    VIO_Real   dx, dy, dz;

    dx = (VIO_Real) Point_x(*p2) - (VIO_Real) Point_x(*p1);
    dy = (VIO_Real) Point_y(*p2) - (VIO_Real) Point_y(*p1);
    dz = (VIO_Real) Point_z(*p2) - (VIO_Real) Point_z(*p1);

    return( sqrt( dx * dx + dy * dy + dz * dz ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : points_within_distance
@INPUT      : p1
              p2
              distance
@OUTPUT     : 
@RETURNS    : TRUE if points within distance
@DESCRIPTION: Determines if two points are within the specified distance.  Tries
              to be efficient by not performing a sqrt.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  VIO_BOOL  points_within_distance(
    VIO_Point  *p1,
    VIO_Point  *p2,
    VIO_Real   distance )
{
    VIO_Real  dx, dy, dz;

    dx = (VIO_Real) Point_x(*p1) - (VIO_Real) Point_x(*p2);
    dy = (VIO_Real) Point_y(*p1) - (VIO_Real) Point_y(*p2);
    dz = (VIO_Real) Point_z(*p1) - (VIO_Real) Point_z(*p2);

    return( dx * dx + dy * dy + dz * dz <= distance * distance );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : apply_point_to_min_and_max
@INPUT      : point
              min_point   - current min position
              max_point   - current max position
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the range of space specified by min_point and max_point
              to the smallest rectangular range containing both the current
              range and the specified point.  This is
              performed by expanding the rectangular range if needed, to
              contain the specified point.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  apply_point_to_min_and_max(
    VIO_Point   *point,
    VIO_Point   *min_point,
    VIO_Point   *max_point )
{
    int    c;

    for_less( c, 0, VIO_N_DIMENSIONS )
    {
        if( Point_coord(*point,c) < Point_coord(*min_point,c) )
            Point_coord(*min_point,c) = Point_coord(*point,c);

        if( Point_coord(*point,c) > Point_coord(*max_point,c) )
            Point_coord(*max_point,c) = Point_coord(*point,c);
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : expand_min_and_max_points
@INPUT      : min_point
              max_point
@OUTPUT     : min_to_check
              max_to_check
@RETURNS    : 
@DESCRIPTION: Sets the rectangular range of space specified by min_point and
              max_point to the rectangular range containing the four argument
              points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  expand_min_and_max_points(
    VIO_Point   *min_point,
    VIO_Point   *max_point,
    VIO_Point   *min_to_check,
    VIO_Point   *max_to_check )
{
    int    c;

    for_less( c, 0, VIO_N_DIMENSIONS )
    {
        if( Point_coord(*min_to_check,c) < Point_coord(*min_point,c) )
            Point_coord(*min_point,c) = Point_coord(*min_to_check,c);
        if( Point_coord(*max_to_check,c) > Point_coord(*max_point,c) )
            Point_coord(*max_point,c) = Point_coord(*max_to_check,c);
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_range_points
@INPUT      : n_points
              points
@OUTPUT     : min_corner
              max_corner
@RETURNS    : 
@DESCRIPTION: Finds the rectangular range containing all the points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_range_points(
    int                n_points,
    VIO_Point              points[],
    VIO_Point              *min_corner,
    VIO_Point              *max_corner )
{
    int     i;

    if( n_points > 0 )
    {
        *min_corner = points[0];
        *max_corner = points[0];

        for_less( i, 1, n_points )
        {
            apply_point_to_min_and_max( &points[i], min_corner, max_corner );
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_points_centroid
@INPUT      : n_points
              points
@OUTPUT     : centroid
@RETURNS    : 
@DESCRIPTION: Finds the centroid of the set of points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_points_centroid(
    int     n_points,
    VIO_Point   points[],
    VIO_Point   *centroid )
{
    int   i;
    VIO_Real  x, y, z;

    x = 0.0;
    y = 0.0;
    z = 0.0;

    for_less( i, 0, n_points )
    {
        x += (VIO_Real) Point_x(points[i]);
        y += (VIO_Real) Point_y(points[i]);
        z += (VIO_Real) Point_z(points[i]);
    }

    fill_Point( *centroid, x / (VIO_Real) n_points,
                           y / (VIO_Real) n_points,
                           z / (VIO_Real) n_points );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : reverse_vectors
@INPUT      : n_vectors
              vectors
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Reverses the direction of each vector in the list.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI   void     reverse_vectors(
    int       n_vectors,
    VIO_Vector    vectors[] )
{
    int     i;

    for_less( i, 0, n_vectors )
        SCALE_VECTOR( vectors[i], vectors[i], -1.0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_angle_between_points
@INPUT      : prev_point
              this_one
              next_point
@OUTPUT     : 
@RETURNS    : angle in radians.
@DESCRIPTION: Returns the angle between the three points, in the range 0 to
              M_PI.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  VIO_Real  get_angle_between_points(
    VIO_Point  *prev_point,
    VIO_Point  *this_point,
    VIO_Point  *next_point )
{
    VIO_Real    angle, c;
    VIO_Vector  v1, v2;

    SUB_POINTS( v1, *prev_point, *this_point );
    SUB_POINTS( v2, *next_point, *this_point );

    NORMALIZE_VECTOR( v1, v1 );
    NORMALIZE_VECTOR( v2, v2 );

    c = DOT_VECTORS( v1, v2 );

    if( c >= 1.0 )
        angle = 0.0;
    else if( c <= -1.0 )
        angle = M_PI;
    else
        angle = acos( c );

    return( angle );
}
