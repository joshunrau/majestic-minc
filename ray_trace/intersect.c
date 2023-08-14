#if HAVE_CONFIG_H
#include "config.h"
#endif
#include  <ray_trace.h>

static  void  interpolate_over_polygon(
    VIO_Point        *point,
    int              n_points,
    VIO_Point        points[],
    VIO_BOOL         flat_shading_flag,
    VIO_Vector       normals[],
    Colour_flags     colour_flag,
    VIO_Colour       colours[],
    VIO_Vector       *normal,
    VIO_Colour       *colour )
{
    int              i;
    VIO_Real             weights[MAX_POINTS_PER_POLYGON];
    VIO_Real             r, g, b, a;
    VIO_Vector           scaled_normal;

    if( colour_flag == PER_VERTEX_COLOURS )
    {
        r = 0.0;
        g = 0.0;
        b = 0.0;
        a = 0.0;
    }

    get_polygon_interpolation_weights( point, n_points, points, weights );

    if( flat_shading_flag )
        find_polygon_normal( n_points, points, normal );
    else
        fill_Vector( *normal, 0.0, 0.0, 0.0 );

    for_less( i, 0, n_points )
    {
        if( colour_flag == PER_VERTEX_COLOURS )
        {
            r += weights[i] * get_Colour_r_0_1(colours[i] );
            g += weights[i] * get_Colour_g_0_1(colours[i] );
            b += weights[i] * get_Colour_b_0_1(colours[i] );
            a += weights[i] * get_Colour_a_0_1(colours[i] );
        }

        if( !flat_shading_flag )
        {
            SCALE_VECTOR( scaled_normal, normals[i], weights[i] );
            ADD_VECTORS( *normal, *normal, scaled_normal );
        }
    }

    NORMALIZE_VECTOR( *normal, *normal );

    if( colour_flag == PER_VERTEX_COLOURS )
        *colour = make_rgba_Colour_0_1( r, g, b, a );
}

  VIO_BOOL  ray_intersects_a_polygon(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    object_struct    *object,
    VIO_BOOL          flat_shading_flag,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Real             *dist )
{
    polygons_struct  *polygons;
    int              size, i, object_index, point_index;
    VIO_Point            points[MAX_POINTS_PER_POLYGON];
    VIO_Vector           normals[MAX_POINTS_PER_POLYGON];
    VIO_Colour           colours[MAX_POINTS_PER_POLYGON];
    VIO_BOOL          found;

    found = intersect_ray_with_object( origin, direction, object,
                                       &object_index, dist,
                                       (VIO_Real **) NULL ) > 0;

    if( found && point != (VIO_Point *) NULL )
    {
        polygons = get_polygons_ptr( object );

        GET_POINT_ON_RAY( *point, *origin, *direction, *dist );

        if( polygons->colour_flag == ONE_COLOUR )
            *colour = polygons->colours[0];
        else if( polygons->colour_flag == PER_ITEM_COLOURS )
            *colour = polygons->colours[object_index];

        size = get_polygon_points( polygons, object_index, points );

        for_less( i, 0, size )
        {
            point_index = polygons->indices[
                            POINT_INDEX(polygons->end_indices,object_index,i)];

            normals[i] = polygons->normals[point_index];

            if( polygons->colour_flag == PER_VERTEX_COLOURS )
                colours[i] = polygons->colours[point_index];
        }

        interpolate_over_polygon( point, size, points, flat_shading_flag,
                                  normals, polygons->colour_flag,
                                  colours, normal, colour );
    }

    return( found );
}

  VIO_BOOL  ray_intersects_a_quadmesh(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    object_struct    *object,
    VIO_BOOL          flat_shading_flag,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Real             *dist )
{
    quadmesh_struct  *quadmesh;
    int              i, j, m, n, p, object_index;
    int              indices[4];
    VIO_Point            points[4];
    VIO_Vector           normals[4];
    VIO_Colour           colours[4];
    VIO_BOOL          found;

    found = intersect_ray_with_object( origin, direction, object,
                                       &object_index, dist,
                                       (VIO_Real **) NULL ) > 0;

    if( found && point != (VIO_Point *) NULL )
    {
        quadmesh = get_quadmesh_ptr( object );

        GET_POINT_ON_RAY( *point, *origin, *direction, *dist );

        if( quadmesh->colour_flag == ONE_COLOUR )
            *colour = quadmesh->colours[0];
        else if( quadmesh->colour_flag == PER_ITEM_COLOURS )
            *colour = quadmesh->colours[object_index];

        get_quadmesh_n_objects( quadmesh, &m, &n );

        i = object_index / n;
        j = object_index % n;

        get_quadmesh_patch_indices( quadmesh, i, j, indices );

        for_less( p, 0, 4 )
        {
            points[p] = quadmesh->points[indices[p]];
            normals[p] = quadmesh->normals[indices[p]];

            if( quadmesh->colour_flag == PER_VERTEX_COLOURS )
                colours[p] = quadmesh->colours[indices[p]];
        }

        interpolate_over_polygon( point, 4, points, flat_shading_flag,
                                  normals, quadmesh->colour_flag,
                                  colours, normal, colour );
    }

    return( found );
}

static  VIO_Real   get_line_alpha(
    VIO_Point     *point,
    VIO_Point     *p1,
    VIO_Point     *p2 )
{
    VIO_Real    alpha, d;
    VIO_Vector  v, offset;

    SUB_POINTS( v, *p2, *p1 );
    SUB_POINTS( offset, *point, *p1 );

    d = DOT_VECTORS( v, v );

    if( d == 0.0 )
        return( 0.0 );

    alpha = DOT_VECTORS( offset, v ) / d;

    if( alpha < 0.0 )
        alpha = 0.0;
    else if( alpha > 1.0 )
        alpha = 1.0;

    return( alpha );
}

  VIO_BOOL  ray_intersects_a_line(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    object_struct    *object,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Real             *dist )
{
    int           line, seg, p1, p2, object_index;
    VIO_Real          alpha, r1, g1, b1, a1, r2, g2, b2, a2;
    VIO_BOOL       found;
    VIO_Point         centre, point1, point2;
    lines_struct  *lines;

    found = intersect_ray_with_object( origin, direction, object,
                                       &object_index, dist,
                                       (VIO_Real **) NULL ) > 0;

    if( found && point != (VIO_Point *) NULL )
    {
        lines = get_lines_ptr( object );

        get_line_segment_index( lines, object_index, &line, &seg );

        GET_POINT_ON_RAY( *point, *origin, *direction, *dist );

        if( lines->colour_flag == ONE_COLOUR )
            *colour = lines->colours[0];
        else if( lines->colour_flag == PER_ITEM_COLOURS )
            *colour = lines->colours[line];

        p1 = lines->indices[POINT_INDEX(lines->end_indices,line,seg)];
        p2 = lines->indices[POINT_INDEX(lines->end_indices,line,seg+1)];
        point1 = lines->points[p1];
        point2 = lines->points[p2];

        alpha = get_line_alpha( point, &point1, &point2 );

        INTERPOLATE_POINTS( centre, point1, point2, alpha );
        SUB_POINTS( *normal, *point, centre );
        NORMALIZE_VECTOR( *normal, *normal );

        if( lines->colour_flag == PER_VERTEX_COLOURS )
        {
            r1 = get_Colour_r_0_1( lines->colours[p1] );
            g1 = get_Colour_g_0_1( lines->colours[p1] );
            b1 = get_Colour_b_0_1( lines->colours[p1] );
            a1 = get_Colour_a_0_1( lines->colours[p1] );
            r2 = get_Colour_r_0_1( lines->colours[p2] );
            g2 = get_Colour_g_0_1( lines->colours[p2] );
            b2 = get_Colour_b_0_1( lines->colours[p2] );
            a2 = get_Colour_a_0_1( lines->colours[p2] );

            *colour = make_rgba_Colour_0_1( r1 + (r2 - r1) * alpha,
                                            g1 + (g2 - g1) * alpha,
                                            b1 + (b2 - b1) * alpha,
                                            a1 + (a2 - a1) * alpha );
        }
    }

    return( found );
}

  VIO_BOOL  ray_intersects_a_marker(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    object_struct    *object,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Real             *dist )
{
    VIO_Real           orig, delta, t_min, t_max, t1, t2, pos1, pos2;
    int            object_index, dim;
    VIO_BOOL        found;
    marker_struct  *marker;

    marker = get_marker_ptr( object );

    t_min = 0.0;
    t_max = 0.0;

    for_less( dim, 0, VIO_N_DIMENSIONS )
    {
        orig = RPoint_coord( *origin, dim );
        delta = RVector_coord( *direction, dim );
        pos1 = RPoint_coord( marker->position, dim ) - marker->size / 2.0;
        pos2 = RPoint_coord( marker->position, dim ) + marker->size / 2.0;
        if( delta < 0.0 )
        {
            t1 = (pos2 - orig) / delta;
            t2 = (pos1 - orig) / delta;
        }
        else if( delta > 0.0 )
        {
            t1 = (pos1 - orig) / delta;
            t2 = (pos2 - orig) / delta;
        }
        else
        {
            if( orig < pos1 || orig > pos2 )
                return( FALSE );
        }

        if( t1 > t_min )
            t_min = t1;
        if( dim == 0 || t2 < t_max )
            t_max = t2;
    }

    if( t1 > t2 )
        return( FALSE );

    found = intersect_ray_with_object( origin, direction, object,
                                       &object_index, dist,
                                       (VIO_Real **) NULL ) > 0;

    if( found && point != (VIO_Point *) NULL )
    {
        marker = get_marker_ptr( object );

        GET_POINT_ON_RAY( *point, *origin, *direction, *dist );

        *colour = marker->colour;

        SUB_POINTS( *normal, *point, marker->position );
        NORMALIZE_VECTOR( *normal, *normal );
    }

    return( found );
}
