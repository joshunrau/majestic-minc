#if HAVE_CONFIG_H
#include "config.h"
#endif
#include  <ray_trace.h>

static  VIO_BOOL  ray_intersects_objects(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    int              n_objects,
    ray_trace_object objects[],
    VIO_BOOL          use_flags[],
    int              *obj_index,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Surfprop         *spr,
    VIO_Real             *closest_dist );

static  void  get_ray(
    view_struct  *view,
    VIO_Real         x_pixel,
    VIO_Real         y_pixel,
    VIO_Point        *origin,
    VIO_Vector       *direction )
{
    VIO_Vector   x_offset, y_offset, z_offset, offset;
    VIO_Real     x_position, y_position;
    VIO_Point    screen_position;

    x_position = -view->window_width / 2.0 + view->window_width *
                 (x_pixel + 0.5) / (VIO_Real) view->x_viewport_size;
    y_position = -view->window_height / 2.0 + view->window_height *
                 (y_pixel + 0.5) / (VIO_Real) view->y_viewport_size;

    if( view->perspective_flag )
    {
        fill_Vector( screen_position,
                     RPoint_x(view->eye) +
                     x_position * RVector_x(view->horizontal) +
                     y_position * RVector_x(view->vertical) +
                     view->window_distance * RVector_x(view->line_of_sight),
                     RPoint_y(view->eye) +
                     x_position * RVector_y(view->horizontal) +
                     y_position * RVector_y(view->vertical) +
                     view->window_distance * RVector_y(view->line_of_sight),
                     RPoint_z(view->eye) +
                     x_position * RVector_z(view->horizontal) +
                     y_position * RVector_z(view->vertical) +
                     view->window_distance * RVector_z(view->line_of_sight) );

        SCALE_VECTOR( offset, view->horizontal, view->stereo_offset );
        ADD_POINT_VECTOR( *origin, view->eye, offset );
        SUB_POINTS( *direction, screen_position, *origin );
        NORMALIZE_VECTOR( *direction, *direction );
    }
    else
    {
        SCALE_VECTOR( x_offset, view->horizontal, x_position );
        SCALE_VECTOR( y_offset, view->vertical, y_position );
        SCALE_VECTOR( z_offset, view->line_of_sight, view->window_distance );
        ADD_POINT_VECTOR( screen_position, view->eye, z_offset );
        ADD_POINT_VECTOR( screen_position, screen_position, x_offset );
        ADD_POINT_VECTOR( screen_position, screen_position, y_offset );

        SUB_POINT_VECTOR( *origin, screen_position, z_offset );
        *direction = z_offset;
        NORMALIZE_VECTOR( *direction, *direction );
    }
}

static  void  get_light_component(
    int               n_objects,
    ray_trace_object  objects[],
    light_struct      *light,
    VIO_BOOL           doing_shadows,
    VIO_Real              shadow_offset,
    VIO_Vector            *view_direction,
    VIO_Point             *point,
    VIO_Vector            *normal,
    VIO_Surfprop          *spr,
    VIO_Real              r,
    VIO_Real              g,
    VIO_Real              b,
    VIO_Real              *r_comp,
    VIO_Real              *g_comp,
    VIO_Real              *b_comp )
{
    VIO_Real            n_dot_l, n_dot_h, attenuation, dist, denominator;
    VIO_Vector          highlight, to_light;

    *r_comp = 0.0;
    *g_comp = 0.0;
    *b_comp = 0.0;

    switch( light->type )
    {
    case  DIRECTIONAL_LIGHT:
    case  POINT_LIGHT:

        if( light->type == DIRECTIONAL_LIGHT )
        {
            SCALE_VECTOR( to_light, light->direction, -1.0 );
        }
        else
        {
            SUB_POINTS( to_light, light->position, *point );
            NORMALIZE_VECTOR( to_light, to_light );
        }

        n_dot_l = DOT_VECTORS( *normal, to_light );

        if( n_dot_l <= 0.0 )
            break;

        if( doing_shadows )
        {
            VIO_Point           origin;
            VIO_Vector          offset;

            SCALE_VECTOR( offset, to_light, shadow_offset );
            ADD_POINT_VECTOR( origin, *point, offset );

            if( ray_intersects_objects( &origin, &to_light, n_objects, objects,
                                        NULL, (int *) NULL,
                                        (VIO_Point *) NULL, (VIO_Vector *) NULL,
                                        (VIO_Colour *) NULL, (VIO_Surfprop *) NULL,
                                        &dist ) &&
                (light->type == DIRECTIONAL_LIGHT ||
                 dist < distance_between_points( point, &light->position )) )
            {
                break;
            }
        }

        if( n_dot_l > 0.0 )
        {
            *r_comp = r * n_dot_l * (VIO_Real) Surfprop_d( *spr );
            *g_comp = g * n_dot_l * (VIO_Real) Surfprop_d( *spr );
            *b_comp = b * n_dot_l * (VIO_Real) Surfprop_d( *spr );
        }

        SUB_VECTORS( highlight, *view_direction, to_light );
        NORMALIZE_VECTOR( highlight, highlight );
        n_dot_h = -DOT_VECTORS( *normal, highlight );

        if( n_dot_h > 0.0 )
        {
            n_dot_h = pow( n_dot_h, (VIO_Real) Surfprop_se( *spr ) );
            *r_comp += r * n_dot_h * (VIO_Real) Surfprop_s( *spr );
            *g_comp += g * n_dot_h * (VIO_Real) Surfprop_s( *spr );
            *b_comp += b * n_dot_h * (VIO_Real) Surfprop_s( *spr );
        }

        if( light->type == POINT_LIGHT &&
            (light->attenuation[0] != 0.0 ||
             light->attenuation[1] != 0.0 ||
             light->attenuation[2] != 0.0) )
        {
            dist = distance_between_points( &light->position, point );
            denominator = light->attenuation[0] + dist * (
                          light->attenuation[1] + dist * light->attenuation[2]);
      
            if( denominator <= 1.0 )
                attenuation = 1.0;
            else
                attenuation = 1.0 / denominator;
        }
        else
            attenuation = 1.0;

        *r_comp *= attenuation * get_Colour_r_0_1( light->colour );
        *g_comp *= attenuation * get_Colour_g_0_1( light->colour );
        *b_comp *= attenuation * get_Colour_b_0_1( light->colour );

        break;
    }
}

static  VIO_BOOL  ray_intersects_regular_object(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    object_struct    *object,
    render_struct    *render,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Surfprop         *spr,
    VIO_Real             *dist )
{
    model_struct      *model;
    VIO_BOOL           found, flat_shading;
    int               i;
    VIO_Point             a_point;
    VIO_Vector            a_normal;
    VIO_Colour            a_colour;
    VIO_Surfprop          a_spr;
    VIO_Real              this_dist;
    static  VIO_Surfprop  default_line_spr = { 0.3f, 0.6f, 0.6f, 40.0f, 1.0f };

    if( render == NULL )
        flat_shading = TRUE;
    else
        flat_shading = render->flat_shading_flag;

    switch( get_object_type( object ) )
    {
    case MODEL:
        model = get_model_ptr( object );
        found = FALSE;
        *dist = 0.0;

        for_less( i, 0, model->n_objects )
        {
            if( ray_intersects_regular_object( origin, direction,
                            model->objects[i], render,
                            &a_point, &a_normal, &a_colour, &a_spr, &this_dist))
            {
                if( !found || this_dist < *dist )
                {
                    found = TRUE;

                    *dist = this_dist;

                    if( point != (VIO_Point *) NULL )
                    {
                        *point = a_point;
                        *normal = a_normal;
                        *colour = a_colour;
                        *spr = a_spr;
                    }
                }
            }
        }
        break;

    case POLYGONS:
        found = ray_intersects_a_polygon( origin, direction, object,
                                          flat_shading, point,
                                          normal, colour, dist );
        if( found )
            *spr = *get_object_surfprop(object);
        break;

    case QUADMESH:
        found = ray_intersects_a_quadmesh( origin, direction, object,
                                           flat_shading, point,
                                           normal, colour, dist );
        if( found )
            *spr = *get_object_surfprop(object);
        break;

    case LINES:
        found = ray_intersects_a_line( origin, direction, object,
                                       point, normal, colour, dist );
        if( found )
            *spr = default_line_spr;
        break;

    case MARKER:
        found = ray_intersects_a_marker( origin, direction, object,
                                         point, normal, colour, dist );
        if( found )
            *spr = default_line_spr;
        break;

    default:
        print_error( "ray_intersects_object():  not implemented\n" );
        break;
    }

    return( found );
}

static  void  sort_reals(
    int    n,
    VIO_Real   values[] )
{
    VIO_Real  swap;
    int   i, j, best;

    for_less( i, 0, n-1 )
    {
        best = i;
        for_less( j, i+1, n )
        {
            if( values[j] < values[best] )
                best = j;
        }

        swap = values[i];
        values[i] = values[best];
        values[best] = swap;
    }
}

#define  TOLERANCE  1.0e-6
#define  DEFAULT_SIZE   100

typedef  struct
{
    VIO_Real   *start_intervals;
    VIO_Real   start_intervals_fixed[DEFAULT_SIZE];
    int    start_size;
    VIO_Real   *end_intervals;
    VIO_Real   end_intervals_fixed[DEFAULT_SIZE];
    int    end_size;

    VIO_Real   *new_start_intervals;
    VIO_Real   new_start_intervals_fixed[DEFAULT_SIZE];
    int    new_start_size;
    VIO_Real   *new_end_intervals;
    VIO_Real   new_end_intervals_fixed[DEFAULT_SIZE];
    int    new_end_size;
} intervals_struct;

static  void  initialize_intervals(
    intervals_struct   *is )
{
    is->start_intervals = is->start_intervals_fixed;
    is->end_intervals = is->end_intervals_fixed;
    is->new_start_intervals = is->new_start_intervals_fixed;
    is->new_end_intervals = is->new_end_intervals_fixed;

    is->start_size = DEFAULT_SIZE;
    is->end_size = DEFAULT_SIZE;
    is->new_start_size = DEFAULT_SIZE;
    is->new_end_size = DEFAULT_SIZE;
}

static  void  delete_intervals(
    intervals_struct   *is )
{
    if( is->start_size > DEFAULT_SIZE )
        FREE( is->start_intervals );
    if( is->end_size > DEFAULT_SIZE )
        FREE( is->end_intervals );
    if( is->new_start_size > DEFAULT_SIZE )
        FREE( is->new_start_intervals );
    if( is->new_end_size > DEFAULT_SIZE )
        FREE( is->new_end_intervals );
}

#define  ADD_ELEMENT( intervals, n_alloced, index, value ) \
    { \
        if( (index) >= (n_alloced) ) \
        { \
            if( (n_alloced) == DEFAULT_SIZE ) \
            { \
                VIO_Real  *_ptr; \
                _ptr = NULL;  \
                SET_ARRAY_SIZE( _ptr, 0, (index)+1, DEFAULT_CHUNK_SIZE ); \
                (void) memcpy( (void *) _ptr, (void *) intervals, \
                               (size_t) (index) * sizeof(*(intervals))); \
                (intervals) = _ptr; \
            } \
            else \
            { \
                SET_ARRAY_SIZE( intervals, n_alloced, (index) + 1, \
                                DEFAULT_CHUNK_SIZE ); \
            } \
            (n_alloced) = (index) + 1; \
        } \
        (intervals)[index] = (value); \
        ++(index); \
    }

static  int  get_ray_clip_intervals(
    VIO_Point             *origin,
    VIO_Vector            *direction,
    clip_struct       *clip,
    intervals_struct  *is )
{
    VIO_Real     dist, *intersections;
    VIO_Real     start_invalid, end_invalid, start, end, *tmp;
    int      n, n_intersections, new_n, object_index, i, j, k, n_new_intervals;
    int      n_to_do, tmp_int;
    VIO_BOOL  clip_inside;

    is->start_intervals[0] = 0.0;
    is->end_intervals[0] = 1.0e60;
    n = 1;

    for_less( i, 0, clip->n_clip_objects )
    {
        n_intersections = intersect_ray_with_object( origin, direction,
                                                     clip->clip_objects[i],
                                                     &object_index, &dist,
                                                     &intersections );

        sort_reals( n_intersections, intersections );

        new_n = 0;
        for_less( k, 0, n_intersections )
        {
            if( k == 0 || !numerically_close( intersections[k],
                                              intersections[k-1], TOLERANCE ) )
            {
                intersections[new_n] = intersections[k];
                ++new_n;
            }
        }
        n_intersections = new_n;

#ifdef DEBUG
        if( n_intersections % 2 == 1 )
        {
            print( "N intersections: %d : ", n_intersections );
            for_less( k, 0, n_intersections )
                print( " %g", intersections[k] );
            print( "\n" );
        }
#endif

        clip_inside = clip->clip_inside_flags[i];

        if( n_intersections == 0 && clip_inside || n_intersections % 2 == 1 )
        {
            if( n_intersections > 0 )
                FREE( intersections );
            continue;
        }

        n_to_do = n_intersections / 2;
        if( !clip_inside )
            ++n_to_do;

        start_invalid = 0.0;
        for_less( j, 0, n_to_do )
        {
            if( clip_inside )
            {
                start_invalid = intersections[2*j];
                end_invalid = intersections[2*j+1];
            }
            else
            {
                if( j == 0 )
                    start_invalid = 0.0;
                else
                    start_invalid = intersections[2*j-1];

                if( j < n_intersections/2 )
                    end_invalid = intersections[2*j];
                else
                    end_invalid = 1.0e60;
            }

            n_new_intervals = 0;
            for_less( k, 0, n )
            {
                start = is->start_intervals[k];
                end = is->end_intervals[k];

                if( end <= start_invalid || start >= end_invalid )
                {
                    ADD_ELEMENT( is->new_start_intervals,
                                 is->new_start_size, n_new_intervals, start );
                    --n_new_intervals;
                    ADD_ELEMENT( is->new_end_intervals,
                                 is->new_end_size, n_new_intervals, end );
                }
                else if( start < start_invalid || end > end_invalid )
                {
                    if( start < start_invalid )
                    {
                        ADD_ELEMENT( is->new_start_intervals,
                                     is->new_start_size,
                                     n_new_intervals, start );
                        --n_new_intervals;
                        ADD_ELEMENT( is->new_end_intervals, is->new_end_size,
                                     n_new_intervals, start_invalid );
                    }
                    if( end > end_invalid )
                    {
                        ADD_ELEMENT( is->new_start_intervals,
                             is->new_start_size, n_new_intervals,
                             end_invalid );
                        --n_new_intervals;
                        ADD_ELEMENT( is->new_end_intervals, is->new_end_size,
                                     n_new_intervals, end );
                    }
                }
            }

            tmp = is->start_intervals;
            is->start_intervals = is->new_start_intervals;
            is->new_start_intervals = tmp;

            tmp = is->end_intervals;
            is->end_intervals = is->new_end_intervals;
            is->new_end_intervals = tmp;

            tmp_int = is->start_size;
            is->start_size = is->new_start_size;
            is->new_start_size = tmp_int;

            tmp_int = is->end_size;
            is->end_size = is->new_end_size;
            is->new_end_size = tmp_int;

            n = n_new_intervals;
        }

        if( n_intersections > 0 )
            FREE( intersections );

        if( n == 0 )
            break;
    }

    return( n );
}

static  VIO_BOOL  ray_intersects_object(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    ray_trace_object *object,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Surfprop         *spr,
    VIO_Real             *dist,
    intervals_struct *is )
{
    int      interval, n_intervals;
    VIO_BOOL  found;
    VIO_Point    clip_origin;

    n_intervals = get_ray_clip_intervals( origin, direction,
                                          &object->clip, is );

    found = FALSE;
    for_less( interval, 0, n_intervals ) 
    {
        GET_POINT_ON_RAY( clip_origin, *origin, *direction,
                          is->start_intervals[interval] );

        if( object->regular_object_flag )
        {
            found = ray_intersects_regular_object( &clip_origin, direction,
                        object->object, &object->render,
                        point, normal, colour, spr, dist );
        }
        else
        {
            found = ray_intersects_a_volume( &clip_origin, direction,
                                             object->render.volumes[
                                                 object->render.n_volumes-1].
                                                            volume,
                                             object->render.volumes[
                                                 object->render.n_volumes-1].
                                                            continuity,
                                             &object->done_bits,
                                             &object->surface_bits,
                                             object->threshold, point, normal,
                                             dist );

            if( found )
            {
                *colour = object->colour;
                *spr = object->spr;
            }
        }

        if( found )
        {
            *dist += is->start_intervals[interval];
            if( *dist <= is->end_intervals[interval] )
                break;
            found = FALSE;
        }
    }

    return( found );
}

static  VIO_BOOL  ray_intersects_objects(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    int              n_objects,
    ray_trace_object objects[],
    VIO_BOOL          use_flags[],
    int              *obj_index,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Surfprop         *spr,
    VIO_Real             *closest_dist )
{
    int              i;
    VIO_Point            a_point;
    VIO_Vector           a_normal;
    VIO_Colour           a_colour;
    VIO_Surfprop         a_spr;
    VIO_Real             dist;
    VIO_BOOL          found;
    intervals_struct is;

    found = FALSE;

    *closest_dist = 0.0;

    initialize_intervals( &is );

    for_less( i, 0, n_objects )
    {
        if( (use_flags == NULL || use_flags[i]) &&
            ray_intersects_object( origin, direction, &objects[i],
                                   &a_point, &a_normal, &a_colour, &a_spr,
                                   &dist, &is ) )
        {
            if( !found || dist < *closest_dist )
            {
                found = TRUE;

                *closest_dist = dist;
                if( obj_index != NULL )
                    *obj_index = i;

                if( point != (VIO_Point *) NULL )
                {
                    *point = a_point;
                    *normal = a_normal;
                    *colour = a_colour;
                    *spr = a_spr;
                }
            }
        }
    }

    delete_intervals( &is );

    return( found );
}

static  VIO_BOOL  recursive_ray_trace_one(
    VIO_Point                 *origin,
    VIO_Vector                *direction,
    lights_struct         *lights,
    int                   n_objects,
    ray_trace_object      objects[],
    VIO_BOOL               use_flags[],
    int                   max_depth,
    VIO_Real                  *r,
    VIO_Real                  *g,
    VIO_Real                  *b,
    VIO_Real                  *a )
{
    int             i, v, object_index;
    VIO_Colour          object_colour, col_code, cum_col;
    VIO_Point           intersection_point, new_origin;
    VIO_Vector          normal, offset;
    VIO_Real            voxel[VIO_MAX_DIMENSIONS], value, fill_value;
    VIO_Real            r_cum, g_cum, b_cum, a_cum;
    VIO_Real            r_cc, g_cc, b_cc, a_cc, weight;
    VIO_Real            r_light, g_light, b_light, dist;
    VIO_Real            r_object, g_object, b_object, a_object, opacity;
    VIO_Real            r_behind, g_behind, b_behind, a_behind;
    VIO_Real            dx, dy, dz;
    VIO_Volume          volume;
    VIO_Surfprop        spr;
    VIO_BOOL         using_volume, hit_something;

    hit_something = ray_intersects_objects( origin, direction,
                                n_objects, objects, use_flags,
                                &object_index, &intersection_point, &normal,
                                &object_colour, &spr, &dist );

    if( hit_something )
    {
        if( objects[object_index].render.hit_only_once_flag )
            use_flags[object_index] = FALSE;

        using_volume = FALSE;

        if( objects[object_index].regular_object_flag &&
            objects[object_index].render.n_volumes > 0 )
        {
            r_cum = 0.0;
            g_cum = 0.0;
            b_cum = 0.0;
            a_cum = 0.0;

            for_less( v, 0, objects[object_index].render.n_volumes )
            {
                volume = objects[object_index].render.volumes[v].volume;
                convert_3D_world_to_voxel( volume,
                                           (VIO_Real) Point_x(intersection_point),
                                           (VIO_Real) Point_y(intersection_point),
                                           (VIO_Real) Point_z(intersection_point),
                                           &voxel[VIO_X], &voxel[VIO_Y], &voxel[VIO_Z] );
                if( voxel_is_within_volume( volume, voxel ) ||
                    objects[object_index].render.volumes[v].extend_volume_flag )
                {
                    if( objects[object_index].render.volumes[v].
                        fill_value_specified )
                    {
                        fill_value = objects[object_index].render.volumes[v].
                                       fill_value;
                    }
                    else
                    {
                        fill_value = get_volume_real_min(volume);
                    }

                    if( objects[object_index].render.volumes[v].gradient_flag )
                    {
                        evaluate_volume_in_world(
                                  volume,
                                  (VIO_Real) Point_x(intersection_point),
                                  (VIO_Real) Point_y(intersection_point),
                                  (VIO_Real) Point_z(intersection_point),
                                  objects[object_index].render.volumes[v].
                                                              continuity,
                                  FALSE,
                                  fill_value,
                                  &value, &dx, &dy, &dz,
                                  NULL, NULL, NULL, NULL, NULL, NULL );

                        value = sqrt( dx * dx + dy * dy + dz * dz );
                    }
                    else
                    {
                        (void) evaluate_volume( volume, voxel, NULL,
                                     objects[object_index].render.volumes[v].
                                                               continuity,
                                     FALSE, fill_value,
                                     &value, NULL, NULL );
                    }

                    if( is_an_rgb_volume( volume ) )
                        col_code = (VIO_Colour) value;
                    else
                    {
                        col_code = get_colour_code( &objects[object_index].
                                     render.volumes[v].colour_coding, value );
                    }

                    opacity = objects[object_index].render.volumes[v].opacity;
                    r_cc = get_Colour_r_0_1( col_code );
                    g_cc = get_Colour_g_0_1( col_code );
                    b_cc = get_Colour_b_0_1( col_code );
                    a_cc = opacity * get_Colour_a_0_1( col_code );

                    weight = (1.0 - a_cc) * a_cum;
                    r_cum = weight * r_cum + a_cc * r_cc;
                    g_cum = weight * g_cum + a_cc * g_cc;
                    b_cum = weight * b_cum + a_cc * b_cc;
                    a_cum = weight + a_cc;
                    if( a_cum > 0.0 )
                    {
                        r_cum /= a_cum;
                        g_cum /= a_cum;
                        b_cum /= a_cum;
                    }

                    using_volume = TRUE;
                }
            }

            if( using_volume )
            {
                if( objects[object_index].render.mult_volume_flag )
                {
                    r_object = get_Colour_r_0_1( object_colour );
                    g_object = get_Colour_g_0_1( object_colour );
                    b_object = get_Colour_b_0_1( object_colour );
                    a_object = get_Colour_a_0_1( object_colour );
                    r_cum *= r_object;
                    g_cum *= g_object;
                    b_cum *= b_object;
                    a_cum *= a_object;
                    object_colour = make_rgba_Colour_0_1(
                                          r_cum, g_cum, b_cum, a_cum );
                }
                else if( objects[object_index].render.composite_volume_flag )
                {
                    cum_col = make_rgba_Colour_0_1( r_cum, g_cum, b_cum, a_cum);
                    COMPOSITE_COLOURS( object_colour, cum_col, object_colour );
                }
                else
                {
                    object_colour = make_rgba_Colour_0_1(
                                          r_cum, g_cum, b_cum, a_cum );
                }
            }
        }

        if( !objects[object_index].render.lighting_flag )
        {
            *r = get_Colour_r_0_1( object_colour );
            *g = get_Colour_g_0_1( object_colour );
            *b = get_Colour_b_0_1( object_colour );
            *a = get_Colour_a_0_1( object_colour );
        }
        else
        {
            if( DOT_VECTORS( normal, *direction ) > 0.0 )
                SCALE_VECTOR( normal, normal, -1.0 );

            r_object = get_Colour_r_0_1( object_colour );
            g_object = get_Colour_g_0_1( object_colour );
            b_object = get_Colour_b_0_1( object_colour );

            *a = (VIO_Real) Surfprop_t(spr) * get_Colour_a_0_1(object_colour);

            *r = r_object * get_Colour_r_0_1( lights->ambient_light ) *
                 (VIO_Real) Surfprop_a(spr);
            *g = g_object * get_Colour_g_0_1( lights->ambient_light ) *
                 (VIO_Real) Surfprop_a(spr);
            *b = b_object * get_Colour_b_0_1( lights->ambient_light ) *
                 (VIO_Real) Surfprop_a(spr);

            for_less( i, 0, lights->n_lights )
            {
                get_light_component( n_objects, objects,
                                 &lights->lights[i],
                                 objects[object_index].render.shadow_state,
                                 objects[object_index].render.shadow_offset,
                                 direction,
                                 &intersection_point, &normal, &spr,
                                 r_object, g_object, b_object,
                                 &r_light, &g_light, &b_light );

                *r += r_light;
                *g += g_light;
                *b += b_light;
            }
        }

        if( *a < 1.0 )
        {
            SCALE_VECTOR( offset, *direction, dist +
                                   objects[object_index].render.shadow_offset );
            ADD_POINT_VECTOR( new_origin, *origin, offset );

            if( recursive_ray_trace_one( &new_origin, direction,
                                         lights, n_objects, objects, use_flags,
                                         max_depth - 1,
                                         &r_behind, &g_behind, &b_behind,
                                         &a_behind ) )
            {
                weight = (1.0 - *a) * a_behind;
                *r = weight * r_behind + *a * (*r);
                *g = weight * g_behind + *a * (*g);
                *b = weight * b_behind + *a * (*b);
                *a = weight + *a;

                if( *a > 0.0 )
                {
                    r_cum /= *a;
                    g_cum /= *a;
                    b_cum /= *a;
                }
            }
        }
    }

    return( hit_something );
}

  VIO_Colour  trace_ray(
    view_struct           *view,
    lights_struct         *lights,
    int                   n_objects,
    ray_trace_object      objects[],
    VIO_BOOL               use_flags[],
    int                   x_pixel,
    int                   y_pixel,
    int                   max_depth,
    int                   x_super_sampling,
    int                   y_super_sampling )
{
    VIO_Real     r_sum, g_sum, b_sum, a_sum, r, g, b, a, x_sub_pixel, y_sub_pixel;
    int      x, y, o;
    VIO_Colour   col;
    VIO_Point    origin;
    VIO_Vector   direction;

    r_sum = 0.0;
    g_sum = 0.0;
    b_sum = 0.0;
    a_sum = 0.0;

    for_less( x, 0, x_super_sampling )
    {
        x_sub_pixel = (VIO_Real) x_pixel - 0.5 + ((VIO_Real) x + 0.5) /
                      (VIO_Real) x_super_sampling;
        for_less( y, 0, y_super_sampling )
        {
            y_sub_pixel = (VIO_Real) y_pixel - 0.5 + ((VIO_Real) y + 0.5) /
                          (VIO_Real) y_super_sampling;

            get_ray( view, x_sub_pixel, y_sub_pixel, &origin, &direction );

            for_less( o, 0, n_objects )
                use_flags[o] = TRUE;

            if( recursive_ray_trace_one( &origin, &direction,
                                         lights, n_objects, objects, use_flags,
                                         max_depth, &r, &g, &b, &a ) )
            {
                r_sum += r * a;
                g_sum += g * a;
                b_sum += b * a;
                a_sum += a;
            }
        }
    }

    if( x_super_sampling != 1 || y_super_sampling != 1 )
    {
        if( a_sum > 0.0 )
        {
            r_sum /= a_sum;
            g_sum /= a_sum;
            b_sum /= a_sum;
            a_sum /= (VIO_Real) (x_super_sampling * y_super_sampling);
        }
    }

    if( r_sum < 0.0 )
        r_sum = 0.0;
    else if( r_sum > 1.0 )
        r_sum = 1.0;
    if( g_sum < 0.0 )
        g_sum = 0.0;
    else if( g_sum > 1.0 )
        g_sum = 1.0;
    if( b_sum < 0.0 )
        b_sum = 0.0;
    else if( b_sum > 1.0 )
        b_sum = 1.0;
    if( a_sum < 0.0 )
        a_sum = 0.0;
    else if( a_sum > 1.0 )
        a_sum = 1.0;

    col = make_rgba_Colour_0_1( r_sum, g_sum, b_sum, a_sum );

    return( col );
}

static  void  ray_trace_scanline(
    view_struct           *view,
    lights_struct         *lights,
    int                   n_objects,
    ray_trace_object      objects[],
    int                   max_depth,
    int                   x_super_sampling,
    int                   y_super_sampling,
    int                   y,
    pixels_struct         *pixels )
{
    int              x;
    VIO_BOOL          *use_flags;

    ALLOC( use_flags, n_objects );

    for_less( x, 0, view->x_viewport_size )
    {
        PIXEL_RGB_COLOUR( *pixels, x, y ) = trace_ray( view, lights,
                                  n_objects, objects, use_flags, x, y,
                                  max_depth,
                                  x_super_sampling, y_super_sampling );

    }

    FREE( use_flags );
}

typedef  struct
{
    view_struct           *view;
    lights_struct         *lights;
    int                   n_objects;
    ray_trace_object      *objects;
    int                   max_depth;
    int                   x_super_sampling;
    int                   y_super_sampling;
    int                   y;
    int                   y_increment;
    int                   y_size;
    pixels_struct         *pixels;
} scan_line_struct;

static  void  ray_trace_a_scanline(
    int     y,
    void    *data)
{
    int               y_start, y_end;
    scan_line_struct  *scan_data;

    scan_data = (scan_line_struct *) data;

    y_start = y;
    y_end = MIN( y + scan_data->y_increment, scan_data->y_size );

    for_less( y, y_start, y_end )
    {
        ray_trace_scanline( scan_data->view,
                            scan_data->lights,
                            scan_data->n_objects,
                            scan_data->objects,
                            scan_data->max_depth,
                            scan_data->x_super_sampling,
                            scan_data->y_super_sampling,
                            y,
                            scan_data->pixels );
    }
}

  void  ray_trace_scene(
    view_struct           *view,
    lights_struct         *lights,
    int                   n_objects,
    ray_trace_object      objects[],
    int                   max_depth,
    int                   x_super_sampling,
    int                   y_super_sampling,
    pixels_struct         *pixels )
{
    int                y;
    VIO_progress_struct    progress;
    scan_line_struct   data;

    initialize_progress_report( &progress, FALSE, view->y_viewport_size,
                                "Ray Tracing" );
 
    data.view = view;
    data.lights = lights;
    data.n_objects = n_objects;
    data.objects = objects;
    data.max_depth = max_depth;
    data.x_super_sampling = x_super_sampling;
    data.y_super_sampling = y_super_sampling;
    data.pixels = pixels;
    data.y_size = view->y_viewport_size;
    data.y_increment = 1;

    for( y = 0;  y < view->y_viewport_size;  y += data.y_increment )
    {
        ray_trace_a_scanline( y, (void *) &data );

        update_progress_report( &progress, y+1 );
    }

    terminate_progress_report( &progress );

#ifdef DEBUG
    {
        int            x, y;
        VIO_Colour         col1, col2;
        pixels_struct  copy;

        copy_pixel_region( pixels, 0, 100000, 0, 100000, &copy );

        data.y_increment = 1;
        data.pixels = &copy;

        for_less( y, 0, view->y_viewport_size )
        {
            ray_trace_a_scanline( y, (void *) &data );

            for_less( x, 0, view->x_viewport_size )
            {
                col1 = PIXEL_RGB_COLOUR( *pixels, x, y );
                col2 = PIXEL_RGB_COLOUR( copy, x, y );
                if( col1 != col2 )
                {
                    print( "%d %d %d   %d %d %d\n",
                           get_Colour_r(col1),
                           get_Colour_g(col1),
                           get_Colour_b(col1),
                           get_Colour_r(col2),
                           get_Colour_g(col2),
                           get_Colour_b(col2) );
                    handle_internal_error( "colours" );
                }
            }
        }
    }
#endif
}
