#include  <bicpl.h>

#ifdef ADJUST
static  void  adjust_tetrahedral_sphere(
    VIO_Point            *centre,
    VIO_Real             rx,
    VIO_Real             ry,
    VIO_Real             rz,
    polygons_struct  *polygons,
    int              n_iters );
#endif

int  main(
    int    argc,
    char   *argv[] )
{
    VIO_Status          status;
    VIO_STR          output_filename;
    int             n_triangles;
#ifdef ADJUST
    int             n_adjustment_iters;
#endif
    VIO_Point           centre;
    object_struct   *object;
    polygons_struct *polygons;
    VIO_Real            cx, cy, cz, rx, ry, rz;

    status = VIO_OK;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( "", &output_filename ) ||
        !get_real_argument( 0.0, &cx ) ||
        !get_real_argument( 0.0, &cy ) ||
        !get_real_argument( 0.0, &cz ) ||
        !get_real_argument( 0.0, &rx ) ||
        !get_real_argument( 0.0, &ry ) ||
        !get_real_argument( 0.0, &rz ) )
    {
        fprintf(stderr, "Usage: %s filename cx cy cz rx ry rz n_triangles.\n",
                     argv[0] );
        return( 1 );
    }

    (void) get_int_argument( 4, &n_triangles );
#ifdef ADJUST
    (void) get_int_argument( 0, &n_adjustment_iters );
#endif

    fill_Point( centre, cx, cy, cz );

    object = create_object( POLYGONS );
    polygons = get_polygons_ptr(object);
    create_tetrahedral_sphere( &centre, rx, ry, rz, n_triangles, polygons );

#ifdef ADJUST
    if( n_adjustment_iters > 0 )
    {
        adjust_tetrahedral_sphere( &centre, rx, ry, rz, polygons,
                                   n_adjustment_iters );
    }
#endif

    compute_polygon_normals( polygons );

    status = output_graphics_file( output_filename, ASCII_FORMAT, 1, &object );

    if( status == VIO_OK )
        delete_object( object );

    if( status == VIO_OK )
        print( "Tetrahedron output.\n" );

    return( status != VIO_OK );
}

#ifdef ADJUST
static  void  convert_uv_to_point(
    VIO_Point            *centre,
    VIO_Real             rx,
    VIO_Real             ry,
    VIO_Real             rz,
    VIO_Real             u,
    VIO_Real             v,
    VIO_Point            *point )
{
    VIO_Real   x, y, z, theta, phi, cos_u, sin_u, cos_v, sin_v;

    theta = u * PI * 2.0;
    phi = v * PI;

    cos_u = cos( theta );
    sin_u = sin( theta );
    cos_v = cos( phi );
    sin_v = sin( phi );

    z = (VIO_Real) Point_z(*centre) - rz * cos_v;
    x = (VIO_Real) Point_x(*centre) + rx * sin_v * cos_u;
    y = (VIO_Real) Point_y(*centre) + ry * sin_v * sin_u;

    fill_Point( *point, x, y, z );
}

static  void  convert_point_to_uv(
    VIO_Point            *centre,
    VIO_Real             rx,
    VIO_Real             ry,
    VIO_Real             rz,
    VIO_Point            *point,
    VIO_Real             *u,
    VIO_Real             *v )
{
    VIO_Real   x, y, z, theta, phi;

    x = ((VIO_Real) Point_x(*point) - (VIO_Real) Point_x(*centre)) / rx;
    y = ((VIO_Real) Point_y(*point) - (VIO_Real) Point_y(*centre)) / ry;
    z = ((VIO_Real) Point_z(*point) - (VIO_Real) Point_z(*centre)) / rz;

    phi = acos( -z );
    theta = 2.0 * PI - compute_clockwise_rotation( x, y );

    *u = theta / PI / 2.0;
    *v = phi / PI;
}

static  VIO_Real  evaluate_rms(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               **neighbours,
    VIO_Real              parameters[] )
{
    int    p, p1;
    VIO_Real   dist, sum_dist;
    VIO_Point  centre = { 0.0f, 0.0f, 0.0f };

    for_less( p, 0, polygons->n_points )
    {
        convert_uv_to_point( &centre, 1.0, 1.0, 1.0,
                             parameters[2*p], parameters[2*p+1],
                             &polygons->points[p] );
    }

    sum_dist = 0.0;

    for_less( p, 0, polygons->n_points-1 )
    {
        for_less( p1, p+1, polygons->n_points )
        {
            dist = distance_between_points( &polygons->points[p1],
                                            &polygons->points[p] );
            if( dist == 0.0 )
                dist = 1.0e30;
            else
                dist = 1.0 / dist / dist;

            sum_dist += dist;
        }
    }

    return( sum_dist );
}
    

#define  TOLERANCE  1.0e-10

typedef  struct
{
    polygons_struct  *polygons;
    int              *n_neighbours;
    int              **neighbours;
} func_data;

static  VIO_Real  evaluate_rms_func(
    VIO_Real   coefs[],
    void   *void_data )
{
    func_data   *data;

    data = (func_data *) void_data;

    return( evaluate_rms( data->polygons,
                          data->n_neighbours, data->neighbours, coefs ) );
}

static  void  adjust_tetrahedral_sphere(
    VIO_Point            *centre,
    VIO_Real             rx,
    VIO_Real             ry,
    VIO_Real             rz,
    polygons_struct  *polygons,
    int              n_iters )
{
    int         p, poly, size, vertex;
    int         *n_neighbours, **neighbours, *total_neighbours;
    int         total_n_neighbours;
    func_data   data;
    VIO_Real        *step_sizes, *improved, *initial, rms;
    VIO_BOOL     interior;

    ALLOC( step_sizes, polygons->n_points * 2 );
    ALLOC( initial, polygons->n_points * 2 );
    ALLOC( improved, polygons->n_points * 2 );

    for_less( p, 0, polygons->n_points )
    {
        step_sizes[2*p] = 1.0e-2;
        step_sizes[2*p+1] = 1.0e-2;
        convert_point_to_uv( centre, rx, ry, rz, &polygons->points[p],
                             &initial[2*p],
                             &initial[2*p+1] );
    }

    check_polygons_neighbours_computed( polygons );

    ALLOC( n_neighbours, polygons->n_points );
    ALLOC( neighbours, polygons->n_points );
    for_less( p, 0, polygons->n_points )
        n_neighbours[p] = 0;

    total_n_neighbours = 0;
    for_less( poly, 0, polygons->n_items )
    {
        size = GET_OBJECT_SIZE( *polygons, poly );
        for_less( vertex, 0, size )
        {
            p = polygons->indices[
                  POINT_INDEX(polygons->end_indices,poly,vertex)];
            if( n_neighbours[p] > 0 )
                continue;

            n_neighbours[p] = get_neighbours_of_point( polygons, poly,
                                                       vertex,
                                                       NULL, 0, &interior );
            total_n_neighbours += n_neighbours[p];
        }
    }

    ALLOC( total_neighbours, total_n_neighbours );
    total_n_neighbours = 0;
    for_less( p, 0, polygons->n_points )
    {
        neighbours[p] = &total_neighbours[total_n_neighbours];
        total_n_neighbours += n_neighbours[p];
    }

    for_less( p, 0, polygons->n_points )
        n_neighbours[p] = 0;

    for_less( poly, 0, polygons->n_items )
    {
        size = GET_OBJECT_SIZE( *polygons, poly );
        for_less( vertex, 0, size )
        {
            p = polygons->indices[
                  POINT_INDEX(polygons->end_indices,poly,vertex)];

            if( n_neighbours[p] > 0 )
                continue;

            n_neighbours[p] = get_neighbours_of_point( polygons, poly, vertex,
                                                       neighbours[p],
                                                       total_n_neighbours,
                                                       &interior );
        }
    }

    rms = evaluate_rms( polygons, n_neighbours, neighbours, initial );
    print( "Initial: %g\n", rms );

    data.polygons = polygons;
    data.n_neighbours = n_neighbours;
    data.neighbours = neighbours;

    (void) gradient_steps_minimize_function( polygons->n_points * 2, initial,
                                             step_sizes,
                                             evaluate_rms_func,
                                             (void *) &data,
                                             2, n_iters,
                                             TOLERANCE, improved );

    rms = evaluate_rms( polygons, n_neighbours, neighbours, improved );
    print( "Final: %g\n", rms );

    for_less( p, 0, polygons->n_points )
    {
        convert_uv_to_point( centre, rx, ry, rz,
                             improved[2*p], improved[2*p+1],
                             &polygons->points[p] );
    }

    FREE( n_neighbours );
    FREE( neighbours );
    FREE( total_neighbours );
}
#endif
