#ifndef  DEF_geom_prototypes
#define  DEF_geom_prototypes

#ifdef __cplusplus
extern "C" {
#endif

BICAPI  double  fast_approx_sqrt(
    double  y );

BICAPI  VIO_Real  fast_approx_distance_between_points(
    VIO_Point  *p1,
    VIO_Point  *p2 );

BICAPI  int  clip_polygon_against_box(
    int     n_points,
    VIO_Point   points[],
    VIO_Real    x_min,
    VIO_Real    x_max,
    VIO_Real    y_min,
    VIO_Real    y_max,
    VIO_Real    z_min,
    VIO_Real    z_max,
    int     n_output_points,
    VIO_Point   output_points[] );

BICAPI  int  clip_polygon_against_plane(
    int     n_points,
    VIO_Point   points[],
    VIO_Real    plane_constant,
    VIO_Vector  *normal,
    VIO_Point   output_points[] );

BICAPI  void  split_polygon_with_plane(
    int     n_points,
    VIO_Point   points[],
    VIO_Real    plane_constant,
    VIO_Vector  *normal,
    int     *n_in,
    VIO_Point   in_points[],
    int     *n_out,
    VIO_Point   out_points[] );

BICAPI  void  get_closest_point_on_line_segment(
    VIO_Point  *point,
    VIO_Point  *p1,
    VIO_Point  *p2,
    VIO_Point  *closest_point );

BICAPI  VIO_Real  get_distance_to_line_segment(
    VIO_Point  *point,
    VIO_Point  *p1,
    VIO_Point  *p2,
    VIO_Real   *alpha );

BICAPI  VIO_Real  get_point_object_distance_sq(
    VIO_Point                 *point,
    object_struct         *object,
    int                   obj_index,
    VIO_Point                 *object_point );

BICAPI  VIO_Real  get_point_object_distance(
    VIO_Point                 *point,
    object_struct         *object,
    int                   obj_index,
    VIO_Point                 *object_point );

BICAPI  VIO_Real  get_point_object_vertex_distance(
    VIO_Point                 *point,
    object_struct         *object,
    int                   obj_index,
    int                   *object_vertex );

BICAPI  VIO_Real  find_closest_point_on_object(
    VIO_Point           *point,
    object_struct   *object,
    int             *obj_index,
    VIO_Point           *point_on_object );

BICAPI  VIO_Real  find_closest_vertex_on_object(
    VIO_Point           *point,
    object_struct   *object,
    int             *vertex_on_object );

BICAPI  void  get_polygon_vertex_curvatures(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    VIO_Real              smoothing_distance,
    VIO_Real              low_threshold,
    VIO_Real              curvatures[] );

BICAPI  void  flatten_around_vertex(
    VIO_Point     *vertex,
    int       n_neighbours,
    VIO_Point     neighbours[],
    VIO_BOOL   closed_flag,
    VIO_Real      x_flat[],
    VIO_Real      y_flat[] );

BICAPI  void  flatten_around_vertex_to_sphere(
    VIO_Real      radius,
    VIO_Point     *vertex,
    int       n_neighbours,
    VIO_Point     neighbours[],
    VIO_Real      x_sphere[],
    VIO_Real      y_sphere[],
    VIO_Real      z_sphere[] );

BICAPI  int  compute_distances_from_point(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    VIO_Point             *point,
    int               poly,
    VIO_Real              max_distance,
    VIO_BOOL           distances_initialized,
    float             distances[],
    int               *list[] ) ;

BICAPI  void  find_polygon_normal_no_normalize(
    int      n_points,
    VIO_Point    points[],
    VIO_Real     *nx,
    VIO_Real     *ny,
    VIO_Real     *nz );

BICAPI  void  find_polygon_normal(
    int      n_points,
    VIO_Point    points[],
    VIO_Vector   *normal );

BICAPI  void   get_plane_through_points(
    int      n_points,
    VIO_Point    points[],
    VIO_Vector   *normal,
    VIO_Real     *plane_constant );

BICAPI  VIO_Real  distance_from_plane(
    VIO_Point    *point,
    VIO_Vector   *plane_normal,
    VIO_Real     plane_constant );

BICAPI  VIO_Real  distance_from_line(
    VIO_Point    *point,
    VIO_Point    *end_point1,
    VIO_Point    *end_point2 );

BICAPI  VIO_BOOL   line_segment_intersects_plane(
    VIO_Point   *p1,
    VIO_Point   *p2,
    VIO_Point   *plane_origin,
    VIO_Vector  *plane_normal,
    VIO_Point   *intersection_point );

BICAPI  VIO_BOOL  get_nearest_point_on_lines(
    VIO_Point   *origin1,
    VIO_Vector  *delta1,
    VIO_Point   *origin2,
    VIO_Vector  *delta2,
    VIO_Point   *nearest_point );

BICAPI  VIO_BOOL  clip_line_to_box(
    VIO_Point    *origin,
    VIO_Vector   *direction,
    VIO_Real     x_min,
    VIO_Real     x_max,
    VIO_Real     y_min,
    VIO_Real     y_max,
    VIO_Real     z_min,
    VIO_Real     z_max,
    VIO_Real     *t_min,
    VIO_Real     *t_max );

BICAPI  void  create_line_circle(
    VIO_Point            *centre,
    int              plane_axis,
    VIO_Real             x_radius,
    VIO_Real             y_radius,
    int              n_points,
    lines_struct     *lines );

BICAPI  void  get_polygon_interpolation_weights(
    VIO_Point       *point,
    int         n_points,
    VIO_Point       points[],
    VIO_Real        weights[] );

BICAPI  void  map_point_between_polygons(
    polygons_struct  *p1,
    int              poly_index,
    VIO_Point            *p1_point,
    polygons_struct  *p2,
    VIO_Point            *p2_point );

BICAPI  VIO_Real  map_point_to_unit_sphere(
    polygons_struct  *p,
    VIO_Point            *point,
    polygons_struct  *unit_sphere,
    VIO_Point            *unit_sphere_point );

BICAPI  void  map_unit_sphere_to_point(
    polygons_struct  *unit_sphere,
    VIO_Point            *unit_sphere_point,
    polygons_struct  *p,
    VIO_Point            *point );

BICAPI  void  polygon_transform_points(
    polygons_struct   *src_polygons,
    polygons_struct   *dest_polygons,
    int               n_points,
    VIO_Point             src_points[],
    VIO_Point             dest_points[] );

BICAPI  void  map_sphere_to_uv(
    VIO_Real    x,
    VIO_Real    y,
    VIO_Real    z,
    VIO_Real    *u,
    VIO_Real    *v );

BICAPI  void  map_uv_to_sphere(
    VIO_Real    u,
    VIO_Real    v,
    VIO_Real    *x,
    VIO_Real    *y,
    VIO_Real    *z );

BICAPI   void   find_path_between_polygons(
    int            polygon1,
    int            polygon2,
    int            n_polygons,
    int            end_indices[],
    VIO_SCHAR   visibilities[],
    int            neighbours[],
    VIO_BOOL        *path_exists,
    int            *path_length,
    int            *path[] );

BICAPI  void  create_unit_tetrahedron(
    polygons_struct  *polygons );

BICAPI  void  create_unit_cube(
    polygons_struct  *polygons );

BICAPI  void  create_unit_octohedron(
    polygons_struct  *polygons );

BICAPI  void  create_unit_icosahedron(
    polygons_struct  *polygons );

BICAPI  void   intersect_planes_with_polygons(
    polygons_struct   *polygons,
    VIO_Point             *plane_origin,
    VIO_Vector            *plane_normal,
    lines_struct      *lines );

BICAPI  void   intersect_planes_with_quadmesh(
    quadmesh_struct   *quadmesh,
    VIO_Point             *plane_origin,
    VIO_Vector            *plane_normal,
    lines_struct      *lines );

BICAPI  VIO_BOOL  null_Point(
    VIO_Point   *p );

BICAPI  VIO_BOOL  null_Vector(
    VIO_Vector   *v );

BICAPI  VIO_Real  distance_between_points(
    VIO_Point  *p1,
    VIO_Point  *p2 );

BICAPI  VIO_BOOL  points_within_distance(
    VIO_Point  *p1,
    VIO_Point  *p2,
    VIO_Real   distance );

BICAPI  void  apply_point_to_min_and_max(
    VIO_Point   *point,
    VIO_Point   *min_point,
    VIO_Point   *max_point );

BICAPI  void  expand_min_and_max_points(
    VIO_Point   *min_point,
    VIO_Point   *max_point,
    VIO_Point   *min_to_check,
    VIO_Point   *max_to_check );

BICAPI  void  get_range_points(
    int                n_points,
    VIO_Point              points[],
    VIO_Point              *min_corner,
    VIO_Point              *max_corner );

BICAPI  void  get_points_centroid(
    int     n_points,
    VIO_Point   points[],
    VIO_Point   *centroid );

BICAPI   void     reverse_vectors(
    int       n_vectors,
    VIO_Vector    vectors[] );

BICAPI  VIO_Real  get_angle_between_points(
    VIO_Point  *prev_point,
    VIO_Point  *this_point,
    VIO_Point  *next_point );

BICAPI  VIO_Real  sq_distance_between_points(
    VIO_Point  *p1,
    VIO_Point  *p2 );

BICAPI  VIO_Real  find_point_polygon_distance_sq(
    VIO_Point     *point,
    int       n_points,
    VIO_Point     poly_points[],
    VIO_Point     *closest_point );

BICAPI  VIO_Real  find_point_polygon_distance(
    VIO_Point     *point,
    int       n_points,
    VIO_Point     poly_points[],
    VIO_Point     *closest_point );

BICAPI  int  find_closest_polygon_point(
    VIO_Point              *point,
    polygons_struct    *polygons,
    VIO_Point              *closest_point );

BICAPI  void  create_polygons_sphere(
    VIO_Point            *centre,
    VIO_Real             x_size,
    VIO_Real             y_size,
    VIO_Real             z_size,
    int              n_up,
    int              n_around,
    VIO_BOOL          subdividing_flag,
    polygons_struct  *polygons );

BICAPI  int  get_sphere_point_index(
    int   up,
    int   around,
    int   n_up,
    int   n_around );

BICAPI  VIO_BOOL  is_this_sphere_topology(
    polygons_struct  *polygons );

BICAPI  VIO_BOOL  get_tessellation_of_polygons_sphere(
    polygons_struct  *polygons,
    int              *tess );

BICAPI  int  get_tessellation_with_n_points(
    int   n_points );

BICAPI  void  half_sample_sphere_tessellation(
    polygons_struct   *polygons,
    polygons_struct   *half );

BICAPI   void  initialize_intersect_directions( void );

BICAPI  VIO_Real  *get_intersect_directions( void );

BICAPI  VIO_BOOL  point_within_polygon(
    VIO_Point   *pt,
    int     n_points,
    VIO_Point   points[],
    VIO_Vector  *polygon_normal );

BICAPI  VIO_BOOL  line_intersects_ellipsoid(
    VIO_Point    *line_origin,
    VIO_Vector   *line_direction,
    VIO_Point    *sphere_centre,
    VIO_Real     x_size,
    VIO_Real     y_size,
    VIO_Real     z_size,
    VIO_Real     *t_min,
    VIO_Real     *t_max );

BICAPI  VIO_BOOL  ray_intersects_sphere(
    VIO_Point       *origin,
    VIO_Vector      *direction,
    VIO_Point       *centre,
    VIO_Real        radius,
    VIO_Real        *dist );

BICAPI  void  intersect_ray_object(
    VIO_Point                 *origin,
    VIO_Vector                *direction,
    object_struct         *object,
    int                   obj_index,
    int                   *closest_obj_index,
    VIO_Real                  *closest_dist,
    int                   *n_intersections,
    VIO_Real                  *distances[] );

BICAPI  int  intersect_ray_with_object(
    VIO_Point           *origin,
    VIO_Vector          *direction,
    object_struct   *object,
    int             *obj_index,
    VIO_Real            *dist,
    VIO_Real            *distances[] );

BICAPI  void   remove_invisible_polygons(
    polygons_struct  *polygons,
    VIO_SCHAR     visibilities[] );

BICAPI  VIO_Real  get_smooth_surface_curvature(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    int               poly,
    int               vertex,
    VIO_BOOL           distances_initialized,
    float             distances[],
    VIO_Real              smoothing_distance );

BICAPI  void  smooth_lines(
    lines_struct  *lines,
    VIO_Real          smooth_length );

BICAPI  void  create_line_spline(
    lines_struct  *lines,
    int           n_curve_segments,
    lines_struct  *new_lines );

BICAPI  void  smooth_polygon(
    polygons_struct  *polygons,
    VIO_Real             max_dist_from_original,
    VIO_Real             fraction_to_move,
    VIO_Real             stop_threshold,
    VIO_Real             normal_ratio,
    VIO_BOOL          range_flag,
    volume_struct    *volume,
    int              min_value,
    int              max_value );

BICAPI  VIO_BOOL  get_interpolation_weights_2d(
    VIO_Real   x,
    VIO_Real   y,
    int    n_points,
    VIO_Real   xs[],
    VIO_Real   ys[],
    VIO_Real   weights[] );

BICAPI  VIO_BOOL  get_prediction_weights_2d(
    VIO_Real   x,
    VIO_Real   y,
    int    n_points,
    VIO_Real   xs[],
    VIO_Real   ys[],
    VIO_Real   *x_weights[2],
    VIO_Real   *x_constant,
    VIO_Real   *y_weights[2],
    VIO_Real   *y_constant );

BICAPI  VIO_BOOL  get_prediction_weights_3d(
    VIO_Real   x,
    VIO_Real   y,
    VIO_Real   z,
    int    n_points,
    VIO_Real   xs[],
    VIO_Real   ys[],
    VIO_Real   zs[],
    VIO_Real   *x_weights[3],
    VIO_Real   *y_weights[3],
    VIO_Real   *z_weights[3] );

BICAPI  VIO_BOOL is_single_closed_curve(
    lines_struct   *lines );

BICAPI  void  subdivide_lines(
    lines_struct  *lines );

BICAPI  void  subdivide_polygons(
    polygons_struct  *polygons );

BICAPI  void  subdivide_polygons_indices(
  polygons_struct  *polygons,
  int              *data_indices[]);

BICAPI  VIO_Real  get_polygon_2d_area(
    int      n_points,
    VIO_Point    points[] );

BICAPI  VIO_Real  get_polygon_surface_area(
    int      n_points,
    VIO_Point    points[] );

BICAPI  VIO_Real  get_polygons_surface_area(
    polygons_struct  *polygons );

BICAPI  VIO_BOOL  is_this_tetrahedral_topology(
    polygons_struct   *polygons );

BICAPI  int  get_tetra_tessellation_with_n_points(
    int   n_points );

BICAPI  void  create_tetrahedral_sphere(
    VIO_Point            *centre,
    VIO_Real             rx,
    VIO_Real             ry,
    VIO_Real             rz,
    int              n_triangles,
    polygons_struct  *polygons );

BICAPI  void  half_sample_tetrahedral_tessellation(
    polygons_struct  *polygons,
    polygons_struct  *half );

BICAPI  int  convert_lines_to_tubes(
    lines_struct     *lines,
    int              n_around,
    VIO_Real             radius,
    quadmesh_struct  *quadmeshes[] );

BICAPI  void   create_slice_quadmesh(
    VIO_Volume           volume,
    int              axis_index,
    VIO_Real             voxel_position,
    int              x_tess,
    int              y_tess,
    VIO_Real             x_min,
    VIO_Real             x_max,
    VIO_Real             y_min,
    VIO_Real             y_max,
    quadmesh_struct  *quadmesh );

BICAPI  void   create_slice_3d(
    VIO_Volume           volume,
    VIO_Point            *origin,
    VIO_Vector           *normal,
    polygons_struct  *polygons );

#ifdef __cplusplus
}
#endif

#endif
