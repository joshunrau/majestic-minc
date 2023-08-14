#ifndef  DEF_RAY_TRACE_PROTOTYPES
#define  DEF_RAY_TRACE_PROTOTYPES

 VIO_BOOL  ray_intersects_a_polygon(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    object_struct    *object,
    VIO_BOOL          flat_shading_flag,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Real             *dist );

 VIO_BOOL  ray_intersects_a_quadmesh(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    object_struct    *object,
    VIO_BOOL          flat_shading_flag,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Real             *dist );

 VIO_BOOL  ray_intersects_a_line(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    object_struct    *object,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Real             *dist );

 VIO_BOOL  ray_intersects_a_marker(
    VIO_Point            *origin,
    VIO_Vector           *direction,
    object_struct    *object,
    VIO_Point            *point,
    VIO_Vector           *normal,
    VIO_Colour           *colour,
    VIO_Real             *dist );

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
    int                   y_super_sampling );

 void  ray_trace_scene(
    view_struct           *view,
    lights_struct         *lights,
    int                   n_objects,
    ray_trace_object      objects[],
    int                   max_depth,
    int                   x_super_sampling,
    int                   y_super_sampling,
    pixels_struct         *pixels );

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
    VIO_Real               *dist );

 void  get_range_of_volume(
    VIO_Volume    volume,
    VIO_Point     *min_point,
    VIO_Point     *max_point );
#endif
