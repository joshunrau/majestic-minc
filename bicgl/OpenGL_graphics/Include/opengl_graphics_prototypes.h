#ifndef  DEF_opengl_graphics_prototypes
#define  DEF_opengl_graphics_prototypes

  VIO_Colour  make_rgba_Colour(
    int    r,
    int    g,
    int    b,
    int    a );

  int  get_Colour_r(
    VIO_Colour   colour );

  int  get_Colour_g(
    VIO_Colour   colour );

  int  get_Colour_b(
    VIO_Colour   colour );

  int  get_Colour_a(
    VIO_Colour   colour );

  void  GS_set_point(
    VIO_Point  *point );

  void  GS_set_normal(
    VIO_Vector  *normal );

  void  GS_set_colour(
    VIO_Colour  colour );

  void  GS_set_colour_index(
    VIO_Colour  colour );

  void  GS_set_ambient_and_diffuse_mode(
    VIO_BOOL  state );

  void  GS_initialize_surface_property(
    GSwindow        window );

  void  GS_set_surface_property(
    GSwindow       window,
    VIO_Colour         col,
    VIO_Surfprop       *surfprop );

  void  GS_set_line_width(
    VIO_Real  line_width );

  void  GS_curve(
    float  geom[4][3] );

  void  GS_begin_point( void );

  void  GS_end_point( void );

  void  GS_begin_line( void );

  void  GS_end_line( void );

  void  GS_begin_closed_line( void );

  void  GS_end_closed_line( void );

  void  GS_begin_polygon( void );

  void  GS_end_polygon( void );

  void  GS_begin_quad_strip( void );

  void  GS_end_quad_strip( void );

  void  GS_set_raster_position(
    VIO_Real  x,
    VIO_Real  y,
    VIO_Real  z );

  void  GS_set_pixel_zoom(
    VIO_Real  x_zoom,
    VIO_Real  y_zoom );

  void  GS_draw_colour_map_pixels(
    int             x_viewport_min,
    int             y_viewport_min,
    pixels_struct   *pixels );

  void  GS_draw_rgb_pixels(
    int             x_viewport_min,
    int             y_viewport_min,
    pixels_struct   *pixels );

  void  GS_read_pixels(
    VIO_BOOL         colour_map_state,
    int             x_min,
    int             x_max,
    int             y_min,
    int             y_max,
    VIO_Colour          pixels[] );

  VIO_Real  GS_get_character_height(
    Font_types       type,
    VIO_Real             size );

  void  GS_draw_text(
    Font_types   type,
    VIO_Real         size,
    VIO_STR       string );

  VIO_Real  GS_get_text_length(
    VIO_STR           str,
    Font_types       type,
    VIO_Real             size );

  void  GS_set_update_function(
    void  (*func)( Window_id ) );

  void  GS_set_update_overlay_function(
    void  (*func)( Window_id ) );

  void  GS_set_resize_function(
    void  (*func)( Window_id, int, int, int, int ) );

  void  GS_set_key_down_function(
    void  (*func)( Window_id, int, int, int, int ) );

  void  GS_set_key_up_function(
    void  (*func)( Window_id, int, int, int, int ) );

  void  GS_set_mouse_movement_function(
    void  (*func)( Window_id, int, int ) );

  void  GS_set_left_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) );

  void  GS_set_left_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) );

  void  GS_set_middle_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) );

  void  GS_set_middle_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) );

  void  GS_set_right_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) );

  void  GS_set_right_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) );

  void  GS_set_scroll_down_function(
    void  (*func)( Window_id, int, int, int ) );

  void  GS_set_scroll_up_function(
    void  (*func)( Window_id, int, int, int ) );

  void  GS_set_iconify_function(
    void  (*func)( Window_id ) );

  void  GS_set_deiconify_function(
    void  (*func)( Window_id ) );

  void  GS_set_enter_function(
    void  (*func)( Window_id ) );

  void  GS_set_leave_function(
    void  (*func)( Window_id ) );

  void  GS_set_quit_function(
    void  (*func)( Window_id ) );

  void  GS_set_update_flag(
    GSwindow   window  );

  void  GS_add_timer_function(
    VIO_Real          seconds,
    void          (*func) ( void * ),
    void          *data );

  void  GS_add_idle_function(
    void  (*func) ( void * ),
    void          *data );

  void  GS_remove_idle_function(
    void  (*func) ( void * ),
    void          *data );

  void  GS_event_loop( void );

void GS_exit_loop(void);

  void  GS_set_mouse_position(
    int       x_screen,
    int       y_screen );

  VIO_BOOL  GS_are_mouse_coordinates_in_screen_space( void );

  void  GS_initialize_lights(
    GSwindow        window );

  void  delete_lights(
    GSwindow        window );

  void  GS_set_ambient_light(
    GSwindow       window,
    VIO_Colour         colour );

  void  GS_define_light(
    GSwindow        window,
    int             light_index,
    Light_types     type,
    VIO_Colour          colour,
    VIO_Vector          *direction,
    VIO_Point           *position,
    VIO_Real            spot_exponent,
    VIO_Real            spot_angle );

  void  GS_set_light_state(
    GSwindow        window,
    int             light_index,
    VIO_BOOL         state );

  void  redefine_lights(
    GSwindow        window );

  void  GS_set_shade_model(
    Shading_types  type );

  void  GS_set_blend_function( Transparency_modes );

  void  GS_set_lighting_state(
    GSwindow     window,
    VIO_BOOL      state );

  void  GS_backface_culling_state(
    VIO_BOOL     state );

  void  GS_set_n_curve_segments(
    int      n_segments );

  void  GS_load_transform(
    VIO_Transform  *transform );

  void  GS_mult_transform(
    VIO_Transform  *transform );

  void  GS_get_transform(
    VIO_Transform  *trans );

  void  GS_ortho(
    VIO_Real   x_min,
    VIO_Real   x_max,
    VIO_Real   y_min,
    VIO_Real   y_max,
    VIO_Real   z_min,
    VIO_Real   z_max );

  void  GS_ortho_2d(
    int   x_min,
    int   x_max,
    int   y_min,
    int   y_max );

  void  GS_frustum(
    VIO_Real   x_min,
    VIO_Real   x_max,
    VIO_Real   y_min,
    VIO_Real   y_max,
    VIO_Real   z_min,
    VIO_Real   z_max );

  void  GS_initialize_window_view(
    GSwindow   window );

  void  GS_push_transform( void );

  void  GS_pop_transform( void );

  void  GS_set_viewport(
    int            x_min,
    int            x_max,
    int            y_min,
    int            y_max );

  void  clear_overlay_planes( void );

  void  GS_set_matrix_mode(
    Matrix_modes   mode );

  Window_id  GS_get_window_id(
    GSwindow  window );

  void  GS_set_current_window(
    GSwindow        window );

  void  GS_initialize( void );

  void  GS_terminate( void );

  VIO_BOOL  GS_can_switch_double_buffering( void );

  VIO_BOOL  GS_can_switch_colour_map_mode( void );

  VIO_Status  GS_create_window(
    GSwindow       window,
    VIO_STR         title,
    int            x_pos,
    int            y_pos,
    int            width,
    int            height,
    VIO_BOOL        colour_map_desired,
    VIO_BOOL        double_buffer_desired,
    VIO_BOOL        depth_buffer_desired,
    int            n_overlay_planes_desired,
    VIO_BOOL        *actual_colour_map_flag,
    VIO_BOOL        *actual_double_buffer_flag,
    VIO_BOOL        *actual_depth_buffer_flag,
    int            *actual_n_overlay_planes,
    GSwindow       parent);

  void  print_info( void );

  void  GS_set_window_title(
    GSwindow   window,
    VIO_STR     title );

  VIO_BOOL  GS_has_transparency_mode( void );

  VIO_BOOL  GS_has_rgb_mode( void );

  VIO_BOOL  GS_set_double_buffer_state(
    GSwindow       window,
    VIO_BOOL        flag );

  VIO_BOOL  GS_set_colour_map_state(
    GSwindow       window,
    VIO_BOOL        flag );

  int  GS_get_n_colour_map_entries(
    GSwindow  window,
    VIO_BOOL   double_buffer_state );

  void  GS_set_colour_map_entry(
    GSwindow        window,
    Bitplane_types  bitplanes,
    int             ind,
    VIO_Colour          colour );

  VIO_BOOL  GS_is_double_buffer_supported( void );

  VIO_BOOL  GS_is_depth_buffer_supported( void );

  void  GS_set_depth_function(
    Depth_functions  func );

  void  GS_set_depth_buffer_state(
    VIO_BOOL         flag );

  VIO_Status  GS_delete_window(
    GSwindow   window );

  int  GS_get_monitor_width( void );

  int  GS_get_monitor_height( void );

  void  GS_get_window_size(
    int        *x_size,
    int        *y_size );

  void  GS_get_window_position(
    int        *x_position,
    int        *y_position );

  void  GS_clear_depth_buffer( void );

  void  GS_clear_overlay( void );

  void  GS_clear_window(
    GSwindow        window,
    int             x_size,
    int             y_size,
    Bitplane_types  current_bitplanes,
    VIO_BOOL         colour_map_state,
    VIO_BOOL         zbuffer_state,
    VIO_Colour          colour );

  void  GS_clear_viewport(
    GSwindow        window,
    int             x_viewport_min,
    int             x_viewport_max,
    int             y_viewport_min,
    int             y_viewport_max,
    int             x_size,
    int             y_size,
    Bitplane_types  current_bitplanes,
    VIO_BOOL         colour_map_state,
    VIO_BOOL         zbuffer_state,
    VIO_Colour          colour );

  void  GS_flush( void );

  void  GS_append_to_last_update(
    GSwindow   window,
    VIO_BOOL    zbuffer_state,
    int        x_size,
    int        y_size );

  int  GS_get_num_overlay_planes( void );

  void  GS_set_bitplanes( 
    GSwindow         window,
    Bitplane_types   bitplanes );

  void  GS_set_overlay_colour_map(
    GSwindow         window,
    int              ind,
    VIO_Colour           colour );

  void  GS_swap_buffers( void );

  void  GS_set_visibility(GSwindow window, VIO_BOOL is_visible);

void GS_set_geometry(GSwindow window, int x, int y, int width, int height);

void GS_draw_polygons_one_colour(GSwindow window, Shading_modes shading_mode,
                                 polygons_struct *polygons);

void GS_draw_polygons_per_vertex_colours(GSwindow window,
                                         Shading_modes shading_mode,
                                         polygons_struct *polygons);

void GS_draw_triangles_one_colour(GSwindow window, Shading_modes shading_mode,
                                  polygons_struct *polygons);

void GS_draw_triangles_per_item_colours(GSwindow window,
                                        Shading_modes shading_mode,
                                        polygons_struct *polygons);

void GS_draw_quadmesh_one_colour(GSwindow window, Shading_modes shading_mode,
                                 quadmesh_struct *quadmesh);

void GS_draw_quadmesh_per_item_colours(GSwindow         window,
                                       Shading_modes    shading_mode,
                                       quadmesh_struct *quadmesh );

void GS_draw_quadmesh_per_vertex_colours(GSwindow window, 
                                         Shading_modes shading_mode,
                                         quadmesh_struct *quadmesh);

void GS_draw_lines(GSwindow window, lines_struct *lines);
#endif
