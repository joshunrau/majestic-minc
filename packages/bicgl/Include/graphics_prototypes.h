#ifndef DEF_graphics_prototypes
#define DEF_graphics_prototypes
 void initialize_display_interrupts(
 Gwindow window );

/**
 * \defgroup Drawing Drawing functions.
 * @{
 */
void G_set_drawing_interrupt_state( Gwindow window, VIO_BOOL state );

VIO_BOOL G_get_drawing_interrupt_state( Gwindow window );

void G_set_drawing_interrupt_check_n_objects( Gwindow window, int interval );

void G_set_drawing_interrupt_time( Gwindow window, VIO_Real interval );

void G_start_interrupt_test( Gwindow window );

VIO_BOOL G_get_interrupt_occurred( Gwindow window );

void G_clear_drawing_interrupt_flag( Gwindow window );

void set_continuation_flag( Gwindow window, VIO_BOOL state );

VIO_BOOL view_is_stereo( Gwindow window );

VIO_BOOL is_triangular( polygons_struct *polygons );
/**
 * Draw the polygons_struct in the selected window.
 *
 * The polygons_struct defines the vertices, colours, surface properties, 
 * and faces of a polygon in 3D. In practice, polygons_structs are generally
 * triangularizations of complex surfaces.
 *
 * This function will draw the object using the current lighting,
 * view mode, and transformation matrices.
 *
 * \param window The window's handle.
 * \param polygons A pointer to the polygons_struct to draw.
 */
void G_draw_polygons( Gwindow window, polygons_struct *polygons );

/**
 * Draw a quadmesh_struct in the selected window.
 *
 * The quadmesh_struct defines the vertices, colours, surface properties, 
 * and faces of a quadmesh in 3D.
 *
 * This function will draw the object using the current lighting,
 * view mode, and transformation matrices.
 *
 * \param window The window's handle.
 * \param quadmesh A pointer to the quadmesh_struct to draw.
 */
void G_draw_quadmesh( Gwindow window, quadmesh_struct *quadmesh );

/**
 * Draw a lines_struct in the selected window.
 *
 * The lines_struct defines the vertices, colours, and segments
 * of a set of lines in 3D.
 *
 * This function will draw the object using the current lighting,
 * view mode, and transformation matrices.
 *
 * \param window The window's handle.
 * \param lines A pointer to the lines_struct to draw.
 */
void G_draw_lines( Gwindow window, lines_struct *lines );

/**
 * Calculate the height in pixels of text using the given font parameters.
 * \param type The font type.
 * \param size The font size.
 * \returns The height of the selected font.
 */
VIO_Real G_get_text_height( Font_types type, VIO_Real size );

/**
 * Calculate the width in pixels of text using the given string and 
 * font parameters.
 * \param str The string whose width we want to calculate.
 * \param type The font type.
 * \param size The font size.
 * \returns The width of the string when drawn in the selected font.
 */
VIO_Real G_get_text_length( VIO_STR str, Font_types type, VIO_Real size );

/**
 * Draw a text object in the given window.
 *
 * The text_struct defines the text, colours, and position of text in
 * 3D. 
 *
 * This function will draw the object using the current lighting,
 * view mode, and transformation matrices.
 *
 * \param window The window's handle.
 * \param text A pointer to the text_struct of the object.
 *
 */
void G_draw_text( Gwindow window, text_struct *text );

void G_draw_marker( Gwindow window, marker_struct *marker );

void G_draw_pixels( Gwindow window, pixels_struct *pixels );

void G_read_pixels( Gwindow window, int x_min, int x_max, int y_min, int y_max,
                    VIO_Colour pixels[] );

/**
 * Draw a graphics object in a window.
 * 
 * Based on the object's type, will call one of G_draw_polygons(), 
 * G_draw_quadmesh(), G_draw_lines(), etc.
 *
 * \param window The window's handle.
 * \param object The object to draw.
 */
void draw_object( Gwindow window, object_struct *object );

/**
 * @}
 */

/**
 * \defgroup system System functions.
 * @{
 */
VIO_BOOL G_get_shift_key_state( void );

VIO_BOOL G_get_ctrl_key_state( void );

VIO_BOOL G_get_alt_key_state( void );
/**
 * @}
 */

/**
 * \defgroup event Event handling functions.
 * @{
 */
/** 
 * The type signature of callback functions for generic window events.
 * \param window The handle of the window.
 * \param data Opaque data handle associated with this callback.
 */
typedef void (*G_generic_callback_t)(Gwindow window, void *data);

void G_set_update_function( Gwindow window, G_generic_callback_t func,
                            void *func_data );

void G_set_update_overlay_function( Gwindow window, G_generic_callback_t func,
                                    void *func_data );

/** 
 * The type signature of callback functions for window resize events.
 * \param window The handle of the window that was resized.
 * \param x The new x position.
 * \param y The new y position.
 * \param width The new window width.
 * \param height The new window height.
 * \param data Opaque data handle associated with this callback.
 */
typedef void (*G_resize_callback_t)( Gwindow window, int x, int y,
                                     int width, int height, void *data );

void G_set_resize_function( Gwindow window, G_resize_callback_t func,
                            void *func_data );

/** 
 * The type signature of callback functions for keyboard events.
 * \param window The handle of the window.
 * \param key The key code associated with the event.
 * \param data Opaque data handle associated with this callback.
 */
typedef void (*G_key_callback_t)( Gwindow window, int key, void *data);

void G_set_key_down_function( Gwindow window, G_key_callback_t func,
                              void *func_data );

void G_set_key_up_function( Gwindow window, G_key_callback_t func,
                            void *func_data );

/** 
 * The type signature of callback functions for mouse events.
 * \param window The handle of the window.
 * \param x The x position of the mouse event within the window.
 * \param y The y position of the mouse event within the window.
 * \param data Opaque data handle associated with this callback.
 */
typedef void (*G_mouse_callback_t)( Gwindow, int, int, void *);

void G_set_mouse_movement_function( Gwindow window, G_mouse_callback_t func,
                                    void *func_data );

void G_set_left_mouse_down_function( Gwindow window, G_mouse_callback_t func,
                                     void *func_data );

void G_set_left_mouse_up_function( Gwindow window, G_mouse_callback_t func,
                                   void *func_data );

void G_set_middle_mouse_down_function( Gwindow window, G_mouse_callback_t func,
                                       void *func_data );

void G_set_middle_mouse_up_function( Gwindow window, G_mouse_callback_t func,
                                     void *func_data );

void G_set_right_mouse_down_function( Gwindow window, G_mouse_callback_t func,
                                      void *func_data );

void G_set_right_mouse_up_function( Gwindow window, G_mouse_callback_t func,
                                    void *func_data );

void G_set_scroll_down_function( Gwindow window, G_mouse_callback_t func,
                                 void *func_data );

void G_set_scroll_up_function( Gwindow window, G_mouse_callback_t func,
                               void *func_data );

void G_set_iconify_function( Gwindow window, G_generic_callback_t func,
                             void *func_data );

void G_set_deiconify_function( Gwindow window, G_generic_callback_t func,
                               void *func_data );

void G_set_window_enter_function( Gwindow window, G_generic_callback_t func,
                                  void *func_data );

void G_set_window_leave_function( Gwindow window, G_generic_callback_t func,
                                  void *func_data );

void G_set_window_quit_function( Gwindow window, G_generic_callback_t func,
                                 void *func_data );

void initialize_callbacks_for_window( Gwindow window );

void G_set_update_flag( Gwindow window );

void G_add_timer_function( VIO_Real seconds, void (*func) ( void * ),
                           void *data );

void G_add_idle_function( void (*func) ( void * ), void *data );

void G_remove_idle_function( void (*func) ( void * ), void *data );

void G_main_loop( void );

void G_exit_loop(void);
/**
 * @}
 */

/**
 * \defgroup mouse Mouse functions
 * @{
 */
void G_set_mouse_position( int x_screen, int y_screen );

VIO_BOOL G_is_mouse_in_window( Gwindow window );

void G_get_mouse_screen_position( int *x_screen_pos, int *y_screen_pos );

VIO_BOOL G_get_mouse_position( Gwindow window, int *x_pixel_pos, 
                               int *y_pixel_pos );

VIO_BOOL G_get_mouse_position_0_to_1( Gwindow window,
                                      VIO_Real *x_pos, VIO_Real *y_pos );
VIO_BOOL G_get_left_mouse_button( void );

VIO_BOOL G_get_middle_mouse_button( void );

VIO_BOOL G_get_right_mouse_button( void );

/**
 * @}
 */
void G_set_default_update_min_interval( VIO_Real seconds );

void G_set_window_update_min_interval( Gwindow window, VIO_Real seconds );

void initialize_graphics_struct(
 graphics_struct *graphics );

void delete_graphics_struct(
 graphics_struct *graphics );

void set_graphics_viewport(
 graphics_struct *graphics,
 int viewport_index,
 int x_min,
 int x_max,
 int y_min,
 int y_max );

void get_graphics_viewport(
 graphics_struct *graphics,
 int viewport_index,
 int *x_min,
 int *x_max,
 int *y_min,
 int *y_max );

void set_graphics_viewport_background(
 graphics_struct *graphics,
 int viewport_index,
 VIO_Colour background_colour,
 int background_colour_index );

void add_object_to_viewport(
 graphics_struct *graphics,
 int viewport_index,
 Bitplane_types bitplane,
 object_struct *object );

void remove_object_from_viewport(
 graphics_struct *graphics,
 int viewport_index,
 Bitplane_types bitplane,
 object_struct *object );

void set_viewport_update_flag( 
 graphics_struct *graphics,
 int viewport_index,
 Bitplane_types bitplane );

VIO_BOOL get_viewport_update_flag( 
 graphics_struct *graphics,
 int viewport_index,
 Bitplane_types bitplane,
 int buffer );

void set_bitplanes_clear_flag( 
 graphics_struct *graphics,
 Bitplane_types bitplane );

void set_viewport_bitplane_objects_visibility(
 graphics_struct *graphics,
 int viewport_index,
 Bitplane_types bitplane,
VIO_BOOL visibility );

void set_viewport_objects_visibility(
 graphics_struct *graphics,
 int viewport_index,
VIO_BOOL visibility );

VIO_BOOL redraw_out_of_date_viewports(
 graphics_struct *graphics,
 Gwindow window,
 int current_buffer );

void initialize_graphics_lights(
 Gwindow window );

void G_set_ambient_light(
 Gwindow window,
 VIO_Colour colour );

void G_define_light(
 Gwindow window,
 int light_index,
 Light_types type,
 VIO_Colour colour,
 VIO_Vector *direction,
 VIO_Point *position,
 VIO_Real spot_exponent,
 VIO_Real spot_angle );

void G_set_light_state(
 Gwindow window,
 int light_index,
VIO_BOOL state );

 Random_mask_type get_random_order_mask(
 int size );

void G_set_shaded_state( Gwindow window, Shading_modes state );

void G_set_shading_type( Gwindow window, Shading_types type );

VIO_BOOL G_get_lighting_state(
 Gwindow window );

void G_set_lighting_state(
 Gwindow window,
VIO_BOOL state );

void G_backface_culling_state(
 Gwindow window,
VIO_BOOL state );

void G_set_render_lines_as_curves_state(
 Gwindow window,
VIO_BOOL state );

void G_set_n_curve_segments(
 Gwindow window,
 int n_segments );

void G_set_markers_labels_visibility(
 Gwindow window,
VIO_BOOL state );

void initialize_window_view(
 Gwindow window );

void set_view_for_eye(
 Gwindow window,
 int which_eye );

void G_set_view_type(
 Gwindow window,
 View_types view_type );

void G_set_3D_view(
 Gwindow window,
 VIO_Point *origin,
 VIO_Vector *line_of_sight,
 VIO_Vector *up_direction,
 VIO_Real front_clip_distance,
 VIO_Real back_clip_distance,
VIO_BOOL perspective_flag,
 VIO_Real perspective_distance,
VIO_BOOL stereo_flag,
 VIO_Real eye_separation,
 VIO_Real window_width,
 VIO_Real window_height );

void G_set_modeling_transform(
 Gwindow window,
 VIO_Transform *transform );

void G_transform_point(
 Gwindow window,
 VIO_Point *point,
 View_types view_type,
 int *x_pixel,
 int *y_pixel );

void update_transforms(
 Gwindow window );

void G_push_transform(
 Gwindow window,
 VIO_Transform *transform );

void G_pop_transform(
 Gwindow window );

void window_was_resized(
 Gwindow window );

void G_set_viewport(
 Gwindow window,
 int x_min,
 int x_max,
 int y_min,
 int y_max );

int get_n_graphics_windows( void );

Gwindow get_nth_graphics_window( int i );

void set_current_window( Gwindow window );

Gwindow get_current_window( void );

/**
 * \defgroup Windows Windowing functions.
 * @{
 */

/**
 * Create a new top-level window.
 * \param title A string to use as the window's title.
 * \param x_pos The initial column position of the window in pixels.
 * \param y_pos The initial row position of the window in pixels.
 * \param width The initial width of the window in pixels.
 * \param height The initial height of the window in pixels.
 * \param colour_map_desired TRUE if this window should use colour map mode.
 * If false, the window will be created in RGB mode. Note that colour map
 * mode is now nearly obsolete.
 * \param double_buffer_desired TRUE if this window should use double 
 * buffering. This is also nearly obsolete.
 * \param depth_buffer_desired TRUE if this window should implement a z-buffer
 * for 3D object placement. This is still useful.
 * \param n_overlay_planes_desired The number of overlay planes to use. Not obsolete.
 * \param window A pointer to the Gwindow handle to be filled with the new
 * window's information.
 * \returns VIO_OK if the operation is successful.
 */
VIO_Status G_create_window( VIO_STR title, int x_pos, int y_pos,
                            int width, int height,
                            VIO_BOOL colour_map_desired,
                            VIO_BOOL double_buffer_desired,
                            VIO_BOOL depth_buffer_desired,
                            int n_overlay_planes_desired,
                            Gwindow *window );

/**
 * Create a new window, which _may_be a child of another window.
 * \param title A string to use as the window's title.
 * \param x_pos The initial column position of the window in pixels.
 * \param y_pos The initial row position of the window in pixels.
 * \param width The initial width of the window in pixels.
 * \param height The initial height of the window in pixels.
 * \param colour_map_desired TRUE if this window should use colour map mode.
 * If false, the window will be created in RGB mode. Note that colour map
 * mode is now nearly obsolete.
 * \param double_buffer_desired TRUE if this window should use double 
 * buffering. This is also nearly obsolete.
 * \param depth_buffer_desired TRUE if this window should implement a z-buffer
 * for 3D object placement. This is still useful.
 * \param n_overlay_planes_desired The number of overlay planes to use. Not obsolete.
 * \param window A pointer to the Gwindow handle to be filled with the new
 * window's information.
 * \param parent A pointer to the parent window for this window. If NULL, the
 * window will be created as a top-level window.
 * \returns VIO_OK if the operation is successful.
 */
VIO_Status G_create_child_window( VIO_STR title, int x_pos, int y_pos,
                                  int width, int height,
                                  VIO_BOOL colour_map_desired,
                                  VIO_BOOL double_buffer_desired,
                                  VIO_BOOL depth_buffer_desired,
                                  int n_overlay_planes_desired,
                                  Gwindow *window, Gwindow parent);

/**
 * Delete the given window.
 * \param window The handle of the window to delete.
 * \returns VIO_OK if the operation was successful.
 */
VIO_Status G_delete_window( Gwindow window );

/**
 * Set the title of the window.
 * \param window The handle of the window.
 * \param title The new title to use.
 */
void G_set_window_title( Gwindow window, VIO_STR title );

/**
 * Delete all windows and exit.
 */
void G_terminate( void );

/**
 * Check whether double buffering is supported by the graphics system.
 * \returns TRUE if double buffering is possible.
 */
VIO_BOOL G_is_double_buffer_supported( void );

/**
 * Returns the current double buffering state of the window.
 * \param window The window's handle.
 * \returns TRUE if this window has double buffering enabled.
 */
VIO_BOOL G_get_double_buffer_state( Gwindow window );

/**
 * Attempts to set the current double buffering state of the window.
 * \param window The window's handle.
 * \param The desired double buffering state.
 */
void G_set_double_buffer_state( Gwindow window, VIO_BOOL flag );

/**
 * Check whether depth buffering is supported by the graphics system.
 * \returns TRUE if depth buffering is possible.
 */
VIO_BOOL G_is_depth_buffer_supported( void );

/**
 * Returns the current depth buffering state of the window.
 * \param window The window's handle.
 * \returns TRUE if this window has depth buffering enabled.
 */
VIO_BOOL G_get_zbuffer_state( Gwindow window );

/**
 * Attempts to set the current depth buffering state of the window.
 * \param window The window's handle.
 * \param The desired depth buffering state.
 */
void G_set_zbuffer_state( Gwindow window, VIO_BOOL flag );

/**
 * Returns the current colour map state of the window.
 * \param window The window's handle.
 * \returns TRUE if this window has colour map mode enabled, FALSE if RGB.
 */
VIO_BOOL G_get_colour_map_state( Gwindow window );

/**
 * Attempts to set the current colour map mode of the window.
 * \param window The window's handle.
 * \param TRUE if we want colour map mode, FALSE for RGB(A).
 */
void G_set_colour_map_state( Gwindow window, VIO_BOOL flag );

/**
 * Get the number of colour map entries for this window.
 * \param window The window's handle.
 * \returns The number of colour map entries.
 */
int G_get_n_colour_map_entries( Gwindow window );

/**
 * Set a colour map entry for this window.
 * \param window The window's handle.
 * \param ind The colour index to set.
 * \param colour The RGBA value to associate with this colour index.
 */
void G_set_colour_map_entry( Gwindow window, int ind, VIO_Colour colour );

void restore_bitplanes( Gwindow window, Bitplane_types bitplane );

/** 
 * Get the width of the current physical display, in pixels.
 * \returns The display width.
 */
int G_get_monitor_width( void );

/** 
 * Get the height of the current physical display, in pixels.
 * \returns The display height.
 */
int G_get_monitor_height( void );

VIO_Real G_get_monitor_widths_to_eye( void );

/**
 * Get the size in pixels of the window.
 * \param window The window's handle.
 * \param x_size A pointer to the variable which will hold the width.
 * \param y_size A pointer to the variable which will hold the height.
 */
void G_get_window_size( Gwindow window, int *x_size, int *y_size );

/**
 * Get the location in pixels of the window.
 * \param window The window's handle.
 * \param x_pos A pointer to the variable which will hold the column position.
 * \param y_pos A pointer to the variable which will hold the row position.
 */
void G_get_window_position( Gwindow window, int *x_pos, int *y_pos );

/**
 * Set the background colour for the window.
 * \param window The window's handle.
 * \param colour The new background colour to use.
 */
void G_set_background_colour( Gwindow window, VIO_Colour colour );

/**
 * Get the background colour of the window.
 * \param window The window's handle
 * \returns The current background colour.
 */
VIO_Colour G_get_background_colour( Gwindow window );

/**
 * Clear the window to the background colour.
 */
void G_clear_window( Gwindow window );

void G_clear_viewport( Gwindow window, VIO_Colour colour );

void G_set_automatic_clear_state( Gwindow window, VIO_BOOL state );

/**
 * Called to indicate that a window redraw is completed - in some cases,
 * this will swap the foreground and background buffers.
 * \param window The window's handle.
 */
void G_update_window( Gwindow window );

void check_window_cleared( Gwindow window );

/**
 * Check whether the graphics system supports overlays.
 * \returns TRUE if overlay planes are supported.
 */
VIO_BOOL G_has_overlay_planes( void );

/**
 * Check whether this window has overlay planes enabled.
 * \param window The window's handle.
 * \returns TRUE if the window has an overlay plane enabled.
 */
VIO_BOOL G_window_has_overlay_planes( Gwindow window );

/**
 * Set the current window. Some operations implicitly operate on 
 * the current window, so we need to be able to force a particular
 * choice for the current window.
 * \param The window's handle.
 */
void G_set_current_window(Gwindow window);

void G_set_bitplanes( Gwindow window, Bitplane_types bitplanes );

Bitplane_types G_get_bitplanes( Gwindow window );

/**
 * Check whether the graphics system allows the double buffering mode
 * to be changed. 
 */
VIO_BOOL G_can_switch_double_buffering( void );

/**
 * Check whether the graphics system allows the colour mapping mode
 * to be changed. 
 */
VIO_BOOL G_can_switch_colour_map_mode( void );

/**
 * Check whether the graphics system supports transparency mode.
 */
VIO_BOOL G_has_transparency_mode( void );

/** 
 * Check whether the graphics system supports RGB(A) mode.
 */
VIO_BOOL G_has_rgb_mode( void );

/**
 * Get the number of overlay planes allowed by the graphics system.
 */
int G_get_n_overlay_planes( void );

/**
 * Set an entry in the colour map of the overlay.
 * \param window The window's handle.
 * \param ind The colour map index to change.
 * \param colour The new colour value to associate with this index.
 */
void G_set_overlay_colour_map( Gwindow window, int ind, VIO_Colour colour );

void G_append_to_last_update( Gwindow window );

void G_continue_last_update( Gwindow window );

/**
 * Set the transparency state of the window.
 * \param window The window's handle.
 * \param state 0 to disable, 1 for normal, 2 for additive.
 */
void G_set_transparency_state( Gwindow window, Transparency_modes state );

/**
 * Make the window visible or invisible.
 * \param window The window's handle.
 * \param is_visible If TRUE, the window will be made visible. If FALSE, it
 * will be made invisible.
 */
void G_set_visibility( Gwindow window, VIO_BOOL is_visible );

/**
 * Set the position and size of the window. If a negative value is
 * passed for a position or size, that aspect of the geometry will
 * not be changed.
 * \param window The window's handle.
 * \param x The new x position.
 * \param y The new y position.
 * \param width The new window width.
 * \param height The new window height.
 */
void G_set_geometry(Gwindow window, int x, int y, int width, int height);
/** 
 * @}
 */
#endif
