#ifndef  DEF_glut_window_prototypes
#define  DEF_glut_window_prototypes

  void  WS_initialize( void );

  VIO_Status  WS_create_window(
    VIO_STR                 title,
    int                    initial_x_pos,
    int                    initial_y_pos,
    int                    initial_x_size,
    int                    initial_y_size,
    VIO_BOOL                colour_map_mode,
    VIO_BOOL                double_buffer_flag,
    VIO_BOOL                depth_buffer_flag,
    int                    n_overlay_planes,
    VIO_BOOL                *actual_colour_map_mode,
    VIO_BOOL                *actual_double_buffer_flag,
    VIO_BOOL                *actual_depth_buffer_flag,
    int                    *actual_n_overlay_planes,
    WSwindow               window,
    WSwindow               parent );

  VIO_BOOL  WS_set_double_buffer_state(
    WSwindow               window,
    VIO_BOOL                double_buffer_flag );

  VIO_BOOL  WS_set_colour_map_state(
    WSwindow               window,
    VIO_BOOL                colour_map_flag );

  void  WS_set_window_title(
    WSwindow   window,
    VIO_STR     title );

  void  WS_delete_window(
    WSwindow  window );

  VIO_BOOL  WS_window_has_overlay_planes(
    WSwindow  window );

  void  WS_set_current_window(
    WSwindow          window );

  Window_id   WS_get_window_id(
    WSwindow  window );

  void  WS_set_bitplanes(
    WSwindow          window,
    Bitplane_types    bitplanes );

  int    WS_get_n_overlay_planes( void );

  void  WS_get_window_position(
    int          *x_pos,
    int          *y_pos );

  void  WS_get_window_size(
    int          *x_size,
    int          *y_size );

  void  glut_set_colour_entry(
    int      ind,
    VIO_Real     r,
    VIO_Real     g,
    VIO_Real     b );

  void  WS_set_colour_map_entry(
    WSwindow          window,
    Bitplane_types    bitplane,
    int               ind,
    VIO_Colour            colour );

  void  WS_set_overlay_colour_map_entry(
    WSwindow          window,
    int               ind,
    VIO_Colour            colour );

  void  WS_swap_buffers( void );

  void  WS_draw_text(
    Font_types  type,
    VIO_Real        size,
    VIO_STR      string );

  VIO_Real  WS_get_character_height(
    Font_types       type,
    VIO_Real             size );

  VIO_Real  WS_get_text_length(
    VIO_STR           str,
    Font_types       type,
    VIO_Real             size );

  void  WS_get_screen_size(
    int   *x_size, 
    int   *y_size  );

  void  WS_set_mouse_position(
    int       x_screen,
    int       y_screen );

  void  WS_set_update_function(
    void  (*func)( Window_id ) );

  void  WS_set_update_overlay_function(
    void  (*func)( Window_id ) );

  void  WS_set_resize_function(
    void  (*func)( Window_id, int, int, int, int ) );

  void  WS_set_key_down_function(
    void  (*func)( Window_id, int, int, int, int ) );

  void  WS_set_key_up_function(
    void  (*func)( Window_id, int, int, int, int ) );

  void  WS_set_mouse_movement_function(
    void  (*func)( Window_id, int, int ) );

  void  WS_set_left_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) );

  void  WS_set_left_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) );

  void  WS_set_middle_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) );

  void  WS_set_middle_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) );

  void  WS_set_right_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) );

  void  WS_set_right_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) );

  void WS_set_scroll_up_function(
    void  (*func)(Window_id, int, int, int) );

  void WS_set_scroll_down_function(
    void  (*func)(Window_id, int, int, int) );

  void  WS_set_iconify_function(
    void  (*func)( Window_id ) );

  void  WS_set_deiconify_function(
    void  (*func)( Window_id ) );

  void  WS_set_enter_function(
    void  (*func)( Window_id ) );

  void  WS_set_leave_function(
    void  (*func)( Window_id ) );

  void  WS_set_quit_function(
    void  (*func)( Window_id ) );

  void  WS_add_timer_function(
    VIO_Real          seconds,
    void          (*func) ( void * ),
    void          *data );

  void  WS_add_idle_function(
    void  (*func) ( void * ),
    void          *data );

  void  WS_remove_idle_function(
    void  (*func) ( void * ),
    void          *data );

  void  WS_event_loop( void );

void WS_exit_loop(void);

  void  WS_set_update_flag(
    WSwindow   window  );

  void  WS_set_visibility(WSwindow window, VIO_BOOL is_visible);

void  WS_set_geometry(WSwindow window, int x, int y, int width, int height);

  void  copy_X_colours(
    int  n_colours_to_copy );
#endif
