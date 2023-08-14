#ifndef  DEF_x_window_prototypes
#define  DEF_x_window_prototypes

  Display   *X_get_display( void );

  int  X_get_screen( void );

  VIO_Status  X_create_window_with_visual(
    VIO_STR           title,
    int              initial_x_pos,
    int              initial_y_pos,
    int              initial_x_size,
    int              initial_y_size,
    VIO_BOOL          colour_map_mode,  
    XVisualInfo      *visual,
    Colormap         cmap,
    X_window_struct  *window );

  VIO_Status  X_create_overlay_window(
    X_window_struct  *window,
    XVisualInfo      *visual,
    X_window_struct  *overlay_window );

  void  X_delete_window(
    X_window_struct  *window );

  int  X_get_n_colours(
    X_window_struct          *window );

  void  X_set_colour_map_entry(
    X_window_struct          *window,
    int                      ind,
    int                      r,
    int                      g,
    int                      b );

  void  X_get_window_geometry(
    X_window_struct          *window,
    int                      *x_position,
    int                      *y_position,
    int                      *x_size,
    int                      *y_size );

  void  X_get_screen_size(
    int    *x_size,
    int    *y_size );

  VIO_BOOL  X_get_event(
    Event_types          *event_type,
    Window_id            *window,
    event_info_struct    *info );

  VIO_BOOL  X_get_font_name(
    Font_types       type,
    int              size,
    VIO_STR           *font_name );

  void  X_set_mouse_position(
    int                      x,
    int                      y );
#endif
