#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H  

#include  <volume_io.h>
#include  <GS_graphics.h>

  void  GS_set_update_function(
    void  (*func)( Window_id ) )
{
    WS_set_update_function( func );
}

  void  GS_set_update_overlay_function(
    void  (*func)( Window_id ) )
{
    WS_set_update_overlay_function( func );
}

  void  GS_set_resize_function(
    void  (*func)( Window_id, int, int, int, int ) )
{
    WS_set_resize_function( func );
}

  void  GS_set_key_down_function(
    void  (*func)( Window_id, int, int, int, int ) )
{
    WS_set_key_down_function( func );
}

  void  GS_set_key_up_function(
    void  (*func)( Window_id, int, int, int, int ) )
{
    WS_set_key_up_function( func );
}

  void  GS_set_mouse_movement_function(
    void  (*func)( Window_id, int, int ) )
{
    WS_set_mouse_movement_function( func );
}

  void  GS_set_left_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) )
{
    WS_set_left_mouse_down_function( func );
}

  void  GS_set_left_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) )
{
    WS_set_left_mouse_up_function( func );
}

  void  GS_set_middle_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) )
{
    WS_set_middle_mouse_down_function( func );
}

  void  GS_set_middle_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) )
{
    WS_set_middle_mouse_up_function( func );
}

  void  GS_set_right_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) )
{
    WS_set_right_mouse_down_function( func );
}

  void  GS_set_right_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) )
{
    WS_set_right_mouse_up_function( func );
}

  void  GS_set_scroll_down_function(
    void  (*func)( Window_id, int, int, int ) )
{
    WS_set_scroll_down_function( func );
}

  void  GS_set_scroll_up_function(
    void  (*func)( Window_id, int, int, int ) )
{
    WS_set_scroll_up_function( func );
}

  void  GS_set_iconify_function(
    void  (*func)( Window_id ) )
{
    WS_set_iconify_function( func );
}

  void  GS_set_deiconify_function(
    void  (*func)( Window_id ) )
{
    WS_set_deiconify_function( func );
}

  void  GS_set_enter_function(
    void  (*func)( Window_id ) )
{
    WS_set_enter_function( func );
}

  void  GS_set_leave_function(
    void  (*func)( Window_id ) )
{
    WS_set_leave_function( func );
}

  void  GS_set_quit_function(
    void  (*func)( Window_id ) )
{
    WS_set_quit_function( func );
}

  void  GS_set_update_flag(
    GSwindow   window  )
{
    WS_set_update_flag( &window->WS_window );
}

  void  GS_add_timer_function(
    VIO_Real          seconds,
    void          (*func) ( void * ),
    void          *data )
{
    WS_add_timer_function( seconds, func, data );
}

  void  GS_add_idle_function(
    void  (*func) ( void * ),
    void          *data )
{
    WS_add_idle_function( func, data );
}

  void  GS_remove_idle_function(
    void  (*func) ( void * ),
    void          *data )
{
    WS_remove_idle_function( func, data );
}

  void  GS_event_loop( void )
{
    WS_event_loop();
}

void GS_exit_loop( void )
{
  WS_exit_loop();
}

  void  GS_set_mouse_position(
    int       x_screen,
    int       y_screen )
{
    WS_set_mouse_position( x_screen, y_screen );
}

  VIO_BOOL  GS_are_mouse_coordinates_in_screen_space( void )
{
    return( FALSE );
}
