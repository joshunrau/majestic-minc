#include  <volume_io.h>
#include  <GS_graphics.h>
#include  <bicpl/queue.h>
#include  <gl/device.h>

#define  SLEEP_WHEN_IDLE_IN_EVENT_LOOP   0.02    /* 20 milliseconds */
 
static  void  (*display_callback) ( Window_id );
static  void  (*display_overlay_callback) ( Window_id );
static  void  (*resize_callback) ( Window_id, int, int, int, int );
static  void  (*key_down_callback) ( Window_id, int, int, int, int );
static  void  (*key_up_callback) ( Window_id, int, int, int, int );
static  void  (*mouse_motion_callback) ( Window_id, int, int );
static  void  (*left_down_callback) ( Window_id, int, int, int );
static  void  (*left_up_callback) ( Window_id, int, int, int );
static  void  (*middle_down_callback) ( Window_id, int, int, int );
static  void  (*middle_up_callback) ( Window_id, int, int, int );
static  void  (*right_down_callback) ( Window_id, int, int, int );
static  void  (*right_up_callback) ( Window_id, int, int, int );
static  void  (*iconify_callback) ( Window_id );
static  void  (*deiconify_callback) ( Window_id );
static  void  (*enter_callback) ( Window_id );
static  void  (*leave_callback) ( Window_id );
static  void  (*quit_callback) ( Window_id );

static  int        current_modifier = 0;
static  int        current_mouse_x = 0;
static  int        current_mouse_y = 0;
static  Window_id  event_window = 0;

static  GSwindow   *windows = NULL;
static  int        n_windows;

  void  GS_set_update_function(
    void  (*func)( Window_id ) )
{
    display_callback = func;
}

  void  GS_set_update_overlay_function(
    void  (*func)( Window_id ) )
{
    display_overlay_callback = func;
}

  void  GS_set_resize_function(
    void  (*func)( Window_id, int, int, int, int ) )
{
    resize_callback = func;
}

  void  GS_set_key_down_function(
    void  (*func)( Window_id, int, int, int, int ) )
{
    key_down_callback = func;
}

  void  GS_set_key_up_function(
    void  (*func)( Window_id, int, int, int, int ) )
{
    key_up_callback = func;
}

  void  GS_set_mouse_movement_function(
    void  (*func)( Window_id, int, int ) )
{
    mouse_motion_callback = func;
}

  void  GS_set_left_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) )
{
    left_down_callback = func;
}

  void  GS_set_left_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) )
{
    left_up_callback = func;
}

  void  GS_set_middle_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) )
{
    middle_down_callback = func;
}

  void  GS_set_middle_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) )
{
    middle_up_callback = func;
}

  void  GS_set_right_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) )
{
    right_down_callback = func;
}

  void  GS_set_right_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) )
{
    right_up_callback = func;
}

  void  GS_set_iconify_function(
    void  (*func)( Window_id ) )
{
    iconify_callback = func;
}

  void  GS_set_deiconify_function(
    void  (*func)( Window_id ) )
{
    deiconify_callback = func;
}

  void  GS_set_enter_function(
    void  (*func)( Window_id ) )
{
    enter_callback = func;
}

  void  GS_set_leave_function(
    void  (*func)( Window_id ) )
{
    leave_callback = func;
}

  void  GS_set_quit_function(
    void  (*func)( Window_id ) )
{
    quit_callback = func;
}

  void  GS_set_update_flag(
    GSwindow   window  )
{
    window->update_required = TRUE;
}

typedef struct
{
    void  (*function) ( void * );
    void              *data;
    VIO_Real              alarm_time;
} timer_info_struct;

static  timer_info_struct   *timers;
static  int                 n_timers = 0;
static  VIO_Real                min_time;

static  void   recompute_min_time( void )
{
    int   i;

    min_time = -1.0;

    for_less( i, 0, n_timers )
    {
        if( min_time < 0.0 || timers[i].alarm_time < min_time )
            min_time = timers[i].alarm_time;
    }
}

  void  GS_add_timer_function(
    VIO_Real          seconds,
    void          (*func) ( void * ),
    void          *data )
{
    timer_info_struct  info;

    info.function = func;
    info.data = data;
    info.alarm_time = current_realtime_seconds() + seconds;

    ADD_ELEMENT_TO_ARRAY( timers, n_timers, info, 1 );

    recompute_min_time();
}

typedef struct
{
    void  (*function) ( void * );
    void              *data;
} idle_info_struct;

static  idle_info_struct   *idles;
static  int                n_idles = 0;

  void  GS_add_idle_function(
    void          (*func) ( void * ),
    void          *data )
{
    idle_info_struct  info;

    info.function = func;
    info.data = data;

    ADD_ELEMENT_TO_ARRAY( idles, n_idles, info, 1 );
}

  void  GS_remove_idle_function(
    void          (*func) ( void * ),
    void          *data )
{
    int   i;

    for_less( i, 0, n_idles )
    {
        if( idles[i].function == func && idles[i].data == data )
            break;
    }

    if( i >= n_idles )
    {
        handle_internal_error( "WS_delete_idle_func" );
        return;
    }

    DELETE_ELEMENT_FROM_ARRAY( idles, n_idles, i, 1 );
}

  void  GS_event_loop( void )
{
    Window_id  save_win;
    long       dev, x_pos, y_pos, x_size, y_size;
    short      val;
    int        i;
    VIO_Real       current_time;
    VIO_BOOL    idle;

    while( TRUE )
    {
        idle = TRUE;

        for_less( i, 0, n_windows )
        {
            if( windows[i]->update_required )
            {
                windows[i]->update_required = FALSE;
                (*display_callback)( GS_get_window_id(windows[i]) );
                idle = FALSE;
            }
        }

        if( min_time >= 0.0 &&
            (current_time = current_realtime_seconds()) >= min_time )
        {
            i = 0;
            while( i < n_timers )
            {
                if( current_time >= timers[i].alarm_time )
                {
                    (*timers[i].function) ( timers[i].data );
                    DELETE_ELEMENT_FROM_ARRAY( timers, n_timers, i, 1 );
                    idle = FALSE;
                }
                else
                    ++i;
            }

            recompute_min_time();
        }

        if( qtest() )
        {
            idle = FALSE;
            dev = qread( &val );

            switch( dev )
            {
            case INPUTCHANGE:
                if( val == 0 )
                {
                    (*leave_callback) ( event_window );
                }
                else
                {
                    event_window = (Window_id) val;
                    (*enter_callback) ( (Window_id) val );
                }
                break;

            case WINQUIT:
            case WINSHUT:
                (*quit_callback) ( (Window_id) val );
                break;

            case DRAWOVERLAY:
                break;

            case DEPTHCHANGE:
                (*display_callback) ( (Window_id) val );
                (*display_overlay_callback) ( (Window_id) val );
                break;

            case WINTHAW:
                (*deiconify_callback) ( (Window_id) val );
                break;

            case WINFREEZE:
                (*iconify_callback) ( (Window_id) val );
                break;

            case REDRAW:
                save_win = winget();
                winset( (Window_id) val );
                getorigin( &x_pos, &y_pos );
                getsize( &x_size, &y_size );
                winset( save_win );
                (*resize_callback) ( (Window_id) val,
                                     (int) x_pos, (int) y_pos,
                                     (int) x_size, (int) y_size );
                break;

            case LEFTARROWKEY:
                if( val )
                {
                    (*key_down_callback)( event_window, LEFT_ARROW_KEY,
                                       current_mouse_x, 
                                       current_mouse_y, current_modifier );
                }
                else
                {
                    (*key_up_callback)( event_window, LEFT_ARROW_KEY,
                                     current_mouse_x, 
                                     current_mouse_y, current_modifier );
                }
                break;

            case RIGHTARROWKEY:
                if( val )
                {
                    (*key_down_callback)( event_window, RIGHT_ARROW_KEY,
                                       current_mouse_x, 
                                       current_mouse_y, current_modifier );
                }
                else
                {
                    (*key_up_callback)( event_window, RIGHT_ARROW_KEY,
                                     current_mouse_x, 
                                     current_mouse_y, current_modifier );
                }
                break;

            case DOWNARROWKEY:
                if( val )
                {
                    (*key_down_callback)( event_window, DOWN_ARROW_KEY,
                                       current_mouse_x, 
                                       current_mouse_y, current_modifier );
                }
                else
                {
                    (*key_up_callback)( event_window, DOWN_ARROW_KEY,
                                     current_mouse_x, 
                                     current_mouse_y, current_modifier );
                }
                break;

            case UPARROWKEY:
                if( val )
                {
                    (*key_down_callback)( event_window, UP_ARROW_KEY,
                                       current_mouse_x, 
                                       current_mouse_y, current_modifier );
                }
                else
                {
                    (*key_up_callback)( event_window, UP_ARROW_KEY,
                                     current_mouse_x, 
                                     current_mouse_y, current_modifier );
                }
                break;

            case LEFTSHIFTKEY:
            case RIGHTSHIFTKEY:
                if( val )
                    current_modifier |= SHIFT_KEY_BIT;
                else
                    current_modifier &= ~SHIFT_KEY_BIT;
                break;

            case LEFTCTRLKEY:
            case RIGHTCTRLKEY:
                if( val )
                    current_modifier |= CTRL_KEY_BIT;
                else
                    current_modifier &= ~CTRL_KEY_BIT;
                break;

            case LEFTALTKEY:
            case RIGHTALTKEY:
                if( val )
                    current_modifier |= ALT_KEY_BIT;
                else
                    current_modifier &= ~ALT_KEY_BIT;
                break;

            case KEYBD:
                (*key_down_callback)( event_window, (int) val,
                                   current_mouse_x, 
                                   current_mouse_y, current_modifier );
                break;

            case MOUSEX:
            case MOUSEY:
                if( dev == MOUSEX )
                    current_mouse_x = (int) val;
                else
                    current_mouse_y = (int) val;

                (*mouse_motion_callback)( event_window, current_mouse_x, 
                                       current_mouse_y );
                break;

            case LEFTMOUSE:
                if( val )
                {
                    (*left_down_callback)( event_window, current_mouse_x, 
                                        current_mouse_y, current_modifier );
                }
                else
                {
                    (*left_up_callback)( event_window, current_mouse_x, 
                                      current_mouse_y, current_modifier );
                }
                break;

            case MIDDLEMOUSE:
                if( val )
                {
                    (*middle_down_callback)( event_window, current_mouse_x, 
                                          current_mouse_y, current_modifier );
                }
                else
                {
                    (*middle_up_callback)( event_window, current_mouse_x, 
                                        current_mouse_y, current_modifier );
                }
                break;

            case RIGHTMOUSE:
                if( val )
                {
                    (*right_down_callback)( event_window, current_mouse_x, 
                                         current_mouse_y, current_modifier );
                }
                else
                {
                    (*right_up_callback)( event_window, current_mouse_x, 
                                       current_mouse_y, current_modifier );
                }
                break;

            case QFULL:
                print_error( "Error: GL event queue full.\n" );
                break;

            default:
                break;
            }
        }

        if( idle )
        {
            if( n_idles > 0 )
            {
                for_less( i, 0, n_idles )
                    (*idles[i].function) ( idles[i].data );
            }
            else
            {
                sleep_program( SLEEP_WHEN_IDLE_IN_EVENT_LOOP );
            }
        }
    }   
}

  void  initialize_window_events(
    GSwindow   window )
{
    window->update_required = TRUE;

    qdevice( KEYBD );

    qdevice( LEFTARROWKEY );
    qdevice( RIGHTARROWKEY );
    qdevice( DOWNARROWKEY );
    qdevice( UPARROWKEY );
    qdevice( LEFTSHIFTKEY );
    qdevice( RIGHTSHIFTKEY );
    qdevice( LEFTALTKEY );
    qdevice( RIGHTALTKEY );
    qdevice( LEFTCTRLKEY );
    qdevice( RIGHTCTRLKEY );

    qdevice( REDRAW );
    qdevice( DEPTHCHANGE );
    qdevice( DRAWOVERLAY );
    qdevice( INPUTCHANGE );

    qdevice( MOUSEX );
    qdevice( MOUSEY );
    qdevice( LEFTMOUSE );
    qdevice( MIDDLEMOUSE );
    qdevice( RIGHTMOUSE );

    qdevice( WINFREEZE );
    qdevice( WINTHAW );
    qdevice( WINQUIT );
    qdevice( WINSHUT );

    qdevice( QFULL );

    ADD_ELEMENT_TO_ARRAY( windows, n_windows, window, 1 );
}

  void  delete_window_events(
    GSwindow   window )
{
    int  i;

    for_less( i, 0, n_windows )
    {
        if( windows[i] == window )
            break;
    }

    if( i >= n_windows )
    {
        handle_internal_error( "delete_window_events" );
        return;
    }

    DELETE_ELEMENT_FROM_ARRAY( windows, n_windows, i, 1 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_set_mouse_position
@INPUT      : x_screen
              y_screen
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the mouse position relative to the full screen.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  void  GS_set_mouse_position(
    int   x_screen,
    int   y_screen )
{
    setvaluator( MOUSEX, (short) x_screen, 0, 10000 );
    setvaluator( MOUSEY, (short) y_screen, 0, 10000 );
}

  VIO_BOOL  GS_are_mouse_coordinates_in_screen_space( void )
{
    return( TRUE );
}
