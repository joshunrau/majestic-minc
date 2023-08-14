 #ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H  

#include  <volume_io.h>
#include  <WS_graphics.h>

#define  DEFAULT_WINDOW_X_SIZE    300
#define  DEFAULT_WINDOW_Y_SIZE    300

static  void  set_event_callbacks_for_current_window( void );


// TO QUERY SCALING FACTOR
#if !defined(GLUT_WINDOW_SCALE)
#	define GLUT_WINDOW_SCALE 199
#endif

#define MAX_GLUT_ID 9
/**
 * This array is used to store pointers to our private state information.
 * It is indexed by the GLUT window ID, so it relies on our knowing that
 * GLUT window ID's are small integers.
 */
static WSwindow windows[MAX_GLUT_ID + 1];

/*----------------------------------------------------- */

void  WS_initialize( void )
{
    static  VIO_BOOL  initialized = FALSE;
    int              argc = 1;
    char             *argv[] = { "no_program_name" };

    if( !initialized )
    {
        int s0;
        initialized = TRUE;
        glutInit( &argc, argv );
        /*ONLY ON APPLE?*/
        /*
        glutInitDisplayString("rgba double hidpi");
        s0 = glutGet(GLUT_WINDOW_SCALE);
        printf("s0=%d\n",s0);
        */
    }
}

static  Window_id  create_GLUT_window(
    VIO_STR                title,
    int                    initial_x_pos,
    int                    initial_y_pos,
    int                    initial_x_size,
    int                    initial_y_size,
    VIO_BOOL               colour_map_mode,
    VIO_BOOL               double_buffer_flag,
    VIO_BOOL               depth_buffer_flag,
    int                    n_overlay_planes,
    VIO_BOOL               *actual_colour_map_mode_ptr,
    VIO_BOOL               *actual_double_buffer_flag_ptr,
    VIO_BOOL               *actual_depth_buffer_flag_ptr,
    int                    *actual_n_overlay_planes_ptr,
    Window_id              parent_id)
{
    unsigned  int      mode;
    Window_id          window_id;
    
    if( initial_x_pos >= 0 && initial_y_pos >= 0 )
        glutInitWindowPosition( initial_x_pos, initial_y_pos );
    else
        glutInitWindowPosition( -1, -1 );

    /*
     * Set the initial window size. If the specified value is 
     * less than or equal to zero, we use a value derived from
     * the screen size.
     */
    if( initial_x_size <= 0 || initial_y_size <= 0 )
    {
        int screen_x = glutGet(GLUT_SCREEN_WIDTH);
        int screen_y = glutGet(GLUT_SCREEN_HEIGHT);
        if (screen_x <= 0 || screen_y <= 0)
        {
            initial_x_size = DEFAULT_WINDOW_X_SIZE;
            initial_y_size = DEFAULT_WINDOW_Y_SIZE;
        }
        else if (screen_x > screen_y)
        {
            initial_x_size = initial_y_size = screen_y / 2;
        }
        else
        {
            initial_x_size = initial_y_size = screen_x / 2;
        }
    }
    glutInitWindowSize( initial_x_size, initial_y_size );

    // Set up the mode.
    mode = GLUT_RGB;
    if (double_buffer_flag)
    {
      mode |= GLUT_DOUBLE;
    }
    if (depth_buffer_flag)
    {
      mode |= GLUT_DEPTH;
    }
    glutInitDisplayMode(mode);

    if (parent_id > 0)
    {
      window_id = glutCreateSubWindow( parent_id, initial_x_pos, initial_y_pos,
                                       initial_x_size, initial_y_size );
    }
    else
    {
      window_id = glutCreateWindow( title );
      glutPopWindow();
    }

    if( window_id < 1 )
    { 
        print_error( "Could not open GLUT window for OpenGL\n" );
        return( -1 );
    }

    /* Freeglut does not support overlays.
     */
    set_event_callbacks_for_current_window();

    /* Some NVIDIA systems lie about RGBA mode, so we have to 
     * ignore the value returned by GLUT_WINDOW_RGBA.
     */
    if (actual_colour_map_mode_ptr != NULL)
      *actual_colour_map_mode_ptr = FALSE;
    if (actual_double_buffer_flag_ptr != NULL)
      *actual_double_buffer_flag_ptr = glutGet(GLUT_WINDOW_DOUBLEBUFFER);
    if (actual_depth_buffer_flag_ptr != NULL)
      *actual_depth_buffer_flag_ptr = glutGet(GLUT_WINDOW_DEPTH_SIZE) > 0;
    if (actual_n_overlay_planes_ptr != NULL)
      *actual_n_overlay_planes_ptr = 0;
    return( window_id );
}

static  void  delete_GLUT_window(
    Window_id   window_id )
{
    glutDestroyWindow( window_id );
}

VIO_Status  WS_create_window(
    VIO_STR                title,
    int                    initial_x_pos,
    int                    initial_y_pos,
    int                    initial_x_size,
    int                    initial_y_size,
    VIO_BOOL               colour_map_mode,
    VIO_BOOL               double_buffer_flag,
    VIO_BOOL               depth_buffer_flag,
    int                    n_overlay_planes,
    VIO_BOOL               *actual_colour_map_mode,
    VIO_BOOL               *actual_double_buffer_flag,
    VIO_BOOL               *actual_depth_buffer_flag,
    int                    *actual_n_overlay_planes,
    WSwindow               window,
    WSwindow               parent )
{
    window->window_id = create_GLUT_window( title, 
                                            initial_x_pos,
                                            initial_y_pos,
                                            initial_x_size,
                                            initial_y_size,
                                            colour_map_mode,
                                            double_buffer_flag,
                                            depth_buffer_flag,
                                            n_overlay_planes,
                                            actual_colour_map_mode,
                                            actual_double_buffer_flag,
                                            actual_depth_buffer_flag,
                                            actual_n_overlay_planes,
                                            (parent != NULL) ? parent->window_id : 0);


    if( window->window_id < 1 )
        return( VIO_ERROR );

    glutSetWindow( window->window_id );

    window->init_x = glutGet(GLUT_WINDOW_X);
    window->init_y = glutGet(GLUT_WINDOW_Y);
    window->border_width = 0;
    window->border_height = 0;
    window->is_visible = TRUE;
    window->is_new = TRUE;

    if (window->window_id >= 0 && window->window_id <= MAX_GLUT_ID)
    {
      windows[window->window_id] = window;
    }
    return( VIO_OK );
}

VIO_BOOL  WS_set_double_buffer_state(
    WSwindow               window,
    VIO_BOOL               double_buffer_flag )
{
    return( double_buffer_flag );
}

VIO_BOOL  WS_set_colour_map_state(
    WSwindow               window,
    VIO_BOOL                colour_map_flag )
{
    return( colour_map_flag );
}

  void  WS_set_window_title(
    WSwindow   window,
    VIO_STR     title )
{
    glutSetWindowTitle( title );
}


  void  WS_delete_window(
    WSwindow  window )
{
    delete_GLUT_window( window->window_id );
}

static  Window_id  get_current_event_window( void )
{
    return( glutGetWindow() );
}

  VIO_BOOL  WS_window_has_overlay_planes(
    WSwindow  window )
{
    return( FALSE );
}

static  void  set_window_normal_planes(
    WSwindow  window )
{
    glutSetWindow( window->window_id );
    if (!window->is_visible) 
        glutHideWindow();
}

static  void  set_window_overlay_planes(
    WSwindow  window )
{
    glutSetWindow( window->window_id );
}

  void  WS_set_current_window(
    WSwindow          window )
{
    WS_set_bitplanes( window, NORMAL_PLANES );
}

  Window_id   WS_get_window_id(
    WSwindow  window )
{
    return( window->window_id );
}

  void  WS_set_bitplanes(
    WSwindow          window,
    Bitplane_types    bitplanes )
{
    if( bitplanes == OVERLAY_PLANES )
        set_window_overlay_planes( window );
    else
        set_window_normal_planes( window );
}


  int    WS_get_n_overlay_planes( void )
{
    return( 0 );
}

/**
 * Get the window position. 
 *
 * Glut is a bit weird in how it sets the position, in that the
 * position we pass to glutInitWindowPosition() will NOT be the final
 * window position reported by GLUT_WINDOW_X and GLUT_WINDOW_Y,
 * because these give the coordinates of the "user" area of the
 * window, minus the border or "decorations". We deduce the size of
 * the decorations and adjust the window position accordingly here, in
 * order to make the values returned from this function identical to
 * the values passed into create_GLUT_window().
 */
void  WS_get_window_position(
    int          *x_pos,
    int          *y_pos )
{
    Window_id window_id = get_current_event_window();
    *x_pos = glutGet( GLUT_WINDOW_X );
    *y_pos = glutGet( GLUT_WINDOW_Y );
    if (window_id >= 0 && window_id <= MAX_GLUT_ID)
    {
      WSwindow window = windows[window_id];
      if (window != NULL)
      {
        *x_pos -= window->border_width;
        *y_pos -= window->border_height;
      }
    }
}

  void  WS_get_window_size(
    int          *x_size,
    int          *y_size )
{
    *x_size = glutGet( (GLenum) GLUT_WINDOW_WIDTH );
    *y_size = glutGet( (GLenum) GLUT_WINDOW_HEIGHT );
}

void  glut_set_colour_entry(
    int      ind,
    VIO_Real     r,
    VIO_Real     g,
    VIO_Real     b )
{
    glutSetColor( ind, (float) r, (float) g, (float) b );
}

void  WS_set_colour_map_entry(
    WSwindow          window,
    Bitplane_types    bitplane,
    int               ind,
    VIO_Colour            colour )
{
    glut_set_colour_entry( ind,
                           get_Colour_r_0_1(colour),
                           get_Colour_g_0_1(colour),
                           get_Colour_b_0_1(colour) );
}

  void  WS_set_overlay_colour_map_entry(
    WSwindow          window,
    int               ind,
    VIO_Colour            colour )
{
}

  void  WS_swap_buffers( void )
{
  if(glutGetWindow()) //VF:a hack
      glutSwapBuffers();
}

static  struct
        {
            int     height;
            void    *font;
        }        known_fonts[] =
              {
                {10, GLUT_BITMAP_HELVETICA_10},
                {12, GLUT_BITMAP_HELVETICA_12},
                {13,  GLUT_BITMAP_8_BY_13},
                {15 , GLUT_BITMAP_9_BY_15},
                {18, GLUT_BITMAP_HELVETICA_18},
                {24, GLUT_BITMAP_TIMES_ROMAN_24},
                {10, GLUT_BITMAP_TIMES_ROMAN_10}
              };


static  void  *lookup_font(
    Font_types       type,
    VIO_Real             size,
    int              *actual_height )
{
    VIO_Real    diff, min_diff;
    int     which, best;
    void    *font;

    if( type == FIXED_FONT )
    {
        font = GLUT_BITMAP_8_BY_13;
        if( actual_height != NULL )
            *actual_height = 13;
    }
    else
    {
        min_diff = 0.0;
        best = 0;
        for_less( which, 0, VIO_SIZEOF_STATIC_ARRAY(known_fonts) )
        {
            diff = VIO_FABS( (VIO_Real) known_fonts[which].height - size );
            if( which == 0 || diff < min_diff )
            {
                best = which;
                min_diff = diff;
            }
        }

        font = known_fonts[best].font;
        if( actual_height != NULL )
            *actual_height = known_fonts[best].height;
    }

    return( font );
}

  void  WS_draw_text(
    Font_types  type,
    VIO_Real        size,
    VIO_STR      string )
{
    int   i;
    void  *font;

    font = lookup_font( type, size, NULL );

    for_less( i, 0, string_length( string ) )
        glutBitmapCharacter( font, (int) string[i] );
}

  VIO_Real  WS_get_character_height(
    Font_types       type,
    VIO_Real             size )
{
    int   height;

    (void) lookup_font( type, size, &height );

    return( (VIO_Real) height );
}

  VIO_Real  WS_get_text_length(
    VIO_STR           str,
    Font_types       type,
    VIO_Real             size )
{
    int    i, len;
    void   *font;

    font = lookup_font( type, size, NULL );

    len = 0;
    for_less( i, 0, (int) strlen( str ) )
        len += glutBitmapWidth( font, (int) str[i] );

    return( (VIO_Real) len );
}

  void  WS_get_screen_size(
    int   *x_size, 
    int   *y_size  )
{
    *x_size = glutGet( (GLenum) GLUT_SCREEN_WIDTH );
    *y_size = glutGet( (GLenum) GLUT_SCREEN_HEIGHT );
}

  void  WS_set_mouse_position(
    int       x_screen,
    int       y_screen )
{
}

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
static  void  (*scroll_up_callback)(Window_id, int, int, int);
static  void  (*scroll_down_callback)(Window_id, int, int, int);

  void  WS_set_update_function(
    void  (*func)( Window_id ) )
{
    display_callback = func;
}

  void  WS_set_update_overlay_function(
    void  (*func)( Window_id ) )
{
    display_overlay_callback = func;
}

  void  WS_set_resize_function(
    void  (*func)( Window_id, int, int, int, int ) )
{
    resize_callback = func;
}

  void  WS_set_key_down_function(
    void  (*func)( Window_id, int, int, int, int ) )
{
    key_down_callback = func;
}

  void  WS_set_key_up_function(
    void  (*func)( Window_id, int, int, int, int ) )
{
    key_up_callback = func;
}

  void  WS_set_mouse_movement_function(
    void  (*func)( Window_id, int, int ) )
{
    mouse_motion_callback = func;
}

  void  WS_set_left_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) )
{
    left_down_callback = func;
}

  void  WS_set_left_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) )
{
    left_up_callback = func;
}

  void  WS_set_middle_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) )
{
    middle_down_callback = func;
}

  void  WS_set_middle_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) )
{
    middle_up_callback = func;
}

  void  WS_set_right_mouse_down_function(
    void  (*func)( Window_id, int, int, int ) )
{
    right_down_callback = func;
}

  void  WS_set_right_mouse_up_function(
    void  (*func)( Window_id, int, int, int ) )
{
    right_up_callback = func;
}

void WS_set_scroll_up_function(
   void (*func)(Window_id, int, int, int) )
{
    scroll_up_callback = func;
}

void WS_set_scroll_down_function(
   void (*func)(Window_id, int, int, int) )
{
    scroll_down_callback = func;
}
                    
  void  WS_set_iconify_function(
    void  (*func)( Window_id ) )
{
    iconify_callback = func;
}

  void  WS_set_deiconify_function(
    void  (*func)( Window_id ) )
{
    deiconify_callback = func;
}

  void  WS_set_enter_function(
    void  (*func)( Window_id ) )
{
    enter_callback = func;
}

  void  WS_set_leave_function(
    void  (*func)( Window_id ) )
{
    leave_callback = func;
}

  void  WS_set_quit_function(
    void  (*func)( Window_id ) )
{
    quit_callback = func;
}

static  int  get_keyboard_modifiers( void )
{
    int  modifier, glut_mod;

    modifier = 0;
    glut_mod = glutGetModifiers();

    if( (glut_mod & GLUT_ACTIVE_SHIFT) != 0 )
        modifier |= SHIFT_KEY_BIT;
    if( (glut_mod & GLUT_ACTIVE_CTRL) != 0 )
        modifier |= CTRL_KEY_BIT;
    if( (glut_mod & GLUT_ACTIVE_ALT) != 0 )
        modifier |= ALT_KEY_BIT;

    return( modifier );
}

static  int  flip_window_y(
    int   window_y )
{
    return( glutGet( (GLenum) GLUT_WINDOW_HEIGHT ) - 1 - window_y );
}

static  void  display_function( void )
{
    (*display_callback) ( get_current_event_window() );
}

static  void  resize_function(
    int   width,
    int   height )
{
    Window_id  window_id;
    int        x, y;

    x = glutGet( (GLenum) GLUT_WINDOW_X );
    y = glutGet( (GLenum) GLUT_WINDOW_Y );

    window_id = get_current_event_window();

    if (window_id >= 0 && window_id <= MAX_GLUT_ID)
    {
      WSwindow window = windows[window_id];
      if (window != NULL)
      {
        if (window->is_new)
        {
          window->is_new = FALSE;
          if (x >= window->init_x)
          {
            window->border_width = x - window->init_x;
          }
          if (y >= window->init_y) 
          {
            window->border_height = y - window->init_y;
          }
        }

        x -= window->border_width;
        y -= window->border_height;
      }
    }

    (*resize_callback) ( window_id, x, y, width, height );
}

static  void  keyboard_function(
    unsigned  char  key,
    int             x,
    int             y )
{
    y = flip_window_y( y );
    (*key_down_callback) ( get_current_event_window(), (int) key, x, y,
                           get_keyboard_modifiers() );
}

static  void  special_keyboard_function(
    int     key,
    int     x,
    int     y )
{
    int   translated = -1;

    y = flip_window_y( y );

    /* Glut's special key codes overlap with normal ASCII codes. We
     * translate them into codes > 128 so that they can be dealt with
     * using a single keyboard handler.
     */
    switch( key )
    {
    case GLUT_KEY_F1:       translated = BICGL_F1_KEY; break;
    case GLUT_KEY_F2:       translated = BICGL_F2_KEY; break;
    case GLUT_KEY_F3:       translated = BICGL_F3_KEY; break;
    case GLUT_KEY_F4:       translated = BICGL_F4_KEY; break;
    case GLUT_KEY_F5:       translated = BICGL_F5_KEY; break;
    case GLUT_KEY_F6:       translated = BICGL_F6_KEY; break;
    case GLUT_KEY_F7:       translated = BICGL_F7_KEY; break;
    case GLUT_KEY_F8:       translated = BICGL_F8_KEY; break;
    case GLUT_KEY_F9:       translated = BICGL_F9_KEY; break;
    case GLUT_KEY_F10:      translated = BICGL_F10_KEY; break;
    case GLUT_KEY_F11:      translated = BICGL_F11_KEY; break;
    case GLUT_KEY_F12:      translated = BICGL_F12_KEY; break;
    case GLUT_KEY_LEFT:     translated = LEFT_ARROW_KEY;  break;
    case GLUT_KEY_UP:       translated = UP_ARROW_KEY;  break;
    case GLUT_KEY_RIGHT:    translated = RIGHT_ARROW_KEY;  break;
    case GLUT_KEY_DOWN:     translated = DOWN_ARROW_KEY;  break;
    case GLUT_KEY_PAGE_UP:  translated = BICGL_PGUP_KEY; break;
    case GLUT_KEY_PAGE_DOWN:translated = BICGL_PGDN_KEY; break;
    case GLUT_KEY_HOME:     translated = BICGL_HOME_KEY; break;
    case GLUT_KEY_END:      translated = BICGL_END_KEY; break;
    case GLUT_KEY_INSERT:   translated = BICGL_INSERT_KEY; break;
    default:
      return;
    }
    (*key_down_callback) ( get_current_event_window(), translated,
                           x, y, get_keyboard_modifiers() );
}

static  void  mouse_button_function(
    int     button,
    int     state,
    int     x,
    int     y )
{
    int         modifiers;
    Window_id   window_id;

    window_id = get_current_event_window();
    modifiers = get_keyboard_modifiers();
    y = flip_window_y( y );

    switch( button )
    {
    case GLUT_LEFT_BUTTON:
        if( state == GLUT_DOWN )
            (*left_down_callback) ( window_id, x, y, modifiers );
        else
            (*left_up_callback) ( window_id, x, y, modifiers );
        break;

    case GLUT_MIDDLE_BUTTON:
        if( state == GLUT_DOWN )
            (*middle_down_callback) ( window_id, x, y, modifiers );
        else
            (*middle_up_callback) ( window_id, x, y, modifiers );
        break;

    case GLUT_RIGHT_BUTTON:
        if( state == GLUT_DOWN )
            (*right_down_callback) ( window_id, x, y, modifiers );
        else
            (*right_up_callback) ( window_id, x, y, modifiers );
        break;

    case 3:
      if (state == GLUT_DOWN)
          (*scroll_up_callback)(window_id, x, y, modifiers);
      break;

    case 4:
      if (state == GLUT_DOWN)
          (*scroll_down_callback)(window_id, x, y, modifiers);
      break;
    }
}

static  void  mouse_motion_function(
    int     x,
    int     y )
{
    y = flip_window_y(y);
    (*mouse_motion_callback) ( get_current_event_window(), x, y );
}

static  void  entry_function(
    int     state )
{
    Window_id   window_id;

    window_id = get_current_event_window();

    if( state == GLUT_LEFT )
        (*leave_callback) ( window_id );
    else
        (*enter_callback) ( window_id );
}


static  void  set_event_callbacks_for_current_window(void)
{
    glutDisplayFunc( display_function );
    glutReshapeFunc( resize_function );
    glutKeyboardFunc( keyboard_function );
    glutSpecialFunc( special_keyboard_function );
    glutMouseFunc( mouse_button_function );
    glutMotionFunc( mouse_motion_function );
    glutPassiveMotionFunc( mouse_motion_function );
    glutEntryFunc( entry_function );
}

typedef struct
{
    void  (*function) ( void * );
    void              *data;
    VIO_BOOL           active;
} callback_info_struct;

static  callback_info_struct   *timers;
static  int                    n_timers = 0;

static  void  global_timer_function(
    int   index )
{
    if( index < 0 || index >= n_timers )
    {
        handle_internal_error( "global_timer_function" );
        return;
    }

    (*timers[index].function)( timers[index].data );

    timers[index].active = FALSE;
}

  void  WS_add_timer_function(
    VIO_Real          seconds,
    void          (*func) ( void * ),
    void          *data )
{
    int                   i;

    for_less( i, 0, n_timers )
    {
        if( !timers[i].active )
            break;
    }

    if( i >= n_timers )
    {
        SET_ARRAY_SIZE( timers, n_timers, n_timers+1, 1 );
        ++n_timers;
    }

    timers[i].function = func;
    timers[i].data = data;
    timers[i].active = TRUE;

    glutTimerFunc( (unsigned int) (1000.0 * seconds + 0.5),
                   global_timer_function, i );
}

static  callback_info_struct   *idles;
static  int                    n_idles = 0;

static  void  global_idle_function( void )
{
    int   i;

    for_less( i, 0, n_idles )
    {
        (*idles[i].function) ( idles[i].data );
    }
}

  void  WS_add_idle_function(
    void  (*func) ( void * ),
    void          *data )
{
    callback_info_struct  info;

    if( n_idles == 0 )
        glutIdleFunc( global_idle_function );

    info.function = func;
    info.data = data;

    ADD_ELEMENT_TO_ARRAY( idles, n_idles, info, 1 );
}

  void  WS_remove_idle_function(
    void  (*func) ( void * ),
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

    if( n_idles == 0 )
        glutIdleFunc( NULL );
}

  void  WS_event_loop( void )
{
#ifdef GLUT_ACTION_GLUTMAINLOOP_RETURNS
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
                  GLUT_ACTION_GLUTMAINLOOP_RETURNS);
#endif /* GLUT_ACTION_ON_WINDOW_CLOSE */
    glutMainLoop();
}

void WS_exit_loop(void)
{
#ifdef GLUT_ACTION_GLUTMAINLOOP_RETURNS
   glutLeaveMainLoop();
#else
   exit( 0 );
#endif
}

  void  WS_set_update_flag(
    WSwindow   window  )
{
    glutPostRedisplay();
}

void WS_set_geometry(WSwindow window, int x, int y, int cx, int cy)
{
    glutSetWindow( window->window_id );
    if (x >= 0 && y >= 0)
      glutPositionWindow(x, y);
    if (cx > 0 && cy > 0)
      glutReshapeWindow(cx, cy);
}

void WS_set_visibility(WSwindow window, VIO_BOOL is_visible)
{
  window->is_visible = is_visible;

  glutSetWindow( window->window_id );
  if (is_visible)
  {
    glutShowWindow();
  }
  else
  {
    glutHideWindow();
  }
}
