#include  <internal_volume_io.h>
#include  <x_windows.h>
#include  <X11/keysym.h>

static  void  bind_special_keys( void );

  Display   *X_get_display( void )
{
    static  VIO_BOOL  first = TRUE;
    static  Display  *display;

    if( first )
    {
        first = FALSE;
        display = XOpenDisplay( 0 );

        if( display == NULL )
        {
            print_error( "Aborting, could not open X display.\n" );
            exit( 1 );
        }
    }

    return( display );
}

  int  X_get_screen( void )
{
    static  VIO_BOOL  first = TRUE;
    static  int      default_screen;

    if( first )
    {
        first = FALSE;
        default_screen = DefaultScreen( X_get_display() );
    }

    return( default_screen );
}

/* ARGSUSED */

static   int   wait_for_window(
    Display *display,
    XEvent  *event,
    char    *arg )
{
    Window   *window;

    window = (void *) arg;

    if( event->type == Expose && event->xexpose.window == *window )
        return( TRUE );
    else
        return( FALSE );
}

  Status  X_create_window_with_visual(
    STRING           title,
    int              initial_x_pos,
    int              initial_y_pos,
    int              initial_x_size,
    int              initial_y_size,
    VIO_BOOL          colour_map_mode,  
    XVisualInfo      *visual,
    Colormap         cmap,
    X_window_struct  *window )
{
    VIO_BOOL                  setting_position;
    long                     event_mask;
    unsigned long            cwa_mask;
    XEvent                   event;
    XSetWindowAttributes     cwa;
    XSizeHints               normal_hints;

    if( initial_x_size <= 0 )
        initial_x_size = 100;

    if( initial_y_size <= 0 )
        initial_y_size = 100;

    cwa.colormap = cmap;
    cwa.border_pixel = 0;
    cwa.event_mask = StructureNotifyMask | ExposureMask;
    cwa_mask = CWColormap|CWBorderPixel;

    if( initial_x_pos < 0 || initial_y_pos < 0 )
    {
        initial_x_pos = 0;
        initial_y_pos = 0;
        setting_position = FALSE;
    }
    else
    {
        int    screen_x_size, screen_y_size;

        X_get_screen_size( &screen_x_size, &screen_y_size );

        if( initial_y_size >= 0 )
            initial_y_pos = screen_y_size-1 - (initial_y_pos+initial_y_size-1);
        else
            initial_y_pos = screen_y_size - 1 - initial_y_pos;

        setting_position = TRUE;
    }

    window->window_id = XCreateWindow( X_get_display(),
                               RootWindow( X_get_display(),
                                           X_get_screen() ),
                               initial_x_pos, initial_y_pos,
                               (unsigned int) initial_x_size,
                               (unsigned int) initial_y_size,
                               0, visual->depth, InputOutput, visual->visual,
                               cwa_mask, &cwa );

    XStoreName( X_get_display(), window->window_id, title );

    if( setting_position )
    {
        normal_hints.flags = USPosition;
        normal_hints.x = initial_x_pos;
        normal_hints.y = initial_y_pos;
        XSetWMProperties( X_get_display(), window->window_id,
                          NULL, NULL, NULL, 0,
                          &normal_hints, NULL, NULL );
    }

    event_mask = KeyPressMask |
                 KeyReleaseMask |
                 PointerMotionMask |
                 ButtonPressMask |
                 ButtonReleaseMask |
                 ExposureMask |
                 StructureNotifyMask |
                 FocusChangeMask;
    XSelectInput( X_get_display(), window->window_id, event_mask );

    bind_special_keys();

    XMapWindow( X_get_display(), window->window_id );

    XIfEvent( X_get_display(), &event, wait_for_window,
              (char *) (&window->window_id) );

    window->visual = visual;
    window->colour_map_mode = colour_map_mode;
    if( colour_map_mode )
    {
        cmap = XCreateColormap( X_get_display(), window->window_id,
                                visual->visual, AllocAll );
        XSetWindowColormap( X_get_display(), window->window_id, cmap );

        window->colour_map = cmap;
    }

    return( VIO_OK );
}

  Status  X_create_overlay_window(
    X_window_struct  *window,
    XVisualInfo      *visual,
    X_window_struct  *overlay_window )
{
    int                      x_size, y_size, x_pos, y_pos;
    XSetWindowAttributes     cwa;
    Colormap                 cmap;
    XEvent                   event;

    window->visual = visual;
    window->colour_map_mode = TRUE;

    cmap = XCreateColormap( X_get_display(),
                            RootWindow( X_get_display(),
                                        X_get_screen() ),
                            visual->visual, AllocAll );

    X_get_window_geometry( window, &x_pos, &y_pos, &x_size, &y_size );

    cwa.colormap = cmap;
    cwa.background_pixmap = None;
    cwa.border_pixel = 0;

    overlay_window->window_id = XCreateWindow( X_get_display(),
                                      window->window_id,
                                      0, 0,
                                      (unsigned int) x_size,
                                      (unsigned int) y_size,
                                      0, visual->depth, InputOutput,
                                      visual->visual,
                                      CWBackPixmap|CWBorderPixel|CWColormap,
                                      &cwa );
    XMapWindow( X_get_display(), overlay_window->window_id );
    XIfEvent( X_get_display(), &event, wait_for_window,
              (char *) (&overlay_window->window_id) );
    XSetWMColormapWindows( X_get_display(), window->window_id,
                           &overlay_window->window_id, 1 );

    return( VIO_OK );
}

  void  X_delete_window(
    X_window_struct  *window )
{
    XFree( (void *) window->visual );

    if( window->colour_map_mode )
        XFreeColormap( X_get_display(), window->colour_map );

    XDestroyWindow( X_get_display(), window->window_id );
}

  int  X_get_n_colours(
    X_window_struct          *window )
{
    return( window->visual->colormap_size );
}

  void  X_set_colour_map_entry(
    X_window_struct          *window,
    int                      ind,
    int                      r,
    int                      g,
    int                      b )
{
    XColor   def;

    if( !window->colour_map_mode )
    {
        return;
    }

    def.red = (unsigned short) ((unsigned short) r << 8);
    def.green = (unsigned short) ((unsigned short) g << 8);
    def.blue = (unsigned short) ((unsigned short) b << 8);
    def.flags = DoRed | DoGreen | DoBlue;
    def.pixel = (unsigned long) ind;

    XStoreColor( X_get_display(), window->colour_map, &def );
}

  void  X_get_window_geometry(
    X_window_struct          *window,
    int                      *x_position,
    int                      *y_position,
    int                      *x_size,
    int                      *y_size )
{
    int           x, y;
    unsigned int  width, height, border_width, depth;
    int           screen_width, screen_height;
    Window        root;

    if( XGetGeometry( X_get_display(), window->window_id,
                      &root, &x, &y, &width, &height, &border_width, &depth ) )
    {
        X_get_screen_size( &screen_width, &screen_height );
        *x_position = x;
        *y_position = screen_height - 1 - y;
        *x_size = (int) width;
        *y_size = (int) height;
    }
    else
    {
        *x_position = 0;
        *y_position = 0;
        *x_size = 0;
        *y_size = 0;
    }
}

  void  X_get_screen_size(
    int    *x_size,
    int    *y_size )
{
    static  unsigned   int  width, height, border_width, depth;
    static  VIO_BOOL    first = TRUE;
    int                x, y;
    Window             root;

    if( first )
    {
        first = FALSE;

        if( !XGetGeometry( X_get_display(), RootWindow(X_get_display(),
                                                       X_get_screen()),
                           &root, &x, &y, &width, &height, &border_width,
                           &depth ) )
        {
            width = 1000;
            height = 1000;
        }
    }

    *x_size = (int) width;
    *y_size = (int) height;
}

/*--------------------------------- events ----------------------------- */

static  VIO_BOOL  translate_key(
    XEvent   *x_event,
    int      *key )
{
    int             char_count;
    char            buffer[10];
    KeySym          keysym;
    XComposeStatus  compose;

    char_count = XLookupString( &x_event->xkey, buffer, 10, &keysym, &compose );

    *key = (int) ((unsigned char *) buffer) [0];

    return( char_count >= 1 );
}

  VIO_BOOL  X_get_event(
    Event_types          *event_type,
    Window_id            *window,
    event_info_struct    *info )
{
    XEvent   x_event;

    *event_type = NO_EVENT;
    *window = 0;
    
    while( *event_type == NO_EVENT && XPending( X_get_display() ) )
    {
        XNextEvent( X_get_display(), &x_event );
    
        switch( x_event.type )
        {
        case ButtonPress:
            *window = x_event.xbutton.window;
            info->x_mouse = x_event.xbutton.x;
            info->y_mouse = x_event.xbutton.y;
            switch( x_event.xbutton.button )
            {
            case Button1:
                *event_type = LEFT_MOUSE_DOWN_EVENT;  break;
            case Button2:
                *event_type = MIDDLE_MOUSE_DOWN_EVENT;  break;
            case Button3:
                *event_type = RIGHT_MOUSE_DOWN_EVENT;  break;
            }
            break;
    
        case ButtonRelease:
            *window = x_event.xbutton.window;
            info->x_mouse = x_event.xbutton.x;
            info->y_mouse = x_event.xbutton.y;
            switch( x_event.xbutton.button )
            {
            case Button1:
                *event_type = LEFT_MOUSE_UP_EVENT;  break;
            case Button2:
                *event_type = MIDDLE_MOUSE_UP_EVENT;  break;
            case Button3:
                *event_type = RIGHT_MOUSE_UP_EVENT;  break;
            }
            break;
    
        case ConfigureNotify:
            *event_type = WINDOW_RESIZE_EVENT;
            *window = x_event.xconfigure.window;
            info->x_position = x_event.xconfigure.x;
            info->y_position = x_event.xconfigure.y;
            info->x_size = x_event.xconfigure.width;
            info->y_size = x_event.xconfigure.height;
            break;
    
        case DestroyNotify:
            *event_type = WINDOW_QUIT_EVENT;
            *window = x_event.xdestroywindow.window;
            break;
    
        case EnterNotify:
            *event_type = WINDOW_ENTER_EVENT;
            *window = x_event.xcrossing.window;
            break;
    
        case LeaveNotify:
            *event_type = WINDOW_LEAVE_EVENT;
            *window = x_event.xcrossing.window;
            break;
    
        case FocusIn:
            *event_type = WINDOW_ENTER_EVENT;
            *window = x_event.xfocus.window;
            break;
    
        case FocusOut:
            *event_type = WINDOW_LEAVE_EVENT;
            *window = x_event.xfocus.window;
            break;
    
        case Expose:
            *event_type = WINDOW_REDRAW_EVENT;
            *window = x_event.xexpose.window;
            break;
    
        case KeyPress:
            if( translate_key( &x_event, &info->key_pressed ) )
            {
                *event_type = KEY_DOWN_EVENT;
                *window = x_event.xkey.window;
            }
            break;
    
        case KeyRelease:
            if( translate_key( &x_event, &info->key_pressed ) )
            {
                *event_type = KEY_UP_EVENT;
                *window = x_event.xkey.window;
            }
            break;
    
        case MotionNotify:
            *event_type = MOUSE_MOVEMENT_EVENT;
            *window = x_event.xmotion.window;
            info->x_mouse = x_event.xmotion.x;
            info->y_mouse = x_event.xmotion.y;
            break;
    
        case MapNotify:
            *event_type = WINDOW_DEICONIZED_EVENT;
            *window = x_event.xmap.window;
            break;
    
        case UnmapNotify:
            *event_type = WINDOW_ICONIZED_EVENT;
            *window = x_event.xunmap.window;
            break;
    
        default:
            *event_type = NO_EVENT;
            break;
        }
    }

    return( *event_type != NO_EVENT );
}

static  void  bind_key(
    KeySym   keysym,
    int      output_key_value )
{
    unsigned char  string[1];

    string[0] = (unsigned char) output_key_value;

    XRebindKeysym( X_get_display(), keysym, NULL,
                   0, string, 1 );

    /* --- in order to get the up-key event for a modifier,
           we need to bind the key with itself as a modifier */

    if( keysym == XK_Shift_L || keysym == XK_Shift_R ||
        keysym == XK_Control_L || keysym == XK_Control_R ||
        keysym == XK_Alt_L || keysym == XK_Alt_R )
    {
        XRebindKeysym( X_get_display(), keysym, &keysym, 1, string, 1 );
    }
}

static  void  bind_special_keys( void )
{
    int      i;
    static   VIO_BOOL  first = TRUE;
    static  struct  { KeySym keysym;  int  character; }  table[] =
         {
             { XK_Left,       LEFT_ARROW_KEY },
             { XK_Right,      RIGHT_ARROW_KEY },
             { XK_Down,       DOWN_ARROW_KEY },
             { XK_Up,         UP_ARROW_KEY },
             { XK_Shift_L,    LEFT_SHIFT_KEY },
             { XK_Shift_R,    RIGHT_SHIFT_KEY },
             { XK_Control_L,  LEFT_CTRL_KEY },
             { XK_Control_R,  RIGHT_CTRL_KEY },
             { XK_Alt_L,      LEFT_ALT_KEY },
             { XK_Alt_R,      RIGHT_ALT_KEY }
         };

    if( !first )
        return;

    first = FALSE;

    for_less( i, 0, VIO_SIZEOF_STATIC_ARRAY(table) )
        bind_key( table[i].keysym, table[i].character );
}

#define  MAX_NAMES   10
#define  MAX_FONT_SIZE_ERROR   5

static  VIO_BOOL  find_font(
    Font_types       type,
    int              size,
    STRING           *font_name )
{
    int      n_returned;
    VIO_BOOL  found;
    char     pattern[EXTREMELY_LARGE_STRING_SIZE];
    char     **names;
    int      dpi, family, width, slant, weight;
#ifdef NOT_NEEDED
    static   char  *dpis[] = { "100", "75" };
    static   char  *families[] = { "helvetica", "times", "courier" };
    static   char  *widths[] = { "normal", "narrow" };
    static   char  *slants[] = { "r", "o", "i" };
    static   char  *weights[] = { "medium", "bold" };
#endif
    static   char  *dpis[] = { "100", "75" };
    static   char  *families[] = { "helvetica" };
    static   char  *widths[] = { "normal" };
    static   char  *slants[] = { "r" };
    static   char  *weights[] = { "medium" };

    found = FALSE;

    for_less( dpi, 0, VIO_SIZEOF_STATIC_ARRAY(dpis) )
    {
        for_less( family, 0, VIO_SIZEOF_STATIC_ARRAY(families) )
        {
            for_less( width, 0, VIO_SIZEOF_STATIC_ARRAY(widths) )
            {
                for_less( slant, 0, VIO_SIZEOF_STATIC_ARRAY(slants) )
                {
                    for_less( weight, 0, VIO_SIZEOF_STATIC_ARRAY(weights) )
                    {
                        if( type == FIXED_FONT )
                        {
                            (void) strcpy( pattern, "fixed" );
                        }
                        else
                        {
                            (void) sprintf( pattern,
                              "*-%s-%s-%s-%s--\?\?-%d-%s-*",
                              families[family],
                              weights[weight],
                              slants[slant],
                              widths[width],
                              size * 10,
                              dpis[dpi] );
                        }

                        names = XListFonts( X_get_display(), pattern,
                                            MAX_NAMES, &n_returned );

                        found = (n_returned > 0 );

                        if( found || type == FIXED_FONT )
                            break;
                    }
                    if( found || type == FIXED_FONT )
                        break;
                }

                if( found || type == FIXED_FONT )
                    break;
            }
            if( found || type == FIXED_FONT )
                break;
        }

        if( found || type == FIXED_FONT )
            break;
    }


    if( found )
    {
        *font_name = create_string( names[0] );
        XFreeFontNames( names );
    }

    return( found );
}

  VIO_BOOL  X_get_font_name(
    Font_types       type,
    int              size,
    STRING           *font_name )
{
    VIO_BOOL  found;
    int      offset;

    offset = 0;

    do
    {
        found = find_font( type, size + offset, font_name );

        if( offset <= 0 )
            offset = -offset + 1;
        else
            offset = -offset;
    }
    while( !found &&
           (type != FIXED_FONT || offset <= MAX_FONT_SIZE_ERROR) );

    return( found );
}

  void  X_set_mouse_position(
    int                      x,
    int                      y )
{
    int    x_size, y_size;

    X_get_screen_size( &x_size, &y_size );

    XWarpPointer( X_get_display(), None, RootWindow(X_get_display(),
                                                    X_get_screen()),
                  0, 0, 0, 0, x, y_size-1-y );
}
