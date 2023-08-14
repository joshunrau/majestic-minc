#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H  

#include  <volume_io.h>
#include  <WS_windows.h>

#if 0
#define  USE_STORED_FONT_ONLY
#endif

static  VIO_BOOL  use_stored_font_only( void )
{
    static  VIO_BOOL  first = TRUE;
#ifdef USE_STORED_FONT_ONLY
    static  VIO_BOOL  use_stored_flag = TRUE;
#else
    static  VIO_BOOL  use_stored_flag = FALSE;
#endif

    if( first )
    {
        first = FALSE;
        if( getenv("USE_STORED_FONT_ONLY") != NULL )
        {
            if( equal_strings( getenv("USE_STORED_FONT_ONLY"), "no" ) )
                use_stored_flag = FALSE;
            else
                use_stored_flag = TRUE;
        }
    }

    return( use_stored_flag );
}

static  VIO_BOOL  GLX_supported( void )
{
    static  VIO_BOOL  first = TRUE;
    static  VIO_BOOL  supported = FALSE;
    int              error, event;

    if( first )
    {
        first = FALSE;
        supported = glXQueryExtension( X_get_display(), &error, &event );
    }

    return( supported );
}

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
    WS_window_struct       *window )
{
    int           attrib[100], n_attrib, flag, min_size, max_size, mid;
    int           colour_req_index;
    VIO_Status        status;
    Colormap      cmap;
    XVisualInfo   *visual;

    if( !GLX_supported() )
    {
        return( VIO_ERROR );
    }

    n_attrib = 0;

    attrib[n_attrib++] = GLX_LEVEL;
    attrib[n_attrib++] = 0;

    if( colour_map_mode )
    {
        attrib[n_attrib++] = GLX_BUFFER_SIZE;
        min_size = 0;
        max_size = 4096;
        colour_req_index = n_attrib;
        attrib[n_attrib++] = max_size;
    }
    else
    {
        attrib[n_attrib++] = GLX_RGBA;
        attrib[n_attrib++] = GLX_RED_SIZE;
        attrib[n_attrib++] = 1;
        attrib[n_attrib++] = GLX_GREEN_SIZE;
        attrib[n_attrib++] = 1;
        attrib[n_attrib++] = GLX_BLUE_SIZE;
        attrib[n_attrib++] = 1;
    }

    if( double_buffer_flag )
    {
        attrib[n_attrib++] = GLX_DOUBLEBUFFER;
    }

    if( depth_buffer_flag )
    {
        attrib[n_attrib++] = GLX_DEPTH_SIZE;
        attrib[n_attrib++] = 1;
    }

    attrib[n_attrib++] = None;

    /*--- if in colour map mode, do a binary search to find the largest colour
          map */

    if( colour_map_mode )
    {
        while( min_size < max_size )
        {
            mid = (min_size + max_size + 1) / 2;

            attrib[colour_req_index] = mid;

            visual = glXChooseVisual( X_get_display(), X_get_screen(), attrib );

            if( visual == NULL )
                max_size = mid-1;
            else
            {
                min_size = mid;
                XFree( (void *) visual );
            }
        }

        attrib[colour_req_index] = max_size;
    }

    visual = glXChooseVisual( X_get_display(), X_get_screen(), attrib );

    if( visual == NULL )
    {
        print_error( "Cannot find matching visual.\n" );
        return( VIO_ERROR );
    }

    glXGetConfig( X_get_display(), visual, GLX_RGBA, &flag );
    *actual_colour_map_mode = (flag == GL_FALSE);

    glXGetConfig( X_get_display(), visual, GLX_DOUBLEBUFFER, &flag );
    *actual_double_buffer_flag = (flag == GL_TRUE);

    glXGetConfig( X_get_display(), visual, GLX_DEPTH_SIZE, &flag );
    *actual_depth_buffer_flag = (flag > 0);

/*
    cmap = DefaultColormap( X_get_display(), X_get_screen() );
*/
    cmap = XCreateColormap( X_get_display(), RootWindow(X_get_display(),
                                                        X_get_screen()),
                            visual->visual, AllocNone );


    status = X_create_window_with_visual( title,
                              initial_x_pos, initial_y_pos,
                              initial_x_size, initial_y_size, 
                              colour_map_mode, visual, cmap,
                              &window->x_window );

    if( status == VIO_OK )
    {
        window->graphics_context = glXCreateContext( X_get_display(),
                                                     visual, NULL, TRUE );
    }

    window->overlay_present = FALSE;

    *actual_n_overlay_planes = 0;

    if( n_overlay_planes > 0 && n_overlay_planes <= WS_get_n_overlay_planes() )
    {
        n_attrib = 0;
        attrib[n_attrib++] = GLX_LEVEL;
        attrib[n_attrib++] = 1;
        attrib[n_attrib++] = GLX_BUFFER_SIZE;
        attrib[n_attrib++] = 1 << n_overlay_planes;
        attrib[n_attrib++] = None;

        visual = glXChooseVisual( X_get_display(), X_get_screen(), attrib );

        if( visual != NULL )
        {
            glXGetConfig( X_get_display(), visual, GLX_BUFFER_SIZE, &flag );

            *actual_n_overlay_planes = 0;
            while( flag > 1 )
            {
                ++(*actual_n_overlay_planes);
                flag >>= 1;
            } 

            if( X_create_overlay_window( &window->x_window,
                                         visual, &window->overlay_window ) ==VIO_OK)
            {
                window->overlay_present = TRUE;
                window->overlay_context = glXCreateContext( X_get_display(),
                                                            visual, None, TRUE);
            }
            else
            {
                XFree( (void *) visual );
                visual = NULL;
            }
        }
    }

    window->n_fonts = 0;

    return( status );
}

  void  WS_delete_window(
    WS_window_struct  *window )
{
    if( window->n_fonts > 0 )
    {
        FREE( window->font_list_bases );
    }

    glXDestroyContext( X_get_display(), window->graphics_context );
    X_delete_window( &window->x_window );
    if( window->overlay_present )
    {
        glXDestroyContext( X_get_display(), window->overlay_context );
        X_delete_window( &window->overlay_window );
    }
}

  VIO_BOOL  WS_window_has_overlay_planes(
    WS_window_struct  *window )
{
    return( window->overlay_present );
}

static  void  set_window_normal_planes(
    WS_window_struct  *window )
{
    if( !glXMakeCurrent( X_get_display(), window->x_window.window_id,
                         window->graphics_context ) )
    {
        print_error( "set_window_normal_planes():  could not make current.\n" );
    }
}

static  void  set_window_overlay_planes(
    WS_window_struct  *window )
{
    if( !glXMakeCurrent( X_get_display(), window->overlay_window.window_id,
                         window->overlay_context ) )
    {
        print_error( "set_window_overlay_planes():  could not make current.\n");
    }
}

  void  WS_set_current_window(
    WS_window_struct  *window,
    Bitplane_types    current_bitplanes )
{
    if( current_bitplanes == OVERLAY_PLANES )
        set_window_overlay_planes( window );
    else
        set_window_normal_planes( window );
}


  int    WS_get_n_overlay_planes( void )
{
    return( 0 );
#ifdef NOT_YET
    int           n_attrib, attrib[10], min_size, mid, max_size, ind, n_bits;
    XVisualInfo   *visual;


    n_attrib = 0;

    attrib[n_attrib++] = GLX_LEVEL;
    attrib[n_attrib++] = 1;
    attrib[n_attrib++] = GLX_BUFFER_SIZE;
    ind = n_attrib;
    attrib[n_attrib++] = 1;
    attrib[n_attrib++] = None;

    min_size = 0;
    max_size = 4096;

    while( min_size < max_size )
    {
        mid = (min_size + max_size + 1) / 2;
        attrib[ind] = mid;

        visual = glXChooseVisual( X_get_display(), X_get_screen(), attrib );

        if( visual == NULL )
            max_size = mid-1;
        else
        {
            XFree( (void *) visual );
            min_size = mid;
        }
    }

    n_bits = 0;

    while( max_size > 1 )
    {
        ++n_bits;
        max_size >>= 1;
    }

    return( n_bits );
#endif
}

  VIO_BOOL  WS_get_event(
    Event_types          *event_type,
    Window_id            *window,
    event_info_struct    *info )
{
    return( X_get_event( event_type, window, info ) );
}

  void  WS_get_window_position(
    WS_window_struct  *window,
    int               *x_pos,
    int               *y_pos )
{
    int    x_size, y_size;

    X_get_window_geometry( &window->x_window, x_pos, y_pos, &x_size, &y_size );
}

  void  WS_get_window_size(
    WS_window_struct  *window,
    int               *x_size,
    int               *y_size )
{
    int    x_pos, y_pos;

    X_get_window_geometry( &window->x_window, &x_pos, &y_pos, x_size, y_size );
}

  void  WS_set_colour_map_entry(
    WS_window_struct  *window,
    int               ind,
    VIO_Colour            colour )
{
    X_set_colour_map_entry( &window->x_window, ind,
                            get_Colour_r(colour),
                            get_Colour_g(colour),
                            get_Colour_b(colour) );
}

  void  WS_set_overlay_colour_map_entry(
    WS_window_struct  *window,
    int               ind,
    VIO_Colour            colour )
{
    X_set_colour_map_entry( &window->overlay_window, ind,
                            get_Colour_r(colour),
                            get_Colour_g(colour),
                            get_Colour_b(colour) );
}

  void  WS_swap_buffers(
    WS_window_struct  *window )
{
    glXSwapBuffers( X_get_display(), window->x_window.window_id );
}

  VIO_BOOL  WS_get_font(
    Font_types       type,
    VIO_Real             size,
    WS_font_info     *font_info )
{
    if( use_stored_font_only() )
    {
        if( type == FIXED_FONT )
        {
            (void) strcpy( font_info->font_name, "stored_fixed_font" );
            font_info->x_font_info = 0;
            return( TRUE );
        }
        else
            return( FALSE );
    }
    else
    {
        Font         x_font;

        if( X_get_font_name( type, (int) size, &font_info->font_name ) )
        {
            x_font = XLoadFont( X_get_display(), font_info->font_name );
            font_info->x_font_info = XQueryFont( X_get_display(), x_font );

            return( TRUE );
        }

        return( FALSE );
    }
}

/* ARGSUSED */

  void  WS_build_font_in_window(
    WS_window_struct  *window,
    int               font_index,
    WS_font_info      *font_info )
{
    int      i, first, last, listBase;
    Font     x_font;

    if( use_stored_font_only() )
    {
        first = 0;
        last = get_fixed_font_n_chars() - 1;
    }
    else
    {
        first = (int) font_info->x_font_info->min_char_or_byte2;
        last = (int) font_info->x_font_info->max_char_or_byte2;
    }

    listBase = (int) glGenLists( last + 1 );

    if( listBase == 0 )
    {
        print_error( "WS_build_font_in_window(): out of display lists.\n" );
    }
    else
    {
        if( use_stored_font_only() )
            create_fixed_font( (GLuint) listBase );
        else
        {
            x_font = XLoadFont( X_get_display(), font_info->font_name );

            glXUseXFont( x_font, first, last-first+1, listBase+first );
        }
    }

    if( font_index >= window->n_fonts )
    {
        SET_ARRAY_SIZE( window->font_list_bases, window->n_fonts, font_index+1,
                        DEFAULT_CHUNK_SIZE );
        for_less( i, window->n_fonts, font_index )
            window->font_list_bases[i] = -1;
        window->n_fonts = font_index + 1;
    }

    window->font_list_bases[font_index] = listBase;
}

/* ARGSUSED */

  void  WS_delete_font_in_window(
    WS_window_struct     *window,
    int                  font_index,
    WS_font_info         *font_info )
{
    int   last, listBase;

    if( use_stored_font_only() )
        last = get_fixed_font_n_chars() - 1;
    else
        last = (int) font_info->x_font_info->max_char_or_byte2;

    listBase = window->font_list_bases[font_index];

    if( listBase > 0 )
        glDeleteLists( (GLuint) listBase, last + 1 );

    window->font_list_bases[font_index] = -1;
}

/* ARGSUSED */

  VIO_BOOL  WS_set_font(
    WS_window_struct     *window,
    int                  font_index )
{
    if( window->font_list_bases[font_index] > 0 )
        glListBase( (GLuint) window->font_list_bases[font_index] );

    return( window->font_list_bases[font_index] > 0 );
}

  void  WS_delete_font(
    WS_font_info  *info )
{
    if( !use_stored_font_only() )
    {
        XFreeFont( X_get_display(), info->x_font_info );
        delete_string( info->font_name );
    }
}

  VIO_Real  WS_get_character_height(
    WS_font_info  *font_info )
{
    if( use_stored_font_only() )
        return( get_fixed_font_height() );
    else
        return( (VIO_Real) font_info->x_font_info->ascent );
}

  VIO_Real  WS_get_text_length(
    WS_font_info     *font_info,
    STRING           str )
{
    int    i, len, min_char, max_char;

    if( use_stored_font_only() )
    {
        min_char = 0;
        max_char = get_fixed_font_n_chars()-1;
    }
    else
    {
        min_char = (int) font_info->x_font_info->min_char_or_byte2;
        max_char = (int) font_info->x_font_info->max_char_or_byte2;
    }

    len = 0;
    for_less( i, 0, (int) strlen( str ) )
    {
        if( min_char <= (int) str[i] && (int) str[i] <= max_char )
        {
            if( use_stored_font_only() )
                len += (int) get_fixed_font_width( str[i] );
            else
            {
                len += font_info->x_font_info->
                         per_char[(int) str[i] - min_char].width;
            }
        }
    }

    return( (VIO_Real) len );
}

  void  WS_get_screen_size(
    int   *x_size, 
    int   *y_size  )
{
    X_get_screen_size( x_size, y_size );
}
