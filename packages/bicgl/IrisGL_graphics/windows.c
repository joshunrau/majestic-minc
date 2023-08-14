 
#include  <volume_io.h>
#include  <GS_graphics.h>
#include  <gl/get.h>

  Window_id  GS_get_window_id( 
    GSwindow  window )
{
    return( window->window_id );
}

  void  GS_set_current_window(
    GSwindow          window )
{
    winset( window->window_id );
}

#define  CALL_GFLUSH_EVERY_NTH_TIME  200

#define  OVERLAY_CLEAR_INDEX       0
#define  DEFAULT_N_CURVE_SEGMENTS  8

#define  CATMULL_ROM_ID            1

static  void  initialize_window(
    GSwindow       window,
    VIO_BOOL        colour_map_desired,
    VIO_BOOL        double_buffer_desired,
    VIO_BOOL        depth_buffer_desired,
    int            n_overlay_planes_desired );

static  long     zdepth_clear;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_initialize
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Initializes any graphics type stuff, is called just before first
              graphics window is opened.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  void  GS_initialize( void )
{
    static  Matrix   catmull_rom_matrix = {
                                              { -0.5f,  1.5f, -1.5f,  0.5f },
                                              {  1.0f, -2.5f,  2.0f, -0.5f },
                                              { -0.5f,  0.0f,  0.5f,  0.0f },
                                              {  0.0f,  1.0f,  0.0f,  0.0f },
                                          };

    noport();
    foreground();
    (void) winopen( "" );

    zdepth_clear = getgdesc( GD_ZMAX );

    defbasis( CATMULL_ROM_ID, catmull_rom_matrix );
    curvebasis( CATMULL_ROM_ID );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_terminate
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Terminates the use of graphics, freeing up memory.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  void  GS_terminate( void )
{
    gexit();
}

  VIO_BOOL  GS_can_switch_double_buffering( void )
{
    return( TRUE );
}

  VIO_BOOL  GS_can_switch_colour_map_mode( void )
{
    return( TRUE );
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_create_window
@INPUT      : title
              x_pos
              y_pos
              width
              height
@OUTPUT     : window
@RETURNS    : 
@DESCRIPTION: Creates a window of the specified size and shape.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  Status  GS_create_window(
    GSwindow       window,
    STRING         title,
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
    int            *actual_n_overlay_planes )
{
    Status     status;
    Window_id  save_window_id;

    save_window_id = winget();

    if( width > 0 && height > 0 )
    {
        if( x_pos >= 0 && y_pos >= 0 )
        {
            prefposition( (long) x_pos, (long) (x_pos + width - 1),
                          (long) y_pos, (long) (y_pos + height - 1) );
        }
        else
        {
            prefsize( (long) width, (long) height );
        }
    }

    foreground();

    window->window_id = winopen( title );

    GS_set_current_window( window );

    if( width > 0 && height > 0 )
    {
        winconstraints();
        winconstraints();
    }

    if( GS_get_window_id(window) >= 0 )
    {
        if( !GS_is_double_buffer_supported() )
            double_buffer_desired = FALSE;

        if( !GS_is_depth_buffer_supported() )
            depth_buffer_desired = FALSE;

        *actual_colour_map_flag = colour_map_desired;
        *actual_double_buffer_flag = double_buffer_desired;
        *actual_depth_buffer_flag = depth_buffer_desired;
        *actual_n_overlay_planes = n_overlay_planes_desired;

        initialize_window( window,
                           colour_map_desired,
                           double_buffer_desired,
                           depth_buffer_desired,
                           n_overlay_planes_desired );
        status = VIO_OK;
    }
    else
        status = VIO_ERROR;

    if( save_window_id >= 0 )
        winset( save_window_id );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_window
@INPUT      : window
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Initializes the window, called when the window is created.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

static  void  initialize_window(
    GSwindow       window,
    VIO_BOOL        colour_map_desired,
    VIO_BOOL        double_buffer_desired,
    VIO_BOOL        depth_buffer_desired,
    int            n_overlay_planes_desired )
{
    (void) GS_set_double_buffer_state( window, double_buffer_desired );

    (void) GS_set_colour_map_state( window, colour_map_desired );

    subpixel( TRUE );

    if( n_overlay_planes_desired > 0 )
    {
        overlay( (long) n_overlay_planes_desired );
        gconfig();

        GS_set_overlay_colour_map( window, 1, RED );
        GS_set_overlay_colour_map( window, 2, GREEN );
        GS_set_overlay_colour_map( window, 3, BLUE );
    }

    initialize_window_events( window );
}

/* ARGSUSED */

  void  GS_set_window_title(
    GSwindow   window,
    STRING     title )
{
    wintitle( title );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_has_transparency_mode
@INPUT      : 
@OUTPUT     : 
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Returns true if transparency mode is supported.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  VIO_BOOL  GS_has_transparency_mode( void )
{
    return( getgdesc( GD_BLEND ) != 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_has_rgb_mode
@INPUT      : 
@OUTPUT     : 
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Returns true if rgb mode is supported.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  VIO_BOOL  GS_has_rgb_mode( void )
{
    return( getgdesc( GD_BITS_NORM_SNG_RED ) > 0 &&
            getgdesc( GD_BITS_NORM_SNG_GREEN ) > 0 &&
            getgdesc( GD_BITS_NORM_SNG_BLUE ) > 0 );
}

  VIO_BOOL  GS_set_double_buffer_state(
    GSwindow       window,
    VIO_BOOL        flag )
{
    long   display_mode;

    if( flag )
        doublebuffer();
    else
        singlebuffer();

    gconfig();

    display_mode = getdisplaymode();
    if( display_mode == DMRGBDOUBLE || display_mode == DMDOUBLE )
        flag = TRUE;
    else
        flag = FALSE;

    return( flag );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_set_colour_map_state
@INPUT      : flat
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the colour map state of the window.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  VIO_BOOL  GS_set_colour_map_state(
    GSwindow       window,
    VIO_BOOL        flag )
{
    if( flag )
        cmode();
    else
        RGBmode();

    gconfig();

    return( flag );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_get_n_colour_map_entries
@INPUT      : window
@OUTPUT     : 
@RETURNS    : n colours
@DESCRIPTION: Returns the number of colour map entries for the window, based
              on its current double buffer state.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  int  GS_get_n_colour_map_entries(
    GSwindow  window,
    VIO_BOOL   double_buffer_state )
{
    int   n_bits, n_colours, query;

    if( double_buffer_state )
    {
        query = GD_BITS_NORM_DBL_CMODE;
    }
    else
    {
        query = GD_BITS_NORM_SNG_CMODE;
    }

    n_bits = (int) getgdesc( (long) query );

    n_colours = (1 << n_bits);

    if( n_colours == 4096 )   /* from man page on mapcolor, should not use */
        n_colours -= 256;     /* upper 256 colours on G systems            */

    return( n_colours );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_set_colour_map_entry
@INPUT      : window
              index
              colour
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the colour of the index'th colour map entry.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  void  GS_set_colour_map_entry(
    GSwindow        window,
    Bitplane_types  bitplanes,
    int             index,
    VIO_Colour          colour )
{
    static   int  count = 0;

    mapcolor( (Colorindex) index,
              (short) get_Colour_r(colour),
              (short) get_Colour_g(colour),
              (short) get_Colour_b(colour) );

    if( bitplanes == NORMAL_PLANES )
    {
        ++count;
        if( count == CALL_GFLUSH_EVERY_NTH_TIME )
        {
            gflush();
            count = 0;
        }
    }
}

  VIO_BOOL  GS_is_double_buffer_supported( void )
{
    VIO_BOOL   available;

    if( getgdesc( GD_BITS_NORM_DBL_RED ) <= 0 ||
        getgdesc( GD_BITS_NORM_DBL_GREEN ) <= 0 ||
        getgdesc( GD_BITS_NORM_DBL_BLUE ) <= 0 )
    {
        available = FALSE;
    }
    else
    {
        available = TRUE;
    }

    return( available );
}

  VIO_BOOL  GS_is_transparency_supported( void )
{
    VIO_BOOL   available;

    available = (getgdesc( GD_BITS_NORM_DBL_ALPHA ) > 0);

    return( available );
}

  VIO_BOOL  GS_is_depth_buffer_supported( void )
{
    VIO_BOOL   available;

    available = (getgdesc( GD_BITS_NORM_ZBUFFER ) > 0);

    return( available );
}

  void  GS_set_depth_buffer_state(
    VIO_BOOL         flag )
{
    zbuffer( (Boolean) flag );
}

  void  GS_set_depth_function(
    Depth_functions   depth_func )
{
    int   gl_depth_func;

    switch( depth_func )
    {
    case LESS_OR_EQUAL:
        gl_depth_func = ZF_LEQUAL;
        break;
    default:
        gl_depth_func = ZF_LEQUAL;
        break;
    }

    zfunction( (long) gl_depth_func  );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_delete_window
@INPUT      : window
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the graphics window, and frees up any associated memory.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  Status  GS_delete_window(
    GSwindow   window )
{
    Status    status;

    if( GS_get_window_id(window) >= 0 )
    {
        delete_window_events( window );
        winclose( GS_get_window_id(window) );

        status = VIO_OK;
    }
    else
    {
        print_error( "Error:  tried to delete invalid window.\n" );
        status = VIO_ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_get_monitor_width
@INPUT      : 
@OUTPUT     : 
@RETURNS    : width of the monitor in pixels
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  int  GS_get_monitor_width( void )
{
    long   width;

    width = getgdesc( GD_XPMAX );

    return( (int) width );
}

  int  GS_get_monitor_height( void )
{
    long   height;

    height = getgdesc( GD_YPMAX );

    return( (int) height );
}

  void  GS_clear_depth_buffer( void )
{
    zclear();
}

static   void  clear_viewport(
    GSwindow        window,
    Bitplane_types  current_bitplanes,
    VIO_BOOL         colour_map_state,
    VIO_BOOL         zbuffer_state,
    VIO_Colour          colour )
{
    if( current_bitplanes == NORMAL_PLANES )
    {
        if( colour_map_state )
        {
            GS_set_colour_index( colour );
            clear();
            if( zbuffer_state )
                zclear();
        }
        else
        {
            if( zbuffer_state )
                czclear( colour, zdepth_clear );
            else
            {
                GS_set_colour( colour );
                clear();
            }
        }
    }
    else
        GS_clear_overlay();
}

  void  GS_clear_overlay( void )
{
    GS_set_colour_index( OVERLAY_CLEAR_INDEX );
    clear();
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_clear_window
@INPUT      : window
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Clears the entire window to the background colour.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  void  GS_clear_window(
    GSwindow        window,
    int             x_size,
    int             y_size,
    Bitplane_types  bitplane,
    VIO_BOOL         colour_map_state,
    VIO_BOOL         zbuffer_state,
    VIO_Colour          colour )
{
    pushviewport();
    viewport( (Screencoord) 0, (Screencoord) (x_size-1),
              (Screencoord) 0, (Screencoord) (y_size-1) );
    clear_viewport( window, bitplane, colour_map_state, zbuffer_state, colour );
    popviewport();
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_clear_viewport
@INPUT      : window
              colour
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Clears the current viewport to the given colour.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

  void  GS_clear_viewport(
    GSwindow        window,
    int             x_viewport_min,
    int             x_viewport_max,
    int             y_viewport_min,
    int             y_viewport_max,
    int             x_size,
    int             y_size,
    Bitplane_types  bitplane,
    VIO_BOOL         colour_map_state,
    VIO_BOOL         zbuffer_state,
    VIO_Colour          colour )
{
    clear_viewport( window, bitplane, colour_map_state, zbuffer_state, colour );
}

  void  GS_flush( void )
{
    gflush();
}

/* ARGSUSED */

  void  GS_swap_buffers( void )
{
    swapbuffers();
}

  void  GS_append_to_last_update(
    GSwindow   window,
    VIO_BOOL    zbuffer_state,
    int        x_size,
    int        y_size )
{
    VIO_BOOL   zbuffer_on;

    zbuffer_on = zbuffer_state;

    if( zbuffer_on )
        GS_set_depth_buffer_state( FALSE );

    readsource( SRC_FRONT );

    rectcopy( (Screencoord) 0,
              (Screencoord) 0,
              (Screencoord) (x_size-1),
              (Screencoord) (y_size-1),
              (Screencoord) 0,
              (Screencoord) 0 );

    if( zbuffer_on )
        GS_set_depth_buffer_state( TRUE );
}

  int  GS_get_num_overlay_planes( void )
{
    static  VIO_BOOL  first = TRUE;
    static  int      n_planes;

    if( first )
    {
        first = FALSE;
        n_planes = (int) getgdesc( GD_BITS_OVER_SNG_CMODE );
    }

    return( n_planes );
}

  void  GS_set_bitplanes( 
    GSwindow         window,
    Bitplane_types   bitplanes )
{
    switch( bitplanes )
    {
    case NORMAL_PLANES:
        drawmode( NORMALDRAW );
        zwritemask( 0xffffffff );
        break;

    case OVERLAY_PLANES:
        drawmode( NORMALDRAW );
        zwritemask( 0 );
        drawmode( OVERDRAW );
        break;
    }
}

  void  GS_set_overlay_colour_map(
    GSwindow        window,
    int             index,
    VIO_Colour          colour )
{
    mapcolor( (Colorindex) index, (short) get_Colour_r(colour),
                                  (short) get_Colour_g(colour),
                                  (short) get_Colour_b(colour) );
}

  void  GS_get_window_position(
    int       *x_pos,
    int       *y_pos )
{
    long    lx_pos, ly_pos;

    getorigin( &lx_pos, &ly_pos );
    *x_pos = (int) lx_pos;
    *y_pos = (int) ly_pos;
}

  void  GS_get_window_size(
    int       *x_size,
    int       *y_size )
{
    long    lx_size, ly_size;

    getsize( &lx_size, &ly_size );
    *x_size = (int) lx_size;
    *y_size = (int) ly_size;
}
