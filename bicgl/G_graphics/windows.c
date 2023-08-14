#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H

#include  <volume_io.h>
#include  <graphics.h>

/*--------------- maintaining list of known windows ------------------------ */

static  Gwindow        *windows = NULL;
static  int            n_windows = 0;
static  Gwindow        current_window = NULL;


static  Gwindow  create_window_struct( void )
{
    Gwindow     new_window;

    ALLOC( new_window, 1 );

    ADD_ELEMENT_TO_ARRAY( windows, n_windows, new_window, DEFAULT_CHUNK_SIZE );

    return( new_window );
}

static  int  lookup_window_index( Gwindow  window )
{
    int  i;

    for_less( i, 0, n_windows )
    {
        if( windows[i] == window )
            break;
    }

    if( i >= n_windows )
        i = -1;

    return( i );
}

  int  get_n_graphics_windows( void )
{
    return( n_windows );
}

  Gwindow  get_nth_graphics_window(
    int   i )
{
    return( windows[i] );
}

static  void  delete_window_struct( Gwindow   window )
{
    int            ind;

    ind = lookup_window_index( window );

    if( ind >= 0 )
    {
        DELETE_ELEMENT_FROM_ARRAY( windows, n_windows, ind,
                                   DEFAULT_CHUNK_SIZE );
    }

    FREE( window );
}

  void  set_current_window( Gwindow   window )
{
    if( window != current_window )
    {
        if( window != (Gwindow) NULL )
        {
            GS_set_current_window( window->GS_window );
            GS_set_bitplanes( window->GS_window, window->current_bitplanes );
        }

        current_window = window;
    }
}

  Gwindow  get_current_window( void )
{
    return( current_window );
}

/*---------------------- window configurations ------------------------ */

static  void  check_graphics_initialized( void )
{
    static  VIO_BOOL  first = TRUE;

    if( first )
    {
        first = FALSE;
        GS_initialize();
        (void) make_rgba_Colour( 0, 0, 0, 0 );
    }
}

#define  DEFAULT_N_CURVE_SEGMENTS  8

static  void  reinitialize_window(
    Gwindow   window )
{
    int            n_segments;
    Shading_types  type;
    Transparency_modes transparency;
    static  VIO_Surfprop  spr = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };

    set_current_window( window );

    GS_set_matrix_mode( VIEWING_MATRIX );

    type = window->shading_type;
    if ( window->shading_type == FLAT_SHADING )
	window->shading_type = GOURAUD_SHADING;
    else
	window->shading_type = FLAT_SHADING;
    G_set_shading_type( window, type );

    window->lighting_state = !window->lighting_state;
    G_set_lighting_state( window, !window->lighting_state );

    transparency = window->transparency_state;
    if ( window->transparency_state == NO_TRANSPARENCY )
      window->transparency_state = NORMAL_TRANSPARENCY;
    else
      window->transparency_state = NO_TRANSPARENCY;
    G_set_transparency_state( window, transparency );

    window->backface_culling_state = !window->backface_culling_state;
    G_backface_culling_state( window, !window->backface_culling_state );

    n_segments = window->n_curve_segments;
    ++window->n_curve_segments;
    G_set_n_curve_segments( window, n_segments );

    window->marker_labels_visibility = !window->marker_labels_visibility;
    G_set_markers_labels_visibility( window, !window->marker_labels_visibility);

    GS_set_depth_function( LESS_OR_EQUAL );

    GS_set_depth_buffer_state( window->zbuffer_state );

    GS_set_surface_property( window->GS_window, WHITE, &spr );

    update_transforms( window );
}

static  void  initialize_window(
    Gwindow   window )
{
    VIO_Colour     default_background_col;

    default_background_col = make_Colour( 0, 0, 0 );

    window->current_bitplanes = NORMAL_PLANES;
    set_current_window( window );

#ifdef DEBUG
    print_info();
#endif

    GS_get_window_size( &window->x_size, &window->y_size );
    GS_get_window_position( &window->x_origin, &window->y_origin );

    G_set_automatic_clear_state( window, TRUE );

    G_set_background_colour( window, default_background_col );

    window->current_bitplanes = NORMAL_PLANES;
    window->bitplanes_cleared[NORMAL_PLANES] = FALSE;
    window->bitplanes_cleared[OVERLAY_PLANES] = FALSE;

    window->shading_type = GOURAUD_SHADING;
    window->lighting_state = FALSE;
    window->transparency_state = NORMAL_TRANSPARENCY;
    window->backface_culling_state = FALSE;
    window->n_curve_segments = DEFAULT_N_CURVE_SEGMENTS;
    window->marker_labels_visibility = TRUE;

    window->shaded_mode_state = FILLED;
    G_set_shaded_state( window, FILLED );

    G_set_render_lines_as_curves_state( window, FALSE );

    GS_initialize_surface_property( window->GS_window );

    initialize_display_interrupts( window );

    initialize_window_view( window );

    reinitialize_window( window );

    initialize_graphics_lights( window );

    initialize_callbacks_for_window( window );
}

VIO_Status  G_create_child_window(
    VIO_STR        title,
    int            x_pos,
    int            y_pos,
    int            width,
    int            height,
    VIO_BOOL       colour_map_desired,
    VIO_BOOL       double_buffer_desired,
    VIO_BOOL       depth_buffer_desired,
    int            n_overlay_planes_desired,
    Gwindow        *window,
    Gwindow        parent)
{
    VIO_Status   status;
    VIO_BOOL  actual_colour_map_flag;
    VIO_BOOL  actual_double_buffer_flag;
    VIO_BOOL  actual_depth_buffer_flag;
    int      actual_n_overlay_planes;

    check_graphics_initialized();

    *window = create_window_struct();

    ALLOC( (*window)->GS_window, 1 );

    status = GS_create_window( (*window)->GS_window,
                               title, x_pos, y_pos, width, height,
                               colour_map_desired,
                               double_buffer_desired,
                               depth_buffer_desired,
                               MIN( n_overlay_planes_desired,
                                    G_get_n_overlay_planes() ),
                               &actual_colour_map_flag,
                               &actual_double_buffer_flag,
                               &actual_depth_buffer_flag,
                               &actual_n_overlay_planes,
                               parent != NULL ? parent->GS_window : NULL );

    if( status != VIO_OK )
    {
        print_error( "G_create_child_window():  cannot open a graphics window.\n" );
    }
    else
    {
        (*window)->double_buffer_available = TRUE;
        (*window)->double_buffer_state = actual_double_buffer_flag;

        (*window)->colour_map_state = actual_colour_map_flag;
        (*window)->n_overlay_planes = actual_n_overlay_planes;
        (*window)->zbuffer_state = actual_depth_buffer_flag;
        (*window)->zbuffer_available = actual_depth_buffer_flag;

        initialize_window( *window );
    }

    return( status );
}

VIO_Status G_create_window(
    VIO_STR         title,
    int            x_pos,
    int            y_pos,
    int            width,
    int            height,
    VIO_BOOL        colour_map_desired,
    VIO_BOOL        double_buffer_desired,
    VIO_BOOL        depth_buffer_desired,
    int            n_overlay_planes_desired,
    Gwindow        *window )
{
    return G_create_child_window(title,
                                 x_pos, y_pos,
                                 width, height,
                                 colour_map_desired,
                                 double_buffer_desired,
                                 depth_buffer_desired,
                                 n_overlay_planes_desired,
                                 window, NULL);
}

  VIO_Status  G_delete_window(
    Gwindow   window )
{
    VIO_Status    status;

    set_current_window( window );

    if( window->n_overlay_planes > 0 )
    {
        GS_set_bitplanes( window->GS_window, OVERLAY_PLANES );
        G_clear_window( window );
    }

    status = GS_delete_window( window->GS_window );

    FREE( window->GS_window );

    set_current_window( (Gwindow) NULL );

    delete_window_struct( window );

    return( status );
}

  void  G_set_window_title(
    Gwindow   window,
    VIO_STR    title )
{
    set_current_window( window );

    GS_set_window_title( window->GS_window, title );
}

  void  G_terminate( void )
{
    while( n_windows > 0 )
    {
        (void) G_delete_window( windows[0] );
    }

    GS_terminate();
}

  VIO_BOOL  G_is_double_buffer_supported( void )
{
    static  VIO_BOOL   first = TRUE;
    static  VIO_BOOL   available;

    if( first )
    {
        available = GS_is_double_buffer_supported();

        first = FALSE;
    }

    return( available );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_get_double_buffer_state
@INPUT      : window
@OUTPUT     :
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Returns true if the window is in double buffered mode.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  VIO_BOOL  G_get_double_buffer_state(
    Gwindow        window )
{
    return( window->double_buffer_state );
}

  void  G_set_double_buffer_state(
    Gwindow        window,
    VIO_BOOL        flag )
{
    Bitplane_types  save_bitplane;

    if( flag && (!G_can_switch_double_buffering() ||
                 !window->double_buffer_available) )
        flag = FALSE;

    if( flag != window->double_buffer_state )
    {
        set_current_window( window );

        save_bitplane = G_get_bitplanes( window );

        G_set_bitplanes( window, NORMAL_PLANES );

        window->double_buffer_state = GS_set_double_buffer_state(
                                                    window->GS_window, flag );

        reinitialize_window( window );

        restore_bitplanes( window, save_bitplane );
    }
}

  VIO_BOOL  G_is_depth_buffer_supported( void )
{
    static  VIO_BOOL   first = TRUE;
    static  VIO_BOOL   available;

    if( first )
    {
        available = GS_is_depth_buffer_supported();

        first = FALSE;
    }

    return( available );
}

  void  G_set_zbuffer_state(
    Gwindow         window,
    VIO_BOOL         flag )
{
    if( flag && (!G_is_depth_buffer_supported() || !window->zbuffer_available) )
        flag = FALSE;

    if( flag != window->zbuffer_state )
    {
        set_current_window( window );

        GS_set_depth_buffer_state( flag );

        window->zbuffer_state = flag;
    }
}

  VIO_BOOL  G_get_zbuffer_state(
    Gwindow         window )
{
    return( window->zbuffer_state );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_get_colour_map_state
@INPUT      : window
@OUTPUT     :
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Returns true if window is in colour map mode.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  VIO_BOOL  G_get_colour_map_state(
    Gwindow        window )
{
    return( window->colour_map_state );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_set_colour_map_state
@INPUT      : window
              flat
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Sets the colour map state of the window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_set_colour_map_state(
    Gwindow        window,
    VIO_BOOL        flag )
{
    Bitplane_types  save_bitplane;

    if( G_can_switch_colour_map_mode() && flag != window->colour_map_state )
    {
        set_current_window( window );

        save_bitplane = G_get_bitplanes( window );

        G_set_bitplanes( window, NORMAL_PLANES );

        window->colour_map_state = GS_set_colour_map_state( window->GS_window,
                                                            flag );

        reinitialize_window( window );

        restore_bitplanes( window, save_bitplane );
    }
}

  int  G_get_n_colour_map_entries(
    Gwindow        window )
{
    set_current_window( window );

    return( GS_get_n_colour_map_entries( window->GS_window,
                                         G_get_double_buffer_state(window) ) );
}

  void  G_set_colour_map_entry(
    Gwindow         window,
    int             ind,
    VIO_Colour          colour )
{
    Bitplane_types  save_bitplane;

    set_current_window( window );

    save_bitplane = G_get_bitplanes( window );

    if( save_bitplane != NORMAL_PLANES )
        G_set_bitplanes( window, NORMAL_PLANES );

    GS_set_colour_map_entry( window->GS_window,
                             window->current_bitplanes, ind, colour );

    if( save_bitplane != NORMAL_PLANES )
        restore_bitplanes( window, save_bitplane );
}

  void  restore_bitplanes(
    Gwindow         window,
    Bitplane_types  bitplane )
{
    if( window == (Gwindow) NULL )
        window = current_window;

    if( window != (Gwindow) NULL )
        G_set_bitplanes( window, bitplane );
}

  int  G_get_monitor_width( void )
{
    return( GS_get_monitor_width() );
}

  int  G_get_monitor_height( void )
{
    return( GS_get_monitor_height() );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_get_monitor_widths_to_eye
@INPUT      :
@OUTPUT     :
@RETURNS    : ratio of monitor width to the eye distance
@DESCRIPTION: Returns the ratio of the monitor width to the distance to the
              users eye.  Used to help define reasonable 3D views.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  VIO_Real  G_get_monitor_widths_to_eye( void )
{
    return( 1.2 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_get_window_size
@INPUT      : window
@OUTPUT     : x_size
              y_size
@RETURNS    :
@DESCRIPTION: Passes back the size of the window in pixels.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_get_window_size(
    Gwindow        window,
    int            *x_size,
    int            *y_size )
{
    *x_size = window->x_size;
    *y_size = window->y_size;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_get_window_position
@INPUT      : window
@OUTPUT     : x_pos
              y_pos
@RETURNS    :
@DESCRIPTION: Passes back the position of the window in pixels from bottom left.@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

void  G_get_window_position(
    Gwindow        window,
    int            *x_pos,
    int            *y_pos )
{
  set_current_window(window);
  GS_get_window_position(x_pos, y_pos);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_set_background_colour
@INPUT      : window
              colour
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Sets the background colour of the window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_set_background_colour(
    Gwindow         window,
    VIO_Colour          colour )
{
    if( window->colour_map_state )
        window->background_colour_index = colour;
    else
        window->background_colour = colour;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_get_background_colour
@INPUT      : window
@OUTPUT     :
@RETURNS    : colour
@DESCRIPTION: Returns the background colour of the window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : Jan. 24, 1995         David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  VIO_Colour  G_get_background_colour(
    Gwindow         window )
{
    if( window->colour_map_state )
        return( window->background_colour_index );
    else
        return( window->background_colour );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_clear_window
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

  void  G_clear_window(
    Gwindow    window )
{
    VIO_Colour    colour;

    set_current_window( window );

    if( G_get_colour_map_state(window) )
        colour = window->background_colour_index;
    else
        colour = window->background_colour;

    GS_clear_window( window->GS_window,
                     window->x_size,
                     window->y_size,
                     window->current_bitplanes,
                     G_get_colour_map_state(window),
                     G_get_zbuffer_state(window), colour );

    window->bitplanes_cleared[window->current_bitplanes] = TRUE;
}

  void  G_clear_viewport(
    Gwindow    window,
    VIO_Colour     colour )
{
    set_current_window( window );

    GS_clear_viewport( window->GS_window,
                       window->x_viewport_min,
                       window->x_viewport_max,
                       window->y_viewport_min,
                       window->y_viewport_max,
                       window->x_size,
                       window->y_size,
                       window->current_bitplanes,
                       G_get_colour_map_state(window),
                       G_get_zbuffer_state(window), colour );
}

  void  G_set_automatic_clear_state(
    Gwindow    window,
    VIO_BOOL    state )
{
    window->automatic_clear_flag = state;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_update_window
@INPUT      : window
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Called at the end of drawing a frame, this basically does just
              a swap front and back buffers.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_update_window( Gwindow   window )
{
    set_current_window( window );

    check_window_cleared( window );

    if( window->current_bitplanes == NORMAL_PLANES )
    {
        if( window->double_buffer_state )
            GS_swap_buffers();
        GS_flush();
    }
    else
        GS_flush();

    set_continuation_flag( window, FALSE );

    window->bitplanes_cleared[window->current_bitplanes] = FALSE;
}

/* ARGSUSED */

  void  check_window_cleared(
     Gwindow    window  )
{
    if( window->automatic_clear_flag &&
        !window->bitplanes_cleared[window->current_bitplanes] )
    {
        G_clear_window( window );
    }
}

  VIO_BOOL  G_has_overlay_planes( void )
{
    return( G_get_n_overlay_planes() > 0 );
}

  VIO_BOOL  G_window_has_overlay_planes(
    Gwindow   window )
{
    return( window->n_overlay_planes > 0 );
}

static  void  update_blend_function(
    Gwindow         window,
    Bitplane_types  bitplane )
{
    set_current_window( window );

    if( bitplane == OVERLAY_PLANES || !window->transparency_state )
        GS_set_blend_function(0);
    else
        GS_set_blend_function( window->transparency_state );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_set_current_window
@INPUT      : window
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Sets the current window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 2015            Robert D. Vincent
@MODIFIED   :
---------------------------------------------------------------------------- */

void G_set_current_window( Gwindow window )
{
    GS_set_current_window( window->GS_window );
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_set_bitplanes
@INPUT      : window
              bitplanes
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Sets the current bitplanes to OVERLAY_PLANES, NORMAL_PLANES, etc.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_set_bitplanes(
    Gwindow          window,
    Bitplane_types   bitplanes )
{
    if( window->current_bitplanes != bitplanes &&
        window->n_overlay_planes > 0 )
    {
        GS_set_bitplanes( window->GS_window, bitplanes );
        window->current_bitplanes = bitplanes;
        update_blend_function( window, bitplanes );
    }
}

  Bitplane_types  G_get_bitplanes(
    Gwindow          window )
{
    return( window->current_bitplanes );
}

  VIO_BOOL  G_can_switch_double_buffering( void )
{
    return( GS_can_switch_double_buffering() );
}

  VIO_BOOL  G_can_switch_colour_map_mode( void )
{
    return( GS_can_switch_colour_map_mode() );
}

  VIO_BOOL  G_has_transparency_mode( void )
{
    return( GS_has_transparency_mode() );
}

  VIO_BOOL  G_has_rgb_mode( void )
{
    return( GS_has_rgb_mode() );
}

  int  G_get_n_overlay_planes( void )
{
    return( GS_get_num_overlay_planes() );
}

  void  G_set_overlay_colour_map(
    Gwindow         window,
    int             ind,
    VIO_Colour          colour )
{
    Bitplane_types   save_bitplane;

    if( window->n_overlay_planes > 0 )
    {
        set_current_window( window );

        save_bitplane = G_get_bitplanes( window );

        GS_set_bitplanes( window->GS_window, OVERLAY_PLANES );

        GS_set_overlay_colour_map( window->GS_window, ind, colour );

        restore_bitplanes( window, save_bitplane );
    }
}

  void  G_append_to_last_update(
     Gwindow   window )
{
    set_current_window( window );
    GS_append_to_last_update( window->GS_window, G_get_zbuffer_state(window),
                              window->x_size, window->y_size );
    set_continuation_flag( window, TRUE );
    window->bitplanes_cleared[NORMAL_PLANES] = TRUE;
}

  void  G_continue_last_update(
     Gwindow   window )
{
    set_current_window( window );
    set_continuation_flag( window, TRUE );
}

  void  G_set_transparency_state(
    Gwindow            window,
    Transparency_modes state )
{
    if( state != window->transparency_state )
    {
        window->transparency_state = state;
        update_blend_function( window, window->current_bitplanes );
    }
}

void G_set_visibility(Gwindow window, VIO_BOOL is_visible)
{
    GS_set_visibility(window->GS_window, is_visible);
}

void G_set_geometry(Gwindow window, int x, int y, int width, int height)
{
    GS_set_geometry(window->GS_window, x, y, width, height);
}
