#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H

#include  <volume_io.h>
#include  <graphics.h>
#include  <random_order.h>

#include <assert.h>

#define DISALLOW_DRAWING_INTERRUPT 1
#define PARANOID 1

/* #define DEBUG 1 */

#define  MAX_LINE_WIDTH  1000.0f

static  void     draw_marker_as_cube( Gwindow, VIO_Colour );
static  void     draw_marker_as_sphere( Gwindow, VIO_Colour );

#define N_VTX_PER_TRIANGLE 3
/**
 * Check whether a polygons_struct is (most likely) a triangular mesh.
 */
VIO_BOOL
is_triangular(polygons_struct *polygons)
{
  int start_index = 0;
  int end_index = 0;
  int i;

  for (i = 0; i < polygons->n_items; i++)
  {
    start_index = end_index;
    end_index = polygons->end_indices[i];
    if ((end_index - start_index) != N_VTX_PER_TRIANGLE)
    {
      return FALSE;
    }
  }
  return TRUE;
}

void  initialize_display_interrupts( Gwindow        window )
{
    window->interrupt_allowed = FALSE;
    window->continuation_flag = FALSE;
    window->interrupt_occurred = FALSE;
    window->interrupt_interval = 100;
    window->interrupt_time_interval = 0.3;
}

void  G_set_drawing_interrupt_state( Gwindow          window,
                                     VIO_BOOL         state )
{
    window->interrupt_allowed = state;
}

VIO_BOOL  G_get_drawing_interrupt_state( Gwindow          window )
{
    return( window->interrupt_allowed );
}

void  G_set_drawing_interrupt_check_n_objects( Gwindow          window,
                                               int             interval )
{
    window->interrupt_interval = interval;
}

void  G_set_drawing_interrupt_time( Gwindow          window,
                                    VIO_Real             interval )
{
    window->interrupt_time_interval = interval;
}

void  G_start_interrupt_test( Gwindow          window )
{
    VIO_Real   current_time;

    current_time = current_realtime_seconds();
    window->interrupt_time = current_time + window->interrupt_time_interval;
}

VIO_BOOL  G_get_interrupt_occurred( Gwindow          window )
{
    return( window->interrupt_occurred );
}

void  G_clear_drawing_interrupt_flag( Gwindow          window )
{
    window->interrupt_occurred = FALSE;
}

void  set_continuation_flag( Gwindow window, VIO_BOOL state )
{
    window->continuation_flag = state;
}

VIO_BOOL  view_is_stereo( Gwindow        window )
{
    return( window->stereo_flag &&
            !window->colour_map_state &&
            window->current_bitplanes == NORMAL_PLANES &&
            (window->current_view_type == WORLD_VIEW ||
             window->current_view_type == MODEL_VIEW) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : about_to_draw_graphics
@INPUT      : window
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Called just before drawing into a window, to set the current
              window, and clear it, if automatic clearing is enabled.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

static  void  about_to_draw_graphics( Gwindow        window )
{
    set_current_window( window );
    check_window_cleared( window );
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : draw_polygons
@INPUT      : windows
              polygons
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws a set of polygons in the window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

static  void  draw_polygons(
    Gwindow         window,
    polygons_struct *polygons )
{
    int n_iters, i;

    if( !window->shaded_mode_state &&
        polygons->line_thickness > 1.0f &&
        polygons->line_thickness < MAX_LINE_WIDTH )
    {
        GS_set_line_width( (VIO_Real) polygons->line_thickness );
    }

    if( view_is_stereo(window) )
    {
        n_iters = 2;
    }
    else
    {
        n_iters = 1;
    }

    for (i = 0; i < n_iters; i++)
    {
      if (n_iters == 2)
      {
          set_view_for_eye( window, i );
          if (i == 0)
            glColorMask(GL_TRUE, GL_FALSE, GL_FALSE, GL_FALSE);
          else
            glColorMask(GL_FALSE, GL_FALSE, GL_TRUE, GL_FALSE);
      }
      switch( polygons->colour_flag )
      {
      case ONE_COLOUR:
        if (!is_triangular(polygons))
        {
          GS_draw_polygons_one_colour( window->GS_window, window->shaded_mode_state, polygons );
        }
        else
        {
          GS_draw_triangles_one_colour( window->GS_window, window->shaded_mode_state, polygons );
        }
        break;

      case PER_ITEM_COLOURS:
        if (!is_triangular(polygons))
        {
          GS_draw_polygons_one_colour( window->GS_window, window->shaded_mode_state, polygons );
        }
        else
        {
          GS_draw_triangles_per_item_colours( window->GS_window, window->shaded_mode_state, polygons );
        }
        break;

      case  PER_VERTEX_COLOURS:
        GS_draw_polygons_per_vertex_colours( window->GS_window, window->shaded_mode_state, polygons );
        break;
      }
    }

    if( view_is_stereo(window) )
    {
        set_view_for_eye( window, -1 );
        glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    }

    if( !window->shaded_mode_state &&
        polygons->line_thickness > 1.0f &&
        polygons->line_thickness < MAX_LINE_WIDTH )
        GS_set_line_width( 1.0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_draw_polygons
@INPUT      : window
              polygons
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Routine to draw a set of polygons in the window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_draw_polygons(
    Gwindow         window,
    polygons_struct *polygons )
{
    if (polygons->n_points == 0 || polygons->n_items == 0)
    {
        return;                 /* no work to do! */
    }

    about_to_draw_graphics( window );

    draw_polygons( window, polygons );
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_draw_quadmesh
@INPUT      : window
              quadmesh
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws a quadmesh in the window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_draw_quadmesh(
    Gwindow         window,
    quadmesh_struct *quadmesh )
{
    about_to_draw_graphics( window );

    switch( quadmesh->colour_flag )
    {
    case  ONE_COLOUR:
        GS_draw_quadmesh_one_colour( window->GS_window,
                                     window->shaded_mode_state,
                                     quadmesh );
        break;

    case  PER_ITEM_COLOURS:
        GS_draw_quadmesh_per_item_colours( window->GS_window,
                                           window->shaded_mode_state,
                                           quadmesh );
        break;

    case  PER_VERTEX_COLOURS:
        GS_draw_quadmesh_per_vertex_colours( window->GS_window, 
                                             window->shaded_mode_state,
                                             quadmesh );
        break;
    }
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_draw_lines
@INPUT      : window
              lines
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws the lines in the given window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */
void
G_draw_lines(Gwindow window, lines_struct *lines )
{
  int save_lights;

  if (lines->n_points == 0 || lines->n_items == 0)
  {
    /* Nothing to do! */
    return;
  }

  about_to_draw_graphics( window );

  save_lights = G_get_lighting_state( window );
  G_set_lighting_state( window, FALSE );

  if( lines->line_thickness > 1.0 && lines->line_thickness < MAX_LINE_WIDTH )
    GS_set_line_width( lines->line_thickness );

  GS_draw_lines( window->GS_window, lines );

  if( lines->line_thickness > 1.0 && lines->line_thickness < MAX_LINE_WIDTH )
    GS_set_line_width( 1.0 );

  G_set_lighting_state( window, save_lights );
}

/* ------------------------------ text and fonts ------------------- */

  VIO_Real  G_get_text_height(
    Font_types       type,
    VIO_Real             size )
{
    return( GS_get_character_height( type, size ) );
}

  VIO_Real  G_get_text_length(
    VIO_STR           str,
    Font_types       type,
    VIO_Real             size )
{
    return( GS_get_text_length( str, type, size ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_draw_text
@INPUT      : window
              text
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws the text in the given window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_draw_text(
    Gwindow         window,
    text_struct     *text )
{
    about_to_draw_graphics( window );

    GS_set_colour( text->colour );

    GS_set_raster_position( (VIO_Real) Point_x(text->origin),
                            (VIO_Real) Point_y(text->origin),
                            (VIO_Real) Point_z(text->origin) );

    GS_draw_text( text->font, text->size, text->string );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_draw_marker
@INPUT      : window
              marker
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws a 3D marker in the window
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_draw_marker(
    Gwindow         window,
    marker_struct   *marker )
{
    VIO_BOOL     save_lights;
    VIO_Transform   transform;

    about_to_draw_graphics( window );

    save_lights = G_get_lighting_state( window );
    G_set_lighting_state( window, FALSE );

    GS_set_colour( marker->colour );

    if( marker->size <= 0.0 )
    {
      GS_begin_point();
      GS_set_point( &marker->position );
      GS_end_point();
    }
    else
    {
      GS_push_transform();

      make_translation_transform( (VIO_Real) Point_x(marker->position),
                                  (VIO_Real) Point_y(marker->position),
                                  (VIO_Real) Point_z(marker->position),
                                  &transform );
      GS_mult_transform( &transform );

      make_scale_transform( marker->size, marker->size, marker->size,
                            &transform );
      GS_mult_transform( &transform );

      switch( marker->type )
      {
      case  BOX_MARKER:
        draw_marker_as_cube( window, marker->colour );
        break;

      case  SPHERE_MARKER:
        draw_marker_as_sphere( window, marker->colour );
        break;

      default:
        handle_internal_error("G_draw_marker(): unknown marker type.");
        break;
      }

      GS_pop_transform();
    }

    if( window->marker_labels_visibility )
    {
      GS_set_raster_position( (VIO_Real) Point_x(marker->position) +
                              1.5 * (VIO_Real) marker->size / 2.0,
                              (VIO_Real) Point_y(marker->position),
                              (VIO_Real) Point_z(marker->position) );

      GS_draw_text( FIXED_FONT, 10.0, marker->label );
    }

    G_set_lighting_state( window, save_lights );
}

#define N_POINTS_CUBE 24
#define N_FACES_CUBE 6

static void
draw_marker_as_cube( Gwindow window, VIO_Colour colour )
{
    static VIO_Colour colours[1];
    static VIO_Point  points[N_POINTS_CUBE] = {
      { { -0.5f, -0.5f, -0.5f } },
      { { -0.5f, -0.5f,  0.5f } },
      { { -0.5f,  0.5f,  0.5f } },
      { { -0.5f,  0.5f, -0.5f } },

      { {  0.5f, -0.5f, -0.5f } },
      { {  0.5f,  0.5f, -0.5f } },
      { {  0.5f,  0.5f,  0.5f } },
      { {  0.5f, -0.5f,  0.5f } },

      { { -0.5f, -0.5f, -0.5f } },
      { {  0.5f, -0.5f, -0.5f } },
      { {  0.5f, -0.5f,  0.5f } },
      { { -0.5f, -0.5f,  0.5f } },

      { { -0.5f,  0.5f, -0.5f } },
      { { -0.5f,  0.5f,  0.5f } },
      { {  0.5f,  0.5f,  0.5f } },
      { {  0.5f,  0.5f, -0.5f } },

      { { -0.5f, -0.5f, -0.5f } },
      { { -0.5f,  0.5f, -0.5f } },
      { {  0.5f,  0.5f, -0.5f } },
      { {  0.5f, -0.5f, -0.5f } },

      { { -0.5f, -0.5f,  0.5f } },
      { {  0.5f, -0.5f,  0.5f } },
      { {  0.5f,  0.5f,  0.5f } },
      { { -0.5f,  0.5f,  0.5f } } };


    static VIO_Vector normals[N_POINTS_CUBE] = {

      { { -1.0f,  0.0f,  0.0f } },
      { { -1.0f,  0.0f,  0.0f } },
      { { -1.0f,  0.0f,  0.0f } },
      { { -1.0f,  0.0f,  0.0f } },

      { {  1.0f,  0.0f,  0.0f } },
      { {  1.0f,  0.0f,  0.0f } },
      { {  1.0f,  0.0f,  0.0f } },
      { {  1.0f,  0.0f,  0.0f } },

      { {  0.0f, -1.0f,  0.0f } },
      { {  0.0f, -1.0f,  0.0f } },
      { {  0.0f, -1.0f,  0.0f } },
      { {  0.0f, -1.0f,  0.0f } },

      { {  0.0f,  1.0f,  0.0f } },
      { {  0.0f,  1.0f,  0.0f } },
      { {  0.0f,  1.0f,  0.0f } },
      { {  0.0f,  1.0f,  0.0f } },

      { {  0.0f,  0.0f, -1.0f } },
      { {  0.0f,  0.0f, -1.0f } },
      { {  0.0f,  0.0f, -1.0f } },
      { {  0.0f,  0.0f, -1.0f } },

      { {  0.0f,  0.0f,  1.0f } },
      { {  0.0f,  0.0f,  1.0f } },
      { {  0.0f,  0.0f,  1.0f } },
      { {  0.0f,  0.0f,  1.0f } }
     };
     static int end_indices[N_FACES_CUBE] = {
       4, 8, 12, 16, 20, 24
     };
     static int indices[N_POINTS_CUBE] = {
       0,  1,  2,  3,  4,  5,  6,  7,
       8,  9, 10, 11, 12, 13, 14, 15,
       16, 17, 18, 19, 20, 21, 22, 23
     };

     static  polygons_struct  polygons = {
       ONE_COLOUR,
       colours,
       { 0.4, 0.6, 0.6, 30.0, 1.0 },
       1.0,
       N_POINTS_CUBE,
       points,
       normals,
       N_FACES_CUBE,
       end_indices,
       indices,
       (VIO_SCHAR *) NULL,
       (int *) NULL,
       (bintree_struct_ptr) NULL
     };

    polygons.colours[0] = colour;
    draw_polygons( window, &polygons );
}

double DEGS_TO_RAD = M_PI / 180.0;

void
create_sphere(VIO_Point center, double radius, int n_latitude, int n_longitude,
              polygons_struct *polygons)
{
  int p;
  int s;
  int n_indices = 0;
  int n = 2;

  const int n_pitch = n_longitude + 1;
  const VIO_Real pitch_inc = (M_PI / (VIO_Real) n_pitch);
  const VIO_Real rot_inc = (M_PI * 2.0 / (VIO_Real) n_latitude);
  const VIO_Real cx = Point_x(center);
  const VIO_Real cy = Point_y(center);
  const VIO_Real cz = Point_z(center);
  const int n_triangles = n_longitude * 4;
  const int n_quads = n_latitude * (n_longitude - 1);
  const int first_offset = 2;
  const int last_offset = first_offset + (n_latitude * (n_longitude - 1));

  polygons->n_items = n_triangles + n_quads;
  polygons->n_points = 2 + (n_latitude * n_longitude);
  ALLOC(polygons->points, polygons->n_points);
  ALLOC(polygons->normals, polygons->n_points);
  ALLOC(polygons->end_indices, polygons->n_items);
  ALLOC(polygons->indices, 3 * n_triangles + 4 * n_quads);

  /* Generate top and bottom vertices. */
  fill_Point(polygons->points[0], cx, cy + radius, cz);
  fill_Point(polygons->points[1], cx, cy - radius, cz);

  fill_Vector(polygons->normals[0], 0, radius, 0);
  fill_Vector(polygons->normals[1], 0, -radius, 0);

  for (p = 1; p < n_pitch; p++)
  {
    VIO_Real out = fabs(radius * sin(p * pitch_inc));
    VIO_Real y = radius * cos(p * pitch_inc);

    for (s = 0; s < n_latitude; s++)
    {
      VIO_Real x = out * cos(s * rot_inc);
      VIO_Real z = out * sin(s * rot_inc);

      fill_Vector(polygons->normals[n], x, y, z);
      fill_Point(polygons->points[n], x + cx, y + cy, z + cz);
      n++;
    }
  }

  n = 0;

  /* Create the quadrilateral faces between intermediate points.
   */
  for (p = 1; p < n_longitude; p++) /* for each longitude line */
  {
    for (s = 0; s < n_latitude; s++) /* for each latitude line */
    {
      int i = p * n_latitude + s;
      int j = (s == n_latitude - 1) ? i - n_latitude : i;

      polygons->indices[n_indices++] = i - n_latitude + first_offset;
      polygons->indices[n_indices++] = j + 1 - n_latitude + first_offset;
      polygons->indices[n_indices++] = j + 1 + first_offset;
      polygons->indices[n_indices++] = i + first_offset;
      polygons->end_indices[n++] = n_indices;
    }
  }

  /* Create the "triangular fan" that surrounds the top and bottom
   * vertices.
   */
  int n_around = n_longitude * 2;
  for (s = 0; s < n_around; s++)
  {
    int j = (s == n_around - 1) ? -1 : s;

    polygons->indices[n_indices++] = 0; /* index of top  */
    polygons->indices[n_indices++] = j + 1 + first_offset;
    polygons->indices[n_indices++] = s + first_offset;
    polygons->end_indices[n++] = n_indices;

    polygons->indices[n_indices++] = 1; /* index of bottom */
    polygons->indices[n_indices++] = s + last_offset;
    polygons->indices[n_indices++] = j + 1 + last_offset;
    polygons->end_indices[n++] = n_indices;
  }
}

static void
draw_marker_as_sphere( Gwindow window, VIO_Colour colour )
{
  object_struct   *object = create_object( POLYGONS );
  polygons_struct *polygons = get_polygons_ptr( object );
  VIO_Point       pt = { { 0, 0, 0 } };
  initialize_polygons( polygons, colour, NULL );
  create_sphere( pt, 0.5, 16, 16/2, polygons );
  draw_polygons( window, polygons );
  delete_object( object );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_draw_pixels
@INPUT      : window
              pixels
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws the pixels in the window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_draw_pixels(
    Gwindow         window,
    pixels_struct   *pixels )
{
    if( pixels->x_size <= 0 || pixels->y_size <= 0 )
        return;

    about_to_draw_graphics( window );

    if( pixels->x_zoom > 0.0 && pixels->y_zoom > 0.0 )
        GS_set_pixel_zoom( pixels->x_zoom, pixels->y_zoom );
    else
        GS_set_pixel_zoom( 1.0, 1.0 );

    if( window->colour_map_state )
    {
        if( pixels->pixel_type == RGB_PIXEL )
        {
            handle_internal_error(
               "G_draw_pixels(): cannot draw rgb pixels in colour map mode.\n");
            return;
        }

        GS_draw_colour_map_pixels( window->x_viewport_min,
                                   window->y_viewport_min, pixels );
    }
    else
    {
        if( pixels->pixel_type != RGB_PIXEL )
        {
            handle_internal_error(
               "G_draw_pixels(): cannot draw colour map pixels in rgb mode.\n");
            return;
        }

        GS_draw_rgb_pixels( window->x_viewport_min,
                            window->y_viewport_min, pixels );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : G_read_pixels
@INPUT      : window
              x_min
              x_max
              y_min
              y_max
@OUTPUT     : pixels
@RETURNS    :
@DESCRIPTION: Reads the pixels from the frame buffer.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

  void  G_read_pixels(
    Gwindow         window,
    int             x_min,
    int             x_max,
    int             y_min,
    int             y_max,
    VIO_Colour          pixels[] )
{
    set_current_window( window );

    GS_read_pixels( G_get_colour_map_state(window),
                    x_min, x_max, y_min, y_max, pixels );
}
