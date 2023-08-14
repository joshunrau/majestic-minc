#define GL_GLEXT_PROTOTYPES 1

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H  

#include  <graphics.h>
#include  <GS_graphics.h>

#if PARANOID
#define GLCHECK {GLint _temp = glGetError(); if (_temp != GL_NO_ERROR) { fprintf(stderr, "OpenGL reporting error 0x%x at line %d.\n", _temp, __LINE__); } }
#else
#define GLCHECK
#endif

  void  GS_set_point(
    VIO_Point  *point )
{
    glVertex3fv( point->coords );
}

  void  GS_set_normal(
    VIO_Vector  *normal )
{
    glNormal3fv( normal->coords );
}

  void  GS_set_colour(
    VIO_Colour  colour )
{
    glColor4ubv( (GLubyte*) &colour );
}

  void  GS_set_colour_index(
    VIO_Colour  colour )
{
    glIndexi( (GLint) colour );
}

  void  GS_set_ambient_and_diffuse_mode(
    VIO_BOOL  state )
{
    if( state )
    {
        glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
        glEnable( GL_COLOR_MATERIAL );
    }
    else
        glDisable( GL_COLOR_MATERIAL );
}

/* ARGSUSED */

  void  GS_initialize_surface_property(
    GSwindow        window )
{
}

/* ARGSUSED */

  void  GS_set_surface_property(
    GSwindow       window,
    VIO_Colour         col,
    VIO_Surfprop       *surfprop )
{
    float  props[4];
    VIO_Real   a, d, s;

    a = (VIO_Real) Surfprop_a( *surfprop );
    d = (VIO_Real) Surfprop_d( *surfprop );
    s = (VIO_Real) Surfprop_s( *surfprop );

    props[0] = (float) (a * get_Colour_r_0_1( col ));
    props[1] = (float) (a * get_Colour_g_0_1( col ));
    props[2] = (float) (a * get_Colour_b_0_1( col ));
    props[3] = (float) get_Colour_a_0_1( col );
    glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, props );

    props[0] = (float) (d * get_Colour_r_0_1( col ));
    props[1] = (float) (d * get_Colour_g_0_1( col ));
    props[2] = (float) (d * get_Colour_b_0_1( col ));
    props[3] = (float) get_Colour_a_0_1( col );
    glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, props );

    props[0] = (float) s;
    props[1] = (float) s;
    props[2] = (float) s;
    props[3] = 1.0f;
    glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, props );

    props[0] = Surfprop_se( *surfprop );
    glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, props );
}

  void  GS_set_line_width(
    VIO_Real  line_width )
{
    glLineWidth( (GLfloat) line_width );
}

/* ARGSUSED */

  void  GS_curve(
    float  geom[4][3] )
{
  fprintf(stderr, "Curves are not implemented in bicgl.\n");
}

  void  GS_begin_point( void )
{
    glBegin( GL_POINTS );
}

  void  GS_end_point( void )
{
    glEnd();
}

  void  GS_begin_line( void )
{
    glBegin( GL_LINE_STRIP );
}

  void  GS_end_line( void )
{
    glEnd();
}

  void  GS_begin_closed_line( void )
{
    glBegin( GL_LINE_LOOP );
}

  void  GS_end_closed_line( void )
{
    glEnd();
}

  void  GS_begin_polygon( void )
{
    glBegin( GL_POLYGON );
}

  void  GS_end_polygon( void )
{
    glEnd();
}

  void  GS_begin_quad_strip( void )
{
    glBegin( GL_QUAD_STRIP );
}

  void  GS_end_quad_strip( void )
{
    glEnd();
}

  void  GS_set_raster_position(
    VIO_Real  x,
    VIO_Real  y,
    VIO_Real  z )
{
    glRasterPos3d( (double) x, (double) y, (double) z );
}

  void  GS_set_pixel_zoom(
    VIO_Real  x_zoom,
    VIO_Real  y_zoom )
{
    glPixelZoom( (GLfloat) x_zoom, (GLfloat) y_zoom );
}

/* ARGSUSED */

  void  GS_draw_colour_map_pixels(
    int             x_viewport_min,
    int             y_viewport_min,
    pixels_struct   *pixels )
{
    int       x, y, x_size, y_size;

    x = pixels->x_position;
    y = pixels->y_position;
    x_size = pixels->x_size;
    y_size = pixels->y_size;

    GS_set_raster_position( (VIO_Real) x, (VIO_Real) y, 0.0 );

    switch( pixels->pixel_type )
    {
    case COLOUR_INDEX_8BIT_PIXEL:
        glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
        glDrawPixels( x_size, y_size, GL_COLOR_INDEX,
                      GL_UNSIGNED_BYTE,
                      (void *) &PIXEL_COLOUR_INDEX_8(*pixels,0,0) );
        break;

    case COLOUR_INDEX_16BIT_PIXEL:
        glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
        glDrawPixels( x_size, y_size, GL_COLOR_INDEX,
                      GL_UNSIGNED_SHORT,
                      (void *) &PIXEL_COLOUR_INDEX_16(*pixels,0,0) );
        break;
    default:
        break;
    }
}

/* ARGSUSED */

  void  GS_draw_rgb_pixels(
    int             x_viewport_min,
    int             y_viewport_min,
    pixels_struct   *pixels )
{
    int       x, y, x_size, y_size;
    void      *void_ptr;

    x = pixels->x_position;
    y = pixels->y_position;
    x_size = pixels->x_size;
    y_size = pixels->y_size;

    GS_set_raster_position( (VIO_Real) x, (VIO_Real) y, 0.0 );

    void_ptr = (void *) &PIXEL_RGB_COLOUR(*pixels,0,0);
    glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
    glDrawPixels( x_size, y_size, GL_RGBA, GL_UNSIGNED_BYTE, void_ptr );
}

  void  GS_read_pixels(
    VIO_BOOL         colour_map_state,
    int             x_min,
    int             x_max,
    int             y_min,
    int             y_max,
    VIO_Colour          pixels[] )
{
    glReadBuffer( GL_FRONT );

    if( colour_map_state )
    {
        glReadPixels( x_min, y_min, x_max - x_min + 1, y_max - y_min + 1,
                      GL_COLOR_INDEX, GL_UNSIGNED_INT, pixels );
    }
    else
    {
        glReadPixels( x_min, y_min, x_max - x_min + 1, y_max - y_min + 1,
                      GL_RGBA, GL_UNSIGNED_BYTE, pixels );
    }
}

  VIO_Real  GS_get_character_height(
    Font_types       type,
    VIO_Real             size )
{
    return( WS_get_character_height( type, size ) );
}

/* ARGSUSED */

  void  GS_draw_text(
    Font_types   type,
    VIO_Real         size,
    VIO_STR       string )
{
    WS_draw_text( type, size, string );
}

  VIO_Real  GS_get_text_length(
    VIO_STR           str,
    Font_types       type,
    VIO_Real             size )
{
    return( WS_get_text_length( str, type, size ) );
}

#define N_DIM 3
#define N_VTX_PER_TRIANGLE 3

/**
 * Simple helper function to set the opacity in the fragment shader.
 */
static void set_program_opacity(GLint program, GLfloat opacity)
{
    glUniform1f( glGetUniformLocation( program, "opacity" ), opacity );
}

/**
 * Simple helper function to set the "model expansion factor" in the vertex
 * shader.
 */
static void set_program_mxfactor(GLint program, GLfloat mxfactor)
{
    glUniform1f( glGetUniformLocation( program, "mxfactor" ), mxfactor );
}

/**
 * Sets the amount of expansion in the model when drawing the wireframe
 * overlay on top of the shaded surface.
 */
#define MX_OVERLAY 0.9999

static void do_wireframe_overlay( GSwindow window,
                                  polygons_struct *polygons,
                                  VIO_Colour colour)
{
    GLint program = window->programs[PROGRAM_SINGLE];
    glUseProgram( program );
    set_program_opacity( program, Surfprop_t( polygons->surfprop ) );
    set_program_mxfactor( program, MX_OVERLAY );
    GS_set_surface_property( window, colour, &polygons->surfprop );
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_draw_polygons_one_colour
@INPUT      : window
              shading_mode
              polygons
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws a set of polygons which has a single colour for all
              polygons.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

/**
 * Draws a polygonal object (non-triangular) using glDrawElements.
 * Supports either single colour or per-item colours. This relies on
 * glDrawElements supporting GL_POLYGON as a mode argument -
 * apparently this is deprecated behavior.
 *
 * \param window The window in which we are drawing.
 * \param polygons The polygons_struct to draw.
 */
void
GS_draw_polygons_one_colour(GSwindow window, Shading_modes shading_mode, 
                            polygons_struct *polygons)
{
  GLuint vbo_points;
  GLuint vbo_normals;
  GLint loc_position;
  GLint loc_normal = -1;
  GLint program;
  int start_index = 0;
  int end_index = 0;
  int i;

#if DEBUG
  printf("dp_one_colour %d %d (normals %p points %d items %d opacity %f)\n",
         shading_mode, 
         polygons->normals, polygons->n_points, polygons->n_items,
         Surfprop_t(polygons->surfprop));
#endif  /* DEBUG */

  if (polygons->normals != NULL)
  {
    /*
     * We have surface normals, so use the standard program for
     * vertex+attribute shading.
     */
    program = window->programs[PROGRAM_SINGLE];

    /*
     * Set up nondefault surface properties.
     */
    GS_set_surface_property( window, polygons->colours[0], &polygons->surfprop );
  }
  else
  {
    program = window->programs[PROGRAM_TRIVIAL];
  }

  GS_set_colour( polygons->colours[0] );
  if( shading_mode )
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL );
  }
  else
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE );
  }
  glUseProgram(program);

  set_program_opacity(program, Surfprop_t(polygons->surfprop));
  set_program_mxfactor( program, 1.0 );

  glGenBuffers(1, &vbo_points);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_points);
  glBufferData(GL_ARRAY_BUFFER, polygons->n_points * sizeof(VIO_Point),
               polygons->points, GL_STATIC_DRAW);

  loc_position = glGetAttribLocation(program, "position");
  glVertexAttribPointer(loc_position, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(loc_position);
  if (polygons->normals != NULL)
  {
    glGenBuffers(1, &vbo_normals);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_normals);
    glBufferData(GL_ARRAY_BUFFER, polygons->n_points * sizeof(VIO_Point),
                 polygons->normals, GL_STATIC_DRAW);

    loc_normal = glGetAttribLocation(program, "normal");
    glVertexAttribPointer(loc_normal, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(loc_normal);
  }

  if (polygons->colour_flag == PER_ITEM_COLOURS)
  {
    for (i = 0; i < polygons->n_items; i++)
    {
      start_index = end_index;
      end_index = polygons->end_indices[i];
      if (polygons->visibilities == NULL || polygons->visibilities[i])
      {
        GS_set_colour( polygons->colours[i] );
        if (polygons->normals != NULL)
        {
          GS_set_surface_property( window, polygons->colours[i],
                                   &polygons->surfprop);
        }
        glDrawElements(GL_POLYGON, end_index - start_index, GL_UNSIGNED_INT,
                       &polygons->indices[start_index]);
      }
    }
  }
  else
  {
    for (i = 0; i < polygons->n_items; i++)
    {
      start_index = end_index;
      end_index = polygons->end_indices[i];
      if (polygons->visibilities == NULL || polygons->visibilities[i])
      {
        glDrawElements(GL_POLYGON, end_index - start_index, GL_UNSIGNED_INT,
                       &polygons->indices[start_index]);
      }
    }
  }

  /* Clean up.
   */
  if (polygons->normals != NULL)
  {
    glDisableVertexAttribArray(loc_normal);
    glDeleteBuffers(1, &vbo_normals);
  }

  glDisableVertexAttribArray(loc_position);
  glDeleteBuffers(1, &vbo_points);

  glUseProgram(0);
}

/**
 * Choose whether to use black or white to get the best contrast
 * with the given colour.
 */
static VIO_Colour get_contrasting_colour(VIO_Colour colour)
{
    int r = get_Colour_r(colour);
    int g = get_Colour_g(colour);
    int b = get_Colour_b(colour);
    int yiq = ((r * 299) + (g * 587) + (b * 114)) / 1000;
    return (yiq >= 128) ? BLACK : WHITE;
}

/**
 * Draw a triangular mesh in a single colour.
 */
void
GS_draw_triangles_one_colour(GSwindow window, Shading_modes shading_mode, 
                             polygons_struct *polygons)
{
  GLuint vbo_points;
  GLuint vbo_normals;
  GLuint ebo;
  GLint loc_position;
  GLint loc_normal = -1;
  int *eboBuffer = NULL;  /* Temporary storage for visible indices. */
  int n_indices = 0;
  GLenum mode = GL_TRIANGLES;
  GLint program;

#if DEBUG
  printf("dt_one_colour %d (normals %p points %d items %d opacity %f)\n",
         shading_mode,
         polygons->normals, polygons->n_points, polygons->n_items,
         Surfprop_t(polygons->surfprop));
#endif  /* DEBUG */

  n_indices = NUMBER_INDICES(*polygons);

  /* If the visibilities field is set, we need to copy the visible
   * indices to a new buffer. I don't know of any other way to control
   * visibility on a per-polygon basis (perhaps we could set the
   * colour to a transparent value in the fragment shader?).
   */
  if (polygons->visibilities != NULL)
  {
    int i, k = 0;

    eboBuffer = malloc(n_indices * sizeof(int));
    if (eboBuffer == NULL)
    {
      return;
    }

    for (i = 0; i < polygons->n_items; i++)
    {
      if (polygons->visibilities[i])
      {
        int start_index = START_INDEX( polygons->end_indices, i );
        int end_index = polygons->end_indices[i];
        int j;
        for (j = start_index; j < end_index; j++)
        {
          eboBuffer[k++] = polygons->indices[j];
        }
      }
    }
    n_indices = k;
  }

  if (polygons->normals != NULL)
  {
    program = window->programs[PROGRAM_SINGLE];

    /*
     * Set up nondefault surface properties.
     */
    GS_set_surface_property( window, polygons->colours[0], &polygons->surfprop );
  }
  else
  {
    program = window->programs[PROGRAM_TRIVIAL];
  }

  GS_set_colour( polygons->colours[0] );
  GLCHECK;
  switch ( shading_mode )
  {
  case WIREFRAME:
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    break;
  case POINT:
    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    break;
  default:
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    break;
  }

  glUseProgram( program );

  set_program_opacity( program, Surfprop_t(polygons->surfprop) );
  set_program_mxfactor( program, 1.0 );

  glGenBuffers(1, &vbo_points);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_points);
  glBufferData(GL_ARRAY_BUFFER, polygons->n_points * sizeof(VIO_Point),
               polygons->points, GL_STATIC_DRAW);

  loc_position = glGetAttribLocation(program, "position");
  glVertexAttribPointer(loc_position, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(loc_position);

  if (polygons->normals != NULL)
  {
    glGenBuffers(1, &vbo_normals);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_normals);
    glBufferData(GL_ARRAY_BUFFER, polygons->n_points * sizeof(VIO_Point),
                 polygons->normals, GL_STATIC_DRAW);

    loc_normal = glGetAttribLocation(program, "normal");
    glVertexAttribPointer(loc_normal, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(loc_normal);
  }

  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               n_indices * sizeof(int),
               (eboBuffer != NULL) ? eboBuffer : polygons->indices,
               GL_STATIC_DRAW);

  glDrawElements(mode, n_indices, GL_UNSIGNED_INT, 0);
  if ( shading_mode == OVERLAY )
  {
    do_wireframe_overlay( window, polygons,
                          get_contrasting_colour(polygons->colours[0]) );
    glDrawElements(mode, n_indices, GL_UNSIGNED_INT, 0);
  }

  glDeleteBuffers(1, &ebo);

  if (polygons->normals != NULL)
  {
    glDisableVertexAttribArray(loc_normal);
    glDeleteBuffers(1, &vbo_normals);
  }

  glDisableVertexAttribArray(loc_position);
  glDeleteBuffers(1, &vbo_points);

  if (eboBuffer != NULL)
  {
    free(eboBuffer);
  }
  glUseProgram(0);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : draw_triangles_per_item_colours
@INPUT      : window
              shading_mode
              polygons
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws a set of polygons which has a separate colour for each
              polygon.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */
void
GS_draw_triangles_per_item_colours(GSwindow window,
                                   Shading_modes shading_mode,
                                   polygons_struct *polygons)
{
  GLuint vbo_points;
  GLuint vbo_normals;
  GLuint vbo_colours;
  GLint loc_position;
  GLint loc_normal;
  GLint loc_colour;
  GLint loc_surfprop;
  float *bufPoints = NULL; /* Temporary storage for visible vertices. */
  float *bufNormals = NULL;
  GLuint *bufColours = NULL;
  int i, k = 0;
  int start_index = 0;
  int end_index = 0;
  int n_items;
  GLuint program;

#if DEBUG
  printf("tri_per_item: %d %d\n", polygons->n_points, polygons->n_items);
#endif /* DEBUG */

  /* To save a little space and time, allocate all of the memory we
   * will need in one go...
   */
  bufPoints = malloc(polygons->n_items * sizeof(float) * N_VTX_PER_TRIANGLE * N_DIM +
                     polygons->n_items * sizeof(float) * N_VTX_PER_TRIANGLE * N_DIM +
                     polygons->n_items * sizeof(GLuint) * N_VTX_PER_TRIANGLE);
  if (bufPoints == NULL)
  {
    return;
  }
  bufNormals = &bufPoints[polygons->n_items * N_VTX_PER_TRIANGLE * N_DIM];
  bufColours = (GLuint *)&bufNormals[polygons->n_items * N_VTX_PER_TRIANGLE * N_DIM];

  for (i = 0; i < polygons->n_items; i++)
  {
    start_index = end_index;
    end_index = polygons->end_indices[i];
    if (polygons->visibilities == NULL || polygons->visibilities[i])
    {
      int j;
      for (j = start_index; j < end_index; j++)
      {
        int c = polygons->indices[j];
        int nk = k * N_DIM;
        bufPoints[nk + VIO_X] = Point_x(polygons->points[c]);
        bufPoints[nk + VIO_Y] = Point_y(polygons->points[c]);
        bufPoints[nk + VIO_Z] = Point_z(polygons->points[c]);
        bufNormals[nk + VIO_X] = Vector_x(polygons->normals[c]);
        bufNormals[nk + VIO_Y] = Vector_y(polygons->normals[c]);
        bufNormals[nk + VIO_Z] = Vector_z(polygons->normals[c]);
        bufColours[k] = polygons->colours[i];
        k++;
      }
    }
  }
  n_items = k;

  program = window->programs[PROGRAM_VERTEX];

  switch ( shading_mode )
  {
  case WIREFRAME:
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    break;
  case POINT:
    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    break;
  default:
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    break;
  }

  glUseProgram(program);

  set_program_opacity(program, Surfprop_t(polygons->surfprop));

  glGenBuffers(1, &vbo_points);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_points);
  glBufferData(GL_ARRAY_BUFFER, n_items * sizeof(float) * N_VTX_PER_TRIANGLE,
               bufPoints, GL_STATIC_DRAW);

  loc_position = glGetAttribLocation(program, "position");
  glVertexAttribPointer(loc_position, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(loc_position);

  loc_surfprop = glGetUniformLocation(program, "surfprop");
  /* sets the crucial 4 parameters for surface properties:
     ambient, diffuse, specular, shininess
  */
  glUniform4f(loc_surfprop, polygons->surfprop.a, polygons->surfprop.d, polygons->surfprop.s, polygons->surfprop.se);

  glGenBuffers(1, &vbo_normals);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_normals);
  glBufferData(GL_ARRAY_BUFFER, n_items * sizeof(float) * N_VTX_PER_TRIANGLE,
               bufNormals, GL_STATIC_DRAW);

  loc_normal = glGetAttribLocation(program, "normal");
  glVertexAttribPointer(loc_normal, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(loc_normal);

  glGenBuffers(1, &vbo_colours);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_colours);
  glBufferData(GL_ARRAY_BUFFER, n_items * sizeof(GLuint),
               bufColours, GL_STATIC_DRAW);

  loc_colour = glGetAttribLocation(program, "colour");
  glVertexAttribPointer(loc_colour, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
  glEnableVertexAttribArray(loc_colour);

  glDrawArrays(GL_TRIANGLES, 0, n_items);

  if (shading_mode == OVERLAY)
  {
    do_wireframe_overlay(window, polygons, WHITE);
    glDrawArrays(GL_TRIANGLES, 0, n_items);
  }

  glDisableVertexAttribArray(loc_position);
  glDisableVertexAttribArray(loc_normal);
  glDisableVertexAttribArray(loc_colour);

  glDeleteBuffers(1, &vbo_colours);
  glDeleteBuffers(1, &vbo_normals);
  glDeleteBuffers(1, &vbo_points);
  free(bufPoints);
  glUseProgram(0);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_draw_polygons_per_vertex_colours
@INPUT      : window
              shading_mode
              polygons
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws a set of polygons which have a separate colour for each
              vertex.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */
void
GS_draw_polygons_per_vertex_colours(GSwindow window, 
                                    Shading_modes shading_mode,
                                    polygons_struct *polygons)
{
  GLuint vbo_points;
  GLuint vbo_normals;
  GLuint vbo_colours;
  GLuint ebo;
  GLint loc_position;
  GLint loc_normal;
  GLint loc_colour;
  GLint loc_surfprop;
  int *eboBuffer = NULL;        /* Temporary storage for visible indices. */
  int n_indices = 0;
  GLenum mode;
  GLuint program;

  if (!is_triangular(polygons))
  {
    mode = GL_POLYGON;
  }
  else
  {
    mode = GL_TRIANGLES;
  }

  n_indices = NUMBER_INDICES(*polygons);

  /* If the visibilities field is set, we need to copy the visible
   * indices to a new buffer. I don't know of any other way to control
   * visibility on a per-polygon basis (perhaps we could set the
   * colour to a transparent value in the fragment shader?).
   */
  if (polygons->visibilities != NULL)
  {
    int i, k = 0;

    eboBuffer = malloc(n_indices * sizeof(int));
    if (eboBuffer == NULL)
    {
      return;
    }

    for (i = 0; i < polygons->n_items; i++)
    {
      if (polygons->visibilities[i])
      {
        int start_index = START_INDEX( polygons->end_indices, i );
        int end_index = polygons->end_indices[i];
        int j;
        for (j = start_index; j < end_index; j++)
        {
          eboBuffer[k++] = polygons->indices[j];
        }
      }
    }
    n_indices = k;
#if DEBUG
    printf("Displaying %d/%d indices\n", n_indices, NUMBER_INDICES(*polygons));
#endif  /* DEBUG */
  }

  program = window->programs[PROGRAM_VERTEX];

  switch ( shading_mode )
  {
  case WIREFRAME:
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    break;
  case POINT:
    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    break;
  default:
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    break;
  }

  glUseProgram(program);

  set_program_opacity(program, Surfprop_t(polygons->surfprop));

  glGenBuffers(1, &vbo_points);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_points);
  glBufferData(GL_ARRAY_BUFFER, polygons->n_points * sizeof(VIO_Point),
               polygons->points, GL_STATIC_DRAW);

  loc_position = glGetAttribLocation(program, "position");
  glVertexAttribPointer(loc_position, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(loc_position);

  loc_surfprop = glGetUniformLocation(program, "surfprop");
  /* sets the crucial 4 parameters for surface properties:
     ambient, diffuse, specular, shininess
  */
  glUniform4f(loc_surfprop, polygons->surfprop.a, polygons->surfprop.d, polygons->surfprop.s, polygons->surfprop.se);

  glGenBuffers(1, &vbo_normals);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_normals);
  glBufferData(GL_ARRAY_BUFFER, polygons->n_points * sizeof(VIO_Point),
               polygons->normals, GL_STATIC_DRAW);

  loc_normal = glGetAttribLocation(program, "normal");
  glVertexAttribPointer(loc_normal, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(loc_normal);

  glGenBuffers(1, &vbo_colours);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_colours);
  glBufferData(GL_ARRAY_BUFFER, polygons->n_points * sizeof(VIO_Colour),
               polygons->colours, GL_STATIC_DRAW);

  loc_colour = glGetAttribLocation(program, "colour");
  glVertexAttribPointer(loc_colour, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
  glEnableVertexAttribArray(loc_colour);

  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               n_indices * sizeof(int),
               (eboBuffer != NULL) ? eboBuffer : polygons->indices,
               GL_STATIC_DRAW);

  glDrawElements( mode, n_indices, GL_UNSIGNED_INT, 0 );
  if ( shading_mode == OVERLAY )
  {
    do_wireframe_overlay( window, polygons, WHITE );
    glDrawElements( mode, n_indices, GL_UNSIGNED_INT, 0 );
  }

  glDisableVertexAttribArray(loc_position);
  glDisableVertexAttribArray(loc_normal);
  glDisableVertexAttribArray(loc_colour);

  glDeleteBuffers(1, &ebo);
  glDeleteBuffers(1, &vbo_colours);
  glDeleteBuffers(1, &vbo_normals);
  glDeleteBuffers(1, &vbo_points);
  if (eboBuffer != NULL)
  {
    free(eboBuffer);
  }
  glUseProgram(0);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_draw_quadmesh_one_colour
@INPUT      : window
              quadmesh
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws a quadmesh of a single colour into the window.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */
void GS_draw_quadmesh_one_colour(GSwindow window, Shading_modes shading_mode,
                                 quadmesh_struct *quadmesh)
{
  GLuint vbo_points;
  GLuint vbo_normals;
  GLuint ebo;
  GLint loc_position;
  GLint loc_normal = -1;
  int n_items = 0;
  GLint program;
  int m_size, n_size;
  int i, j;
  GLuint *indices;

#if DEBUG
  printf("GS_draw_quadmesh_one_colour %d %d %d %d %d %p (NEW)\n", quadmesh->m_closed, quadmesh->n_closed, quadmesh->m, quadmesh->n, quadmesh->colour_flag, quadmesh->normals);
#endif  /* DEBUG */

  get_quadmesh_n_objects(quadmesh, &m_size, &n_size);
  n_items = quadmesh->m * quadmesh->n;

  if (quadmesh->normals != NULL)
  {
    program = window->programs[PROGRAM_SINGLE];
    /*
     * Set up nondefault surface properties.
     */
    GS_set_surface_property( window, quadmesh->colours[0], &quadmesh->surfprop );
  }
  else
  {
    program = window->programs[PROGRAM_TRIVIAL];
  }

  GS_set_colour( quadmesh->colours[0] );
  if( shading_mode )
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL );
  }
  else
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE );
  }
  glUseProgram(program);

  set_program_opacity(program, Surfprop_t(quadmesh->surfprop));
  set_program_mxfactor( program, 1.0 );

  glGenBuffers(1, &vbo_points);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_points);
  glBufferData(GL_ARRAY_BUFFER, n_items * sizeof(VIO_Point),
               quadmesh->points, GL_STATIC_DRAW);

  loc_position = glGetAttribLocation(program, "position");
  glVertexAttribPointer(loc_position, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(loc_position);
  GLCHECK;

  if (quadmesh->normals != NULL)
  {
    glGenBuffers(1, &vbo_normals);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_normals);
    glBufferData(GL_ARRAY_BUFFER, n_items * sizeof(VIO_Vector),
                 quadmesh->normals, GL_STATIC_DRAW);
    loc_normal = glGetAttribLocation(program, "normal");
    glVertexAttribPointer(loc_normal, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(loc_normal);
    GLCHECK;
  }

  indices = malloc((n_size + 1) * 2 * sizeof(GLuint));

  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

  /* Draw these in m_size separate strips, otherwise we get unwanted
   * extra faces "closing" the surface.
   */
  for (i = 0; i < m_size; i++)
  {
    int k = 0;
    for (j = 0; j <= n_size; j++)
    {
      int used_j = j % quadmesh->n;
      indices[k++] = VIO_IJ(i, used_j, quadmesh->n);
      indices[k++] = VIO_IJ((i + 1) % quadmesh->m, used_j, quadmesh->n);
    }
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, k * sizeof(GLuint),
                 indices, GL_STATIC_DRAW);
    glDrawElements(GL_TRIANGLE_STRIP, k, GL_UNSIGNED_INT, 0);
    GLCHECK;
  }

  free(indices);

  if (quadmesh->normals != NULL)
  {
    glDisableVertexAttribArray(loc_normal);
    glDeleteBuffers(1, &vbo_normals);
    GLCHECK;
  }

  glDeleteBuffers(1, &ebo);
  glDisableVertexAttribArray(loc_position);
  glDeleteBuffers(1, &vbo_points);
  GLCHECK;

  glUseProgram(0);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GSdraw_quadmesh_per_item_colours
@INPUT      : window
              shading_mode
              quadmesh
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws a quadmesh with a separate colour per quadrilateral.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

void  GS_draw_quadmesh_per_item_colours(
    GSwindow         window,
    Shading_modes    shading_mode,
    quadmesh_struct *quadmesh )
{
  printf("draw_quadmesh_per_item_colours\n");
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : GS_draw_quadmesh_per_vertex_colours
@INPUT      : window
              shading_mode
              quadmesh
@OUTPUT     :
@RETURNS    :
@DESCRIPTION: Draws a quadmesh with a separate colour per vertex.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : 1993            David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */
void GS_draw_quadmesh_per_vertex_colours(GSwindow window, 
                                         Shading_modes shading_mode,
                                         quadmesh_struct *quadmesh)
{
  GLuint vbo_points;
  GLuint vbo_normals;
  GLuint vbo_colours;
  GLuint ebo;
  GLint loc_position;
  GLint loc_normal;
  GLint loc_colour;
  GLint loc_surfprop;
  int n_items = 0;
  GLint program;
  int m_size, n_size;
  int i, j;
  GLuint *indices;

#if DEBUG
  printf("GS_draw_quadmesh_per_vertex_colours %d %d %d %d %d %p (NEW)\n", quadmesh->m_closed, quadmesh->n_closed, quadmesh->m, quadmesh->n, quadmesh->colour_flag, quadmesh->normals);
  printf("  surfprop %f %f %f %f %f\n", 
         Surfprop_a(quadmesh->surfprop),
         Surfprop_d(quadmesh->surfprop),
         Surfprop_s(quadmesh->surfprop),
         Surfprop_se(quadmesh->surfprop),
         Surfprop_t(quadmesh->surfprop));
#endif /* DEBUG */

  get_quadmesh_n_objects(quadmesh, &m_size, &n_size);
  n_items = quadmesh->m * quadmesh->n;

  program = window->programs[PROGRAM_VERTEX];
  glUseProgram(program);
  set_program_opacity(program, Surfprop_t(quadmesh->surfprop));

  if( shading_mode!=WIREFRAME )  
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL );
  }
  else
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE );
  }

  glGenBuffers(1, &vbo_points);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_points);
  glBufferData(GL_ARRAY_BUFFER, n_items * sizeof(VIO_Point),
               quadmesh->points, GL_STATIC_DRAW);

  loc_position = glGetAttribLocation(program, "position");
  glVertexAttribPointer(loc_position, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(loc_position);
  GLCHECK;

  loc_surfprop = glGetUniformLocation(program, "surfprop"); /*VF:maybe don't set shinenes to the same value?*/
  glUniform4f(loc_surfprop, quadmesh->surfprop.a, quadmesh->surfprop.d, quadmesh->surfprop.s, quadmesh->surfprop.se);

  glGenBuffers(1, &vbo_normals);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_normals);
  glBufferData(GL_ARRAY_BUFFER, n_items * sizeof(VIO_Vector),
               quadmesh->normals, GL_STATIC_DRAW);
  
  loc_normal = glGetAttribLocation(program, "normal");
  glVertexAttribPointer(loc_normal, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(loc_normal);

  glGenBuffers(1, &vbo_colours);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_colours);
  glBufferData(GL_ARRAY_BUFFER, n_items * sizeof(VIO_Colour),
               quadmesh->colours, GL_STATIC_DRAW);

  loc_colour = glGetAttribLocation(program, "colour");
  glVertexAttribPointer(loc_colour, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
  glEnableVertexAttribArray(loc_colour);

  indices = malloc((n_size + 1) * 2 * sizeof(GLuint));

  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

  /* Draw these in m_size separate strips, otherwise we get unwanted
   * extra faces "closing" the surface.
   */
  for (i = 0; i < m_size; i++)
  {
    int k = 0;
    for (j = 0; j <= n_size; j++)
    {
      int used_j = j % quadmesh->n;
      indices[k++] = VIO_IJ(i, used_j, quadmesh->n);
      indices[k++] = VIO_IJ((i + 1) % quadmesh->m, used_j, quadmesh->n);
    }
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, k * sizeof(GLuint),
                 indices, GL_STATIC_DRAW);
    glDrawElements(GL_TRIANGLE_STRIP, k, GL_UNSIGNED_INT, 0);
    GLCHECK;
  }

  free(indices);

  glDisableVertexAttribArray(loc_position);
  glDisableVertexAttribArray(loc_colour);
  glDisableVertexAttribArray(loc_normal);
  glDeleteBuffers(1, &ebo);
  glDeleteBuffers(1, &vbo_points);
  glDeleteBuffers(1, &vbo_colours);
  glDeleteBuffers(1, &vbo_normals);

  glUseProgram(0);
}

void
GS_draw_lines(GSwindow window, lines_struct *lines)
{
  GLuint vbo_points;            /* OpenGL vertex buffer object */
  GLuint vbo_colours;
  GLint loc_position;           /* Location of the "position" attribute. */
  GLint loc_colour;
  GLint loc_surfprop;
  GLint program;
  int start_index = 0;
  int end_index = 0;
  int i;

  if (lines->colour_flag == PER_VERTEX_COLOURS)
  {
    program = window->programs[PROGRAM_VERTEX];
  }
  else
  {
    GS_set_colour( lines->colours[0] );
    program = window->programs[PROGRAM_TRIVIAL];
  }

  glUseProgram(program);

  set_program_opacity(program, 1.0);

  /* Set up the buffer for our vertex data.
   */
  glGenBuffers(1, &vbo_points);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_points);
  glBufferData(GL_ARRAY_BUFFER, lines->n_points * sizeof(VIO_Point),
               lines->points, GL_STATIC_DRAW);

  /* Assign our vertex data to the "position" value in the shader.
   */
  loc_position = glGetAttribLocation(program, "position");
  glVertexAttribPointer(loc_position, N_DIM, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(loc_position);

  if (lines->colour_flag == PER_VERTEX_COLOURS)
  {
    VIO_Surfprop spr = {1, 0, 0, 0, 1};
    loc_surfprop = glGetUniformLocation(program, "surfprop");
    /* sets the crucial 4 parameters for surface properties:
       ambient, diffuse, specular, shininess
    */
    glUniform4f(loc_surfprop,spr.a, spr.d, spr.s, spr.se);

    /* Set up the buffer for the colour data.
     */
    glGenBuffers(1, &vbo_colours);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_colours);
    glBufferData(GL_ARRAY_BUFFER, lines->n_points * sizeof(GLuint),
                 lines->colours, GL_STATIC_DRAW);

    /* Assign the colour data to the "colour" value in the shader.
     */
    loc_colour = glGetAttribLocation(program, "colour");
    glVertexAttribPointer(loc_colour, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
    glEnableVertexAttribArray(loc_colour);
  }

  /* Now start sending the line segments.
   */
  for (i = 0; i < lines->n_items; i++)
  {
    start_index = end_index;
    end_index = lines->end_indices[i];
    if (lines->colour_flag == PER_ITEM_COLOURS)
    {
      GS_set_colour( lines->colours[i] );
    }
    glDrawElements(GL_LINE_STRIP, end_index - start_index, GL_UNSIGNED_INT,
                   &lines->indices[start_index]);
  }

  /* Clean up after ourselves.
   */
  if (lines->colour_flag == PER_VERTEX_COLOURS)
  {
    glDisableVertexAttribArray(loc_colour);
    glDeleteBuffers(1, &vbo_colours);
  }
  glDisableVertexAttribArray(loc_position);
  glDeleteBuffers(1, &vbo_points);
  glUseProgram(0);
}
