#ifndef  DEF_GRAPHICS_SYSTEM_SPECIFIC
#define  DEF_GRAPHICS_SYSTEM_SPECIFIC

#include <bicpl.h>
#include <WS_graphics.h>

#ifdef HAVE_APPLE_OPENGL_FRAMEWORK
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

typedef  struct
{
    int             light_index;
    Light_types     type;
    VIO_Colour          colour;
    VIO_Vector          direction;
    VIO_Point           position;
    VIO_Real            spot_exponent;
    VIO_Real            spot_angle;
} light_info_struct;

enum opengl_programs
{
  PROGRAM_TRIVIAL,              /* Trivial */
  PROGRAM_SINGLE,               /* Single colour */
  PROGRAM_VERTEX,               /* Per-vertex colours */
  N_PROGRAMS
};

typedef  struct
{
    WS_window_struct   WS_window;

    VIO_BOOL           ambient_set;
    VIO_Real           ambient[3];

    int                n_lights_defined;
    light_info_struct  *lights;

    int                n_light_states;
    VIO_BOOL           *light_indices;
    VIO_BOOL           *light_states;

    GLuint             programs[N_PROGRAMS];
}
GS_window_struct;

typedef  GS_window_struct   *GSwindow;

#include  <opengl_graphics_prototypes.h>

#endif
