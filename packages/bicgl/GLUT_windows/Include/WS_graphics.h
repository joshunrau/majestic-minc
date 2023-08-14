#ifndef  DEF_GLUT_WINDOWS
#define  DEF_GLUT_WINDOWS

#include  <graphics_base.h>

#ifdef HAVE_APPLE_OPENGL_FRAMEWORK
#include <GLUT/glut.h>
#else
#ifdef HAVE_FREEGLUT_H
#include <GL/freeglut.h>
#else
#include  <GL/glut.h>
#endif
#endif

typedef  int  Window_id;

typedef  struct
{
    Window_id    window_id;
    int          init_x, init_y;
    int          border_width, border_height;
    VIO_BOOL     is_new;
    VIO_BOOL     is_visible;
}
WS_window_struct;

typedef  WS_window_struct *  WSwindow;

#include  <glut_window_prototypes.h>

#endif
