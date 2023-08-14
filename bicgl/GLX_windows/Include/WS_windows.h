#ifndef  DEF_GLX_WINDOWS
#define  DEF_GLX_WINDOWS

#include  <graphics_base.h>

#include  <volume_io.h>
#include  <x_windows.h>

#ifdef HAVE_APPLE_OPENGL_FRAMEWORK
#include <OpenGL/glx.h>
#else
#include  <GL/glx.h>
#endif

typedef  struct
{
    X_window_struct   x_window;
    GLXContext        graphics_context;
    VIO_BOOL           overlay_present;
    X_window_struct   overlay_window;
    GLXContext        overlay_context;
    int               n_fonts;
    GLint             *font_list_bases;
}
WS_window_struct;

typedef  struct
{
    VIO_STR        font_name;
    XFontStruct   *x_font_info;
}  WS_font_info;


#include  <glx_window_prototypes.h>

#ifdef  public_was_defined_here
#undef       
#undef       public_was_defined_here
#endif

#endif
