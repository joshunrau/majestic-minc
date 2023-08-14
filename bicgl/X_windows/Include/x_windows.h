#ifndef  DEF_X_WINDOWS
#define  DEF_X_WINDOWS

#include  <bicpl.h>

#include  <X11/Xlib.h>
#include  <X11/Xutil.h>

#ifdef  Status
#undef  Status
#endif

//#include  <window_events.h>

typedef  Window   Window_id;

typedef  struct
{
    Window_id        window_id;
    VIO_BOOL          colour_map_mode;
    Colormap         colour_map;
    int              n_colours;
    unsigned long    *colour_indices;
    XVisualInfo      *visual;
    GC               graphics_context;
}
X_window_struct;

#include  <x_window_prototypes.h>

#ifdef  public_was_defined_here
#undef       public_was_defined_here
#endif


#endif
