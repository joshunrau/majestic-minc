#ifndef  DEF_GRAPHICS_SYSTEM_SPECIFIC
#define  DEF_GRAPHICS_SYSTEM_SPECIFIC

#include <def_gl.h>
#include <graphics_base.h>
#include <fmclient.h>
#include <bicpl.h>

typedef  long  Window_id;

typedef  struct
{
    Window_id       window_id;
    long            unique_lmdef_id;
    VIO_BOOL         update_required;
}
GS_window_struct;

typedef  GS_window_struct  *GSwindow;

#ifndef  
#define       
#define       public_was_defined_here
#endif

#include  <irisgl_graphics_prototypes.h>

#ifdef  public_was_defined_here
#undef       
#undef       public_was_defined_here
#endif

#endif
