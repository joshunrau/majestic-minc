#ifndef  DEF_GRAPHICS_STRUCT
#define  DEF_GRAPHICS_STRUCT

#include  <volume_io/basic.h>
#include  <graphics.h>

typedef  struct
{
    int             n_objects;
    object_struct   **objects;
    VIO_BOOL         update_flag[2];    /* DOUBLE BUFFER */
} bitplane_graphics_struct;

typedef struct
{
    VIO_Colour                      background_colour;
    int                         background_colour_index;
    int                         x_min;
    int                         x_max;
    int                         y_min;
    int                         y_max;
    bitplane_graphics_struct    bitplanes[N_BITPLANE_TYPES];
} viewport_struct;

typedef  struct
{
    VIO_BOOL           clear_bitplane_flags[N_BITPLANE_TYPES][2]; /* doub buff */
    int               n_viewports;
    viewport_struct   *viewports;
} graphics_struct;

#endif
