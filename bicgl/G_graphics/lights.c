#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H  

#include  <volume_io.h>
#include  <graphics.h>

  void   initialize_graphics_lights(
    Gwindow  window )
{
    set_current_window( window );

    GS_initialize_lights( window->GS_window );
}

  void  G_set_ambient_light(
    Gwindow        window,
    VIO_Colour         colour )
{
    set_current_window( window );

    GS_set_ambient_light( window->GS_window, colour );
}

  void  G_define_light(
    Gwindow         window,
    int             light_index,
    Light_types     type,
    VIO_Colour          colour,
    VIO_Vector          *direction,
    VIO_Point           *position,
    VIO_Real            spot_exponent,
    VIO_Real            spot_angle )
{
    set_current_window( window );

    GS_define_light( window->GS_window, light_index, type, colour,
                     direction, position, spot_exponent, spot_angle );
}

  void  G_set_light_state(
    Gwindow         window,
    int             light_index,
    VIO_BOOL         state )
{
    set_current_window( window );

    if( state )
        G_set_view_type( window, WORLD_VIEW );

    GS_set_light_state( window->GS_window, light_index, state );
}
