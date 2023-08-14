 
#include  <volume_io.h>
#include  <GS_graphics.h>

  void   GS_initialize_lights(
    GSwindow  window )
{
    static   float    data[] = { (float) AMBIENT, 1.0f, 1.0f, 1.0f,
                                 (float) LOCALVIEWER, 0.0f,
                                 (float) TWOSIDE, 1.0f,
                                 (float) LMNULL };

    lmdef( DEFLMODEL, (short) window->unique_lmdef_id,
           VIO_SIZEOF_STATIC_ARRAY(data), data );

    lmbind( LMODEL, (short) window->unique_lmdef_id );

    lmbind( LIGHT0, 0 );
    lmbind( LIGHT1, 0 );
    lmbind( LIGHT2, 0 );
    lmbind( LIGHT3, 0 );
    lmbind( LIGHT4, 0 );
    lmbind( LIGHT5, 0 );
    lmbind( LIGHT6, 0 );
    lmbind( LIGHT7, 0 );
}

  void  GS_set_ambient_light(
    GSwindow       window,
    VIO_Colour         colour )
{
    static   float    data[] = { (float) AMBIENT, 0.0f, 0.0f, 0.0f,
                                 (float) LMNULL };

    data[1] = (float) get_Colour_r_0_1(colour);
    data[2] = (float) get_Colour_g_0_1(colour);
    data[3] = (float) get_Colour_b_0_1(colour);

    lmdef( DEFLMODEL, (short) window->unique_lmdef_id,
           VIO_SIZEOF_STATIC_ARRAY(data), data );
}

  void  GS_define_light(
    GSwindow        window,
    int             light_index,
    Light_types     type,
    VIO_Colour          colour,
    VIO_Vector          *direction,
    VIO_Point           *position,
    VIO_Real            spot_exponent,
    VIO_Real            spot_angle )
{
    VIO_Vector  unit_direction;
    int     n_data;
    float   data[25];

    n_data = 0;

    switch( type )
    {
    case DIRECTIONAL_LIGHT:
        NORMALIZE_VECTOR( unit_direction, *direction );
        data[n_data++] = (float) AMBIENT;
        data[n_data++] = 0.0f;
        data[n_data++] = 0.0f;
        data[n_data++] = 0.0f;
        data[n_data++] = (float) LCOLOR;
        data[n_data++] = (float) get_Colour_r_0_1(colour);
        data[n_data++] = (float) get_Colour_g_0_1(colour);
        data[n_data++] = (float) get_Colour_b_0_1(colour);
        data[n_data++] = (float) POSITION;
        data[n_data++] = -Vector_x(unit_direction);
        data[n_data++] = -Vector_y(unit_direction);
        data[n_data++] = -Vector_z(unit_direction);
        data[n_data++] = 0.0f;
        data[n_data++] = (float) LMNULL;
        break;

    case POINT_LIGHT:
        data[n_data++] = (float) AMBIENT;
        data[n_data++] = 0.0f;
        data[n_data++] = 0.0f;
        data[n_data++] = 0.0f;
        data[n_data++] = (float) LCOLOR;
        data[n_data++] = (float) get_Colour_r_0_1(colour);
        data[n_data++] = (float) get_Colour_g_0_1(colour);
        data[n_data++] = (float) get_Colour_b_0_1(colour);
        data[n_data++] = (float) POSITION;
        data[n_data++] = (float) Point_x(*position);
        data[n_data++] = (float) Point_y(*position);
        data[n_data++] = (float) Point_z(*position);
        data[n_data++] = 1.0f;
        data[n_data++] = (float) LMNULL;
        break;

    case SPOT_LIGHT:
        NORMALIZE_VECTOR( unit_direction, *direction );
        data[n_data++] = (float) AMBIENT;
        data[n_data++] = 0.0f;
        data[n_data++] = 0.0f;
        data[n_data++] = 0.0f;
        data[n_data++] = (float) LCOLOR;
        data[n_data++] = (float) get_Colour_r_0_1(colour);
        data[n_data++] = (float) get_Colour_g_0_1(colour);
        data[n_data++] = (float) get_Colour_b_0_1(colour);
        data[n_data++] = (float) POSITION;
        data[n_data++] = Point_x(*position);
        data[n_data++] = Point_y(*position);
        data[n_data++] = Point_z(*position);
        data[n_data++] = (float) SPOTLIGHT;
        data[n_data++] = (float) spot_exponent;
        data[n_data++] = (float) spot_angle;
        data[n_data++] = (float) SPOTDIRECTION;
        data[n_data++] = Vector_x(unit_direction);
        data[n_data++] = Vector_y(unit_direction);
        data[n_data++] = Vector_z(unit_direction);
        data[n_data++] = (float) LMNULL;
        break;
    }

    lmdef( DEFLIGHT, (short) ((int) window->unique_lmdef_id +
                              light_index),
           (short) n_data, data );
}

  void  GS_set_light_state(
    GSwindow        window,
    int             light_index,
    VIO_BOOL         state )
{
    short   gl_light_index;

    switch( light_index )
    {
    case 0:  gl_light_index = LIGHT0;    break;
    case 1:  gl_light_index = LIGHT1;    break;
    case 2:  gl_light_index = LIGHT2;    break;
    case 3:  gl_light_index = LIGHT3;    break;
    case 4:  gl_light_index = LIGHT4;    break;
    case 5:  gl_light_index = LIGHT5;    break;
    case 6:  gl_light_index = LIGHT6;    break;
    case 7:  gl_light_index = LIGHT7;    break;
    }

    if( state )
        lmbind( gl_light_index,
                (short) ((int)window->unique_lmdef_id+light_index) );
    else
        lmbind( gl_light_index, 0 );
}
