 
#include  <volume_io.h>
#include  <GS_graphics.h>

  void  GS_set_shade_model(
    Shading_types  type )
{
    if( type == GOURAUD_SHADING )
        shademodel( GOURAUD );
    else
        shademodel( FLAT );
}

  void  GS_turn_off_blend_function( void )
{
    blendfunction( BF_ONE, BF_ZERO );
}

  void  GS_turn_on_blend_function( void )
{
    blendfunction( BF_SA, BF_MSA );
}

  void  GS_set_lighting_state(
    GSwindow     window,
    VIO_BOOL      state )
{
    if( state )
        lmbind( MATERIAL, (short) window->unique_lmdef_id );
    else
        lmbind( MATERIAL, 0 );
}

  void  GS_backface_culling_state(
    VIO_BOOL     state )
{
    backface( (Boolean) state );
}

  void  GS_set_n_curve_segments(
    int      n_segments )
{
    curveprecision( (short) n_segments );
}
