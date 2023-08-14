
#include  <volume_io.h>
#include  <GS_graphics.h>

  void  GS_set_point(
    VIO_Point  *point )
{
    v3f( point->coords );
}

  void  GS_set_normal(
    VIO_Vector  *normal )
{
    n3f( normal->coords );
}

  void  GS_set_colour(
    VIO_Colour  colour )
{
    cpack( colour );
}

  void  GS_set_colour_index(
    VIO_Colour  colour )
{
    color( (Colorindex) colour );
}

  void  GS_set_ambient_and_diffuse_mode(
    VIO_BOOL  state )
{
    if( state )
        lmcolor( LMC_AD );
    else
        lmcolor( LMC_COLOR );
}

static  long  get_unique_lmdef_id( void )
{
    static  long   unique_lmdef_id = 1;

    unique_lmdef_id += 8;
    if( unique_lmdef_id > 60000 )
        unique_lmdef_id = 1;

    return( unique_lmdef_id );
}

  void  GS_initialize_surface_property(
    GSwindow        window )
{
    static  Surfprop  spr = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };

    window->unique_lmdef_id = get_unique_lmdef_id();
    GS_set_surface_property( window, WHITE, &spr );
}

  void  GS_set_surface_property(
    GSwindow       window,
    VIO_Colour         col,
    Surfprop       *surfprop )
{
    static  float  props[] = { (float) AMBIENT, 0.0f, 0.0f, 0.0f,
                               (float) DIFFUSE, 0.0f, 0.0f, 0.0f,
                               (float) SPECULAR, 0.0f, 0.0f, 0.0f,
                               (float) SHININESS, 0.0f,
                               (float) ALPHA, 0.0f,
                               (float) LMNULL };

    props[1] = (float) ((VIO_Real) Surfprop_a( *surfprop )*get_Colour_r_0_1( col ));
    props[2] = (float) ((VIO_Real) Surfprop_a( *surfprop )*get_Colour_g_0_1( col ));
    props[3] = (float) ((VIO_Real) Surfprop_a( *surfprop )*get_Colour_b_0_1( col ));

    props[5] = (float) ((VIO_Real) Surfprop_d( *surfprop )*get_Colour_r_0_1( col ));
    props[6] = (float) ((VIO_Real) Surfprop_d( *surfprop )*get_Colour_g_0_1( col ));
    props[7] = (float) ((VIO_Real) Surfprop_d( *surfprop )*get_Colour_b_0_1( col ));

    props[9] = Surfprop_s( *surfprop );
    props[10] = Surfprop_s( *surfprop );
    props[11] = Surfprop_s( *surfprop );

    props[13] = Surfprop_se( *surfprop );

    props[15] = Surfprop_t( *surfprop );

    lmdef( DEFMATERIAL, (short) window->unique_lmdef_id,
           VIO_SIZEOF_STATIC_ARRAY(props), props );
}

  void  GS_set_line_width(
    VIO_Real  line_width )
{
    linewidth( (short) line_width );
}

  void  GS_curve(
    float  geom[4][3] )
{
    crv( geom );
}

  void  GS_begin_point( void )
{
    bgnpoint();
}

  void  GS_end_point( void )
{
    endpoint();
}

  void  GS_begin_line( void )
{
    bgnline();
}

  void  GS_end_line( void )
{
    endline();
}

  void  GS_begin_closed_line( void )
{
    bgnclosedline();
}

  void  GS_end_closed_line( void )
{
    endclosedline();
}

  void  GS_begin_polygon( void )
{
    bgnpolygon();
}

  void  GS_end_polygon( void )
{
    endpolygon();
}
  void  GS_begin_quad_strip( void )
{
    bgnqstrip();
}

  void  GS_end_quad_strip( void )
{
    endqstrip();
}

  void  GS_set_raster_position(
    VIO_Real  x,
    VIO_Real  y,
    VIO_Real  z )
{
    cmov( (Coord) x, (Coord) y, (Coord) z );
}

  void  GS_set_pixel_zoom(
    VIO_Real  x_zoom,
    VIO_Real  y_zoom )
{
    rectzoom( (float) x_zoom, (float) y_zoom );
}

  void  GS_draw_colour_map_pixels(
    int             x_viewport_min,
    int             y_viewport_min,
    pixels_struct   *pixels )
{
    int   x, y, x_size, y_size;

    x = pixels->x_position + x_viewport_min;
    y = pixels->y_position + y_viewport_min;
    x_size = pixels->x_size;
    y_size = pixels->y_size;

    switch( pixels->pixel_type )
    {
    case COLOUR_INDEX_8BIT_PIXEL:
        pixmode( PM_SIZE, 8 );
        lrectwrite( (Screencoord) x, (Screencoord) y,
                    (Screencoord) (x + x_size - 1),
                    (Screencoord) (y + y_size - 1),
                    (unsigned long *) ((void *)
                     (pixels->data.pixels_8bit_colour_index)) );
        break;

    case COLOUR_INDEX_16BIT_PIXEL:
        rectwrite( (Screencoord) x, (Screencoord) y,
                   (Screencoord) (x + x_size - 1),
                   (Screencoord) (y + y_size - 1),
                   (Colorindex *)
                    (pixels->data.pixels_16bit_colour_index) );
        break;
    default:
        break;
    }
}

  void  GS_draw_rgb_pixels(
    int             x_viewport_min,
    int             y_viewport_min,
    pixels_struct   *pixels )
{
    int   x, y, x_size, y_size;

    x = pixels->x_position + x_viewport_min;
    y = pixels->y_position + y_viewport_min;
    x_size = pixels->x_size;
    y_size = pixels->y_size;

    lrectwrite( (Screencoord) x, (Screencoord) y,
                (Screencoord) (x + x_size - 1),
                (Screencoord) (y + y_size - 1),
                pixels->data.pixels_rgb );
}

/* ARGSUSED */

  void  GS_read_pixels(
    VIO_BOOL         colour_map_state,
    int             x_min,
    int             x_max,
    int             y_min,
    int             y_max,
    VIO_Colour          pixels[] )
{
    int   i, n_pixels;

    readsource( SRC_FRONT );

    if( colour_map_state )
    {
        (void) rectread( (Screencoord) x_min, (Screencoord) y_min,
                         (Screencoord) x_max, (Screencoord) y_max,
                         (Colorindex *) pixels );
    }
    else
    {
        (void) lrectread( (Screencoord) x_min, (Screencoord) y_min,
                          (Screencoord) x_max, (Screencoord) y_max, pixels );

        n_pixels = (x_max - x_min + 1) * (y_max - y_min + 1);
        for_less( i, 0, n_pixels )
            pixels[i] |= BLACK;
    }
}

static  const  STRING  font_name = "Helvetica";

static  fmfonthandle   base_font_handle = 0;

static  struct
{
    VIO_Real           size;
    fmfonthandle   font;
}
*fonts = NULL;

static  int  n_fonts = 0;
static  int  current_font_index = 0;

static  VIO_BOOL  fonts_present( void )
{
    return( base_font_handle != (fmfonthandle) 0 );
}

static  fmfonthandle  lookup_font(
    Font_types    font_type,
    VIO_Real          size )
{
    static  VIO_BOOL  first = TRUE;

    if( first )
    {
        first = FALSE;
        fminit();

        base_font_handle = fmfindfont( font_name );
    }

    if( font_type == FIXED_FONT || !fonts_present() )
        return( (fmfonthandle) 0 );

    if( current_font_index <= 0 || current_font_index >= n_fonts ||
        fonts[current_font_index].size != size )
    {
        for_less( current_font_index, 0, n_fonts )
        {
            if( fonts[current_font_index].size == size )
                break;
        }
    }

    if( current_font_index >= n_fonts )
    {
        SET_ARRAY_SIZE( fonts, n_fonts, n_fonts+1, 1 );
        fonts[current_font_index].size = size;
        fonts[current_font_index].font = fmscalefont( base_font_handle, size );
    }

    return( fonts[current_font_index].font );
}

  VIO_Real  GS_get_text_length(
    STRING          str,
    Font_types      font_type,
    VIO_Real            font_size )
{
    fmfonthandle   font;
    VIO_Real           width;

    font = lookup_font( font_type, font_size );

    if( font != 0 )
        width = (VIO_Real) fmgetstrwidth( font, str );
    else
        width = (VIO_Real) strwidth( str );

    return( width );
}

  VIO_Real  GS_get_character_height(
    Font_types      font_type,
    VIO_Real            font_size )
{
    fmfonthandle   font;
    VIO_Real           height;

    font = lookup_font( font_type, font_size );

    if( font != 0 )
        height = font_size;
    else
        height = 18.0;

    return( height );
}

  void   GS_draw_text(
    Font_types   font_type,
    VIO_Real         font_size,
    STRING       string )
{
    fmfonthandle   font;

    if( font_type == SIZED_FONT )
    {
        font = lookup_font( font_type, font_size );
        fmsetfont( font );
        fmprstr( string );
    }
    else
        charstr( string );
}

