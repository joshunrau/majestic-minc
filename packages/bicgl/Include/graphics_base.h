#ifndef  DEF_GRAPHICS_BASE
#define  DEF_GRAPHICS_BASE

#include  <bicpl.h>

#define  RETURN_KEY    ((char) 13)
#define  DELETE_KEY    ((char) 127)
#define  BACKSPACE_KEY ((char) 8)
#define  ESCAPE_KEY     ((char) 27)

#define  BICGL_INSERT_KEY 228
#define  BICGL_HOME_KEY   229
#define  BICGL_END_KEY    230
#define  BICGL_F1_KEY     231
#define  BICGL_F2_KEY     232
#define  BICGL_F3_KEY     233
#define  BICGL_F4_KEY     234
#define  BICGL_F5_KEY     235
#define  BICGL_F6_KEY     236
#define  BICGL_F7_KEY     237
#define  BICGL_F8_KEY     238
#define  BICGL_F9_KEY     239
#define  BICGL_F10_KEY    240
#define  BICGL_F11_KEY    241
#define  BICGL_F12_KEY    242
#define  BICGL_PGUP_KEY   243
#define  BICGL_PGDN_KEY   244
#define  LEFT_ARROW_KEY   245
#define  RIGHT_ARROW_KEY  246
#define  DOWN_ARROW_KEY   247
#define  UP_ARROW_KEY     248
#define  LEFT_SHIFT_KEY   249
#define  RIGHT_SHIFT_KEY  250
#define  LEFT_CTRL_KEY    251
#define  RIGHT_CTRL_KEY   252
#define  LEFT_ALT_KEY     253
#define  RIGHT_ALT_KEY    254

#define  SHIFT_KEY_BIT    (1 << 0)
#define  CTRL_KEY_BIT     (1 << 1)
#define  ALT_KEY_BIT      (1 << 2)

typedef  enum  {
                 NORMAL_PLANES,
                 OVERLAY_PLANES,
                 N_BITPLANE_TYPES } Bitplane_types;

typedef  enum  { PIXEL_VIEW, SCREEN_VIEW, WORLD_VIEW, MODEL_VIEW,
                 N_VIEW_TYPES }     View_types;

typedef  enum  { WIREFRAME, FILLED, POINT, OVERLAY } Shading_modes;

typedef  enum  { FLAT_SHADING, GOURAUD_SHADING }  Shading_types;

typedef enum {
  NO_TRANSPARENCY, NORMAL_TRANSPARENCY, ADDITIVE_TRANSPARENCY
} Transparency_modes;

#define  N_LIGHTS  8

typedef  enum { 
                DIRECTIONAL_LIGHT,
                POINT_LIGHT,
                SPOT_LIGHT
              } Light_types;

typedef  enum { VIEWING_MATRIX, PROJECTION_MATRIX } Matrix_modes;

typedef  enum { LESS_OR_EQUAL } Depth_functions;



#define  NO_BELL       0
#define  SHORT_BELL    1
#define  LONG_BELL     2

#endif
