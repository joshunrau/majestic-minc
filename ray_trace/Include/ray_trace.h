#ifndef  DEF_RAY_TRACE_H
#define  DEF_RAY_TRACE_H

#include  <volume_io.h>
#include  <bicpl.h>
#include  <graphics.h>

#define  N_ATTENUATION   3

typedef struct
{
    Light_types   type;
    VIO_Colour        colour;
    VIO_Point         position;
    VIO_Vector        direction;
    VIO_Real          attenuation[N_ATTENUATION];
} light_struct;

typedef struct
{
    int            n_lights;
    light_struct   *lights;
    VIO_Colour         ambient_light;
} lights_struct;

typedef struct
{
    VIO_Point           eye;
    VIO_Real            stereo_offset;
    VIO_Vector          line_of_sight;
    VIO_Vector          horizontal;
    VIO_Vector          vertical;
    VIO_Vector          perpendicular;
    VIO_Real            window_width;
    VIO_Real            window_height;
    VIO_Real            window_distance;
    VIO_BOOL         perspective_flag;
    int             x_viewport_size;
    int             y_viewport_size;
} view_struct;

typedef struct
{
    VIO_Real                   opacity;
    colour_coding_struct   colour_coding;
    VIO_Volume                 volume;
    int                    continuity;
    VIO_BOOL                fill_value_specified;
    VIO_Real                   fill_value;
    VIO_BOOL                gradient_flag;
    VIO_BOOL                extend_volume_flag;
} volume_info;

typedef  struct
{
    VIO_BOOL                shadow_state;
    VIO_Real                   shadow_offset;
    VIO_BOOL                flat_shading_flag;
    VIO_BOOL                lighting_flag;
    VIO_BOOL                mult_volume_flag;
    VIO_BOOL                composite_volume_flag;
    int                    n_volumes;
    volume_info            *volumes;
    VIO_BOOL                hit_only_once_flag;
}
render_struct;

typedef  struct
{
    int             n_clip_objects;
    object_struct   **clip_objects;
    VIO_BOOL         *clip_inside_flags;
}
clip_struct;

typedef  struct
{
    render_struct      render;
    VIO_BOOL            regular_object_flag;
    object_struct      *object;
    clip_struct        clip;
    VIO_Real               threshold;
    VIO_Colour             colour;
    VIO_Surfprop           spr;
    bitlist_3d_struct  done_bits;
    bitlist_3d_struct  surface_bits;
}
ray_trace_object;

#include  <ray_trace_prototypes.h>

#endif
