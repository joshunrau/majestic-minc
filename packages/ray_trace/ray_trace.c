#if HAVE_CONFIG_H
#include "config.h"
#endif
#include  <ray_trace.h>
#include  <bicpl/images.h>

#define  GRAY_STRING                "-gray"
#define  HOT_STRING                 "-hot"
#define  HOT_NEG_STRING             "-hot_inv"
#define  COLD_METAL_STRING          "-cold_metal"
#define  COLD_METAL_NEG_STRING      "-cold_metal_inv"
#define  GREEN_METAL_STRING         "-green_metal"
#define  GREEN_METAL_NEG_STRING     "-green_metal_inv"
#define  LIME_METAL_STRING          "-lime_metal"
#define  LIME_METAL_NEG_STRING      "-lime_metal_inv"
#define  RED_METAL_STRING           "-red_metal"
#define  RED_METAL_NEG_STRING       "-red_metal_inv"
#define  PURPLE_METAL_STRING        "-purple_metal"
#define  PURPLE_METAL_NEG_STRING    "-purple_metal_inv"
#define  SPECTRAL_STRING            "-spectral"
#define  RED_STRING                 "-red"
#define  GREEN_STRING               "-green"
#define  BLUE_STRING                "-blue"
#define  SINGLE_COLOUR_STRING       "-colour"
#define  USER_DEF_STRING            "-usercc"

#define  BINTREE_FACTOR  0.4

#define  FIT_FACTOR      1.0

#define  CROP_BORDER     1

#define  MAX_DEPTH       4

#define  DEFAULT_SHADOW_OFFSET  0.05

#define  DEFAULT_EYE_SEPARATION  7.0

static  VIO_Status  load_file(
    VIO_STR              filename,
    VIO_General_transform   *transform,
    VIO_Real                bintree_factor,
    VIO_Real                line_width_override,
    VIO_Colour              current_marker_colour,
    VIO_Real                current_marker_size,
    Marker_types        current_marker_type,
    int                 *n_objects,
    object_struct       ***objects );

static  void  define_view(
    view_struct   *view,
    VIO_BOOL       perspective_flag,
    VIO_BOOL       eye_specified,
    VIO_Point         *eye,
    VIO_BOOL       view_specified,
    VIO_Vector        *line_of_sight,
    VIO_Vector        *up_dir,
    VIO_Real          specified_window_width,
    VIO_Real          specified_perspective_distance,
    int           which_eye,
    VIO_Real          stereo_eye_offset,
    int           x_size,
    int           y_size,
    VIO_Real          pixel_aspect,
    VIO_Point         *min_range,
    VIO_Point         *max_range,
    VIO_Real          fit_factor );
static  void  define_default_lights(
    lights_struct   *lights );
static  void  orient_lights(
    lights_struct   *lights,
    view_struct     *view );
  VIO_Status  output_rgb_file(
    VIO_STR          filename,
    pixels_struct   *pixels );
static  void  prepare_objects(
    int                 n_objects,
    object_struct       *objects[],
    VIO_General_transform   *transform,
    VIO_Real                bintree_factor,
    VIO_Real                line_width_override );

static  void   get_preset_view(
    VIO_STR   string,
    VIO_Vector   *line_of_sight,
    VIO_Vector   *up_dir );

static  VIO_BOOL   string_is_preset_view(
    VIO_STR    string );

static  VIO_Real  get_Colour_intensity(
    VIO_Colour   col )
{
    static   VIO_Real  coefs[] = { 0.3, 0.59, 0.11 };

    return( coefs[0] * get_Colour_r_0_1(col) +
            coefs[1] * get_Colour_g_0_1(col) +
            coefs[2] * get_Colour_b_0_1(col) );
}

static  void  get_current_transform(
    VIO_BOOL            camera_space_flag,
    VIO_BOOL            eye_specified,
    VIO_Point              *eye,
    VIO_BOOL            view_specified,
    VIO_Vector             *line_of_sight,
    VIO_Vector             *up_dir,
    VIO_General_transform  *current_transform,
    VIO_General_transform  *used_transform );

typedef  struct
{
    VIO_STR              filename;
    VIO_General_transform   transform;
    VIO_Real                bintree_factor;
    VIO_Real                line_width_override;
    int                 n_objects;
    object_struct       **objects;
} file_lookup;

static  void  print_usage( char   executable_name[] ) {

    VIO_STR usage = "\n\
Usage: %s   ...  <object1.obj> ... <object2.obj> ... \n\
\n\
-help | --help | -h | help \n\
            : print this message and quit \n\
\n\
Options are:\n\
\n\
-size <x_size> <y_size> \n\
            : image size (default 300 300) \n\
-window_width <width> \n\
            : width of camera view (default 3) \n\
-aspect <aspect> \n\
            : aspect of image (default 1.0) \n\
-scale <factor> \n\
            : scaling factor (default 1.0) \n\
-output <filename> \n\
            : output filename for image (default image.rgb) \n\
\n\
-light | -nolight \n\
            : turn on and off lighting effects for subsequent objects (default -light) \n\
-ambient <r> <g> <b> \n\
            : ambient light (default 0.3 0.3 0.3) \n\
-directional <x> <y> <z> <r> <g> <b> \n\
            : define the direction x y z of the light and its rgb components \n\
              (default ?) \n\
-point <nx> <ny> <nz> <r> <g> <b> <ax> <ay> <az> \n\
            : define the direction x y z of the light and its rgb components, \n\
              with attenuation factors ax ay az (default ?) \n\
-persp | -ortho \n\
            : perspective or orthographic view (default -persp) \n\
-perspective_distance <dist> \n\
            : distance from viewpoint in perspective view (default ?) \n\
-camera | -model \n\
            : ? (default -model) \n\
-left_eye | -right_eye \n\
            : look from left or right eye (default is center) \n\
-eye_separation <offset> \n\
            : separation distance between eyes (default 7.0) \n\
-eye <nx> <ny> <nz> \n\
            : eyepoint oriention (default ?) \n\
-shadows | -noshadows \n\
            : enable shadows or not (default -noshadows) \n\
-shadow_offset <offset> \n\
            : offset for shadows (default 0.05) \n\
-normals | -nonormals \n\
            : enable calculation of normal vectors or not (default -normals) \n\
-smooth_normals <n_iters> <ratio> \n\
            : smooth the surface normals (default 0 0.8) \n\
-flat | -smooth \n\
            : use flat or smooth shading for polygons (default -smooth) \n\
\n\
-gray              <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-hot               <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-hot_inv           <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-cold_metal        <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-cold_metal_inv    <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-green_metal       <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-green_metal_inv   <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-lime_metal        <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-lime_metal_inv    <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-red_metal         <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-red_metal_inv     <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-purple_metal      <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-purple_metal_inv  <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-spectral          <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-red               <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-green             <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-blue              <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-colour            <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
-usercc            <user_map.rgb>  <min>  <max>  [grad] <volume.mnc>  <interpolation> <opacity> \n\
 \n\
            : colour subsequent objects on the command line with this \n\
              colour coding method.  The values min and max define the \n\
              range of the colour bar.  The interpolation is either -1 \n\
              for nearest neighbour, 0 for trilinear, or 2 for tri-cubic. \n\
              The opacity is a value between 0 and 1 where 1 is fully \n\
              opaque.  If -colour is used, then the colour coding used is \n\
              a scale from black through to the current over_colour. (default none) \n\
 \n\
-reverse_order_colouring | -noreverse_order_colouring \n\
            : colour the volumes in input or reverse input order \n\
              (default -noreverse_order_colouring) \n\
-volume <threshold> { <colour_name> | '1 0 0' | '1 0 0 1' } \n\
            : threshold and colour for volume (default none) \n\
-fill_value <val> | -nofill_value \n\
            : fill value for volume (default none) \n\
-under { <colour_name> | '1 0 0' | '1 0 0 1' } | -over { <colour_name> | '1 0 0' | '1 0 0 1' } \n\
            : set the colour for volume values above or below the min,max \n\
              colour range.  VIO_Colour names include red, green, blue, \n\
              yellow, purple, brown, etc. or transparent (default none) \n\
-bg { <colour_name> | '1 0 0' | '1 0 0 1' } \n\
            : set the background colour (default 0.3 0.3 0.3) \n\
-behind { <colour_name> | '1 0 0' | '1 0 0 1' } \n\
            : (default none) \n\
\n\
-crop | -nocrop \n\
            : crop or do not crop the image (default -nocrop) \n\
\n\
-marker_colour { <colour_name> | '1 0 0' | '1 0 0 1' } \n\
            : colour for the markers (tag points) (default green) \n\
-marker_sphere | -marker_box \n\
            : symbol for the markers (default -marker_sphere) \n\
-marker_size <size> \n\
            : size of the marker symbol (default 1) \n\
-line_width <width> \n\
            : line thickness for drawing line objects (default ?) \n\
\n\
-top | -bottom | -left | -right | -back | -front \n\
            : view orientation for surfaces (default - see -view)\n\
-view <nx> <ny> <nz> <x_up> <y_up> <z_up> \n\
            : user-defined view orientation (default 0.77 -0.18 -0.6 0.55 0.6 0.55) \n\
-clip off | { { inside | outside } <file.obj> } \n\
            : no clip or clip to inside or outside a given object (default off) \n\
\n\
-sup <sampling> \n\
            : super sampling factor (default 1) \n\
-delete_volume <n> \n\
            : delete n_th volume (default none) \n\
-centre <x> <y> <z> \n\
            : centre of rotation for transformation (default none) \n\
-transform [[inverse] <trans.xfm> | identity | xrot <angle> | yrot <angle> | zrot <angle> | \n\
            rot <nx> <ny> <nz> <angle> ] \n\
            : apply a transformation file or a simple transformation (default none)\n\
\n\
Other (undocumented) parameters: \n\
\n\
-one_hit | -multi_hit \n\
            : (default -multi_hit) \n\
-bintree <factor> (default 0.4)\n\
-mult | -nomult \n\
            : (default -nomult) \n\
-composite | -nocomposite \n\
            : (default -nocomposite) \n\
-extend_volume | -no_extend_volume \n\
            : (default -no_extend_volume) \n\
-text | -grey | -dither4 | -dither8 | -rgb \n\
            : text = use 0 bit representation for colours \n\
            : grey = use 1 bit representation for colours \n\
            : dither4 = use 4 bits representation for colours \n\
            : dither8 = use 8 bits representation for colours \n\
            : rgb = use 24 bits representation for colours (default)\n\
\n\
(See http://www.bic.mni.mcgill.ca/~david/Ray_trace/ray_trace_tutorial.html for\n\
a tutorial introduction.) \n\n";

    print_error( usage, executable_name );
}

int  main(
    int   argc,
    char  *argv[] )
{
    FILE                 *file;
    VIO_BOOL              crop_flag, clip_inside_flag;
    VIO_BOOL              over_colour_set, compute_normals;
    VIO_BOOL              reverse_order_colouring_flag;
    VIO_BOOL              perspective_flag, extend_volume_flag;
    VIO_STR               dimension_names[] = { ANY_SPATIAL_DIMENSION,
                                               ANY_SPATIAL_DIMENSION,
                                               ANY_SPATIAL_DIMENSION };
    VIO_STR               volume_filename, output_file;
    VIO_STR               argument, user_def_filename;
    VIO_STR               transform_filename;
    int                  i, n_objects, x, y, x_size, y_size, n_colour_comps;
    int                  n_objects_in_file, super_sampling;
    int                  sizes[VIO_MAX_DIMENSIONS], n_lights;
    int                  axis, which;
    VIO_Point                min_pos, max_pos, min_range, max_range;
    VIO_Real                 intensity, aspect, fit_factor, threshold;
    VIO_Real                 fill_value;
    VIO_BOOL              fill_value_specified;
    VIO_Real                 angle, shadow_offset;
    VIO_Transform            rotation, trans;
    VIO_Real                 x_centre, y_centre, z_centre;
    VIO_Real                 window_width, perspective_distance;
    VIO_Point                centre_of_rotation;
    VIO_General_transform    transform, new_transform, specified_transform;
    VIO_General_transform    used_transform;
    object_struct        **file_objects, *object;
    ray_trace_object     *object_list, ray_object;
    object_traverse_struct  object_traverse;
    view_struct          view;
    lights_struct        lights;
    pixels_struct        pixels, pixels_tmp;
    Colour_coding_types  coding_type;
    VIO_Real                 r, g, b, a, xv, yv, zv, a1, a2, a3;
    VIO_Real                 colour_coding_min, colour_coding_max;
    VIO_Real                 stereo_eye_offset;
    int                  which_eye;
    VIO_BOOL              behind_specified;
    VIO_Colour               background_colour, colour, under_colour, over_colour;
    VIO_Colour               behind_transparency_colour, alternate_behind, col;
    VIO_Colour               this_over_colour;
    VIO_STR               colour_name;
    VIO_Real                 x_line_of_sight, y_line_of_sight, z_line_of_sight;
    VIO_Real                 x_up, y_up, z_up, bintree_factor;
    VIO_Real                 line_width_override;
    VIO_Real                 normal_ratio;
    int                  n_normal_iters, rc, bc, gc, ac;
    VIO_Vector               line_of_sight, up_dir, rot_axis;
    VIO_Real                 x_rot, y_rot, z_rot;
    VIO_BOOL              view_specified, eye_specified, invert_flag;
    VIO_Point                eye;
    VIO_Colour               current_marker_colour;
    VIO_Real                 current_marker_size;
    Marker_types         current_marker_type;
    render_struct        render;
    minc_input_options   options;
    clip_struct          clip;
    static int           colour_comps16[16][3] = {
                                           0,   0,   0,
                                           0,   0, 255,
                                           0, 255,   0,
                                           0, 255, 255,
                                         255,   0,   0,
                                         255,   0, 255,
                                         255, 255,   0,
                                         255, 255, 255,
                                           0, 128, 128,
                                         255, 128, 128,
                                         128,   0, 128,
                                         128, 255, 128,
                                         128, 128,   0,
                                         128, 128, 255,
                                          85,  85,  85,
                                         170, 170, 170
                                    };

    VIO_Colour              colours16[16];
    volume_info         vol_info;
    VIO_BOOL             camera_space_flag;

    initialize_argument_processing( argc, argv );

    for_less( i, 0, 16 )
    {
        colours16[i] = make_Colour( colour_comps16[i][0],
                                    colour_comps16[i][1],
                                    colour_comps16[i][2] );
    }

    view_specified = FALSE;
    eye_specified = FALSE;
    fit_factor = FIT_FACTOR;

    fill_value_specified = FALSE;
    fill_value = 0.0;

    x_size = 300;
    y_size = 300;
    n_objects = 0;
    object_list = NULL;
    fill_Point( min_range, 0.0, 0.0, 0.0 );
    fill_Point( max_range, 0.0, 0.0, 0.0 );
    crop_flag = FALSE;
    perspective_flag = TRUE;
    aspect = 1.0;
    n_colour_comps = 24;
    super_sampling = 1;
    clip.n_clip_objects = 0;
    render.shadow_state = FALSE;
    render.shadow_offset = DEFAULT_SHADOW_OFFSET;
    render.flat_shading_flag = FALSE;
    render.lighting_flag = TRUE;
    render.n_volumes = 0;
    render.volumes = NULL;
    render.mult_volume_flag = FALSE;
    render.composite_volume_flag = FALSE;
    render.hit_only_once_flag = FALSE;

    reverse_order_colouring_flag = FALSE;

    n_normal_iters = 0;
    normal_ratio = 0.8;

    extend_volume_flag = FALSE;

    current_marker_colour = GREEN;
    current_marker_size = 1.0;
    current_marker_type = SPHERE_MARKER;

    output_file = "image.rgb";
    background_colour = make_Colour_0_1( 0.3, 0.3, 0.3 );
    behind_specified = FALSE;
    bintree_factor = BINTREE_FACTOR;
    invert_flag = FALSE;
    create_linear_transform( &transform, NULL );
    over_colour_set = FALSE;
    under_colour = BLACK;
    over_colour = WHITE;
    perspective_distance = -1.0;
    window_width = -1.0;
    line_width_override = -1.0;
    compute_normals = TRUE;
    camera_space_flag = FALSE;
    stereo_eye_offset = DEFAULT_EYE_SEPARATION;
    which_eye = 0;

    set_default_minc_input_options( &options );
    set_minc_input_vector_to_colour_flag( &options, TRUE );

    define_default_lights( &lights );
    n_lights = 0;

    while( get_string_argument( NULL, &argument ) )
    {
        if( equal_strings( argument, "-help" )  ||
            equal_strings( argument, "-h" )  ||
            equal_strings( argument, "--help" )  ||
            equal_strings( argument, "help" ) )
        {
            print_usage( argv[0] );
            return( 0 );
        }
        else if( equal_strings( argument, GRAY_STRING )  ||
            equal_strings( argument, HOT_STRING )  ||
            equal_strings( argument, HOT_NEG_STRING ) ||
            equal_strings( argument, COLD_METAL_STRING ) ||
            equal_strings( argument, COLD_METAL_NEG_STRING ) ||
            equal_strings( argument, GREEN_METAL_STRING ) ||
            equal_strings( argument, GREEN_METAL_NEG_STRING ) ||
            equal_strings( argument, LIME_METAL_STRING ) ||
            equal_strings( argument, LIME_METAL_NEG_STRING ) ||
            equal_strings( argument, RED_METAL_STRING ) ||
            equal_strings( argument, RED_METAL_NEG_STRING ) ||
            equal_strings( argument, PURPLE_METAL_STRING ) ||
            equal_strings( argument, PURPLE_METAL_NEG_STRING ) ||
            equal_strings( argument, SPECTRAL_STRING )  ||
            equal_strings( argument, RED_STRING )  ||
            equal_strings( argument, GREEN_STRING )  ||
            equal_strings( argument, BLUE_STRING ) ||
            equal_strings( argument, SINGLE_COLOUR_STRING ) ||
            equal_strings( argument, USER_DEF_STRING ) )
        {
            this_over_colour = WHITE;
            if( equal_strings( argument, GRAY_STRING ) )
                coding_type = GRAY_SCALE;
            else if( equal_strings( argument, SINGLE_COLOUR_STRING ) )
                coding_type = SINGLE_COLOUR_SCALE;
            else if( equal_strings( argument, HOT_STRING ) )
                coding_type = HOT_METAL;
            else if( equal_strings( argument, HOT_NEG_STRING ) )
                coding_type = HOT_METAL_NEG;
            else if( equal_strings( argument, COLD_METAL_STRING ) )
                coding_type = COLD_METAL;
            else if( equal_strings( argument, COLD_METAL_NEG_STRING ) )
                coding_type = COLD_METAL_NEG;
            else if( equal_strings( argument, GREEN_METAL_STRING ) )
                coding_type = GREEN_METAL;
            else if( equal_strings( argument, GREEN_METAL_NEG_STRING ) )
                coding_type = GREEN_METAL_NEG;
            else if( equal_strings( argument, LIME_METAL_STRING ) )
                coding_type = LIME_METAL;
            else if( equal_strings( argument, LIME_METAL_NEG_STRING ) )
                coding_type = LIME_METAL_NEG;
            else if( equal_strings( argument, RED_METAL_STRING ) )
                coding_type = RED_METAL;
            else if( equal_strings( argument, RED_METAL_NEG_STRING ) )
                coding_type = RED_METAL_NEG;
            else if( equal_strings( argument, PURPLE_METAL_STRING ) )
                coding_type = PURPLE_METAL;
            else if( equal_strings( argument, PURPLE_METAL_NEG_STRING ) )
                coding_type = PURPLE_METAL_NEG;
            else if( equal_strings( argument, SPECTRAL_STRING ) )
                coding_type = SPECTRAL;
            else if( equal_strings( argument, USER_DEF_STRING ) )
            {
                coding_type = USER_DEFINED_COLOUR_MAP;
                if( !get_string_argument( NULL, &user_def_filename ) )
                {
                    print_error( "Error in user def colour coding argument.\n");
                    print_usage( argv[0] );
                    return( 1 );
                }
            }
            else
            {
                if( equal_strings( argument, RED_STRING ) )
                {
                    coding_type = RED_COLOUR_MAP;
                    this_over_colour = RED;
                }
                else if( equal_strings( argument, GREEN_STRING ) )
                {
                    coding_type = GREEN_COLOUR_MAP;
                    this_over_colour = GREEN;
                }
                else if( equal_strings( argument, BLUE_STRING ) )
                {
                    coding_type = BLUE_COLOUR_MAP;
                    this_over_colour = BLUE;
                }
            }

            if( !get_real_argument( 0.0, &colour_coding_min ) ||
                !get_real_argument( 0.0, &colour_coding_max ) ||
                !get_string_argument( NULL, &volume_filename ) )
            {
                print_error( "Error in colour coding argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            if( equal_strings( volume_filename, "grad" ) )
            {
                if( !get_string_argument( "", &volume_filename ) )
                {
                    print_error( "Error in colour coding argument.\n" );
                    print_usage( argv[0] );
                    return( 1 );
                }
                vol_info.gradient_flag = TRUE;
            }
            else
                vol_info.gradient_flag = FALSE;

            vol_info.extend_volume_flag = extend_volume_flag;
            
            if( !get_int_argument( 0, &vol_info.continuity ) ||
                !get_real_argument( 0.0, &vol_info.opacity ) )
            {
                print_error( "Error in colour coding argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            if( input_volume( volume_filename, 3, dimension_names,
                              NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                              TRUE, &vol_info.volume,
                              &options ) != VIO_OK )
                return( 1 );

            get_current_transform( camera_space_flag, eye_specified, &eye,
                                   view_specified, &line_of_sight, &up_dir,
                                   &transform, &used_transform );

            concat_general_transforms(
                           get_voxel_to_world_transform(vol_info.volume),
                           &used_transform, &new_transform );

            delete_general_transform( &used_transform );

            set_voxel_to_world_transform( vol_info.volume, &new_transform );

            if( !over_colour_set )
                over_colour = this_over_colour;

            initialize_colour_coding( &vol_info.colour_coding,
                                      coding_type, under_colour, over_colour,
                                      colour_coding_min, colour_coding_max );

            if( coding_type == USER_DEFINED_COLOUR_MAP )
            {
                if( input_user_defined_colour_coding( &vol_info.colour_coding,
                                                      user_def_filename ) != VIO_OK)
                {
                    print_error( "Error in user defined colour map: %s\n",
                                 user_def_filename );
                }
            }

            vol_info.fill_value_specified = fill_value_specified;
            vol_info.fill_value = fill_value;

            if( reverse_order_colouring_flag )
            {
                SET_ARRAY_SIZE( render.volumes, render.n_volumes,
                                render.n_volumes+1, 1 );
                for_down( i, render.n_volumes, 1 )
                {
                    render.volumes[i] = render.volumes[i-1];
                }
                render.volumes[0] = vol_info;
                ++render.n_volumes;
            }
            else
            {
                ADD_ELEMENT_TO_ARRAY( render.volumes, render.n_volumes,
                                      vol_info, 1 );
            }
        }
        else if( equal_strings( argument, "-delete_volume" ) )
        {
            if( !get_int_argument( 0, &which ) ||
                which < 0 || which >= render.n_volumes )
            {
                print_error( "Error in -delete_volume arguments.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
            DELETE_ELEMENT_FROM_ARRAY( render.volumes, render.n_volumes,
                                       which, 1 );
        }
        else if( equal_strings( argument, "-centre" ) )
        {
            if( !get_real_argument( 0.0, &x_centre ) ||
                !get_real_argument( 0.0, &y_centre ) ||
                !get_real_argument( 0.0, &z_centre ) )
            {
                print_error( "Error in centre arguments.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
            fill_Point( centre_of_rotation, x_centre, y_centre, z_centre );
        }
        else if( equal_strings( argument, "-fill_value" ) )
        {
            if( !get_real_argument( 0.0, &fill_value ) )
            {
                print_error( "Error in -fill_value arguments.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
            fill_value_specified = TRUE;
        }
        else if( equal_strings( argument, "-nofill_value" ) )
        {
            fill_value_specified = FALSE;
        }
        else if( equal_strings( argument, "-transform" ) )
        {
            if( !get_string_argument( "", &transform_filename ) )
            {
                print_error( "Missing transform filename argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            if( equal_strings( transform_filename, "inverse" ) )
            {
                if( !get_string_argument( "", &transform_filename ) )
                {
                    print_error( "Missing transform filename argument.\n" );
                    print_usage( argv[0] );
                    return( 1 );
                }

                invert_flag = TRUE;
            }
            else
                invert_flag = FALSE;

            if( equal_strings( transform_filename, "identity" ) )
            {
                delete_general_transform( &transform );
                create_linear_transform( &transform, NULL );
            }
            else
            {
                if( equal_strings( transform_filename, "xrot" ) ||
                    equal_strings( transform_filename, "yrot" ) ||
                    equal_strings( transform_filename, "zrot" ) )
                {
                    if( !get_real_argument( 0.0, &angle ) )
                    {
                        print_error( "Missing rotation angle argument.\n" );
                        print_usage( argv[0] );
                        return( 1 );
                    }

                    axis = transform_filename[0] - 'x';

                    make_rotation_transform( angle * VIO_DEG_TO_RAD, axis,
                                             &rotation );
                    make_transform_relative_to_point( &centre_of_rotation,
                                                      &rotation, &trans );

                    create_linear_transform( &specified_transform, &trans );
                }
                else if( equal_strings( transform_filename, "rot" ) )
                {
                    if( !get_real_argument( 0.0, &x_rot ) ||
                        !get_real_argument( 0.0, &y_rot ) ||
                        !get_real_argument( 0.0, &z_rot ) ||
                        !get_real_argument( 0.0, &angle ) )
                    {
                        print_error(
                           "Missing rotation axis and angle arguments.\n" );
                        print_usage( argv[0] );
                        return( 1 );
                    }

                    fill_Vector( rot_axis, x_rot, y_rot, z_rot );
                    make_rotation_about_axis( &rot_axis, angle * VIO_DEG_TO_RAD,
                                              &rotation );

                    make_transform_relative_to_point( &centre_of_rotation,
                                                      &rotation, &trans );

                    create_linear_transform( &specified_transform, &trans );
                }
                else
                {
                    if( input_transform_file( transform_filename,
                                              &specified_transform ) !=VIO_OK )
                        return( 1 );
                }

                if( invert_flag )
                {
                    create_inverse_general_transform( &specified_transform,
                                                      &new_transform );
                    delete_general_transform( &specified_transform );
                    specified_transform = new_transform;
                }

                concat_general_transforms( &transform, &specified_transform,
                                           &new_transform );

                delete_general_transform( &specified_transform );
                delete_general_transform( &transform );
                transform = new_transform;
            }
        }
        else if( equal_strings( argument, "-size" ) )
        {
            if( !get_int_argument( 0, &x_size ) ||
                !get_int_argument( 0, &y_size ) ||
                x_size < 1 || y_size < 1 )
            {
                print_error( "Invalid size arguments.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-left_eye" ) )
            which_eye = -1;
        else if( equal_strings( argument, "-right_eye" ) )
            which_eye = 1;
        else if( equal_strings( argument, "-eye_separation" ) )
        {
            if( !get_real_argument( 0.0, &stereo_eye_offset ) )
            {
                print_error( "Error in eye separation argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-marker_sphere" ) )
            current_marker_type = SPHERE_MARKER;
        else if( equal_strings( argument, "-marker_box" ) )
            current_marker_type = BOX_MARKER;
        else if( equal_strings( argument, "-one_hit" ) )
            render.hit_only_once_flag = TRUE;
        else if( equal_strings( argument, "-multi_hit" ) )
            render.hit_only_once_flag = FALSE;
        else if( equal_strings( argument, "-marker_size" ) )
        {
            if( !get_real_argument( 0.0, &current_marker_size ) )
            {
                print_error( "Missing marker size argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-marker_colour" ) )
        {
            VIO_STR  colour_name;

            if( !get_string_argument( NULL, &colour_name ) )
            {
                print_error( "Missing marker colour.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            current_marker_colour = convert_string_to_colour( colour_name );
        }
        else if( equal_strings( argument, "-aspect" ) )
        {
            if( !get_real_argument( 0.0, &aspect ) )
            {
                print_error( "Missing aspect argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-line_width" ) )
        {
            if( !get_real_argument( 0.0, &line_width_override ) )
            {
                print_error( "Missing line width override argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-bintree" ) )
        {
            if( !get_real_argument( 0.0, &bintree_factor ) )
            {
                print_error( "Missing bintree factor argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-sup" ) )
        {
            if( !get_int_argument( 0, &super_sampling ) )
            {
                print_error( "Missing super sampling argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-output" ) )
        {
            if( !get_string_argument( NULL, &output_file ) )
            {
                print_error( "Missing output filename.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-under" ) )
        {
            if( !get_string_argument( NULL, &colour_name ) )
            {
                print_error( "Missing under colour.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
            under_colour = convert_string_to_colour( colour_name );
        }
        else if( equal_strings( argument, "-over" ) )
        {
            if( !get_string_argument( NULL, &colour_name ) )
            {
                print_error( "Missing over colour.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
            over_colour = convert_string_to_colour( colour_name );
            over_colour_set = TRUE;
        }
        else if( equal_strings( argument, "-reverse_order_colouring" ) )
            reverse_order_colouring_flag = TRUE;
        else if( equal_strings( argument, "-noreverse_order_colouring" ) )
            reverse_order_colouring_flag = FALSE;
        else if( equal_strings( argument, "-mult" ) )
            render.mult_volume_flag = TRUE;
        else if( equal_strings( argument, "-nomult" ) )
            render.mult_volume_flag = FALSE;
        else if( equal_strings( argument, "-composite" ) )
            render.composite_volume_flag = TRUE;
        else if( equal_strings( argument, "-nocomposite" ) )
            render.composite_volume_flag = FALSE;
        else if( equal_strings( argument, "-bg" ) )
        {
            VIO_STR  colour_name;

            if( !get_string_argument( NULL, &colour_name ) )
            {
                print_error( "Missing background colour.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            background_colour = convert_string_to_colour( colour_name );
        }
        else if( equal_strings( argument, "-behind" ) )
        {
            VIO_STR  colour_name;

            if( !get_string_argument( NULL, &colour_name ) )
            {
                print_error( "Missing behind_transparency colour.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            behind_transparency_colour = convert_string_to_colour(colour_name);
            behind_specified = TRUE;
        }
        else if( equal_strings( argument, "-flat" ) )
            render.flat_shading_flag = TRUE;
        else if( equal_strings( argument, "-smooth" ) )
            render.flat_shading_flag = FALSE;
        else if( equal_strings( argument, "-shadows" ) )
            render.shadow_state = TRUE;
        else if( equal_strings( argument, "-noshadows" ) )
            render.shadow_state = FALSE;
        else if( equal_strings( argument, "-nonormals" ) )
            compute_normals = FALSE;
        else if( equal_strings( argument, "-normals" ) )
            compute_normals = TRUE;
        else if( equal_strings( argument, "-smooth_normals" ) )
        {
            if( !get_int_argument( 0, &n_normal_iters ) ||
                !get_real_argument( 0.0, &normal_ratio ) )
            {
                print_error( "Error in smooth normals argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-shadow_offset" ) )
        {
            if( !get_real_argument( 0.0, &shadow_offset ) )
            {
                print_error( "Missing shadow offset argument.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            render.shadow_offset = shadow_offset;
        }
        else if( equal_strings( argument, "-text" ) )
            n_colour_comps = 0;
        else if( equal_strings( argument, "-grey" ) )
            n_colour_comps = 1;
        else if( equal_strings( argument, "-dither4" ) )
            n_colour_comps = 4;
        else if( equal_strings( argument, "-dither8" ) )
            n_colour_comps = 8;
        else if( equal_strings( argument, "-rgb" ) )
            n_colour_comps = 24;
        else if( equal_strings( argument, "-light" ) )
            render.lighting_flag = TRUE;
        else if( equal_strings( argument, "-nolight" ) )
            render.lighting_flag = FALSE;
        else if( equal_strings( argument, "-extend_volume" ) )
            extend_volume_flag = TRUE;
        else if( equal_strings( argument, "-no_extend_volume" ) )
            extend_volume_flag = FALSE;
        else if( equal_strings( argument, "-crop" ) )
            crop_flag = TRUE;
        else if( equal_strings( argument, "-nocrop" ) )
            crop_flag = FALSE;
        else if( equal_strings( argument, "-persp" ) )
            perspective_flag = TRUE;
        else if( equal_strings( argument, "-ortho" ) )
            perspective_flag = FALSE;
        else if( equal_strings( argument, "-view" ) )
        {
            if( !get_real_argument( 0.0, &x_line_of_sight ) ||
                !get_real_argument( 0.0, &y_line_of_sight ) ||
                !get_real_argument( 0.0, &z_line_of_sight ) ||
                !get_real_argument( 0.0, &x_up ) ||
                !get_real_argument( 0.0, &y_up ) ||
                !get_real_argument( 0.0, &z_up ) )
            {
                print_error( "Incorrect view specification.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            view_specified = TRUE;
            fill_Point( line_of_sight, x_line_of_sight, y_line_of_sight,
                                       z_line_of_sight );
            fill_Point( up_dir, x_up, y_up, z_up );
        }
        else if( equal_strings( argument, "-camera" ) )
        {
            camera_space_flag = TRUE;
        }
        else if( equal_strings( argument, "-model" ) )
        {
            camera_space_flag = FALSE;
        }
        else if( string_is_preset_view( argument ) )
        {
            view_specified = TRUE;
            get_preset_view( argument, &line_of_sight, &up_dir );
        }
        else if( equal_strings( argument, "-eye" ) )
        {
            if( !get_real_argument( 0.0, &xv ) ||
                !get_real_argument( 0.0, &yv ) ||
                !get_real_argument( 0.0, &zv ) )
            {
                print_error( "Incorrect eyepoint specification.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            eye_specified = TRUE;
            fill_Point( eye, xv, yv, zv );
        }
        else if( equal_strings( argument, "-window_width" ) )
        {
            if( !get_real_argument( 0.0, &window_width ) ||
                window_width <= 0.0 )
            {
                print_error( "Incorrect window width specification.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-perspective_distance" ) )
        {
            if( !get_real_argument( 0.0, &perspective_distance ) ||
                perspective_distance <= 0.0 )
            {
                print_error( "Incorrect perspective distance specification.\n");
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-scale" ) )
        {
            if( !get_real_argument( 0.0, &fit_factor ) || fit_factor <= 0.0 )
            {
                print_error( "Incorrect fit factor specification.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( argument, "-ambient" ) )
        {
            if( !get_real_argument( 0.0, &r ) ||
                !get_real_argument( 0.0, &g ) ||
                !get_real_argument( 0.0, &b ) )
            {
                print_error( "Missing ambient arguments.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
            lights.ambient_light = make_Colour_0_1( r, g, b );
        }
        else if( equal_strings( argument, "-directional" ) )
        {
            if( n_lights == 0 && lights.n_lights > 0 )
            {
                FREE( lights.lights );
                lights.n_lights = 0;
            }

            if( !get_real_argument( 0.0, &xv ) ||
                !get_real_argument( 0.0, &yv ) ||
                !get_real_argument( 0.0, &zv ) ||
                !get_real_argument( 0.0, &r ) ||
                !get_real_argument( 0.0, &g ) ||
                !get_real_argument( 0.0, &b ) )
            {
                print_error( "Missing directional arguments.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
            if( n_lights >= lights.n_lights )
            {
                SET_ARRAY_SIZE( lights.lights, lights.n_lights,
                                lights.n_lights+1, 1 );
                ++lights.n_lights;
            }
            lights.lights[n_lights].type = DIRECTIONAL_LIGHT;
            lights.lights[n_lights].colour = make_Colour_0_1( r, g, b );
            fill_Vector( lights.lights[n_lights].direction, xv, yv, zv );
            ++n_lights;
        }
        else if( equal_strings( argument, "-point" ) )
        {
            if( n_lights == 0 && lights.n_lights > 0 )
            {
                FREE( lights.lights );
                lights.n_lights = 0;
            }

            if( !get_real_argument( 0.0, &xv ) ||
                !get_real_argument( 0.0, &yv ) ||
                !get_real_argument( 0.0, &zv ) ||
                !get_real_argument( 0.0, &r ) ||
                !get_real_argument( 0.0, &g ) ||
                !get_real_argument( 0.0, &b ) ||
                !get_real_argument( 0.0, &a1 ) ||
                !get_real_argument( 0.0, &a2 ) ||
                !get_real_argument( 0.0, &a3 ) )
            {
                print_error( "Missing point light arguments.\n" );
                print_usage( argv[0] );
                return( 1 );
            }
            if( n_lights >= lights.n_lights )
            {
                SET_ARRAY_SIZE( lights.lights, lights.n_lights,
                                lights.n_lights+1, 1 );
                ++lights.n_lights;
            }
            lights.lights[n_lights].type = POINT_LIGHT;
            lights.lights[n_lights].colour = make_Colour_0_1( r, g, b );
            fill_Point( lights.lights[n_lights].position, xv, yv, zv );
            lights.lights[n_lights].attenuation[0] = a1;
            lights.lights[n_lights].attenuation[1] = a2;
            lights.lights[n_lights].attenuation[2] = a3;
            ++n_lights;
        }
        else if( equal_strings( argument, "-volume" ) )
        {
            if( render.n_volumes == 0 )
            {
                print_error( "Must specify a volume before -volume.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            if( !get_real_argument( 0.0, &threshold ) ||
                !get_string_argument( NULL, &colour_name ) )
            {
                print_error( "Incorrect volume specification.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            colour = convert_string_to_colour( colour_name );
            ray_object.render = render;
            if( render.n_volumes > 0 )
            {
                ALLOC( ray_object.render.volumes, render.n_volumes );
                for_less( i, 0, render.n_volumes )
                    ray_object.render.volumes[i] = render.volumes[i];
            }
            ray_object.regular_object_flag = FALSE;
            ray_object.threshold = threshold;
            ray_object.colour = colour;
            get_default_surfprop( &ray_object.spr );
            get_volume_sizes( render.volumes[render.n_volumes-1].volume, sizes);
            create_bitlist_3d( sizes[VIO_X], sizes[VIO_Y], sizes[VIO_Z],
                               &ray_object.done_bits );
            create_bitlist_3d( sizes[VIO_X], sizes[VIO_Y], sizes[VIO_Z],
                               &ray_object.surface_bits );
            ray_object.clip.n_clip_objects = clip.n_clip_objects;
            if( clip.n_clip_objects > 0 )
            {
                ALLOC( ray_object.clip.clip_objects,
                       clip.n_clip_objects );
                ALLOC( ray_object.clip.clip_inside_flags,
                       clip.n_clip_objects );
                for_less( i, 0, clip.n_clip_objects )
                {
                    ray_object.clip.clip_objects[i] =
                                             clip.clip_objects[i];
                    ray_object.clip.clip_inside_flags[i] =
                                             clip.clip_inside_flags[i];
                }
            }

            ADD_ELEMENT_TO_ARRAY( object_list, n_objects, ray_object,
                                  DEFAULT_CHUNK_SIZE );

            get_range_of_volume( render.volumes[render.n_volumes-1].volume,
                                 &min_pos, &max_pos );

            if( n_objects == 1 )
            {
                min_range = min_pos;
                max_range = max_pos;
            }
            else
            {
                apply_point_to_min_and_max( &min_pos,
                                            &min_range, &max_range );
                apply_point_to_min_and_max( &max_pos,
                                            &min_range, &max_range );
            }
        }
        else if( equal_strings( argument, "-clip" ) )
        {
            if( !get_string_argument( NULL, &argument ) )
            {
                print_error( "Error in -clip arguments.\n" );
                print_usage( argv[0] );
                return( 1 );
            }

            if( equal_strings( argument, "off" ) )
                clip.n_clip_objects = 0;
            else
            {
                if( equal_strings( argument, "inside" ) )
                    clip_inside_flag = TRUE;
                else if( equal_strings( argument, "outside" ) )
                    clip_inside_flag = FALSE;
                else
                {
                    print_error( "Error in -clip arguments.\n" );
                    print_usage( argv[0] );
                    return( 1 );
                }

                if( !get_string_argument( NULL, &argument ) )
                {
                    print_error( "Error in -clip arguments.\n" );
                    print_usage( argv[0] );
                    return( 1 );
                }

                get_current_transform( camera_space_flag, eye_specified, &eye,
                                       view_specified, &line_of_sight, &up_dir,
                                       &transform, &used_transform );
                if( load_file( argument,
                               &used_transform, bintree_factor,
                               line_width_override,
                               current_marker_colour, current_marker_size,
                               current_marker_type,
                               &n_objects_in_file, &file_objects ) != VIO_OK )
                    return( 1 );

                delete_general_transform( &used_transform );

                for_less( i, 0, n_objects_in_file )
                {
                    ADD_ELEMENT_TO_ARRAY( clip.clip_objects,
                                          clip.n_clip_objects,
                                          file_objects[i], 1 );
                    --clip.n_clip_objects;
                    ADD_ELEMENT_TO_ARRAY( clip.clip_inside_flags,
                                          clip.n_clip_objects,
                                          clip_inside_flag, 1 );
                }
            }
        }
        else
        {
            get_current_transform( camera_space_flag, eye_specified, &eye,
                                   view_specified, &line_of_sight, &up_dir,
                                   &transform, &used_transform );

            if( load_file( argument,
                           &used_transform, bintree_factor, line_width_override,
                           current_marker_colour, current_marker_size,
                           current_marker_type,
                           &n_objects_in_file, &file_objects ) != VIO_OK )
                return( 1 );

            delete_general_transform( &used_transform );

            initialize_object_traverse( &object_traverse, FALSE,
                                        n_objects_in_file, file_objects );

            while( get_next_object_traverse( &object_traverse, &object ) )
            {
                if( get_object_type( object ) == POLYGONS )
                {
                    if( compute_normals )
                    {
                        average_polygon_normals( get_polygons_ptr(object),
                                                 n_normal_iters, normal_ratio );
                    }
                }
                else if( get_object_type( object ) == QUADMESH )
                {
                    if( compute_normals )
                        compute_quadmesh_normals( get_quadmesh_ptr(object) );
                }
                 
                ray_object.regular_object_flag = TRUE;
                ray_object.object = object;
                ray_object.render = render;

                if( render.n_volumes > 0 )
                {
                    ALLOC( ray_object.render.volumes, render.n_volumes );
                    for_less( i, 0, render.n_volumes )
                        ray_object.render.volumes[i] = render.volumes[i];
                }

                ray_object.clip.n_clip_objects = clip.n_clip_objects;
                if( clip.n_clip_objects > 0 )
                {
                    ALLOC( ray_object.clip.clip_objects,
                           clip.n_clip_objects );
                    ALLOC( ray_object.clip.clip_inside_flags,
                           clip.n_clip_objects );
                    for_less( i, 0, clip.n_clip_objects )
                    {
                        ray_object.clip.clip_objects[i] =
                                                 clip.clip_objects[i];
                        ray_object.clip.clip_inside_flags[i] =
                                                 clip.clip_inside_flags[i];
                    }
                }

                ADD_ELEMENT_TO_ARRAY( object_list, n_objects, ray_object,
                                      DEFAULT_CHUNK_SIZE );

                (void) get_range_of_object( object,
                                            FALSE, &min_pos, &max_pos );

                if( get_object_type(object) == LINES )
                {
                    int   c;
                    VIO_Real  line_width;

                    line_width = (VIO_Real) get_lines_ptr(object)->line_thickness;
                    for_less( c, 0, VIO_N_DIMENSIONS )
                    {
                        Point_coord(min_pos,c) -= (VIO_Point_coord_type) line_width;
                        Point_coord(max_pos,c) += (VIO_Point_coord_type) line_width;
                    }
                }

                if( n_objects == 1 )
                {
                    min_range = min_pos;
                    max_range = max_pos;
                }
                else
                {
                    apply_point_to_min_and_max( &min_pos,
                                                &min_range, &max_range );
                    apply_point_to_min_and_max( &max_pos,
                                                &min_range, &max_range );
                }
            }
        }
    }

    define_view( &view, perspective_flag, eye_specified, &eye,
                 view_specified, &line_of_sight, &up_dir,
                 window_width, perspective_distance, which_eye,
                 stereo_eye_offset,
                 x_size, y_size, aspect, &min_range, &max_range, fit_factor );

    orient_lights( &lights, &view );

    initialize_pixels( &pixels, 0, 0,
                       view.x_viewport_size, view.y_viewport_size, 1.0, 1.0,
                       RGB_PIXEL );

    ray_trace_scene( &view, &lights,
                     n_objects, object_list,
                     MAX_DEPTH,
                     super_sampling,
                     super_sampling, &pixels );

    if( crop_flag )
    {
        pixels_struct  tmp;

        for_less( x, 0, pixels.x_size )
        for_less( y, 0, pixels.y_size )
        {
            col = PIXEL_RGB_COLOUR( pixels, x, y );
            ac = get_Colour_a( col );
            if( ac == 0 )
                PIXEL_RGB_COLOUR( pixels, x, y ) = make_rgba_Colour(0,0,0,0);
        }

        crop_pixels( &pixels, make_rgba_Colour(0,0,0,0), CROP_BORDER, &tmp );
        delete_pixels( &pixels );
        pixels = tmp;
    }

    if( behind_specified )
    {
        rc = get_Colour_r( behind_transparency_colour );
        gc = get_Colour_g( behind_transparency_colour );
        bc = get_Colour_b( behind_transparency_colour );

        behind_transparency_colour = make_Colour( rc, gc, bc );

        if( bc < 255 )
            ++bc;
        else
            --bc;
        alternate_behind = make_Colour( rc, gc, bc );
    }

    for_less( x, 0, pixels.x_size )
    for_less( y, 0, pixels.y_size )
    {
        col = PIXEL_RGB_COLOUR( pixels, x, y );
        a = get_Colour_a_0_1( col );

        if( behind_specified && a == 0.0 )
            col = behind_transparency_colour;
        else
        {
            COMPOSITE_COLOURS( col, col, background_colour );

            r = get_Colour_r_0_1( col );
            g = get_Colour_g_0_1( col );
            b = get_Colour_b_0_1( col );
            a = get_Colour_a_0_1( col );
            col = make_Colour_0_1( r*a, g*a, b*a );

            if( behind_specified && col == behind_transparency_colour )
                col = alternate_behind;
            else
                col = make_rgba_Colour_0_1( r, g, b, a );
        }

        PIXEL_RGB_COLOUR( pixels, x, y ) = col;
    }

    if( n_colour_comps == 0 )
    {
#define  N_LEVELS 6
        static    char  gray_scale[N_LEVELS+1] = " .:oOX";

        for( y = pixels.y_size-1;  y >= 0;  --y )
        {
            for_less( x, 0, pixels.x_size )
            {
                intensity = (VIO_Real) get_Colour_intensity(
                               PIXEL_RGB_COLOUR( pixels, x, y ) );
                print( "%c", gray_scale[ (int) (intensity * (VIO_Real) N_LEVELS)] );
            }
            print( "\n" );
        }
    }
    else
    {
        if( filename_extension_matches( output_file, "rgb" ) )
            n_colour_comps = 24;

        if( n_colour_comps == 1 )
        {
            pixels_tmp = pixels;
            convert_pixels24_to_gray_scale( &pixels_tmp, &pixels );
            delete_pixels( &pixels_tmp );
        }
        else if( n_colour_comps == 4 )
        {
            pixels_tmp = pixels;
            convert_pixels24_to_dithered( &pixels_tmp, &pixels, 16, colours16 );
            delete_pixels( &pixels_tmp );
        }
        else if( n_colour_comps == 8 )
        {
            pixels_tmp = pixels;
            convert_pixels24_to_dithered( &pixels_tmp, &pixels, 256,
                                          get_8bit_rgb_pixel_lookup() );
            delete_pixels( &pixels_tmp );
        }

        if( filename_extension_matches( output_file, "rgb" ) )
            (void) output_rgb_file( output_file, &pixels );
        else
        {
            if( open_file( output_file, WRITE_FILE, BINARY_FORMAT, &file )!= VIO_OK)
                return( 1 );

            if( io_pixels( file, WRITE_FILE, BINARY_FORMAT, &pixels ) != VIO_OK )
                return( 1 );

            (void) close_file( file );
        }
    }

    delete_pixels( &pixels );

    return( 0 );
}

#define  WINDOW_WIDTH_RATIO   0.3

#define  TOLERANCE 1.0e-4

static  void  define_view(
    view_struct   *view,
    VIO_BOOL       perspective_flag,
    VIO_BOOL       eye_specified,
    VIO_Point         *eye,
    VIO_BOOL       view_specified,
    VIO_Vector        *line_of_sight,
    VIO_Vector        *up_dir,
    VIO_Real          specified_window_width,
    VIO_Real          specified_perspective_distance,
    int           which_eye,
    VIO_Real          stereo_eye_offset,
    int           x_size,
    int           y_size,
    VIO_Real          pixel_aspect,
    VIO_Point         *min_range,
    VIO_Point         *max_range,
    VIO_Real          fit_factor )
{
    static  VIO_Vector  default_line_of_sight = { .77f, -.18f, -.6f };
    static  VIO_Vector  default_up_dir = { .55f, .6f, .55f };
    int             c;
    VIO_Vector          used_line_of_sight, used_up, offset;
    VIO_Real            radius, aspect, window_size;
    VIO_Point           centre;

    if( view_specified )
    {
        used_line_of_sight = *line_of_sight;
        used_up = *up_dir;
    }
    else
    {
        for_less( c, 0, VIO_N_DIMENSIONS )
        {
            if( (VIO_Real) Point_coord(*min_range,c) >=
                (VIO_Real) Point_coord(*max_range,c) - TOLERANCE )
                break;
        }

        if( c < VIO_N_DIMENSIONS )
        {
            fill_Vector( used_line_of_sight, 0.0f, 0.0f, 0.0f );
            fill_Vector( used_up, 0.0f, 0.0f, 0.0f );
            Vector_coord( used_line_of_sight, c ) = -1.0f;
            switch( c )
            {
            case VIO_X:
            case VIO_Y:
                Vector_coord( used_up, VIO_Z ) = 1.0f;
                break;
            case VIO_Z:
                Vector_coord( used_up, VIO_Y ) = 1.0f;
                break;
            }
        }
        else
        {
            used_line_of_sight = default_line_of_sight;
            used_up = default_up_dir;
        }
    }

    radius = distance_between_points( min_range, max_range ) / 2.0;
    INTERPOLATE_POINTS( centre, *min_range, *max_range, 0.5 );

    if( specified_perspective_distance > 0.0 &&
        specified_window_width <= 0.0 )
    {
        specified_window_width = WINDOW_WIDTH_RATIO * fit_factor *
                                 specified_perspective_distance;
    }
    else if( specified_perspective_distance <= 0.0 &&
             specified_window_width > 0.0 )
    {
        specified_perspective_distance = specified_window_width /
                                         WINDOW_WIDTH_RATIO / fit_factor;
    }
    else if( specified_perspective_distance <= 0.0 &&
             specified_window_width <= 0.0 )
    {
        specified_window_width = 2.0 * radius;
        specified_perspective_distance = specified_window_width /
                                         WINDOW_WIDTH_RATIO;
        specified_window_width *= fit_factor;
    }

    view->window_distance = specified_perspective_distance;
    window_size = specified_window_width;

    if( eye_specified )
    {
        view->eye = *eye;
        if( !view_specified )
        {
            SUB_POINTS( used_line_of_sight, centre, *eye );
        }
        NORMALIZE_VECTOR( used_line_of_sight, used_line_of_sight );
    }
    else
    {
        NORMALIZE_VECTOR( used_line_of_sight, used_line_of_sight );
        SCALE_VECTOR( offset, used_line_of_sight, -view->window_distance );
        ADD_POINT_VECTOR( view->eye, centre, offset );
    }

    view->line_of_sight = used_line_of_sight;
    SCALE_VECTOR( view->perpendicular, view->line_of_sight, -1.0 );
    CROSS_VECTORS( view->horizontal, used_up, view->perpendicular );
    CROSS_VECTORS( view->vertical, view->perpendicular, view->horizontal );
    NORMALIZE_VECTOR( view->horizontal, view->horizontal );
    NORMALIZE_VECTOR( view->vertical, view->vertical );

    view->x_viewport_size = x_size;
    view->y_viewport_size = y_size;

    aspect = pixel_aspect *
             (VIO_Real) view->y_viewport_size / (VIO_Real) view->x_viewport_size;

    view->window_width = window_size;
    view->window_height = window_size * aspect;

    view->perspective_flag = perspective_flag;

    if( perspective_flag && which_eye != 0 && stereo_eye_offset != 0.0 )
    {
        view->stereo_offset = (VIO_Real) which_eye / 2.0 * stereo_eye_offset /
                              100.0 * view->window_distance;
    }
    else
        view->stereo_offset = 0.0;
}

static  void  orient_lights(
    lights_struct   *lights,
    view_struct     *view )
{
    int      i;
    VIO_Vector   x_offset, y_offset, z_offset;

    for_less( i, 0, lights->n_lights )
    {
        if( lights->lights[i].type == DIRECTIONAL_LIGHT )
        {
            SCALE_VECTOR( x_offset, view->horizontal,
                          Vector_x(lights->lights[i].direction) );
            SCALE_VECTOR( y_offset, view->vertical,
                          Vector_y(lights->lights[i].direction) );
            SCALE_VECTOR( z_offset, view->perpendicular,
                          Vector_z(lights->lights[i].direction) );

            ADD_VECTORS( lights->lights[i].direction, x_offset, y_offset );
            ADD_VECTORS( lights->lights[i].direction,
                         lights->lights[i].direction, z_offset );
            NORMALIZE_VECTOR( lights->lights[i].direction,
                              lights->lights[i].direction );
        }
    }
}

static  void  define_default_lights(
    lights_struct   *lights )
{
    static  VIO_Vector  light_dir[] = { {  1.0f,  1.0f, -1.0f },
                                    { -1.0f,  1.0f, -1.0f },
                                    {  1.0f, -1.0f, -1.0f },
                                    { -1.0f, -1.0f, -1.0f } };
    int      i;

    lights->ambient_light = make_Colour_0_1( 0.3, 0.3, 0.3 );
    lights->n_lights = VIO_SIZEOF_STATIC_ARRAY( light_dir );
    ALLOC( lights->lights, lights->n_lights );
    for_less( i, 0, lights->n_lights )
    {
        lights->lights[i].colour = make_Colour_0_1( 0.6, 0.6, 0.6 );

        lights->lights[i].type = DIRECTIONAL_LIGHT;
        lights->lights[i].direction = light_dir[i];
    }
}

static  void  prepare_objects(
    int                 n_objects,
    object_struct       *objects[],
    VIO_General_transform   *transform,
    VIO_Real                bintree_factor,
    VIO_Real                line_width_override )
{
    int                     i, n_points;
    object_struct           *object;
    object_traverse_struct  object_traverse;
    VIO_Real                    xt, yt, zt;
    VIO_Point                   *points;

    initialize_object_traverse( &object_traverse, FALSE, n_objects, objects );

    while( get_next_object_traverse( &object_traverse, &object ) )
    {
        n_points = get_object_points( object, &points );
        for_less( i, 0, n_points )
        {
             general_transform_point( transform,
                                      (VIO_Real) Point_x(points[i]),
                                      (VIO_Real) Point_y(points[i]),
                                      (VIO_Real) Point_z(points[i]),
                                      &xt, &yt, &zt );

            fill_Point( points[i], xt, yt, zt );
        }

        if( get_object_type( object ) == LINES )
        {
            if( line_width_override > 0.0 )
                get_lines_ptr(object)->line_thickness =
                                           (float) line_width_override;

            if( bintree_factor > 0.0 )
                create_lines_bintree( get_lines_ptr(object),
                   (int) ((VIO_Real) get_lines_ptr(object)->n_points * 5.0 *
                          bintree_factor) );
        }
        else if( get_object_type( object ) == POLYGONS )
        {
            if( bintree_factor > 0.0 )
                create_polygons_bintree( get_polygons_ptr(object),
                   (int) ((VIO_Real) get_polygons_ptr(object)->n_items *
                          bintree_factor) );
        }
        else if( get_object_type( object ) == QUADMESH )
        {
            int  m, n;

            get_quadmesh_n_objects( get_quadmesh_ptr(object),
                                    &m, &n );

            if( bintree_factor > 0.0 )
                create_quadmesh_bintree( get_quadmesh_ptr(object),
                                         (int) ((VIO_Real) m * (VIO_Real) n *
                                                 bintree_factor ) );
        }
    }
}

static  VIO_BOOL  equal_transforms(
    VIO_General_transform  *t1,
    VIO_General_transform  *t2 )
{
    int        i, j;
    VIO_Transform  *l1, *l2;

    if( get_transform_type(t1) != LINEAR ||
        get_transform_type(t2) != LINEAR )
        return( FALSE );

    l1 = get_linear_transform_ptr( t1 );
    l2 = get_linear_transform_ptr( t2 );

    for_less( i, 0, 4 )
    for_less( j, 0, 4 )
        if( Transform_elem(*l1,i,j) != Transform_elem(*l2,i,j) )
            return( FALSE );

    return( TRUE );
}

static  VIO_Status  load_file(
    VIO_STR              filename,
    VIO_General_transform   *transform,
    VIO_Real                bintree_factor,
    VIO_Real                line_width_override,
    VIO_Colour              current_marker_colour,
    VIO_Real                current_marker_size,
    Marker_types        current_marker_type,
    int                 *n_objects,
    object_struct       ***objects )
{
    static  int            n_files_read = 0;
    static  file_lookup    *files_read = NULL;
    file_lookup            entry;
    int                    i;

    for_less( i, 0, n_files_read )
    {
        if( equal_strings( files_read[i].filename, filename ) &&
            equal_transforms( transform, &files_read[i].transform ) &&
            bintree_factor == files_read[i].bintree_factor &&
            line_width_override == files_read[i].line_width_override )
        {
            break;
        }
    }

    if( i >= n_files_read )
    {
        entry.filename = create_string( filename );
        copy_general_transform( transform, &entry.transform );
        entry.bintree_factor = bintree_factor;
        entry.line_width_override = line_width_override;

        if( input_objects_any_format( NULL, filename, current_marker_colour,
                                      current_marker_size, current_marker_type,
                                      &entry.n_objects, &entry.objects ) != VIO_OK )
            return( VIO_ERROR );

        ADD_ELEMENT_TO_ARRAY( files_read, n_files_read, entry, 1 );

        prepare_objects( entry.n_objects, entry.objects, transform,
                         bintree_factor, line_width_override );

    }

    *n_objects = files_read[i].n_objects;
    *objects = files_read[i].objects;

    return( VIO_OK );
}

static  struct
{
    VIO_STR    name;
    VIO_Vector    line_of_sight;
    VIO_Vector    up_dir;
}  preset_views[] = {
              { "top",    {  0.0f,  0.0f, -1.0f }, {  0.0f,  1.0f,  0.0f} },
              { "bottom", {  0.0f,  0.0f,  1.0f }, {  0.0f,  1.0f,  0.0f} },
              { "left",   {  1.0f,  0.0f,  0.0f }, {  0.0f,  0.0f,  1.0f} },
              { "right",  { -1.0f,  0.0f,  0.0f }, {  0.0f,  0.0f,  1.0f} },
              { "back",   {  0.0f,  1.0f,  0.0f }, {  0.0f,  0.0f,  1.0f} },
              { "front",  {  0.0f, -1.0f,  0.0f }, {  0.0f,  0.0f,  1.0f} }
                    };

static  int  convert_string_to_view_index(
    VIO_STR    string )
{
    int   i, found_index;

    found_index = -1;
    for_less( i, 0, VIO_SIZEOF_STATIC_ARRAY( preset_views ) )
    {
        if( equal_strings( &string[1], preset_views[i].name ) )
        {
            found_index = i;
            break;
        }
    }

    return( found_index );
}

static  VIO_BOOL   string_is_preset_view(
    VIO_STR    string )
{
    int   ind;

    ind = convert_string_to_view_index( string );

    return( ind >= 0 );
}

static  void   get_preset_view(
    VIO_STR   string,
    VIO_Vector   *line_of_sight,
    VIO_Vector   *up_dir )
{
    int   ind;

    ind = convert_string_to_view_index( string );

    *line_of_sight = preset_views[ind].line_of_sight;
    *up_dir = preset_views[ind].up_dir;
}

static  void  get_current_transform(
    VIO_BOOL            camera_space_flag,
    VIO_BOOL            eye_specified,
    VIO_Point              *eye,
    VIO_BOOL            view_specified,
    VIO_Vector             *line_of_sight,
    VIO_Vector             *up_dir,
    VIO_General_transform  *current_transform,
    VIO_General_transform  *used_transform )
{
    VIO_Vector               x_axis, y_axis, z_axis;
    VIO_Transform            to_bases, from_bases;
    VIO_General_transform    from_bases_gen, to_bases_gen, new_transform;

    if( !camera_space_flag )
    {
        copy_general_transform( current_transform, used_transform );
        return;
    }

    if( !view_specified || !eye_specified )
    {
        print_error( "Error: -camera must come after -view and -eye\n");
        copy_general_transform( current_transform, used_transform );
        return;
    }

    SCALE_VECTOR( z_axis, *line_of_sight, -1.0 );
    NORMALIZE_VECTOR( z_axis, z_axis );
    CROSS_VECTORS( x_axis, *up_dir, z_axis );
    NORMALIZE_VECTOR( x_axis, x_axis );
    CROSS_VECTORS( y_axis, z_axis, x_axis );
    NORMALIZE_VECTOR( y_axis, y_axis );

    make_change_to_bases_transform( eye, &x_axis, &y_axis, &z_axis, &to_bases );
    make_change_from_bases_transform( eye, &x_axis, &y_axis, &z_axis,
                                      &from_bases );

    create_linear_transform( &to_bases_gen, &to_bases );
    create_linear_transform( &from_bases_gen, &from_bases );

    concat_general_transforms( &from_bases_gen,
                               current_transform,
                               &new_transform );
    delete_general_transform( &from_bases_gen );

    concat_general_transforms( &new_transform,
                               &to_bases_gen,
                               used_transform );

    delete_general_transform( &new_transform );
    delete_general_transform( &to_bases_gen );
}
