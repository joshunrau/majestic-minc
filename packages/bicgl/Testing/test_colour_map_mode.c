
#include  <internal_volume_io.h>
#include  <graphics.h>

#define  COLOUR_INDEX1       200
#define  COLOUR_INDEX2       201
#define  COLOUR1             CYAN
#define  COLOUR2             MAGENTA

#define  LIGHT_INDEX         0

#define  PIXELS_X_ZOOM       1.0
#define  PIXELS_Y_ZOOM       1.0

#define  N_FONTS             10
#define  N_FONTS_TO_DRAW      N_FONTS

#define  TIMER_INCREMENT  0.1

typedef struct
{
    Gwindow          window;
    int              prev_rotation_mouse_x;
    VIO_BOOL          in_rotation_mode;
    pixels_struct    pixels;
    pixels_struct    ball;
    VIO_Real             ball_x;
    VIO_Real             ball_y;
    VIO_Real             ball_x_dir;
    VIO_Real             ball_y_dir;
    lines_struct     lines;
    lines_struct     lines_2d;
    lines_struct     single_point;
    text_struct      text;
    text_struct      font_examples[N_FONTS_TO_DRAW];
    polygons_struct  polygons;
    Transform        modeling_transform;
    VIO_Real             last_message;
    VIO_BOOL          double_buffer_flag;

}main_struct;

static  void  update(
    Gwindow    window,
    void       *update_data )
{
    int           i;
    main_struct   *info;

    info = (main_struct *) update_data;

    G_clear_window( window );

    G_set_zbuffer_state( window, FALSE );
    G_set_lighting_state( window, FALSE );
    G_set_view_type( window, PIXEL_VIEW );

    G_set_colour_map_entry( window, COLOUR_INDEX1, COLOUR1 );
    G_set_colour_map_entry( window, COLOUR_INDEX2, COLOUR2 );

    G_draw_pixels( window, &info->pixels ); 
    G_draw_text( window, &info->text );

    G_update_window( window );
}

static  void  resize_window(
    Gwindow    window,
    int        x,
    int        y,
    int        x_size,
    int        y_size,
    void       *update_data )
{
    main_struct   *info;

    info = (main_struct *) update_data;

    print( "Window resized/moved: %d %d    %d by %d\n", x, y, x_size, y_size );

    fill_Point( info->lines_2d.points[0], 5.0, 5.0, 0.0 );
    fill_Point( info->lines_2d.points[1], (VIO_Real) x_size - 5.0, 5.0, 0.0 );
    fill_Point( info->lines_2d.points[2], (VIO_Real) x_size - 5.0,
                                          (VIO_Real) y_size - 5.0, 0.0 );
    fill_Point( info->lines_2d.points[3], 5.0, (VIO_Real) y_size - 5.0, 0.0 );
    G_set_update_flag( window );
}

static  void  left_mouse_down(
    Gwindow    window,
    int        x,
    int        y,
    void       *update_data )
{
    main_struct   *info;

    info = (main_struct *) update_data;

    info->prev_rotation_mouse_x = x;
    info->in_rotation_mode = TRUE;
}

static  void  update_rotation(
    Gwindow       window,
    main_struct   *info,
    int           x )
{
    VIO_Real        angle_in_degrees;
    Transform   rotation_transform;
    VIO_Point       centre_of_rotation;

    angle_in_degrees = (VIO_Real) (info->prev_rotation_mouse_x - x);

    make_rotation_transform( angle_in_degrees * DEG_TO_RAD, Y,
                             &rotation_transform );

    fill_Point( centre_of_rotation, 0.3, 0.0, 0.0 );
    make_transform_relative_to_point( &centre_of_rotation,
                                      &rotation_transform,
                                      &rotation_transform );
    concat_transforms( &info->modeling_transform, &info->modeling_transform,
                       &rotation_transform );
    G_set_modeling_transform( window, &info->modeling_transform );

    info->prev_rotation_mouse_x = x;

    G_set_update_flag( window );
}

static  void  left_mouse_up(
    Gwindow    window,
    int        x,
    int        y,
    void       *update_data )
{
    main_struct   *info;

    info = (main_struct *) update_data;

    update_rotation( window, info, x );
    info->in_rotation_mode = FALSE;
}

static  void  middle_mouse_down(
    Gwindow    window,
    int        x,
    int        y,
    void       *update_data )
{
    main_struct   *info;

    info = (main_struct *) update_data;

    print( "MIDDLE mouse down\n" );
}

static  void  middle_mouse_up(
    Gwindow    window,
    int        x,
    int        y,
    void       *update_data )
{
    main_struct   *info;

    info = (main_struct *) update_data;

    print( "MIDDLE mouse up\n" );
}

static  void  right_mouse_down(
    Gwindow    window,
    int        x,
    int        y,
    void       *update_data )
{
    main_struct   *info;

    info = (main_struct *) update_data;

    info->double_buffer_flag = !info->double_buffer_flag;

    if( info->double_buffer_flag )
        print( "Double buffering: On\n" );
    else
        print( "Double buffering: Off\n" );

    G_set_double_buffer_state( window, info->double_buffer_flag );
    G_set_update_flag( window );
}

static  void  right_mouse_up(
    Gwindow    window,
    int        x,
    int        y,
    void       *update_data )
{
    main_struct   *info;

    info = (main_struct *) update_data;

    print( "RIGHT mouse up\n" );
}

static  void  mouse_movement(
    Gwindow    window,
    int        x,
    int        y,
    void       *update_data )
{
    main_struct   *info;
    int           x_size, y_size;

    info = (main_struct *) update_data;

    if( info->in_rotation_mode )
        update_rotation( window, info, x );

    G_get_window_size( window, &x_size, &y_size );

    (void) sprintf( info->text.string,
                    "Mouse: %4d,%4d pixels   %4.2f,%4.2f window",
                    x, y, ((VIO_Real) x + 0.5) / (VIO_Real) x_size,
                          ((VIO_Real) y + 0.5) / (VIO_Real) y_size );

    G_set_update_flag( window );
}

static  void  enter_window(
    Gwindow    window,
    void       *update_data )
{
}

static  void  leave_window(
    Gwindow    window,
    void       *update_data )
{
    main_struct   *info;

    info = (main_struct *) update_data;

    (void) sprintf( info->text.string,
                    "Mouse is out of the window" );

    G_set_update_flag( window );
}

static  void  key_down(
    Gwindow    window,
    int        key,
    void       *update_data )
{
    int           i;
    main_struct   *info;

    info = (main_struct *) update_data;

    print( "Key pressed down: \"%c\" (%d)\n", key, key );

    if( key == '' )
    {
        print( "Quitting\n" );

        delete_lines( &info->lines );
        delete_lines( &info->lines_2d );
        delete_lines( &info->single_point );

        delete_polygons( &info->polygons );

        delete_pixels( &info->pixels );

        delete_text( &info->text );

        for_less( i, 0, N_FONTS )
            delete_text( &info->font_examples[i] );

        (void) G_delete_window( info->window );

        G_terminate();

        output_alloc_to_file( NULL );

        exit( 0 );
    }
}

/* ARGSUSED */

static  void  key_up(
    Gwindow    window,
    int        key,
    void       *update_data )
{
    main_struct   *info;

    info = (main_struct *) update_data;

    print( "Key released: \"%c\" (%d)\n", key, key );
}

static  void  timer_function(
    void       *update_data )
{
    int            window_x_size, window_y_size;
    VIO_Real           ball_size;
    main_struct    *info;

    info = (main_struct *) update_data;

    G_get_window_size( info->window, &window_x_size, &window_y_size );

    ball_size = (VIO_Real) info->ball.x_size;

    if( info->ball_x - ball_size < 0.0 ||
        info->ball_x + ball_size >= (VIO_Real) window_x_size ||
        info->ball_y - ball_size < 0.0 ||
        info->ball_y + ball_size >= (VIO_Real) window_y_size )
    {
        info->ball_x = (VIO_Real) window_x_size / 2.0;
        info->ball_y = (VIO_Real) window_y_size / 2.0;
    }

    info->ball_x += info->ball_x_dir;
    if( info->ball_x - ball_size < 0.0 )
    {
        info->ball_x_dir = FABS( info->ball_x_dir );
        info->ball_x += info->ball_x_dir;
    }
    else if( info->ball_x + ball_size >= (VIO_Real) window_x_size )
    {
        info->ball_x_dir = - FABS( info->ball_x_dir );
        info->ball_x += info->ball_x_dir;
    }

    info->ball_y += info->ball_y_dir;
    if( info->ball_y - ball_size < 0.0 )
    {
        info->ball_y_dir = FABS( info->ball_y_dir );
        info->ball_y += info->ball_y_dir;
    }
    else if( info->ball_y + ball_size >= (VIO_Real) window_y_size )
    {
        info->ball_y_dir = - FABS( info->ball_y_dir );
        info->ball_y += info->ball_y_dir;
    }

    info->ball.x_position = ROUND( info->ball_x );
    info->ball.y_position = ROUND( info->ball_y );

    G_set_update_flag( info->window );
    G_add_timer_function( TIMER_INCREMENT, timer_function, (void *) info );
}

static  void  idle_function(
    void       *update_data )
{
    VIO_Real           current_time;
    main_struct    *info;

    info = (main_struct *) update_data;

    current_time = current_realtime_seconds();

    if( current_time >= info->last_message + 1.0 )
    {
        print( "Idling.\n" );
        info->last_message = current_time + 1.0;
    }

    sleep_program( 0.001 );
}

int main(
    int    argc,
    char   *argv[] )
{
    main_struct       info;
    Status            status;
    VIO_BOOL           stereo_flag;
    static Font_types font_types[N_FONTS] = {
                          FIXED_FONT,
                          SIZED_FONT,
                          SIZED_FONT,
                          SIZED_FONT,
                          SIZED_FONT,
                          SIZED_FONT,
                          SIZED_FONT,
                          SIZED_FONT,
                          SIZED_FONT,
                          SIZED_FONT };
    static VIO_Real       font_sizes[N_FONTS] = { 10.0,
                           6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0 };
    static Surfprop   spr = { 0.2f, 0.5f, 0.5f, 20.0f, 1.0f };
    VIO_Point             point;
    VIO_Vector            normal, light_direction;
    VIO_Real              dx, dy;
    int               x_position, y_position, ball_size;
    int               x_pixel, y_pixel;
    int               i, j, pixels_x_size, pixels_y_size, x_size, y_size;
    VIO_Real              eye_separation;
    static VIO_Point      origin = { 0.0f, 0.0f, 2.0f };
    static VIO_Vector     up_direction = { 0.0f, 1.0f, 0.0f };
    static VIO_Vector     line_of_sight = { 0.0f, 0.0f, -1.0f };

    info.in_rotation_mode = FALSE;
    info.double_buffer_flag = TRUE;

    stereo_flag = (argc > 1);

    if( !stereo_flag || sscanf( argv[1], "%lf", &eye_separation ) != 1 )
        eye_separation = 0.0;

    status = G_create_window( "Test Window",
                              100, 600, 300, 300,
                              TRUE, info.double_buffer_flag, TRUE, 0,
                              &info.window );

    G_set_transparency_state( info.window, TRUE );

    if( status != OK )
        return( 1 );

    if( stereo_flag )
        G_set_background_colour( info.window, WHITE );

    G_set_3D_view( info.window, &origin, &line_of_sight, &up_direction,
                   0.01, 4.0, TRUE, 2.0, stereo_flag, eye_separation, 2.0, 2.0 );

    fill_Point( point, -0.3, 0.3, 0.0 );
    G_transform_point( info.window, &point, MODEL_VIEW, &x_pixel, &y_pixel );

    print( "(%g,%g,%g) maps to %d %d in pixels\n",
            Point_x(point), Point_y(point), Point_z(point), x_pixel, y_pixel );

    /* ------- define text to be drawn (text.string filled in later ----- */

    fill_Point( point, 10.0, 10.0, 0.0 );
    initialize_text( &info.text, &point, 1, SIZED_FONT, 14.0 );
    info.text.string = alloc_string( 2000 );
    (void) sprintf( info.text.string,
                                "Mouse: %4d,%4d pixels   %4.2f,%4.2f window",
                                0, 0, 0.0, 0.0 );

    for_less( i, 0, N_FONTS )
    {
        fill_Point( point, 10.0, 30.0 + (VIO_Real) i * 20.0, 0.0 );
        initialize_text( &info.font_examples[i], &point,
                         make_Colour(255,0,255),
                         font_types[i], font_sizes[i] );

        if( font_types[i] == FIXED_FONT )
            info.font_examples[i].string = create_string( "Fixed font: " );
        else
        {
            info.font_examples[i].string = alloc_string( 1000 );
            (void) sprintf( info.font_examples[i].string, "Sized font %g points: ",
                            font_sizes[i] );
        }

        concat_to_string( &info.font_examples[i].string,
                          "abcdefghijklmnopqrstuvwxyz" );
    }

    /* ------------ define line to be drawn  ------------- */

    initialize_lines( &info.lines, make_Colour(0,255,0) );

    info.lines.n_points = 4;
    ALLOC( info.lines.points, 4 );
    fill_Point( info.lines.points[0], 0.0, 0.0, 0.0 );
    fill_Point( info.lines.points[1], 1.0, 1.0, 0.0 );
    fill_Point( info.lines.points[2], -0.3, 1.0, -1.0 );
    fill_Point( info.lines.points[3], 0.5, 0.4, 0.0 );

    info.lines.n_items = 3;
    ALLOC( info.lines.end_indices, info.lines.n_items );
    info.lines.end_indices[0] = 2;
    info.lines.end_indices[1] = 4;
    info.lines.end_indices[2] = 5;

    ALLOC( info.lines.indices, info.lines.end_indices[info.lines.n_items-1] );
    info.lines.indices[0] = 0;
    info.lines.indices[1] = 1;

    info.lines.indices[2] = 0;
    info.lines.indices[3] = 2;

    info.lines.indices[4] = 3;

    /* --- initialize single point */

    initialize_lines( &info.single_point, make_Colour(255,255,255) );

    info.single_point.n_points = 1;
    ALLOC( info.single_point.points, 4 );
    fill_Point( info.single_point.points[0], 0.0, 0.0, 0.0 );
    info.single_point.n_items = 1;
    ALLOC( info.single_point.end_indices, info.single_point.n_items );
    info.single_point.end_indices[0] = 1;
    ALLOC( info.single_point.indices,
           info.single_point.end_indices[info.single_point.n_items-1] );
    info.single_point.indices[0] = 0;

    /* ------------ define 2d line to be drawn  ------------- */

    initialize_lines( &info.lines_2d, make_Colour(0,255,150) );

    G_get_window_size( info.window, &x_size, &y_size );

    info.lines_2d.n_points = 4;
    ALLOC( info.lines_2d.points, 4 );
    fill_Point( info.lines_2d.points[0], 5.0, 5.0, 0.0 );
    fill_Point( info.lines_2d.points[1], (VIO_Real) x_size - 5.0, 5.0, 0.0 );
    fill_Point( info.lines_2d.points[2], (VIO_Real) x_size - 5.0, (VIO_Real) y_size - 5.0,
                                    0.0 );
    fill_Point( info.lines_2d.points[3], 5.0, (VIO_Real) y_size - 5.0, 0.0 );

    info.lines_2d.n_items = 1;
    ALLOC( info.lines_2d.end_indices, info.lines_2d.n_items );
    info.lines_2d.end_indices[0] = 5;

    ALLOC( info.lines_2d.indices, info.lines_2d.end_indices[info.lines_2d.n_items-1] );
    info.lines_2d.indices[0] = 0;
    info.lines_2d.indices[1] = 1;
    info.lines_2d.indices[2] = 2;
    info.lines_2d.indices[3] = 3;
    info.lines_2d.indices[4] = 0;

    /* ------------ define pixels to be drawn  ------------- */

    pixels_x_size = 250;
    pixels_y_size = 250;

    x_position = 10;
    y_position = 10;
    initialize_pixels( &info.pixels, x_position, y_position,
                       pixels_x_size, pixels_y_size,
                       PIXELS_X_ZOOM, PIXELS_Y_ZOOM, COLOUR_INDEX_8BIT_PIXEL );

    for_less( i, 0, pixels_x_size )
    {
        for_less( j, 0, pixels_y_size )
        {
            PIXEL_COLOUR_INDEX_8(info.pixels,i,j) = COLOUR_INDEX1 + (i/10 + j/10) % 2;
        }
    }

    /* ------------ define bouncing ball ------------- */

    ball_size = 50;

    info.ball_x = 10.0;
    info.ball_y = 10.0;
    info.ball_x_dir = 2.324423;
    info.ball_y_dir = 1.7123123;

    initialize_pixels( &info.ball, ROUND(info.ball_x), ROUND(info.ball_y),
                       ball_size, ball_size,
                       1.0, 1.0, RGB_PIXEL );

    for_less( i, 0, ball_size )
    {
        dx = (VIO_Real) i - (VIO_Real) (ball_size-1) / 2.0;
        for_less( j, 0, ball_size )
        {
            dy = (VIO_Real) j - (VIO_Real) (ball_size-1) / 2.0;
            if( dx * dx + dy * dy < (VIO_Real) ((ball_size-1) * (ball_size-1))/4.0)
                PIXEL_RGB_COLOUR(info.ball,i,j) = make_rgba_Colour(0,0,255,255);
            else
                PIXEL_RGB_COLOUR(info.ball,i,j) = make_rgba_Colour(0,0,0,0);
        }
    }

    /* ------------ define polygons to be drawn  ------------- */

    initialize_polygons( &info.polygons, make_Colour(0,255,255), &spr );

    start_new_polygon( &info.polygons );

    fill_Point( point, -0.3, -0.3, 0.0 );
    fill_Vector( normal, 0.0, 0.0, 1.0 );
    add_point_to_polygon( &info.polygons, &point, &normal );
    fill_Point( point, 0.3, -0.3, 0.0 );
    add_point_to_polygon( &info.polygons, &point, &normal );
    fill_Point( point, 0.3, 0.3, 0.0 );
    add_point_to_polygon( &info.polygons, &point, &normal );
    fill_Point( point, -0.3, 0.3, 0.0 );
    add_point_to_polygon( &info.polygons, &point, &normal );


    /* ------------ define lights ----------------- */

    fill_Vector( light_direction, 1.0, 1.0, -1.0 );/* from over left shoulder */

    G_define_light( info.window, LIGHT_INDEX, DIRECTIONAL_LIGHT,
                    make_Colour(255,255,255),
                    &light_direction, (VIO_Point *) 0, 0.0, 0.0 );
    G_set_light_state( info.window, LIGHT_INDEX, TRUE );

    /* --------------------------------------- */
    /* ------------ do main loop ------------- */
    /* --------------------------------------- */

    make_identity_transform( &info.modeling_transform );

    print( "Hold down left button and move mouse to rotate\n" );
    print( "Hit ESC key to exit\n" );

    G_set_update_function( info.window, update, (void *) &info );
    G_set_resize_function( info.window, resize_window, (void *) &info );
    G_set_left_mouse_down_function( info.window, left_mouse_down, (void *) &info );
    G_set_left_mouse_up_function( info.window, left_mouse_up, (void *) &info );
    G_set_middle_mouse_down_function( info.window, middle_mouse_down, (void *)&info);
    G_set_middle_mouse_up_function( info.window, middle_mouse_up, (void *) &info );
    G_set_right_mouse_down_function( info.window, right_mouse_down, (void *) &info );
    G_set_right_mouse_up_function( info.window, right_mouse_up, (void *) &info );
    G_set_mouse_movement_function( info.window, mouse_movement, (void *) &info );
    G_set_window_enter_function( info.window, enter_window, (void *) &info );
    G_set_window_leave_function( info.window, leave_window, (void *) &info );
    G_set_key_down_function( info.window, key_down, (void *) &info );
    G_set_key_up_function( info.window, key_up, (void *) &info );

    G_add_timer_function( TIMER_INCREMENT, timer_function, (void *) &info );

    info.last_message = current_realtime_seconds();
    G_add_idle_function( idle_function, (void *) &info );

    G_main_loop();

    return( 0 );
}
