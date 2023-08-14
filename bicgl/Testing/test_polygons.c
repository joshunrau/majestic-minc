#include  <internal_volume_io.h>
#include  <graphics.h>

static  void  process_polygons( polygons_struct * );
static  void  display_polygons( polygons_struct * );

int  main(
    int  argc,
    char *argv[] )
{
    Status          status;
    char            *input_filename;
    File_formats    format;
    int             n_objects;
    object_struct   **object_list;
    polygons_struct *polygons;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( "", &input_filename ) )
    {
        print( "Need a filename.\n" );
        return( 1 );
    }

    status = input_graphics_file( input_filename, &format, &n_objects,
                                  &object_list );

    if( status != OK )
        return( 1 );

    if( n_objects >= 1 && object_list[0]->object_type == POLYGONS )
    {
        polygons = get_polygons_ptr( object_list[0] );

        process_polygons( polygons );
    }

    display_polygons( polygons );

    delete_object_list( n_objects, object_list );

    return( status );
}

static  void  process_polygons( polygons_struct *polygons )
{
    int    poly_index, polygon_size, vertex, point_index;
    VIO_Point  point;

    print( "Number of points = %d\n", polygons->n_points );

    print( "Number of polygons = %d\n", polygons->n_items );

    for_less( poly_index, 0, polygons->n_items )
    {
        polygon_size = GET_OBJECT_SIZE( *polygons, poly_index );

        print( "Polygon %d:  size %d\n", poly_index, polygon_size );

        for_less( vertex, 0, polygon_size )
        {
            point_index = polygons->indices[
                        POINT_INDEX(polygons->end_indices,poly_index,vertex)];
            point = polygons->points[point_index];
            print( "   points[%d] = %g %g %g\n", point_index,
                   Point_x(point), Point_y(point), Point_z(point) );
        }
        print( "\n" );
    }
}

static  void  display_polygons( polygons_struct *polygons )
{
    window_struct     *window;
    static VIO_Point      origin = { 0.0, 0.0, 3.0 };
    static VIO_Vector     up_direction = { 0.0, 1.0, 0.0 };
    static VIO_Vector     line_of_sight = { 0.0, 0.0, -1.0 };
    VIO_Vector            light_direction;

    (void) G_create_window( "Polygons", -1, -1, -1, -1, &window );

    G_set_3D_view( window, &origin, &line_of_sight, &up_direction,
                   0.01, 6.0, TRUE, 2.0, FALSE, 0.0, 2.0, 2.0 );

    G_set_background_colour( window, PINK );

    fill_Vector( light_direction, 1.0, -1.0, -1.0 );/* from over left shoulder */

    G_define_light( window, 1, DIRECTIONAL_LIGHT, WHITE,
                    &light_direction, (VIO_Point *) 0, 0.0, 0.0 );
    G_set_light_state( window, 1, TRUE );

    G_set_lighting_state( window, TRUE );

    G_set_view_type( window, MODEL_VIEW );
    G_draw_polygons( window, polygons );

    G_update_window( window );

    print( "Hit return to continue: " );

    (void) getchar();
}
