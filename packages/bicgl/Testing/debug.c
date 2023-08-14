#include  <stdio.h>
#include  "../../GLUT/include/glut.h"

static int  test_glut_window(
    int   rgb_mode,
    int   double_buffer_flag,
    int   depth_buffer_flag )
{
    int                rgba, doub, depth;
    unsigned  int      mode;
    int                window_id;
    char               title[1000];

    mode = 0;

    if( rgb_mode )
        mode |= GLUT_RGB;
    else
        mode |= GLUT_INDEX;

    if( double_buffer_flag )
        mode |= GLUT_DOUBLE;
    else
        mode |= GLUT_SINGLE;

    if( depth_buffer_flag )
        mode |= GLUT_DEPTH;
    
    glutInitDisplayMode( mode );

    if( !glutGet( (GLenum) GLUT_DISPLAY_MODE_POSSIBLE ) )
    {
        (void) printf(
              "Test window: GLUT display mode not possible (%d,%d,%d)\n",
              rgb_mode, double_buffer_flag, depth_buffer_flag );
        return( -1 );
    }

    (void) sprintf( title, "Title: %d %d %d",
                    rgb_mode, double_buffer_flag, depth_buffer_flag );

    window_id = glutCreateWindow( title );

    if( window_id < 1 )
    {
        (void) printf( "Test window: Could not open GLUT window\n" );
        return( -1 );
    }

    rgba = glutGet((GLenum) GLUT_WINDOW_RGBA);
    doub = glutGet((GLenum) GLUT_WINDOW_DOUBLEBUFFER);
    depth = glutGet((GLenum) GLUT_WINDOW_DEPTH_SIZE);

    if( rgba != rgb_mode ||
        doub != double_buffer_flag ||
        (depth > 0) != depth_buffer_flag )
    {
        (void) printf(
            "Test window: Window was opened in incorrect mode ( %d %d %d ), got( %d %d %d )\n",
            rgb_mode, double_buffer_flag, depth_buffer_flag,
            rgba, doub, (depth > 0 ) );
    }

    return( window_id );
}

static  void  display_function( void )
{
}

static  void  test_glut_windows( void )
{
    int  i, window_ids[8];

    window_ids[0] = test_glut_window( 0, 0, 0 );
    window_ids[1] = test_glut_window( 0, 0, 1 );
    window_ids[2] = test_glut_window( 0, 1, 0 );
    window_ids[3] = test_glut_window( 0, 1, 1 );
    window_ids[4] = test_glut_window( 1, 0, 0 );
    window_ids[5] = test_glut_window( 1, 0, 1 );
    window_ids[6] = test_glut_window( 1, 1, 0 );
    window_ids[7] = test_glut_window( 1, 1, 1 );

    for( i = 0;  i < 8;  ++i )
    {
        glutSetWindow( window_ids[i] );
        glutDisplayFunc( display_function );
    }
}

int main(
    int    argc,
    char   *argv[] )
{
    glutInit( &argc, argv );

    (void) printf( "\n1:\n" );
    test_glut_windows();

    (void) printf( "\n2:\n" );
    test_glut_windows();

    /*--- glutMainLoop(); ---*/

    return( 0 );
}
