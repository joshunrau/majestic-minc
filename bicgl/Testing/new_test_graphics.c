
#include  <glut.h>
#include  <stdio.h>

#define  ENABLE_BUG

static  void  display_function( void )
{
    (void) printf( "In Display callback of window %d:   rgb mode: %d\n",
                    glutGetWindow(),
                    glutGet( (GLenum) GLUT_WINDOW_RGBA ) );
}

static  void  create_GLUT_window(
    char                   title[],
    int                    initial_x_pos,
    int                    initial_y_pos,
    int                    initial_x_size,
    int                    initial_y_size,
    int                    rgb_mode )
{
    int                rgba;
    unsigned  int      mode;
    int                window_id;

    mode = 0;

    if( rgb_mode )
        mode |= GLUT_RGB;
    else
        mode |= GLUT_INDEX;

    mode |= GLUT_DOUBLE;
    mode |= GLUT_DEPTH;
    
    glutInitWindowPosition( initial_x_pos, initial_y_pos );
    glutInitWindowSize( initial_x_size, initial_y_size );
    glutInitDisplayMode( mode );

#ifdef ENABLE_BUG
    if( !glutGet( (GLenum) GLUT_DISPLAY_MODE_POSSIBLE ) )
        (void) printf( "Display mode NOT possible\n" );
    if( !glutGet( (GLenum) GLUT_DISPLAY_MODE_POSSIBLE ) )
        (void) printf( "Display mode NOT possible\n" );
#endif

    window_id = glutCreateWindow( title );

    if( window_id < 1 )
    {
        (void) printf( "Could not open GLUT window for OpenGL\n" );
        return;
    }

    glutDisplayFunc( display_function );

    glutUseLayer( GLUT_NORMAL );

    rgba = glutGet((GLenum) GLUT_WINDOW_RGBA);

    (void) printf( "Requested rgb mode: %d   ", rgb_mode );
    (void) printf( "Acquired rgb mode: %d\n",
                   glutGet( (GLenum) GLUT_WINDOW_RGBA ) );
}

int main(
    int    argc,
    char   *argv[] )
{
    glutInit( &argc, argv );

    create_GLUT_window( "First Window",  100, 100, 300, 300, 1 );
    create_GLUT_window( "Second Window", 500, 100, 200, 200, 0 );

    glutMainLoop();

    return( 0 );
}
