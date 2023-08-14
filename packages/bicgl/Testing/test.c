
#include  <stdio.h>
#include  <glut.h>

int main(
    int    argc,
    char   *argv[] )
{
    int      mode;
    int      first, second;

    mode = GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH;

    glutInitWindowPosition( 100, 100 );
    glutInitWindowSize( 300, 300 );

    glutInitDisplayMode( mode );

    first = glutCreateWindow( "first" );

    (void) printf( "%d %d %d\n",
                   glutGet((GLenum) GLUT_WINDOW_RGBA),
                   glutGet((GLenum) GLUT_WINDOW_DOUBLEBUFFER),
                   glutGet((GLenum) GLUT_WINDOW_DEPTH_SIZE) );

    mode = GLUT_INDEX | GLUT_DOUBLE | GLUT_DEPTH;

    glutInitWindowPosition( 600, 100 );
    glutInitWindowSize( 300, 300 );

    glutInitDisplayMode( mode );

    second = glutCreateWindow( "second" );

    (void) printf( "%d %d %d\n",
                   glutGet((GLenum) GLUT_WINDOW_RGBA),
                   glutGet((GLenum) GLUT_WINDOW_DOUBLEBUFFER),
                   glutGet((GLenum) GLUT_WINDOW_DEPTH_SIZE) );

    glutMainLoop();

    return( 0 );
}
