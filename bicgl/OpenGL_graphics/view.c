#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H  

#include  <volume_io.h>
#include  <GS_graphics.h>

typedef  double  Matrix[16];

static  void  make_matrix(
    VIO_Transform  *trans,
    Matrix  gl_trans )
{
    int      i, j;

    for_less( i, 0, 4 )
    {
        for_less( j, 0, 4 )
            gl_trans[VIO_IJ(j,i,4)] = Transform_elem( *trans, i, j );
    }
}

  void  GS_load_transform(
    VIO_Transform  *transform )
{
    Matrix  gl_trans;

    make_matrix( transform, gl_trans );

    glLoadMatrixd( gl_trans );
}

  void  GS_mult_transform(
    VIO_Transform  *transform )
{
    Matrix  gl_trans;

    make_matrix( transform, gl_trans );

    glMultMatrixd( gl_trans );
}

  void  GS_get_transform(
    VIO_Transform  *trans )
{
    int         i, j, type;
    Matrix   gl_trans;

    glGetIntegerv( GL_MATRIX_MODE, &type );

    if( type == GL_MODELVIEW )
        glGetDoublev( GL_MODELVIEW_MATRIX, gl_trans );
    else if( type == GL_PROJECTION )
        glGetDoublev( GL_PROJECTION_MATRIX, gl_trans );
    else
        handle_internal_error( "GS_get_transform" );

    for_less( i, 0, 4 )
    {
        for_less( j, 0, 4 )
            Transform_elem( *trans, i, j ) = gl_trans[VIO_IJ(j,i,4)];
    }
}

  void  GS_ortho(
    VIO_Real   x_min,
    VIO_Real   x_max,
    VIO_Real   y_min,
    VIO_Real   y_max,
    VIO_Real   z_min,
    VIO_Real   z_max )
{
    glLoadIdentity();
    glOrtho( x_min, x_max, y_min, y_max, z_min, z_max );
}

  void  GS_ortho_2d(
    int   x_min,
    int   x_max,
    int   y_min,
    int   y_max )
{

    glLoadIdentity();

#ifdef MESA_OPENGL

#ifdef NOT_NEEDED
{
    static  VIO_BOOL  first = TRUE;
    static  VIO_Real     start, end;

    if( first )
    {
        first = FALSE;
        if( getenv( "START" ) == NULL ||
            sscanf( getenv("START"), "%lf", &start ) != 1 )
            start = 0.0;
        if( getenv( "END" ) == NULL ||
            sscanf( getenv("END"), "%lf", &end ) != 1 )
            end = 1.0;
    }

    glOrtho( (VIO_Real) x_min + start, (VIO_Real) x_max + end,
             (VIO_Real) y_min + start, (VIO_Real) y_max + end, -1.0, 1.0 );
}
#endif

    glOrtho( (VIO_Real) x_min, (VIO_Real) x_max + 1.0,
             (VIO_Real) y_min, (VIO_Real) y_max + 1.0, -1.0, 1.0 );
#else

#ifdef __alpha
    glOrtho( (VIO_Real) x_min - 0.5, (VIO_Real) x_max + 1.0,
             (VIO_Real) y_min - 0.5, (VIO_Real) y_max + 1.0, -1.0, 1.0 );
#else
    glOrtho( (VIO_Real) x_min, (VIO_Real) x_max + 1.0,
             (VIO_Real) y_min, (VIO_Real) y_max + 1.0, -1.0, 1.0 );
#endif

#endif
}

  void  GS_frustum(
    VIO_Real   x_min,
    VIO_Real   x_max,
    VIO_Real   y_min,
    VIO_Real   y_max,
    VIO_Real   z_min,
    VIO_Real   z_max )
{
    glLoadIdentity();
    glFrustum( x_min, x_max, y_min, y_max, z_min, z_max );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_window_view
@INPUT      : window
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Initializes the window view, called on creation of the window.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ARGSUSED */

  void  GS_initialize_window_view(
    GSwindow   window )
{
}

  void  GS_push_transform( void )
{
    glPushMatrix();
}

  void  GS_pop_transform( void )
{
    glPopMatrix();
}

  void  GS_set_viewport(
    int            x_min,
    int            x_max,
    int            y_min,
    int            y_max )
{
    glViewport( x_min, y_min, x_max - x_min + 1, y_max - y_min + 1 );
}

  void  clear_overlay_planes( void )
{
#ifdef TO_DO
#endif
}

  void  GS_set_matrix_mode(
    Matrix_modes   mode )
{
    int  gl_mode;

    switch( mode )
    {
    case VIEWING_MATRIX:
        gl_mode = GL_MODELVIEW;
        break;
    case PROJECTION_MATRIX:
        gl_mode = GL_PROJECTION;
        break;
    }

    glMatrixMode( (GLenum) gl_mode );
}
