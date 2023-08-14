 
#include  <volume_io.h>
#include  <GS_graphics.h>

static  void  make_matrix(
    Transform   *trans,
    Matrix      gl_trans )
{
    int      i, j;

    for_less( i, 0, 4 )
    {
        for_less( j, 0, 4 )
        {
            gl_trans[i][j] = (float) Transform_elem( *trans, j, i );
        }
    }
}

  void  GS_mult_transform(
    Transform  *transform )
{
    Matrix  gl_trans;

    make_matrix( transform, gl_trans );

    multmatrix( gl_trans );
}


  void  GS_load_transform(
    Transform  *transform )
{
    Matrix  gl_trans;

    make_matrix( transform, gl_trans );

    loadmatrix( gl_trans );
}

  void  GS_get_transform(
    Transform  *trans )
{
    int      i, j;
    Matrix   gl_trans;

    getmatrix( gl_trans );

    for_less( i, 0, 4 )
    {
        for_less( j, 0, 4 )
        {
            Transform_elem( *trans, j, i ) = (VIO_Real) gl_trans[i][j];
        }
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
    ortho( (Coord) x_min, (Coord) x_max, (Coord) y_min, (Coord) y_max,
           (Coord) z_min, (Coord) z_max );
}

  void  GS_ortho_2d(
    int   x_min,
    int   x_max,
    int   y_min,
    int   y_max )
{
    ortho( (Coord) x_min - (Coord) 0.5, (Coord) x_max + (Coord) 0.5,
           (Coord) y_min - (Coord) 0.5, (Coord) y_max + (Coord) 0.5,
           (Coord) -1.0,                (Coord) 1.0 );
}

  void  GS_frustum(
    VIO_Real   x_min,
    VIO_Real   x_max,
    VIO_Real   y_min,
    VIO_Real   y_max,
    VIO_Real   z_min,
    VIO_Real   z_max )
{
    window( (Coord) x_min, (Coord) x_max, (Coord) y_min, (Coord) y_max,
            (Coord) z_min, (Coord) z_max );
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
    pushmatrix();
}

  void  GS_pop_transform( void )
{
    popmatrix();
}

  void  GS_set_viewport(
    int            x_min,
    int            x_max,
    int            y_min,
    int            y_max )
{
    viewport( (Screencoord) x_min, (Screencoord) x_max,
              (Screencoord) y_min, (Screencoord) y_max );
}

  void  clear_overlay_planes( void )
{
    Matrix          save_projection;

    mmode( MPROJECTION );
    getmatrix( save_projection );
    mmode( MVIEWING );

    pushmatrix();
    pushviewport();

    fullscrn();
/*
    drawmode( NORMALDRAW );
    blendfunction( BF_ONE, BF_ZERO );
*/
    drawmode( OVERDRAW );
    color( 0 );
    clear();
    endfullscrn();

    popviewport();
    popmatrix();

    mmode( MPROJECTION );
    loadmatrix( save_projection );
    mmode( MVIEWING );
}

  void  GS_set_matrix_mode(
    Matrix_modes   mode )
{
    int  gl_mode;

    switch( mode )
    {
    case VIEWING_MATRIX:
        gl_mode = MVIEWING;
        break;
    case PROJECTION_MATRIX:
        gl_mode = MPROJECTION;
        break;
    }

    mmode( (short) gl_mode );
}
