
/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_tps.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: library of routines for warping/mapping transformations.
@METHOD     : The original source code for these routines are taken from
              program VOI, written by Weiqian Dai.
@GLOBALS    : 
@CALLS      : 
@CREATED    : Dec 2, 1991 LC
@MODIFIED   : Mon Apr  5 09:00:54 EST 1993 louis 
                - building new library routines, with prototypes
@MODIFIED   : Jul 6 1995,   David MacDonald removed recipes type code, rewrote
---------------------------------------------------------------------------- */

#include "bicpl_internal.h"

/* prototype definitions: */

static  void  calculate_weights(
    VIO_Real    **values,
    VIO_Real    **INVML,
    VIO_Real    **INVMLY,
    int     n_points,
    int     n_dims,
    int     n_values );

static  void  makeL(
    VIO_Real   **positions,
    VIO_Real   **ML,
    int    n_points,
    int    n_dims );

/* This function will get coefficients of the warping function. */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_nonlinear_warp
@INPUT      : positions
              values
              n_points
              n_dims
              n_values
@OUTPUT     : INVMLY
@RETURNS    : 
@DESCRIPTION: Computes the thin plate splines weights required to interpolate
              the given values at the positions.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_nonlinear_warp(
   VIO_Real     **positions,   /* n_points x n_dims */
   VIO_Real     **values,   /* n_points x n_values */
   VIO_Real     **INVMLY,   /* n_points+1+n_dims x n_values */
   int      n_points,
   int      n_dims,
   int      n_values )
{
   VIO_Real    **ML,**INVML;

   VIO_ALLOC2D( ML, n_points+n_dims+1, n_points+n_dims+1 );
   VIO_ALLOC2D( INVML, n_points+n_dims+1, n_points+n_dims+1 );

   /* This function will build the L matrix */

   makeL( positions, ML, n_points, n_dims );

   (void) invert_square_matrix( n_points+n_dims+1, ML, INVML );

   /*  build the array of deformation vectors   */

   calculate_weights( values, INVML, INVMLY, n_points, n_dims, n_values );

   VIO_FREE2D( ML );
   VIO_FREE2D( INVML );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : makeL
@INPUT      : positions
              n_points
              n_dims
@OUTPUT     : 
              ML
@RETURNS    : 
@DESCRIPTION: Creates the L matrix full of distance weights for the positions,
              used for computing the weights.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  makeL(
    VIO_Real   **positions,
    VIO_Real   **ML,
    int    n_points,
    int    n_dims )
{
    int    i,j;
    VIO_Real   fu;
 
    /* initialize matrix to zero */

    for_less( i, 0, n_points+n_dims+1 )
    {
        for_less( j, 0, n_points+n_dims+1 )
            ML[i][j] = 0.0;
    }
    
    /* set rest of the K matrix as follows */

    for_less( i, 0, n_points )
    {
        for_less( j, i+1, n_points )
        {
            fu = thin_plate_spline_U( positions[i], positions[j], n_dims );
            ML[j][i] = fu;
            ML[i][j] = fu;
        }
    }
 
    /* set the rest of the L matrix */

    for_less( i, 0, n_points )
    {
        ML[n_points][i] = 1.0;
        ML[i][n_points] = 1.0;
        for_less( j, 0, n_dims )
        {
            ML[n_points+1+j][i] = positions[i][j];
            ML[i][n_points+1+j] = positions[i][j];
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : calculate_weights
@INPUT      : YM
              INVML
              n_points
              n_dims
              n_values
@OUTPUT     : INVMLY
@RETURNS    : 
@DESCRIPTION: Computes the weights by multiplying the VIO_Y matrix by the
              inverse of the L matrix.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  calculate_weights(
    VIO_Real    **YM,
    VIO_Real    **INVML,
    VIO_Real    **INVMLY,
    int     n_points,
    int     n_dims,
    int     n_values )
{
    /* L_{-1} VIO_Y = (W|a_{1}, a_{x}, a_{y})^T. */

    matrix_multiply( n_points + n_dims + 1, n_points, n_values,
                     INVML, YM, INVMLY );
}
