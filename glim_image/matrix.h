/* ----------------------------- MNI Header -----------------------------------
@NAME       : matrix.h
@INPUT      : 
@OUTPUT     : (nothing)
@RETURNS    : 
@DESCRIPTION: Header file for general purpose matrix functions
@CREATED    : Sept 11, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */

#ifndef MATRIX_H
#define MATRIX_H

#define COMPILING_MATRIX 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#ifndef lint
static char matrix_h_rcsid[] = "$Header: /private-cvsroot/statistics/glim_image/matrix.h,v 1.2 2005-04-07 07:38:21 rotor Exp $ MINC (MNI)";
#endif

#ifndef INVALID_DATA
   #define INVALID_DATA (-DBL_MAX)
#endif

typedef struct {
   int num_rows;
   int num_columns;
   double **values;
} Matrix;

typedef enum {MATRIX_NEITHER, MATRIX_LEFT, MATRIX_RIGHT, MATRIX_BOTH} Trans_Code;

Matrix *create_matrix(int num_rows, int num_columns); 

void delete_matrix(Matrix *matrix);

void multiply_matrix (Matrix *matrix_result, Matrix *matrix_left,
                      Matrix *matrix_right, Trans_Code transpose);

void multiply_matrix_sym(Matrix *matrix_result, Matrix *matrix_left, 
                         Matrix *matrix_right , Trans_Code transpose);

void multiply_diagonal_left (Matrix *matrix_result, Matrix *matrix_left,
                             Matrix *matrix_right , Trans_Code transpose);

void multiply_diagonal_right (Matrix *matrix_result, Matrix *matrix_left,
                              Matrix *matrix_right , Trans_Code transpose);

int invert_matrix_sym (Matrix *matrix_result, Matrix *matrix);

void transpose_matrix (Matrix *matrix_result, Matrix *matrix);

void copy_matrix (Matrix*matrix_result, Matrix *matrix, int num_rows, 
                  int num_columns, double scalar , Trans_Code transpose);

void add_matrix (Matrix *matrix_result, Matrix *matrix_a, Matrix *matrix_b , Trans_Code transpose);

void sub_matrix (Matrix *matrix_result, Matrix *matrix_pos,
                  Matrix *matrix_neg , Trans_Code transpose); 
       
int print_matrix(FILE *strm, Matrix *matrix);

void sprint_matrix(char **str, Matrix *matrix);

double normalize_vector(Matrix *vector);

double determinant(Matrix *matrix);

void svd_matrix(Matrix **singval, Matrix **buffer, 
          Matrix *matrix, int free_all);

#endif
