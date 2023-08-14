/* ----------------------------- MNI Header -----------------------------------
@NAME       : matrix.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: source code for the functions in matrix.h, general
              purpose matrix functions, mainly dealing with symmetric matrices
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sept. 11, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */

#include "matrix.h"
#include "glim.h"

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_matrix
@INPUT      : matrix_ptr - pointer to Matrix
              num_rows
              num_columns
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: creates a matrix with nrows and ncolumns
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
Matrix *create_matrix(int num_rows, int num_columns)
{
   Matrix *matrix_result;
   double **values;
   int i, j;
   int error_flag = FALSE;
   
   matrix_result = GI_MALLOC(sizeof(Matrix));
   matrix_result->num_rows = num_rows;
   matrix_result->num_columns = num_columns; 

   values = GI_MALLOC(sizeof(*values) * num_rows);
   
   for(i=0; i < num_rows; i++) {
      values[i] = GI_MALLOC(sizeof(double) * num_columns);
      for(j=0; j < num_columns; j++) {
         values[i][j] = 0.0;
      }
   }

   matrix_result->values = values; 

   return matrix_result;

}
	
/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_matrix
@INPUT      : matrix - pointer to Matrix
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: deletes a Matrix
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void delete_matrix(Matrix *matrix)
{
   int i;

   if(matrix != NULL) {
      for(i=0 ; i<matrix->num_rows; i++) {
         GI_FREE(matrix->values[i]);
      }
      GI_FREE(matrix->values);
   }
   else {
      fprintf(stderr,"\nWarning: tried to delete NULL matrix.\n");
      return;
   }

   GI_FREE(matrix);

   return;

}
/* ----------------------------- MNI Header -----------------------------------
@NAME       : multiply_matrix
@INPUT      : matrix_left - pointer to Matrix
              matrix_right - pointer to Matrix
              matrix_result - pointer to Matrix
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: multiplies two matrices together
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void multiply_matrix(Matrix *matrix_result, Matrix *matrix_left,
                     Matrix *matrix_right, Trans_Code transpose) 
{
   int i,j,k;
   double tmp_value;
   double value;

   switch(transpose) {
   case MATRIX_NEITHER:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_right->num_columns) ||
         (matrix_result->num_rows != matrix_left->num_rows)) {
         fprintf(stderr,"\nError in multiply_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_left and matrix_right.\n");
         exit(EXIT_FAILURE); 
      } 

      if(matrix_left->num_columns != matrix_right->num_rows) {
         fprintf(stderr,"\nError in multiply_matrix:\n");
         fprintf(stderr,"Num_rows of matrix_right and num_columns of matrix_left don't match. Tranpsose code: %d.\n", transpose);
         exit(EXIT_FAILURE); 
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in multiply_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");
         exit(EXIT_FAILURE);
      }

      for(i=0 ; i<matrix_result->num_rows ; i++) {
         for(j=0 ; j<matrix_result->num_columns ; j++){
            matrix_result->values[i][j] = 0.0;
            for(k=0 ; k<matrix_right->num_rows ; k++) {
               matrix_result->values[i][j] += (matrix_left->values[i][k] *
                                               matrix_right->values[k][j]);
            }
         }
      } 
	 break;
   
   case MATRIX_LEFT:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_right->num_columns) ||
         (matrix_result->num_rows != matrix_left->num_columns)) {
         fprintf(stderr,"\nError in multiply_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_left and matrix_right.\n");
         exit(EXIT_FAILURE); 
      }

      if(matrix_left->num_rows != matrix_right->num_rows) {
         fprintf(stderr,"\nError in multiply_matrix:\n");
         fprintf(stderr,"Num_rows of matrix_right and num_columns of matrix_left don't match. Tranpsose code: %d.\n", transpose);
         exit(EXIT_FAILURE); 
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in multiply_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");
         exit(EXIT_FAILURE);
      }

      for(i=0 ; i<matrix_result->num_rows ; i++) {
         for(j=0 ; j<matrix_result->num_columns ; j++){
            matrix_result->values[i][j] = 0.0;
            for(k=0 ; k<matrix_right->num_rows ; k++) {
               matrix_result->values[i][j] += (matrix_left->values[k][i] *
                                               matrix_right->values[k][j]);
            }
         }
      } 
      break;
   
   case MATRIX_RIGHT:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_right->num_rows) ||
         (matrix_result->num_rows != matrix_left->num_rows)) {
         fprintf(stderr,"\nError in multiply_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_left and matrix_right.\n");
         exit(EXIT_FAILURE); 
      }

      if(matrix_left->num_columns != matrix_right->num_columns) {
         fprintf(stderr,"\nError in multiply_matrix:\n");
         fprintf(stderr,"Num_columns of matrix_right and num_columns of matrix_left don't match. Transpose code: %d.\n", transpose);
         exit(EXIT_FAILURE); 
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in multiply_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");
         exit(EXIT_FAILURE);
      }

      for(i=0 ; i<matrix_result->num_rows ; i++) {
         for(j=0 ; j<matrix_result->num_columns ; j++){
            matrix_result->values[i][j] = 0.0;
            for(k=0 ; k<matrix_right->num_columns ; k++) {
               matrix_result->values[i][j] += (matrix_left->values[i][k] *
                                               matrix_right->values[j][k]);
            }
         }
      } 
      break;
      
   case MATRIX_BOTH:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_left->num_columns) ||
         (matrix_result->num_rows != matrix_right->num_rows)) {
         fprintf(stderr,"\nError in multiply_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_left and matrix_right.\n");
         exit(EXIT_FAILURE); 
      }

      if(matrix_left->num_rows != matrix_right->num_columns) {
         fprintf(stderr,"\nError in multiply_matrix:\n");
         fprintf(stderr,"Num_rows of matrix_left and num_columns of matrix_right don't match. Transpose code: %d.\n", transpose);
         exit(EXIT_FAILURE); 
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in multiply_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");
         exit(EXIT_FAILURE);
      }

      for(i=0 ; i<matrix_result->num_rows ; i++) {
         for(j=0 ; j<matrix_result->num_columns ; j++){
            matrix_result->values[i][j] = 0.0;
            for(k=0 ; k<matrix_right->num_rows ; k++) {
               matrix_result->values[i][j] += (matrix_left->values[k][i] *
                                               matrix_right->values[j][k]);
            }
         }
      } 
      break;

   default:
      fprintf(stderr,"\nError: %d is not valid for transpose argument of multiply_matrix.\n", transpose);
      fprintf(stderr,"\nValid codes are:\n%d   transpose neither \n%d   transpose left matrix \n%d transpose right matrix \n%d transpose both.\n",
              MATRIX_NEITHER, MATRIX_LEFT, MATRIX_RIGHT, MATRIX_BOTH);
      exit(EXIT_FAILURE);
      break;
   }
   
   return; 

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : multiply_matrix_sym
@INPUT      : matrix_left - pointer to Matrix
              matrix_right - pointer to Matrix
              matrix_result - pointer to Matrix
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: multiplies two matrices together whose result is symmetric
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void multiply_matrix_sym(Matrix *matrix_result, Matrix *matrix_left,
                         Matrix *matrix_right, Trans_Code transpose) 
{
   int i,j,k;
   double tmp_value;
   double value;

   switch(transpose) {
   case MATRIX_NEITHER:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_right->num_columns) ||
         (matrix_result->num_rows != matrix_left->num_rows)) {
         fprintf(stderr,"\nError in matrix_matrix_sym:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_left and matrix_right.\n");
         exit(EXIT_FAILURE); 
      } 

      if(matrix_left->num_columns != matrix_right->num_rows) {
         fprintf(stderr,"\nError in matrix_matrix_sym:\n");
         fprintf(stderr,"Num_rows of matrix_right and num_columns of matrix_left don't match. Tranpsose code: %d.\n", transpose);
         exit(EXIT_FAILURE); 
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in matrix_matrix_sym:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");
         exit(EXIT_FAILURE);
      }

      for(i=0 ; i<matrix_result->num_rows ; i++) {
         for(j=0 ; j<=i ; j++){
            matrix_result->values[i][j] = 0.0;
            for(k=0 ; k<matrix_right->num_rows ; k++) {
               matrix_result->values[i][j] += (matrix_left->values[i][k] *
                                               matrix_right->values[k][j]);
            }
            matrix_result->values[j][i] = matrix_result->values[i][j];
         }
      } 
	 break;
   
   case MATRIX_LEFT:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_right->num_columns) ||
         (matrix_result->num_rows != matrix_left->num_columns)) {
         fprintf(stderr,"\nError in matrix_matrix_sym:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_left and matrix_right.\n");
         exit(EXIT_FAILURE); 
      }

      if(matrix_left->num_rows != matrix_right->num_rows) {
         fprintf(stderr,"\nError in matrix_matrix_sym:\n");
         fprintf(stderr,"Num_rows of matrix_right and num_columns of matrix_left don't match. Tranpsose code: %d.\n", transpose);
         exit(EXIT_FAILURE); 
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in matrix_matrix_sym:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");
         exit(EXIT_FAILURE);
      }

      for(i=0 ; i<matrix_result->num_rows ; i++) {
         for(j=0 ; j<=i ; j++){
            matrix_result->values[i][j] = 0.0;
            for(k=0 ; k<matrix_right->num_rows ; k++) {
               matrix_result->values[i][j] += (matrix_left->values[k][i] *
                                               matrix_right->values[k][j]);
            }
            matrix_result->values[j][i] = matrix_result->values[i][j];
         }
      } 
      break;
   
   case MATRIX_RIGHT:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_right->num_rows) ||
         (matrix_result->num_rows != matrix_left->num_rows)) {
         fprintf(stderr,"\nError in matrix_matrix_sym:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_left and matrix_right.\n");
         exit(EXIT_FAILURE); 
      }

      if(matrix_left->num_columns != matrix_right->num_columns) {
         fprintf(stderr,"\nError in matrix_matrix_sym:\n");
         fprintf(stderr,"Num_columns of matrix_right and num_columns of matrix_left don't match. Transpose code: %d.\n", transpose);
         exit(EXIT_FAILURE); 
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in matrix_matrix_sym:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");
         exit(EXIT_FAILURE);
      }

      for(i=0 ; i<matrix_result->num_rows ; i++) {
         for(j=0 ; j<=i ; j++){
            matrix_result->values[i][j] = 0.0;
            for(k=0 ; k<matrix_right->num_rows ; k++) {
               matrix_result->values[i][j] += (matrix_left->values[i][k] *
                                               matrix_right->values[j][k]);
            }
            matrix_result->values[j][i] = matrix_result->values[i][j];
         }
      } 
      break;
      
   case MATRIX_BOTH:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_left->num_columns) ||
         (matrix_result->num_rows != matrix_right->num_rows)) {
         fprintf(stderr,"\nError in matrix_matrix_sym:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_left and matrix_right.\n");
         exit(EXIT_FAILURE); 
      }

      if(matrix_left->num_rows != matrix_right->num_columns) {
         fprintf(stderr,"\nError in matrix_matrix_sym:\n");
         fprintf(stderr,"Num_rows of matrix_left and num_columns of matrix_right don't match. Transpose code: %d.\n", transpose);
         exit(EXIT_FAILURE); 
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in matrix_matrix_sym:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");
         exit(EXIT_FAILURE);
      }

      for(i=0 ; i<matrix_result->num_rows ; i++) {
         for(j=0 ; j<=i ; j++){
            matrix_result->values[i][j] = 0.0;
            for(k=0 ; k<matrix_right->num_rows ; k++) {
               matrix_result->values[i][j] += (matrix_left->values[k][i] *
                                               matrix_right->values[j][k]);
            }
            matrix_result->values[j][i] = matrix_result->values[i][j];
         }
      } 
      break;

   default:
      fprintf(stderr,"\nError: %d is not valid for transpose argument of matrix_matrix_sym.\n", transpose);
      fprintf(stderr,"\nValid codes are:\n%d   transpose neither \n%d   transpose left matrix \n%d transpose right matrix \n%d transpose both.\n",
              MATRIX_NEITHER, MATRIX_LEFT, MATRIX_RIGHT, MATRIX_BOTH);
      fprintf(stderr,"matrix_left:\n");
      print_matrix(stderr, matrix_left);
      fprintf(stderr, "matrix_right:\n");
      print_matrix(stderr, matrix_right);
      exit(EXIT_FAILURE);
      break;
   }
   
   return; 
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : multiply_diagonal_left
@INPUT      : matrix_left - pointer to Matrix
              matrix_right - pointer to Matrix
              matrix_result - pointer to Matrix
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: multiplies the rows of matrix_right term by term by the
              values of matrix_left (matrix_left is an n by 1 matrix)
              identical to multiplying by a diagonal square matrix with 
              elements diagonal->values[i][i] = matrix_left->values[i][0]
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void multiply_diagonal_left(Matrix *matrix_result, Matrix *matrix_left,
                             Matrix *matrix_right, Trans_Code transpose) 
{
   int i,j;
   double value;

   switch(transpose) {
   case MATRIX_NEITHER:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns < matrix_right->num_columns) ||
         (matrix_result->num_rows < matrix_right->num_rows)) {
         fprintf(stderr,"\nError in multiply_diagonal_left\n");
         fprintf(stderr,"Size of matrix_result doesn't match size of matrix_right.\n");
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }

      if(matrix_left->num_rows != matrix_right->num_rows) {
         fprintf(stderr,"\nError in multiply_diagonal_left\n");
         fprintf(stderr,"Num_rows of matrix_left and num_rows of matrix_right don't match\n");
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in multiply_diagonal_left:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");        
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }  

      for(i=0 ; i<matrix_right->num_rows ; i++) {
         for(j=0 ; j<matrix_right->num_columns ; j++){
            matrix_result->values[i][j] = (matrix_right->values[i][j] *
                                           matrix_left->values[i][0]);
         }
      }

   case MATRIX_RIGHT:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns < matrix_right->num_columns) ||
         (matrix_result->num_rows < matrix_right->num_rows)) {
         fprintf(stderr,"\nError in multiply_diagonal_left\n");
         fprintf(stderr,"Size of matrix_result doesn't match size of matrix_right.\n");        
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }

      if(matrix_left->num_rows != matrix_right->num_columns) {
         fprintf(stderr,"\nError in multiply_diagonal_left\n");
         fprintf(stderr,"Num_rows of matrix_left and num_columns of matrix_right don't match\n");      
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in multiply_diagonal_left:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }  

      for(i=0 ; i<matrix_right->num_rows ; i++) {
         for(j=0 ; j<matrix_right->num_columns ; j++){
            matrix_result->values[i][j] = (matrix_right->values[j][i] *
                                           matrix_left->values[i][0]);
         }
      }

   default:
      fprintf(stderr,"\nError: %d is not valid for transpose argument of multiply_diagonal_left.\n", transpose);
      fprintf(stderr,"\nValid codes are:\n%d   transpose neither \n%d transpose right matrix.\n", MATRIX_NEITHER, MATRIX_RIGHT);
      fprintf(stderr,"matrix_left:\n");
      print_matrix(stderr, matrix_left);
      fprintf(stderr, "matrix_right:\n");
      print_matrix(stderr, matrix_right);
      exit(EXIT_FAILURE);
      break;
  }
		

   return;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : multiply_diagonal_right
@INPUT      : matrix_left - pointer to Matrix
              matrix_right - pointer to Matrix
              matrix_result - pointer to Matrix
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: multiplies the columns of matrix_left term by term by the
              values of matrix_right (matrix_right is an n by 1 matrix)
              identical to multiplying by a diagonal square matrix with 
              elements diagonal->values[i][i] = matrix_right->values[i][0]
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void multiply_diagonal_right (Matrix *matrix_result, Matrix *matrix_left,
                              Matrix *matrix_right, Trans_Code transpose) 
{
   int i,j;
   double value;

   switch(transpose) {
   case MATRIX_NEITHER:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns < matrix_right->num_columns) ||
         (matrix_result->num_rows < matrix_right->num_rows)) {
         fprintf(stderr,"\nError in multiply_diagonal_right\n");
         fprintf(stderr,"Size of matrix_result doesn't match size of matrix_right.\n");   
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }

      if(matrix_left->num_rows != matrix_right->num_rows) {
         fprintf(stderr,"\nError in multiply_diagonal_right\n");
         fprintf(stderr,"Num_rows of matrix_left and num_rows of matrix_right don't match\n");   
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in multiply_diagonal_right:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");   
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }  

      for(i=0 ; i<matrix_right->num_rows ; i++) {
         for(j=0 ; j<matrix_right->num_columns ; j++){
            matrix_result->values[i][j] = (matrix_right->values[i][j] *
                                           matrix_left->values[i][0]);
         }
      }

   case MATRIX_LEFT:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns < matrix_right->num_columns) ||
         (matrix_result->num_rows < matrix_right->num_rows)) {
         fprintf(stderr,"\nError in multiply_diagonal_right\n");
         fprintf(stderr,"Size of matrix_result doesn't match size of matrix_right.\n");   
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }

      if(matrix_left->num_rows != matrix_right->num_columns) {
         fprintf(stderr,"\nError in multiply_diagonal_right\n");
         fprintf(stderr,"Num_rows of matrix_left and num_columns of matrix_right don't match\n");       
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }

      if((matrix_result == matrix_left) || (matrix_result == matrix_right)) {
         fprintf(stderr, "\nError in multiply_diagonal_right:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");       
         fprintf(stderr,"matrix_left:\n");
         print_matrix(stderr, matrix_left);
         fprintf(stderr, "matrix_right:\n");
         print_matrix(stderr, matrix_right);
         exit(EXIT_FAILURE);
      }  

      for(i=0 ; i<matrix_right->num_rows ; i++) {
         for(j=0 ; j<matrix_right->num_columns ; j++){
            matrix_result->values[i][j] = (matrix_right->values[j][i] *
                                           matrix_left->values[i][0]);
         }
      }

   default:
      fprintf(stderr,"\nError: %d is not valid for transpose argument of multiply_diagonal_right.\n", transpose);
      fprintf(stderr,"\nValid codes are:\n%d   transpose neither \n%d transpose left matrix.\n", MATRIX_NEITHER, MATRIX_LEFT);   
      fprintf(stderr,"matrix_left:\n");
      print_matrix(stderr, matrix_left);
      fprintf(stderr, "matrix_right:\n");
      print_matrix(stderr, matrix_right);
      exit(EXIT_FAILURE);
      break;
  }
		

   return;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : invert_matrix_sym
@INPUT      : matrix - pointer to Matrix
              matrix
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: inverts positive definite symmetric square matrix using Gauss-Jordan
              symmetric row sweeper
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int invert_matrix_sym(Matrix *matrix_result, Matrix *matrix)
{
   int i,j,k,p;
   int i_tmp1, i_tmp2, j_tmp1, j_tmp2, tmp;
   double a_kk, abs_a_kk;
   double tol = 0.000001;
   double det; /* for 2X2 matrix */
   /* to be removed */

   /* Check matrix is square first */
   if(matrix->num_rows != matrix->num_columns) {
      fprintf(stderr,"\nError in invert_matrix_sym:\n");
      fprintf(stderr,"Matrix isn't square.\n");
      return FALSE;
   }
      
   if(matrix_result->num_rows != matrix_result->num_columns) {
      fprintf(stderr,"\nError in invert_matrix_sym:\n");
      fprintf(stderr,"Matrix_result isn't square.\n");
      return FALSE;
   }

   if(matrix->num_rows == 1) {
      matrix_result->values[0][0] = 1.0 / matrix->values[0][0];
      return TRUE;
   }
   
   if(matrix->num_rows == 2) {
      det = (matrix->values[0][0] * matrix->values[1][1] -
             matrix->values[1][0] * matrix->values[1][0]);
      matrix_result->values[0][0] = matrix->values[1][1] / det;
      matrix_result->values[1][1] = matrix->values[0][0] / det;
      matrix_result->values[1][0] = - matrix->values[1][0] / det;
      matrix_result->values[0][1] = matrix_result->values[1][0];
      return TRUE;
   }
   
   copy_matrix(matrix_result, matrix, matrix->num_rows, matrix->num_rows, 1.0, MATRIX_NEITHER);

   for(k = 0; k < matrix->num_rows; k++) {

      /* Check for collinearity */

      if (k >= 1) {
         for(i = 0; i< k; i++) {
            if (1.0 /(((matrix_result->values[i][k] * 
                        matrix_result->values[i][k]) / 
                       matrix_result->values[k][k])  -
                      matrix_result->values[i][i]) <= 
                matrix->values[i][i] * tol) {
               fprintf(stderr,"\nWarning: Matrix  to be inverted may be singular.\nColumn %d of matrix may be linearly dependent.\n",k+1);
               return FALSE;
            }
         }
      }

      a_kk = matrix_result->values[k][k];
      abs_a_kk = fabs(a_kk);
      for(i=0; i < matrix_result->num_rows; i++) {
         for(j=0; j <= i; j++) {
            if((i > k) && (j == k))
              matrix_result->values[i][j]=(matrix_result->values[j][i]/
                                           abs_a_kk);
            else if ((i == k) && (j != k))
              matrix_result->values[i][j]=(matrix_result->values[j][i]/
                                           abs_a_kk);
            else if ((i == k) && (j == k))
               matrix_result->values[i][j] = -1.0/a_kk;
            else {
               if (i < k) {
                  i_tmp1 = i;
                  j_tmp1 = k;
               }
               else {
                  i_tmp1 = k;
                  j_tmp1 = i;
               }
               
               if (j < k) {
                  i_tmp2 = j;
                  j_tmp2 = k;
               }
               else {
                  i_tmp2 = k;
                  j_tmp2 = j;
               }
                            
              matrix_result->values[i][j]=(matrix_result->values[j][i]-
                                     matrix_result->values[i_tmp1][j_tmp1]*
                                     matrix_result->values[i_tmp2][j_tmp2]/
                                     a_kk);
            }
         }
      }
      for(i=0; i<matrix_result->num_rows; i++) {
         for(j=0; j<i; j++) {
            matrix_result->values[j][i] = matrix_result->values[i][j];
         }
      }
   }

   for(i=0; i<matrix_result->num_rows; i++) {
      for(j=0; j<= i; j++) {
         matrix_result->values[j][i] = -matrix_result->values[i][j];
         if (j != i)
            matrix_result->values[i][j] = matrix_result->values[j][i];
      }
   }

   return TRUE;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : transpose_matrix
@INPUT      : matrix - pointer to Matrix
              matrix_result - pointer to Matrix
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: transposes matrix
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void transpose_matrix (Matrix *matrix_result, Matrix *matrix)
{
   int i,j;
   double value;

    /* Check if num_rows and num_columns are right */
   if((matrix_result->num_columns != matrix->num_rows) || 
      (matrix_result->num_rows != matrix->num_columns)) {
      fprintf(stderr,"\nError in transpose_matrix:\n");
      fprintf(stderr,"Size of matrix_result doesn't match matrix.\n");
      exit(EXIT_FAILURE);
   }

   if(matrix_result == matrix) {
      fprintf(stderr, "\nError in transpose_matrix:\n");
      fprintf(stderr, "Overwriting matrices is not allowed.\n");
      exit(EXIT_FAILURE);
   }

   for(i=0; i<matrix_result->num_rows; i++) {
      for(j=0; j<matrix_result->num_columns; j++) {
         matrix_result->values[i][j] = matrix->values[j][i];
      }
   } 

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : copy_matrix
@INPUT      : matrix - pointer to Matrix
              matrix_result - pointer to Matrix
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: copies the first num_rows rows and num_columns columns of
              a Matrix
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void copy_matrix (Matrix *matrix_result, Matrix *matrix, int num_rows, 
                  int num_columns, double scalar, Trans_Code transpose)
{
   int i,j;
   double value;

   if(matrix_result == matrix) {
      fprintf(stderr, "\nError in copy_matrix:\n");
      fprintf(stderr, "Overwriting matrices is not allowed.\n");   
      fprintf(stderr,"matrix_result:\n");
      print_matrix(stderr, matrix_result);
      fprintf(stderr, "matrix:\n");
      print_matrix(stderr, matrix);
      exit(EXIT_FAILURE);
   }

   switch(transpose) {
   case MATRIX_NEITHER:
      if((num_columns > matrix_result->num_columns) || (num_rows > matrix_result->num_rows)) {
         fprintf(stderr,"\nError in copy_matrix: num_rows or num_columns of matrix_result not big enough.\nnum_rows_in: %d, num_columns_in: %d\nnum_rows_out: %d num_columns_out: %d\n",
            num_rows, num_columns, matrix_result->num_rows, matrix_result->num_columns);
         fprintf(stderr,"matrix_result:\n");
         print_matrix(stderr, matrix_result);
         fprintf(stderr, "matrix:\n");
         print_matrix(stderr, matrix);
         exit(EXIT_FAILURE);
      }
            
      for(i=0; i<num_rows; i++) {
         for(j=0; j<num_columns; j++) {
            matrix_result->values[i][j] = scalar * matrix->values[i][j];
         }
      }
      break;
   
   case MATRIX_LEFT:
         if((num_columns > matrix_result->num_rows) || (num_rows > matrix_result->num_columns)) {
         fprintf(stderr,"\nError in copy_matrix: num_rows or num_columns of matrix_result not big enough.\nnum_rows_in: %d, num_columns_in: %d\nnum_rows_out: %d num_columns_out: %d\n",
            num_rows, num_columns, matrix_result->num_columns, matrix_result->num_rows);
         fprintf(stderr,"matrix_result:\n");
         print_matrix(stderr, matrix_result);
         fprintf(stderr, "matrix:\n");
         print_matrix(stderr, matrix);
         exit(EXIT_FAILURE);
      }
            
      for(i=0; i<num_rows; i++) {
         for(j=0; j<num_columns; j++) {
            matrix_result->values[j][i] = scalar * matrix->values[i][j];
         }
      }
      break;

   case MATRIX_RIGHT:
         if((num_columns > matrix_result->num_rows) || (num_rows > matrix_result->num_columns)) {
         fprintf(stderr,"\nError in copy_matrix: num_rows or num_columns of matrix_result not big enough.\nnum_rows_in: %d, num_columns_in: %d\nnum_rows_out: %d num_columns_out: %d\n",
            num_rows, num_columns, matrix_result->num_columns, matrix_result->num_rows);
         fprintf(stderr,"matrix_result:\n");
         print_matrix(stderr, matrix_result);
         fprintf(stderr, "matrix:\n");
         print_matrix(stderr, matrix);
         exit(EXIT_FAILURE);
      }
            
      for(i=0; i<num_rows; i++) {
         for(j=0; j<num_columns; j++) {
            matrix_result->values[j][i] = scalar * matrix->values[i][j];
         }
      }
      break;

   case MATRIX_BOTH:
         if((num_columns > matrix_result->num_rows) || (num_rows > matrix_result->num_columns)) {
         fprintf(stderr,"\nError in copy_matrix: num_rows or num_columns of matrix_result not big enough.\nnum_rows_in: %d, num_columns_in: %d\nnum_rows_out: %d num_columns_out: %d\n",
            num_rows, num_columns, matrix_result->num_columns, matrix_result->num_rows);
         fprintf(stderr,"matrix_result:\n");
         print_matrix(stderr, matrix_result);
         fprintf(stderr, "matrix:\n");
         print_matrix(stderr, matrix);
         exit(EXIT_FAILURE);
      }
            
      for(i=0; i<num_rows; i++) {
         for(j=0; j<num_columns; j++) {
            matrix_result->values[j][i] = scalar * matrix->values[i][j];
         }
      }
      break;

   default:
      fprintf(stderr,"\nError: %d is not valid for transpose argument of copy_matrix.\n", transpose);
      fprintf(stderr,"\nValid codes are:\n%d   transpose neither \n%d   transpose left matrix \n%d transpose right matrix \n%d transpose both.\n",
              MATRIX_NEITHER, MATRIX_LEFT, MATRIX_RIGHT, MATRIX_BOTH);
      fprintf(stderr,"\nNote that anything other than MATRIX_NEITHER will have copy_matrix transpose the matrix being copied.\n");
      fprintf(stderr,"matrix_result:\n");
      print_matrix(stderr, matrix_result);
      fprintf(stderr, "matrix:\n");
      print_matrix(stderr, matrix);
      exit(EXIT_FAILURE);
      break;
   }

   
   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : add_matrix
@INPUT      : matrix_a - pointer to Matrix
              matrix_b - pointer to Matrix
              matrix_result - pointer to Matrix
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: adds two matrices together
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void add_matrix (Matrix *matrix_result, Matrix *matrix_a, Matrix *matrix_b,
   Trans_Code transpose)
{
   int i,j;
   double value;

   switch(transpose) {    
   case MATRIX_NEITHER:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_a->num_columns) || 
         (matrix_result->num_rows != matrix_a->num_rows)) {
         fprintf(stderr,"\nError in add_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_a and matrix_b.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }
       
      if((matrix_a->num_columns != matrix_b->num_columns) ||
         (matrix_a->num_rows != matrix_b->num_rows)) {
         fprintf(stderr,"\nError in add_matrix:\n");
         fprintf(stderr,"Size of matrix_a doesn't match matrix_b.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }

      if((matrix_result == matrix_a) || (matrix_result == matrix_b)) {
         fprintf(stderr, "\nError in add_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }
   
      for(i=0; i<matrix_a->num_rows; i++) {
         for(j=0; j<matrix_a->num_columns; j++) {
            matrix_result->values[i][j] = matrix_a->values[i][j] +
               matrix_b->values[i][j];
         }
      }
      break;
      
   case MATRIX_LEFT:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_a->num_rows) || 
         (matrix_result->num_rows != matrix_a->num_columns)) {
         fprintf(stderr,"\nError in add_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_a and matrix_b.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }

      if((matrix_a->num_rows != matrix_b->num_columns) ||
         (matrix_a->num_columns != matrix_b->num_rows)) {
         fprintf(stderr,"\nError in add_matrix:\n");
         fprintf(stderr,"Size of matrix_a doesn't match matrix_b.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }
       
      if((matrix_result == matrix_a) || (matrix_result == matrix_b)) {
         fprintf(stderr, "\nError in add_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }
   
      for(i=0; i<matrix_a->num_rows; i++) {
         for(j=0; j<matrix_a->num_columns; j++) {
            matrix_result->values[i][j] = matrix_a->values[j][i] +
               matrix_b->values[i][j];
         }
      }
      break;

   case MATRIX_RIGHT:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_a->num_columns) || 
         (matrix_result->num_rows != matrix_a->num_rows)) {
         fprintf(stderr,"\nError in add_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_a and matrix_b.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }

      if((matrix_a->num_rows != matrix_b->num_columns) ||
         (matrix_a->num_columns != matrix_b->num_rows)) {
         fprintf(stderr,"\nError in add_matrix:\n");
         fprintf(stderr,"Size of matrix_a doesn't match matrix_b.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }
       
      if((matrix_result == matrix_a) || (matrix_result == matrix_b)) {
         fprintf(stderr, "\nError in add_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }
   
      for(i=0; i<matrix_a->num_rows; i++) {
         for(j=0; j<matrix_a->num_columns; j++) {
            matrix_result->values[i][j] = matrix_a->values[i][j] +
               matrix_b->values[j][i];
         }
      }
      break;

   case MATRIX_BOTH:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_a->num_columns) || 
         (matrix_result->num_rows != matrix_a->num_rows)) {
         fprintf(stderr,"\nError in add_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_a and matrix_b.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }

      if((matrix_a->num_columns != matrix_b->num_columns) ||
         (matrix_a->num_rows != matrix_b->num_rows)) {
         fprintf(stderr,"\nError in add_matrix:\n");
         fprintf(stderr,"Size of matrix_a doesn't match matrix_b.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }
       
      if((matrix_result == matrix_a) || (matrix_result == matrix_b)) {
         fprintf(stderr, "\nError in add_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");     
         fprintf(stderr,"matrix_a:\n");
         print_matrix(stderr, matrix_a);
         fprintf(stderr, "matrix_b:\n");
         print_matrix(stderr, matrix_b);
         exit(EXIT_FAILURE);
      }
   
      for(i=0; i<matrix_a->num_rows; i++) {
         for(j=0; j<matrix_a->num_columns; j++) {
            matrix_result->values[i][j] = matrix_a->values[j][i] +
               matrix_b->values[j][i];
         }
      }
      break;

   default:
      fprintf(stderr,"\nError: %d is not valid for transpose argument of add_matrix.\n", transpose);
      fprintf(stderr,"\nValid codes are:\n%d   transpose neither \n%d   transpose left matrix \n%d transpose right matrix \n%d transpose both.\n",
              MATRIX_NEITHER, MATRIX_LEFT, MATRIX_RIGHT, MATRIX_BOTH);
      fprintf(stderr,"matrix_a:\n");
      print_matrix(stderr, matrix_a);
      fprintf(stderr, "matrix_b:\n");
      print_matrix(stderr, matrix_b);
      exit(EXIT_FAILURE);
      break;
   }

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : sub_matrix
@INPUT      : matrix_pos - pointer to Matrix
              matrix_neg - pointer to Matrix
              matrix_result - pointer to Matrix
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: subtracts matrix_neg from matrix_pos
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void sub_matrix(Matrix *matrix_result, Matrix *matrix_pos, Matrix *matrix_neg,
                Trans_Code transpose)
{
   int i,j;
   double value;

   switch(transpose) {    
   case MATRIX_NEITHER:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_pos->num_columns) || 
         (matrix_result->num_rows != matrix_pos->num_rows)) {
         fprintf(stderr,"\nError in sub_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_pos and matrix_neg.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }
       
      if((matrix_pos->num_columns != matrix_neg->num_columns) ||
         (matrix_pos->num_rows != matrix_neg->num_rows)) {
         fprintf(stderr,"\nError in sub_matrix:\n");
         fprintf(stderr,"Size of matrix_pos doesn't match matrix_neg.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }

      if((matrix_result == matrix_pos) || (matrix_result == matrix_neg)) {
         fprintf(stderr, "\nError in sub_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }
   
      for(i=0; i<matrix_pos->num_rows; i++) {
         for(j=0; j<matrix_pos->num_columns; j++) {
            matrix_result->values[i][j] = matrix_pos->values[i][j] -
               matrix_neg->values[i][j];
         }
      }
      break;
      
   case MATRIX_LEFT:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_pos->num_rows) || 
         (matrix_result->num_rows != matrix_pos->num_columns)) {
         fprintf(stderr,"\nError in sub_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_pos and matrix_neg.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }

      if((matrix_pos->num_rows != matrix_neg->num_columns) ||
         (matrix_pos->num_columns != matrix_neg->num_rows)) {
         fprintf(stderr,"\nError in sub_matrix:\n");
         fprintf(stderr,"Size of matrix_pos doesn't match matrix_neg.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }
       
      if((matrix_result == matrix_pos) || (matrix_result == matrix_neg)) {
         fprintf(stderr, "\nError in sub_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }
   
      for(i=0; i<matrix_pos->num_rows; i++) {
         for(j=0; j<matrix_pos->num_columns; j++) {
            matrix_result->values[i][j] = matrix_pos->values[j][i] -
               matrix_neg->values[i][j];
         }
      }
      break;

   case MATRIX_RIGHT:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_pos->num_columns) || 
         (matrix_result->num_rows != matrix_pos->num_rows)) {
         fprintf(stderr,"\nError in sub_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_pos and matrix_neg.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }

      if((matrix_pos->num_rows != matrix_neg->num_columns) ||
         (matrix_pos->num_columns != matrix_neg->num_rows)) {
         fprintf(stderr,"\nError in sub_matrix:\n");
         fprintf(stderr,"Size of matrix_pos doesn't match matrix_neg.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }
       
      if((matrix_result == matrix_pos) || (matrix_result == matrix_neg)) {
         fprintf(stderr, "\nError in sub_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }
   
      for(i=0; i<matrix_pos->num_rows; i++) {
         for(j=0; j<matrix_pos->num_columns; j++) {
            matrix_result->values[i][j] = matrix_pos->values[i][j] -
               matrix_neg->values[j][i];
         }
      }
      break;

   case MATRIX_BOTH:
      /* Check if num_rows and num_columns are right */
      if((matrix_result->num_columns != matrix_pos->num_columns) || 
         (matrix_result->num_rows != matrix_pos->num_rows)) {
         fprintf(stderr,"\nError in sub_matrix:\n");
         fprintf(stderr,"Size of matrix_result doesn't match matrix_pos and matrix_neg.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }

      if((matrix_pos->num_columns != matrix_neg->num_columns) ||
         (matrix_pos->num_rows != matrix_neg->num_rows)) {
         fprintf(stderr,"\nError in sub_matrix:\n");
         fprintf(stderr,"Size of matrix_pos doesn't match matrix_neg.\n");   
         fprintf(stderr,"matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg:\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }
       
      if((matrix_result == matrix_pos) || (matrix_result == matrix_neg)) {
         fprintf(stderr, "\nError in sub_matrix:\n");
         fprintf(stderr, "Overwriting matrices is not allowed.\n");   
         fprintf(stderr,"Matrix_pos:\n");
         print_matrix(stderr, matrix_pos);
         fprintf(stderr, "matrix_neg\n");
         print_matrix(stderr, matrix_neg);
         exit(EXIT_FAILURE);
      }
   
      for(i=0; i<matrix_pos->num_rows; i++) {
         for(j=0; j<matrix_pos->num_columns; j++) {
            matrix_result->values[i][j] = matrix_pos->values[j][i] -
               matrix_neg->values[j][i];
         }
      }
      break;

   default:
      fprintf(stderr,"\nError: %d is not valid for transpose argument of sub_matrix.\n", transpose);
      fprintf(stderr,"\nValid codes are:\n%d   transpose neither \n%d   transpose left matrix \n%d transpose right matrix \n%d transpose both.\n",
              MATRIX_NEITHER, MATRIX_LEFT, MATRIX_RIGHT, MATRIX_BOTH);
      fprintf(stderr,"matrix_pos:\n");
      print_matrix(stderr, matrix_pos);
      fprintf(stderr, "matrix_neg:\n");
      print_matrix(stderr, matrix_neg);
      exit(EXIT_FAILURE);
      break;
   }

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : print_matrix
@INPUT      : pointer to a Matrix
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: routine for printing out contents of a Matrix
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */
int print_matrix(FILE *strm, Matrix *matrix) 
{ 
   int i,j;

   fprintf(strm,"\n");

   for(i=0; i<matrix->num_rows; i++) {
      fprintf(strm,"Row %d: ",i+1);
      for(j=0; j<matrix->num_columns; j++) {
         fprintf(strm,"%3.3f  ",matrix->values[i][j]);
      }
      fprintf(strm,"\n");
   }

   fprintf(strm,"\n");

   return TRUE;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : sprint_matrix
@INPUT      : pointer to a Matrix
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: routine for printing out contents of a Matrix to a string
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */
void sprint_matrix(char **str, Matrix *matrix) 
{ 
   int i,j;
   char *tmp;
   char *tmp2;
   tmp = GI_MALLOC(sizeof(*tmp) * 64);

   tmp2 = GI_MALLOC(sizeof(*tmp2) * ((matrix->num_rows * matrix->num_columns) *
                                  20 + matrix->num_columns * 10 + 1));

   sprintf(tmp2, "\n");

   for(i=0; i<matrix->num_rows; i++) {
      for(j=0; j<matrix->num_columns; j++) {
         sprintf(tmp, "%4.4e  ",matrix->values[i][j]);
         strcat(tmp2, tmp);
      }
      strcat(tmp2,"\n");
   }

   strcat(tmp2,"\n");

   *str = GI_MALLOC(sizeof(**str) * (strlen(tmp2) + 1));

   strcpy(*str, tmp2);

   GI_FREE(tmp);
   GI_FREE(tmp2);

   return;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : normalize_vector
@INPUT      : pointer to a Matrix
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: routine to normalize a 1 by N or an N by 1 vector,
              returns norm of vector
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */
double normalize_vector(Matrix *vector)
{
   int i;
   double norm, tmp_sum;
   
   tmp_sum = 0.0;
   
   if(vector->num_columns == 1) {
      for(i=0; i<vector->num_rows; i++) {
         if (vector->values[i][0] != INVALID_DATA)
            tmp_sum += vector->values[i][0] * vector->values[i][0];
      }
   
      norm = sqrt(tmp_sum);
   
      for(i=0; i<vector->num_rows; i++) {
         vector->values[i][0] /= norm;
      }
   }
   else if(vector->num_rows == 1) {
      for(i=0; i<vector->num_columns; i++) {
         if (vector->values[0][i] != INVALID_DATA)
            tmp_sum += vector->values[0][i] * vector->values[0][i];
      }
   
      norm = sqrt(tmp_sum);
   
      for(i=0; i<vector->num_rows; i++) {
         if (vector->values[i][0] != INVALID_DATA)
            vector->values[0][i] /= norm;
      }
   }
   else {
      fprintf(stderr,"\nError in normalize_vector, matrix is not a vector.\nmatrix:\n");
      print_matrix(stderr, vector);
      exit(EXIT_FAILURE);
   }
   
   return norm;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : determinant
@INPUT      : pointer to a Matrix
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: routine to return determinant of a matrix
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */
double determinant(Matrix *matrix) 
{
   int i, j, k;
   double tmp_sum, determinant, sign;
   Matrix *matrix_tmp;
   
   /* Check if matrix is square */
   if(matrix->num_rows != matrix->num_columns) {
      fprintf(stderr,"\nError in determinant, matrix isn't square.\nmatrix:\n");
      print_matrix(stderr, matrix);
      exit(EXIT_FAILURE);
   }
   
   sign = 1.0;
   
   if(matrix->num_rows > 3) {
      fprintf(stderr,"\nError in determinant, higher than 3D not ready.\nmatrix:\n");
      print_matrix(stderr, matrix);
      exit(EXIT_FAILURE);
   }
   
   if(matrix->num_rows == 3) {
      determinant = (matrix->values[0][0] *
      	             (matrix->values[1][1] * matrix->values[2][2] - 
      	              matrix->values[1][2] * matrix->values[2][1]) -
      	             (matrix->values[0][1] *
      	             (matrix->values[1][0] * matrix->values[2][2] - 
      	              matrix->values[1][2] * matrix->values[2][0])) +
      	             (matrix->values[0][2] *
      	             (matrix->values[1][0] * matrix->values[2][1] -
      	              matrix->values[1][1] * matrix->values[2][0])));
   }
   
   if(matrix->num_rows == 2) {
      determinant = ((matrix->values[0][0] *
      	             matrix->values[1][1]) - 
      	             (matrix->values[1][0] * 
      	             matrix->values[0][1]));
   }

   return determinant;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : svd_matrix
@INPUT      : matrix - pointer to a Matrix
@OUTPUT     : singval - singular values of matrix
              buffer - matrix containing matrices U S and V of svd decomposition
@RETURNS    : 
@DESCRIPTION: routine to calculate svd of a matrix
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 16, 1998 (J. Taylor), source code taken off of:
              usc.edu/pub/C-numanal/svd.c.Z, modified to fit Matrix routines.
@MODIFIED   : 
---------------------------------------------------------------------------- */
 void svd_matrix(Matrix **singval, Matrix **buffer, 
          Matrix *matrix, int free_all)
{
  int i, j, k, est_col_rank, rot_count, sweep_count, sweep_limit;
  double eps, e2, tol, vt, p, h2, x0, y0, q, r, c0, s0, c2, d1, d2;
  int num_rows, num_columns;
  eps = .0000001;  /* some small tolerance value */

  num_rows = matrix->num_rows;
  num_columns = matrix->num_columns;
  
  *buffer = create_matrix(2 * (num_rows + num_columns), num_columns);
  
  *singval = create_matrix(num_rows, 1);
  
  sweep_limit = num_columns/4;
  if (sweep_limit < 6.0)
    sweep_limit = 6;

  sweep_count = 0;

  e2 = 10.0*num_rows*eps*eps;
  tol = eps*.1;
  est_col_rank = num_columns;
  
  for(i=0; i<num_columns; i++) {
     for(j=0; j<num_columns; j++) {
        (*buffer)->values[i][j] = matrix->values[i][j];
     }
  }
    
  for (i=0; i<num_columns; i++) {
    for (j=0; j<num_columns; j++) {
	  (*buffer)->values[num_rows+i][j] = 0.0;
	  (*buffer)->values[num_rows+i][i] = 1.0;
    }
  }
  
  rot_count = est_col_rank*(est_col_rank-1)/2;
 
  while ((rot_count != 0) && (sweep_count <= sweep_limit))
    {
      rot_count = est_col_rank*(est_col_rank-1)/2;
      sweep_count++;
      for (j=0; j<est_col_rank-1; j++)
	{
	  for (k=j+1; k<est_col_rank; k++)
	    {
	      p = 0.0;
	      q = 0.0;
	      r = 0.0;
	      for (i=0; i<num_rows; i++) {
		  x0 = (*buffer)->values[i][j];
		  y0 = (*buffer)->values[i][k];
		  p += x0*y0;
            q += x0*x0; 
            r += y0*y0;
		}
	     (*singval)->values[j][0] = q; 
	     (*singval)->values[k][0] = r;
	     if (q >= r) {
		  if ((q<=e2*(*singval)->values[0][0]) || (fabs(p)<=tol*q))
		    rot_count--;
		  else {
		      p /= q; r = 1 - r/q;
		      vt = sqrt(4*p*p+r*r);
		      c0 = sqrt(.5*(1+r/vt)); 
		      s0 = p/(vt*c0);
		      for (i=0; i<num_rows+num_columns; i++) {
			  d1 = (*buffer)->values[i][j]; 
			  d2 = (*buffer)->values[i][k];
			  (*buffer)->values[i][j] = d1*c0+d2*s0;
			  (*buffer)->values[i][k] = -d1*s0+d2*c0;
			}
		   }
		}
	     else {
		  p /= r;
		  q /= (r-1);
		  vt = sqrt(4*p*p+q*q);
		  s0 = sqrt(.5*(1-q/vt));
		  if (p<0)
		     s0 = -s0;
		  c0 = p/(vt*s0);
		  for (i=0; i<num_rows+num_columns; i++) {
		      d1 = (*buffer)->values[i][j];
		      d2 = (*buffer)->values[i][k];
		      (*buffer)->values[i][j] = d1*c0+d2*s0;
		      (*buffer)->values[i][k] = -d1*s0+d2*c0;
	          }
		  }
	     }
     }
      fprintf(stderr, "Sweep = %d  # of rotations performed = %d\n", sweep_count, rot_count);
/*      while (est_col_rank >= num_columns) {
         while ((*singval)->values[est_col_rank - 1][0] <= (*singval)->values[0][0]*tol+tol*tol) {
	       est_col_rank = est_col_rank - 1;
         }
      } */
   }
   
  for(i=0; i<num_rows; i++) {
      for(j=0; j<num_columns; j++) {
         (*buffer)->values[i][j] /= sqrt((*singval)->values[j][0]);
      }
  } 
     
  if (sweep_count > sweep_limit)
    fprintf(stderr, "Sweep Limit exceeded\n");

  print_matrix(stderr, *singval);
  print_matrix(stderr, *buffer);
  return;
}

