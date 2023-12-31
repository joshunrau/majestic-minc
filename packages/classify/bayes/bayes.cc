/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1997, John Sled
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: bayes.cc,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2002-03-20 22:16:35 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : bayes.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: bayesian classifier
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 21, 1996 (John Sled)
@MODIFIED   : $Log: bayes.cc,v $
@MODIFIED   : Revision 1.1  2002-03-20 22:16:35  jason
@MODIFIED   : Initial revision
@MODIFIED   :
@MODIFIED   : Revision 1.2  2000/08/16 18:18:28  jgsled
@MODIFIED   : Fixed problem in which the normalization factors for the exponentials
@MODIFIED   : were not being initialized when the training data was read from a file.
@MODIFIED   :
@MODIFIED   : Revision 1.1.1.1  1997/02/11 00:06:42  alex
@MODIFIED   : Sources for classify, copied from Vasken Kollokian
@MODIFIED   :
 * Revision 1.1  1996/08/29  03:54:52  vasco
 * Initial revision
 *
---------------------------------------------------------------------------- */
extern "C" {
#include <volume_io.h>
#include <limits.h>
#include <math.h>
}
#include "../class_globals.h" 


void bayesian_allocate_memory(void);
VIO_Real matrix_determinant(int dimension, VIO_Real **matrix);
void scale_matrices(VIO_Real ***matrix, VIO_Real scale);


/* locally defined global variables */
static VIO_Real   **mean_feature_matrix;   /* matrix to reflect mean features */
static VIO_Real   ***covariance_matrix;    /* covariance matrices for each class */
static VIO_Real   ***inv_covariance_matrix;   /* inverse of covariance matrices 
                                      for each class */
static int    *mean_feature_class_vector; 
                                /* vector to reflect mean feature classes */
static VIO_Real   *normalize_class_vector;    
                               /* normalization factor for each class */
static VIO_Real   **distance_vector;          /* temporary storage for x - u */
static VIO_Real   *fuzzy_bayes_vector;
static VIO_Real   stdev_scale_factor;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : bayesian_init_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 21, 1996 (John Sled)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void bayesian_init_training(char *param_filename /* parameter is ignored */)
{
  FILE *param_file;

  /* if load_training is not used, then allocate space for struct */
  if ( !load_train_filename ) {

    bayesian_allocate_memory();
  }

  /* check to see if the filename is there */
  if ( param_filename && !file_exists(param_filename)  ) {

    (void) fprintf(stderr,"File `%s' doesn't exist !\n ", param_filename);
    exit(EXIT_FAILURE);

  }

  if( !param_filename) {
    stdev_scale_factor = 1.0;
  }
  else {
    if (verbose) 
      fprintf(stdout, "Loading the parameter file %s\n", param_filename);      
    
    /* open the parameter file, and read the values  */
    param_file = fopen(param_filename, "r");
    
    if ( param_file == NULL) {

      fprintf(stderr, "Cannot open %s\n", param_filename);
      exit(EXIT_FAILURE);
    }
    
   /* scan for standard deviation scale factor */
    fscanf( param_file, "scale=%lf\n", &stdev_scale_factor);
     
    fclose(param_file);

    if ( stdev_scale_factor <= 0 ) {
      fprintf(stderr, "Scale parameter cannot be zero or negative.\n");
      exit(EXIT_FAILURE);
    }
  }

  if ( debug > 2) {

    fprintf(stderr, "apriori flag = %s\n", 
            (apriori)? "TRUE": "FALSE");      
    fprintf(stderr, "scale factor for deviation = %g\n", stdev_scale_factor);
  }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : bayesian_allocate_memory
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 21, 1996 (John Sled)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void bayesian_allocate_memory(void)
{
  /* reserve area for the mean feature matrix */
    VIO_ALLOC2D(mean_feature_matrix, num_classes, num_features);

    /* reserve area for the mean_feature_class_vector */
    ALLOC(mean_feature_class_vector, num_classes);

    /* reserve area for the normalize_class_vector */
    ALLOC(normalize_class_vector, num_classes);

    /* reserve area for the fuzzy_bayes_vector */
    ALLOC( fuzzy_bayes_vector, num_classes );

    /* reserve area for the covariance matrices */
    VIO_ALLOC3D(covariance_matrix, num_classes, num_features, num_features);
    VIO_ALLOC3D(inv_covariance_matrix, num_classes, num_features, num_features);
    
    /* reserve working space for distance vectors */
    VIO_ALLOC2D(distance_vector, num_classes, num_features);

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : bayesian_train_samples
@INPUT      : 
@OUTPUT     : 
@RETURNS    : ?
@DESCRIPTION: takes a feature matrix and trains a classifier on it.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 21, 1996 (John Sled)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void bayesian_train_samples(void)
{
  
  int      i, j, k, l;         /* counters - samples, features, classes */  

  if (verbose)
    (void) fprintf(stderr, "Training samples ... \n");


  /* initialize mean and covariance feature matrices,
     feature class vector & num of class samples */
  for_less( i, 0, num_classes) {

    mean_feature_class_vector[i] = INT_MAX;  /* big number to denote vacancy */

    /* check that there are enough samples */
    if(class_count[i] < num_features + 1) {
      (void) fprintf(stderr, "Not enough samples to train classifier\n");
      exit(EXIT_FAILURE);
    }


    for_less( j, 0, num_features) { 
      mean_feature_matrix[i][j] = 0.0;
      for_less(k, 0, num_features)
        covariance_matrix[i][j][k] = 0.0;
    }
  }

  /* Compute sample means */
   
  /* repeat for the total number of samples */
  for_less( i, 0, num_samples) {  

    /* repeat for the total number of classes */   
    for_less( l, 0, num_classes ) { 
    
      if ( mean_feature_class_vector[l] == INT_MAX ||    /* unoccupied spot */
	   mean_feature_class_vector[l] == class_column[i] ) { 
        
	for_less( j, 0, num_features ) 
	  mean_feature_matrix[l][j] += feature_matrix[i][j];

	/* if unoccupied, then assign class */
	if ( mean_feature_class_vector[l] == INT_MAX ) 
	  mean_feature_class_vector[l] = class_column[i];

	break;
      }

    } /* for l */

  } /* for i */


  if (verbose)
    (void) fprintf(stderr, "Generating mean feature matrix ...\n\n");

  for_less( i, 0, num_classes) 
    for_less( j, 0, num_features) 
      if (class_count[i] != 0)
	mean_feature_matrix[i][j] /=  (VIO_Real) class_count[i];

  if (debug > 2 ) {

    fprintf( stderr, "Printing mean_feature_matrix ...\n");

    for_less( i, 0, num_classes) {
      for_less( j, 0, num_features) 
	fprintf( stderr, "%lf ", mean_feature_matrix[i][j]);
      fprintf( stderr, "%d\n", mean_feature_class_vector[i]);      
    }

    fprintf( stderr, "-----\n");

  }

  /* Compute sample covariances */

  if (verbose)
    (void) fprintf(stderr, "Generating covariance feature matrix ...\n\n");
   
  /* repeat for the total number of samples */
  for_less( i, 0, num_samples) {  

    /* repeat for the total number of classes */   
    for_less( l, 0, num_classes ) { 
    
      if ( mean_feature_class_vector[l] == class_column[i] ) { 
	
        for_less( j, 0, num_features ) 
          for_less( k, j, num_features ) 
            covariance_matrix[l][j][k] += 
            (feature_matrix[i][j] - mean_feature_matrix[l][j]) * 
              (feature_matrix[i][k] - mean_feature_matrix[l][k]);
      }

    } /* for l */

  } /* for i */
  
  /* fill in lower half of covariance matrixes */
  for_less( l, 0, num_classes ) 
    for_less( j, 0, num_features ) 
      for_less( k, j+1, num_features )
        covariance_matrix[l][k][j] = covariance_matrix[l][j][k]; 
  
  /* normalize for number of samples taking in account uncertainty in mean */
  for_less( l, 0, num_classes )  
    for_less( j, 0, num_features )  
      for_less( k, 0, num_features ) 
        covariance_matrix[l][j][k] /= ((VIO_Real) class_count[l] - 1);

  if (debug > 2 ) {

    fprintf( stderr, "Printing covariance_matrix ...\n");

    for_less( i, 0, num_classes) {
      fprintf( stderr, "\nCovariance for class %d\n", 
               mean_feature_class_vector[i]);
      for_less( j, 0, num_features) { 
        for_less( k, 0, num_features) 
          fprintf( stderr, "%lf ", covariance_matrix[i][j][k]);
        fprintf(stderr, "\n");
      }
      fprintf( stderr, "-----\n");
    }
  }

  if (verbose)
    (void) fprintf(stderr, "Inverting covariance feature matrices ...\n\n");

  for_less( l, 0, num_classes )  
    if(!invert_square_matrix(num_features, covariance_matrix[l], 
                             inv_covariance_matrix[l])) {
      (void) fprintf(stderr, "Covariance matrix for class %d is singular,"
                     " Training of classifier failed.\n", l);
              
      exit(EXIT_FAILURE);
    }
     
  /* compute normalization factors */
  for_less( l, 0, num_classes ) {
    normalize_class_vector[l] = 
      sqrt(matrix_determinant(num_features, inv_covariance_matrix[l]) /
           pow(2.0 * M_PI, num_features)); 
  }
  if (debug > 6 ) {

    fprintf( stderr, "normalize_class_vector[] = ");

    for_less( l, 0, num_classes ) {

      fprintf( stderr, "%g ", normalize_class_vector[l]);
    }
    fprintf( stderr, "\n");

  }

  /* put in scale factor for standard deviation */
  scale_matrices(inv_covariance_matrix, 
                 1.0/(stdev_scale_factor*stdev_scale_factor));
  if (debug > 2 ) {
    fprintf( stderr, "Scaling inverse covariance by %g\n",
             1.0/(stdev_scale_factor*stdev_scale_factor));
  }

}
 
/* ----------------------------- MNI Header -----------------------------------
@NAME       : scale_matrices
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: scale matrices by given factor
@METHOD     : 
@GLOBALS    : num_classes, num_features
@CALLS      : 
@CREATED    : August 21, 1996 (John Sled)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void scale_matrices(VIO_Real ***matrix, VIO_Real scale)
{
  int i, j, k;

  if(matrix) {
    for_less(i, 0, num_classes)
      for_less(j, 0, num_features)
        for_less(k, 0, num_features)
          matrix[i][j][k] *= scale;
  }
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : bayesian_classify_sample
@INPUT      : 
@OUTPUT     : 
@RETURNS    : sample class
@DESCRIPTION: given a feature vector and its size, return a class
@METHOD     : compute   Wi  where  Wi = Pi / sum(Pi)
                 Pi = (det(Si)/(2 pi)^d)^0.5 exp(-0.5 (x-u) Si (x -u)) p(i)
                 where
                    Si is the inverse of the ith covariance matrix
                    d  is the number of features
                    u  is the class mean vector
                    x  is the observed feature vector
                    p(i) is the apriori probability that tissue is ith class
@GLOBALS    : feature_vector, apriori
@CALLS      : 
@CREATED    : August 21, 1996 (John Sled)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void bayesian_classify_sample(int *class_num, double *class_prob, 
                              int *class_labels)
{

  int           i, j, k;               /* counters */
  VIO_Real          exponent, sum, max_prob;
  int           class_index;

  /* calculate distance vectors */
  for_less ( i, 0, num_classes ) 
    for_less ( j, 0, num_features) 
      distance_vector[i][j] = feature_vector[j] - mean_feature_matrix[i][j];

  /* caluculate weights for each class */
  sum = 0.0;  /* sum of weights */
  for_less (i, 0, num_classes ) {
    
    exponent = 0.0;
    
    /* do diagonal */
    for_less(j, 0, num_features) {
      exponent += distance_vector[i][j] * distance_vector[i][j] *
                    inv_covariance_matrix[i][j][j];
    }
    /* do off diagonal once and count it twice */
    for_less(j, 0, num_features)
      for_less(k, j+1, num_features) {
        exponent += 2.0 * distance_vector[i][j] * 
                       inv_covariance_matrix[i][j][k] * distance_vector[i][k];
      } 
    
    fuzzy_bayes_vector[i] = normalize_class_vector[i] * exp(-0.5*exponent);

    if (apriori) { /* modify probability using spatial priors */

      /* fuzzy_bayes_vector[i] *= apriori_vector[i]; */
      fuzzy_bayes_vector[i] *= apriori_vector[mean_feature_class_vector[i] - 1]; //change from Daniel
    }

    sum += fuzzy_bayes_vector[i];
  }

  /* normalize weights and select maximum */
  class_index = 0;
  max_prob = 0.0;
  if(sum > 0) {
    for_less(i, 0, num_classes) {
      fuzzy_bayes_vector[i] /= sum;
      /* check for new most probably class */
      if(fuzzy_bayes_vector[i] > max_prob) {
        max_prob = fuzzy_bayes_vector[i];
        class_index = i;
      }
    }
  }

  /* if fuzzy values are expected, copy them over */
  if ( class_prob )
    for_less( i, 0, num_classes ) {

      class_prob[i] = fuzzy_bayes_vector[i];
      class_labels[i] = mean_feature_class_vector[i]; 
    }


  /* set discrete classification */
  *class_num = mean_feature_class_vector[class_index];
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : bayesian_load_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 21, 1996 (John Sled)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void bayesian_load_training(char *load_train_filename)
{
  int i, j, k;
  VIO_Real feature_value;

  FILE *learn_file;

  /* open the input training  file */
  learn_file = fopen(load_train_filename, "r");

  if ( learn_file == NULL) {
    fprintf(stderr, "Cannot open %s\n", load_train_filename);
    exit(EXIT_FAILURE);
  }

  /* scan for the number of features */
  fscanf( learn_file, "num_of_features = %d\n", &num_features);

  /* scan for the max number of classes M\n*/
  fscanf( learn_file, "num_of_classes = %d\n", &num_classes);

  if (debug) {
    fprintf( stderr, "num_of_features = %d\n", num_features);
    fprintf( stderr, "num_of_classes = %d\n", num_classes);
  }

  bayesian_allocate_memory();

  /* reserve a character pointer ( char *) for each class name */
  ALLOC( class_name, num_classes ); 
 
  /* load feature matrix from learn file */
  if (verbose)
    fprintf(stderr, "Loading the training file...\n");

  for_less( k, 0, num_classes) {
    
    fscanf(learn_file, "class = %d\n", &mean_feature_class_vector[k]);    

    /* allocate some space for each class_name[k], 5 for ex and
       simulate itoa, this could be replaced by a more elegant looking
       code, when I have time */
    ALLOC( class_name[k], 5);
    sprintf( class_name[k], "%d", mean_feature_class_vector[k]); 

    for_less( i, 0, num_features) {
      fscanf(learn_file, "%lf\n", &feature_value);
      mean_feature_matrix[k][i] = feature_value ;
    }

    for_less( i, 0, num_features) 
      for_less( j, 0, num_features) {
        fscanf(learn_file, "%lf\n", &feature_value);
        inv_covariance_matrix[k][i][j] = feature_value ;
    }

  }

  if (debug > 2 ) {

    fprintf( stderr, "Printing mean_feature_matrix ...\n");
    for_less( i, 0, num_classes) {
      for_less( j, 0, num_features) 
	fprintf( stderr, "%f ", mean_feature_matrix[i][j]);
      fprintf( stderr, "%d\n", mean_feature_class_vector[i]);      
    }

    fprintf( stderr, "-----\n");

    fprintf( stderr, "Printing inverse covariance_matrix ...\n");

    for_less( i, 0, num_classes) {
      fprintf( stderr, "\nCovariance for class %d\n", 
               mean_feature_class_vector[i]);
      for_less( j, 0, num_features) { 
        for_less( k, 0, num_features) 
          fprintf( stderr, "%lf ", inv_covariance_matrix[i][j][k]);
        fprintf(stderr, "\n");
      }
      fprintf( stderr, "-----\n");
    }
  }

  fclose(learn_file);

  /* put in scale factor for standard deviation */
  scale_matrices(inv_covariance_matrix, 
                 1.0/(stdev_scale_factor*stdev_scale_factor));
  if (debug > 2 ) {
    fprintf( stderr, "Scaling inverse covariance by %g\n",
             1.0/(stdev_scale_factor*stdev_scale_factor));
  }

  /* compute normalization factors */
  int l;
  for_less( l, 0, num_classes ) {
    normalize_class_vector[l] = 
      sqrt(matrix_determinant(num_features, inv_covariance_matrix[l]) /
           pow(2.0 * M_PI, num_features)); 
  }
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : bayesian_save_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 21, 1996 (John Sled)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void bayesian_save_training(char *save_train_filename)
{

  int i, j, k;                 /* counters */
  FILE *train_file;         /* save filename */



  /* remove scale factor for standard deviation */
  scale_matrices(inv_covariance_matrix, 
                 (stdev_scale_factor*stdev_scale_factor));
  if (debug > 2 ) {
    fprintf( stderr, "Scaling inverse covariance by %g\n",
             1.0/(stdev_scale_factor*stdev_scale_factor));
  }

  train_file = fopen(save_train_filename, "w");


  /* see if the file could be opened */
  if ( train_file == NULL) {

    printf("Cannot open %s\n", save_train_filename);
    exit(EXIT_FAILURE);

  }

  /* store in the file the number of features trained on */
  fprintf(train_file, "num_of_features = %d\n", num_features);

  /* store in the file the max number of classes trained on */
  fprintf(train_file, "num_of_classes = %d\n", num_classes);
  
  for_less ( i, 0, num_classes ) {
    fprintf(train_file, "class = %d\n", mean_feature_class_vector[i]); 
    for_less ( j, 0, num_features )
      fprintf(train_file, "%lf ", mean_feature_matrix[i][j]); 
    fprintf(train_file, "\n\n");
    for_less ( j, 0, num_features )
      for_less( k, 0, num_features )
        fprintf(train_file, "%lf ", inv_covariance_matrix[i][j][k]); 
    fprintf(train_file, "\n");
  }

  fclose(train_file);


  /* put in scale factor for standard deviation */
  scale_matrices(inv_covariance_matrix, 
                 (stdev_scale_factor*stdev_scale_factor));
  if (debug > 2 ) {
    fprintf( stderr, "Scaling inverse covariance by %g\n",
             1.0/(stdev_scale_factor*stdev_scale_factor));
  }
}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : matrix_determinant
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : the simple slow method that you would do by hand
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 21, 1996 (John Sled)
@MODIFIED   : 
---------------------------------------------------------------------------- */
VIO_Real matrix_determinant(int dimension, VIO_Real **matrix)
{
  VIO_Real determinant = 0;
  unsigned col;
  int factor, i, j;
  VIO_Real **sub_matrix;

  if (dimension > 1) {
    VIO_ALLOC2D(sub_matrix, dimension-1, dimension-1);
    for(col = 0, factor = 1; col < dimension; col++, factor *= -1) {
      /* create minor */
      for_less(i, 1, dimension) {
        for_less(j, 0, col) 
          sub_matrix[i-1][j] = matrix[i][j];
        for_less(j, col+1, dimension)
          sub_matrix[i-1][j-1] = matrix[i][j];
      }
      /* compute recursively */
      determinant += matrix[0][col] * 
        matrix_determinant(dimension-1, sub_matrix) * factor;
    }
    VIO_FREE2D(sub_matrix);
  }
  else 
    determinant = matrix[0][0];

  return determinant; 
}







