/* ----------------------------- MNI Header -----------------------------------
@NAME       : smoothness.h
@INPUT      : (nothing)
@OUTPUT     : (nothing)
@RETURNS    : 
@DESCRIPTION: routines to estimate smoothness of random field for glim_image
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

#ifndef INVALID_DATA
   #define INVALID_DATA (-DBL_MAX)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <minc.h>
#include <volume_io.h>
#ifndef COMPILING_MATRIX
   #include "matrix.h"
#endif

/* Typedefs */

typedef struct {
   Matrix *lambda;
   Matrix *lambda_pool;
   double num_calc;
   int out_id;
   int is_gaussian;
   int ring_size;
   double **data;
   double *deriv;
   double constant;
   double fwhm;
   double smoothness;
   double int_lambda;
   double search_vol;
   double smoothness_int;
   double fwhm_int;
   double threshold;
   double fwhm_pool;
   double smoothness_pool;
   double threshold_pool;
   char *outfile;
} Fwhm_Info;

typedef struct {
   int is_gaussian;
   int num_buffers;
   int num_dim;
   int num_test;
   int is_avg;
   int *index;
   int num_check;
   double test_voxel;
   double constant;
   double num_calc;
   int current;
   int *check;
   int ring_size;
   int total;
   int *element_sizes;
   int *coord;
   int *sizes;
   double *step;
   double *test;
   Fwhm_Info *fwhm_avg;
   Fwhm_Info *fwhm_gaussian;
   Fwhm_Info **fwhm_general;
   Fwhm_Info *fwhm_simple;
   char **beta_tmpfiles;
   char *scale_tmpfile;
} Lambda;

/* Structure typedefs */

Lambda *create_lambda_buffer(int num_buffers, int num_test, int num_dim,
                             int *sizes, int is_gaussian, int is_avg_fwhm,
                             char *fwhm_gaussian_file);

void calculate_smoothness(Lambda *lambda_buffer);

int delete_lambda_buffer(Lambda *lambda_buffer);

int convert_coord_to_index(int *coord, Lambda *lambda_buffer);

void convert_index_to_coord(int index, int *coord, Lambda *lambda_buffer);

Fwhm_Info *create_general_buffer(Lambda *lambda_buffer, int is_gaussian);

void calculate_smoothness_ind(Fwhm_Info *fwhm_info, Lambda *lambda_buffer);

int average_smoothness(char *tmp_lambda_file, char *lambda_file, 
                        double scalar);

double pvalue_t(double t, double int_fwhm, int df);

double thresh_t(double pvalue, double int_fwhm, int df);
   



