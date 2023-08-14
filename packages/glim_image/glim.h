/* ----------------------------- MNI Header -----------------------------------
@NAME       : glim.h
@INPUT      : 
@OUTPUT     : (nothing)
@RETURNS    : 
@DESCRIPTION: Header file for generalized linear models functions
@CREATED    : Sept 11, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */

#ifndef GLIM_H
#define GLIM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "matrix.h"

#define MAX_ROW_LENGTH 1000
#define MAX_NUM_CONTRASTS 50
#define MAX_ITERATIONS 10
#define MAX_PATH_LENGTH 100
#define MAX_FILENAME_LENGTH 100
#ifndef INVALID_DATA
   #define INVALID_DATA (-DBL_MAX)
#endif

#define ROW_SEPARATOR ';'
#define ENTRY_SEPARATOR ','

#ifndef FALSE
#  define FALSE 0
#endif
#ifndef TRUE
#  define TRUE 1
#endif

#define GI_MALLOC(size) ((void *) malloc(size))
#define GI_FREE(ptr) free(ptr)
#define GI_REALLOC(ptr, size) ((void *) realloc(ptr, size))
#define GI_CALLOC(nelem, elsize) ((void *) calloc(nelem, elsize))


/* Define families, link and variance functions as integers help compare */

/* Families */

typedef enum {GAUSSIAN,GAMMA,POISSON,BINOMIAL,INV_GAUSS,QUASI,CHI_SQ,
                 EXPONENTIAL} Family;

/* Links */

typedef enum {LOGIT, PROBIT, CLOGLOG, IDENTITY, INVERSE, INVERSE_SQ,
                 LOG, SQRT, LINK_DEFAULT} Link;

typedef enum {NORMAL, DERIV, INVER} Link_Mode;    

/* Variance Functions */

typedef enum {CONST, MU, MU_SQUARED, MU_ONE_MINUS_MU, 
                 MU_CUBED, VARIANCE_DEFAULT} Variance_Function;

/* Working Correlation Structure */

typedef enum {DIAGONAL_CORR, AUTO_REG, CORR_DEFAULT} Corr_Struct;

/* Convergence Error Codes */
typedef enum {CONT, CONVERGE, MAX_ITER, SING_MATRIX, ZERO_DEV, MISSING_DATA,
                 IS_SAME, NOT_IN_SEARCH, ZERO_DEG_FREE} Conv_Status;

/* Output Options for Contrast Matrices */

typedef enum {BETA_HAT, STDEV_BETA, F_STAT, T_STAT, P_CORR} Output_Mode;

typedef enum {POOLED_SD, VOXEL_SD} Stdev_Mode;

typedef struct {
   Output_Mode out_mode;
   Stdev_Mode stdev_mode;
   char *outfile;
   char *tmpfile;
   int out_id;
   double pooled_sd;
   double pooled_dev;
   double deg_free;
   double df_offset;
   double maxVoxelSdBasedValue;
} Avg_Info;

typedef struct {
   Avg_Info **avg_info;
   int num_pool;
   int num_test;
   char *deviance_file;
} Avg_Array;

typedef struct {
   Matrix *contrast_matrix;
   Matrix *deriv;
   Output_Mode out_mode;
   Stdev_Mode stdev_mode;
   char *outfile;
   char *tmpfile;
   int column_num;
   double deg_free;
   double pooled_dev;
   double pooled_sd;
   int out_id;
   double int_lambda;
   char *fwhm_outfile;
   double search_vol;
   double num_fwhm;
   double fwhm_int;
   double maxVoxelSdBasedValue;
} Contrast_Matrix;

typedef struct {
   int num_columns;
   int num_contrasts;
   int num_avg;
   int num_test;
   Avg_Array *avg_array;
   Contrast_Matrix **contrast_matrix_array;
   int *test_map;
} Contrast_Mat_Array;

typedef struct {
   char *in_type;
   char *out_type;
   char *outfile_name;
   char *raw_contrast;
   char *stdev_type;
   int is_avg;
} Contrast_Raw;

typedef struct {
   int num_output;
   Contrast_Raw **contrast_array;
   int num_avg;
} Contrast_Raw_Array;

typedef struct {
   int num_values;
   double *values;
} Double_Array;
   
typedef struct {
   int num_iter;
   int max_iter;
   double tol;
   double eps_sing;
   double maxVoxelSdBasedValue;
   int binomial_n;
   int do_fit;
   int is_scale;
   int est_scale;
   int is_canonical;
   int in_mask;
   int use_mu;
   int use_mu_switch;
   int dev_id;
   int avg_dev_id;
   int scale_id;
   int err_id;
   int pool_id;
   int count_id;
   int search_id;
   int deg_free_id;
   int is_fwhm;
   int fit_all;
   int scylla;
   int num_missing;
   Conv_Status status;
   Conv_Status last_status;
} Glm_Control;

typedef struct {
   Matrix *beta_hat;
   Matrix *info_inv;
   Matrix *response;
   Matrix *mu;
   Matrix *z_i;
   Matrix *x_beta;
   Matrix *weight;
   Matrix *sqrt_weight;
   Matrix *g_prime;
   Matrix *output;
   Matrix ***tmp_con;
   Matrix *tmp_weight;
   Matrix *design_matrix;
   Matrix *lm_matrix;
   Matrix *info_mat;
   Matrix *sigma;
   Matrix *resid;
   Matrix *avg_resid;
   Matrix *missing;
   Family family;
   Link link;
   Variance_Function variance_function;
   Corr_Struct corr_struct;
   Glm_Control *control;
   Contrast_Mat_Array *contrast;
   double scale_est;
   double pearson;
   double scale;
   double distance;
   double deviance;
   double num_pool;
   double avg_num_pool;
   double pooled_dev;
   double avg_pooled_dev;
   double deg_free;
   double deg_free_tmp;
   int num_missing;
   int num_var;
   char **family_names;
   char **link_names;
   char **variance_function_names;
   char **corr_struct_names;
   char *deviance_file;
   char *scale_file;
   char *error_file;
   char *search_mask;
} Glm_Object;

typedef struct {
   Contrast_Mat_Array *contrast;
   Matrix *design_matrix;
   Family family;
   Link link;
   Variance_Function variance_function;
   Corr_Struct corr_struct;
   Glm_Control *control;
   double scale;
   double scale_est;
} Initial_Data;

typedef struct {
   char **files;
   double **time;
} Time_Map;

/* Function prototypes */

void delete_tmpfiles(FILE **tmpfile_list);

int get_contrast_from_raw( Contrast_Raw_Array *contrast_list,
                                      Contrast_Mat_Array **contrast_matrices);

int get_contrast_from_matrix(char *contrast_string, int num_columns,
                             Contrast_Matrix **cur_matrix);

int get_contrast_from_file(char *contrast_string, int num_columns,
                             Contrast_Matrix **cur_matrix);

int get_contrast_from_column(char *contrast_string, int num_columns,
                             Contrast_Matrix **cur_matrix);

int get_double_list(char *double_string, Double_Array **double_array);

void delete_glm_matrices(Glm_Object *glm_obj, int is_last);

void calculate_weight(Matrix *matrix_result, Glm_Object *glm_obj);

void calculate_sigma_inv(Matrix *matrix_result, Glm_Object *glm_obj);

void calculate_link(Matrix *matrix_result, Glm_Object *glm_obj, 
                    Link_Mode mode);

void calculate_contrasts(Glm_Object *glm_obj);

void general_least_sq(Glm_Object *glm_obj);

void initialize_info_matrix(Glm_Object *glm_obj);

double calculate_deviance(Glm_Object *glm_obj);

void fit_glm(Glm_Object *glm_obj);

void initialize_mu(Glm_Object *glm_obj);

void initialize_glm(Glm_Object *glm_obj, Initial_Data *initial_data);

void create_contrast(Glm_Object *glm_obj);

int get_design_matrix(char *design_filename, char ***infiles, 
                      Matrix **design_matrix, char *path, int num_response,
                      int max_num_row, int max_num_col);

int get_design_matrix_file(char *design_filename, char ***infiles, 
                      Matrix **design_matrix, char *path, int num_response,
                      int max_num_row, int max_num_col);

int get_design_matrix_stdin(char ***infiles, 
                      Matrix **design_matrix, char *path, int num_response,
                      int max_num_row, int max_num_col);

void create_family_names(char ***family_names, char ***link_names, 
                         char ***variance_function_names,
                         char ***corr_struct_names);

Family get_family(char *family_str, char **family_names);

Link get_link(char *link_str, char **link_names);

Variance_Function get_variance_function(char *variance_function_str,
                                        char **variance_function_names);

Corr_Struct get_corr_struct(char *corr_struct_str, 
                           char **corr_struct_names);

#endif
