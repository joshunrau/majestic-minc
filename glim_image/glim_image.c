/* ----------------------------- MNI Header -----------------------------------
@NAME       : glim_image
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (nothing)
@RETURNS    : status
@DESCRIPTION: Program to fit a generalized linear statistical model to a
              set of minc images and output various diagnostic images,
@CREATED    : June 18, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <voxel_loop.h>
#include <time_stamp.h>
#include <ParseArgv.h>
#include <volume_io.h>
#include "glim.h"
#include "smoothness.h"

#define VIO_BOOL_DEFAULT (-1)
#define INVALID_DATA (-DBL_MAX)

/* typedefs */

typedef struct {
   Glm_Object *glm_obj;
   Lambda *lambda_buffer;
   char **tmp_files;
} Program_Data;

/* Function prototypes */
int get_t_stat_raw(char *dst, char *key, int argc,
                     char **argv);

int get_f_stat_raw(char *dst, char *key, int argc,
                     char **argv);

int get_avg_raw(char *dst, char *key, int argc,
                     char **argv);

void voxel_function_glm(void *caller_data, long num_voxels, 
                    int input_num_buffers, int input_vector_length, 
                    double *input_data[],
                    int output_num_buffers, int output_vector_length, 
                    double *output_data[],
                    Loop_Info *loop_info);

void voxel_function_pool_con(void *caller_data, long num_voxels, 
                    int input_num_buffers, int input_vector_length, 
                    double *input_data[],
                    int output_num_buffers, int output_vector_length, 
                    double *output_data[],
                    Loop_Info *loop_info);

void voxel_function_pool_avg(void *caller_data, long num_voxels, 
                    int input_num_buffers, int input_vector_length, 
                    double *input_data[],
                    int output_num_buffers, int output_vector_length, 
                    double *output_data[],
                    Loop_Info *loop_info);

void start_function(void *caller_data, long num_voxels, 
                    int output_num_buffers, int output_vector_length, 
                    double *output_data[], Loop_Info *loop_info);

void end_function(void *caller_data, long num_voxels, 
                  int output_num_buffers, int output_vector_length, 
                  double *output_data[], Loop_Info *loop_info);

void create_stat_variable(Program_Data *program_data, char **input_files);

void create_glim_parameters(Glm_Object *glm_obj, char **input_files, 
                            char **param_str);

void get_fwhm_filenames(Lambda *lambda_buffer, Glm_Object *glm_obj);

int update_lambda_buffer(Lambda *lambda_buffer, Glm_Object *glm_obj, 
                         double check);

/* Argument variables */

/* Glim_image specific variables */

char *family_str = NULL;
char *link_str = NULL;
char *variance_function_str = NULL;
char *corr_struct_str = NULL;
int use_mu_switch = FALSE;
int binomial_n = 1;
int num_response = 1;
int max_num_row = 300;
int max_num_col = 100;
static Contrast_Raw_Array contrast_raw_list = {0, NULL};
double tol = 0.0001;
double eps_sing = 0.0001;
double max_iter= MAX_ITERATIONS;
double scale_par = -1;
double maxVoxelSdBasedValue = 99.0;
int est_scale = FALSE;
int chi_sq_nu = 0;
int is_fwhm_flag = FALSE;
int fit_all = FALSE;
int scylla = FALSE;
char *deg_free_file = NULL;
char *deviance_file = NULL;
char *avg_deviance_file = NULL;
char *scale_file = NULL;
char *pool_mask = NULL;
char *search_mask = NULL;
char *fwhm_gaussian_file = NULL;
char *fwhm_simple_file = NULL;
char *fwhm_avg_file = NULL;
char *error_file = NULL;
char *path_tmp = NULL;
FILE *tmpfile_list = NULL;
char *tmpfile_name = NULL;
int num_tmpfiles = 0;
static long num_fitted = 0;
static long num_skipped = 0;
static long num_singular = 0;

/* Voxel loop specific variables */

int clobber = FALSE;
int verbose = TRUE;
nc_type datatype = NC_UNSPECIFIED;
int is_signed = TRUE;
double valid_range[2] = {0,0};
int copy_all_header = FALSE;
int max_buffer_size_in_kb = 10 * 1024;
int check_dim_info = TRUE;

/* Argument table */
ArgvInfo argTable[] = {
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
       "Stats options:"},
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
       "Common options:"},
   {"-family", ARGV_STRING, (char *) NULL, (char *) &family_str,
       "Specify exponential family for GLM."},
   {"-link", ARGV_STRING, (char *) NULL, (char *) &link_str,
       "Specify link function for GLM."},
   {"-var_func", ARGV_STRING, (char *) NULL, (char *) &variance_function_str,
       "Specify variance function for GLM."},
   {"-t_stat", ARGV_GENFUNC, (char *) get_t_stat_raw, 
       (char *) &contrast_raw_list, 
       "Defines t-statistics wanted for output"},
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
       "T-statistics are defined as follows: -t_stat <out.mnc> type (one of t_stat, \nstdev, corr, beta) stdev_type (one of pool or voxel) intype (one of matrix, \ncolumn or file) string (string defining contrast depends on intype)"},
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
       "For example, the following: \n    -t_stat out.mnc corr pool column 4 \ndefines an output file that will be a partial correlation testing column 4 \nof the design matrix"},
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
       "    -t_stat out.mnc beta voxel matrix '0 1.5 0 -2.3'  \ndefines an output file that will have the of 1.5*beta2 - 2.3*beta4 \n    where the betas are the parameters estimated in the regression."},
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
       "    -t_stat out.mnc stdev pool file contrast1.txt   \ndefines an output file that will be the stdev of the t_stat defined in the \ntextfile contrast1.txt with a pooled stdev. This doesn't make sense because \nwe would not use a pooled stdev to look at the stdev as it varies over the\nimage, so glim_image spits this out as a warning."},
   {"-f_stat", ARGV_GENFUNC, (char *) get_f_stat_raw, 
       (char *) &contrast_raw_list, 
       "Defines F-statistics wanted for output"},
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
       "F_statistics are defined in the same way as t_statistics except there is no \ntype, this is because the variance of an F-stat is a matrix, as is the \nnumerator. For this reason, column is also not an intype for F_statistics \nbecause a matrix with more than two rows is needed for an F_stat."},
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
       "For example: \n    -f_stat out.mnc pool matrix '0 1 2 0 ; 1 2 3 0' \ndefines an output file that is an F-stat testing whether \n    beta2 + 2*beta3 = 0 and beta1 + beta2 + beta3 = 0 \nsimultaneously. This comes into play for anova designs when there are factors\nwith multiple levels."},
   {"-est_scale", ARGV_CONSTANT, (char *) TRUE, (char *) &est_scale,
      "Use estimated scale (for binomial and poisson families)"},
   {"-scale_par", ARGV_FLOAT, (char *) 1, (char *) &scale_par,
       "Set scale parameter, (default: estimate it)"},
   {"-max_iter", ARGV_FLOAT, (char *) 1, (char *) &max_iter, 
       "Maximum iterations for Fisher scoring."},
   {"-tolerance", ARGV_FLOAT, (char *) 1, (char *) &tol,
       "Tolerance for convergence of fitted parameters."},
   {"-eps_sing", ARGV_FLOAT, (char *) 1, (char *) &eps_sing,
       "Epsilon added to repsonse values that would cause \nsingularities."},
   {"-maxVoxelSdBasedValue", ARGV_FLOAT, (char *) 1,
       (char *) &maxVoxelSdBasedValue, "Maximum absolute value of t_statistic image."},
   {"-avg", ARGV_GENFUNC, (char *) get_avg_raw, (char *) &contrast_raw_list,
       "File name for means (missing values not skipped. Needs sd_type and out_type."},
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
       "Avg is used for comparison of means at the same time as fitting a design\n file, the output types are t_stat, beta, corr and stdev; stdev types are \nvoxel or pool."},
   {"-deviance", ARGV_STRING, (char *) NULL, (char *) &deviance_file,
       "Output an image of deviances (default: FALSE)"},
   {"-avg_deviance", ARGV_STRING, (char *) NULL, (char *) &avg_deviance_file,
       "Output image of deviances for avg."},
   {"-scale", ARGV_STRING, (char *) NULL, (char *) &scale_file,
       "Output image of scale estimate (default: FALSE)"},
   {"-errfile", ARGV_STRING, (char *) NULL, (char *) &error_file,
       "Output convergence error codes."},
   {"-use_mu", ARGV_CONSTANT, (char *) TRUE, (char *) &use_mu_switch,
       "Use previous fitted values instead of response for \ninitial values."},
   {"-pool_sd", ARGV_STRING, (char *) NULL, (char *) &pool_mask,
       "Mask file for region over which to pool s.d."},
   {"-fwhm", ARGV_CONSTANT, (char *) TRUE, (char *) &is_fwhm_flag,
       "Calculate FWHM estimate to be used for assessing significance."},
   {"-fwhm_gauss", ARGV_STRING, (char *) NULL, (char *) &fwhm_gaussian_file,
       "Output fwhm image for gaussian fit."},
   {"-fwhm_simple", ARGV_STRING, (char *) NULL, (char *) &fwhm_simple_file,
       "Output non-weight corrected FWHM estimate for non-Gaussian data."},
   {"-fwhm_avg", ARGV_STRING, (char *) NULL, (char *) &fwhm_avg_file,
       "Output fwhm image for comparison of means fit."},
   {"-search", ARGV_STRING, (char *) NULL, (char *) &search_mask,
       "Mask file for search region."},
   {"-fit_all", ARGV_CONSTANT, (char *) TRUE, (char *) &fit_all,
       "Fit GLM at all voxels with no missing values, not just search region."},
   {"-fit_missing", ARGV_STRING, (char *) TRUE, (char *) &deg_free_file,
       "Fit GLM at all voxels, even with missing values. Must provide name for\ndegrees of freedom volume."},
   {"-max_num_row", ARGV_INT, (char *) 1, (char *) &max_num_row,
       "Maximum number of rows for design_file. Should be \nincreased for large files."},
   {"-max_num_col", ARGV_INT, (char *) 1, (char *) &max_num_col,
       "Maximum number of rows for design_file. Should be \nincreased for large files."},  
   {"-scylla", ARGV_CONSTANT, (char *) TRUE, (char *) &scylla,
       "Flag to see if computer is at neuro or math dept."},
    {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
       "General options:"},
   {"-path", ARGV_STRING, (char *) NULL, (char *) &path_tmp,
       "Path for input files and design_matrix file."},
   {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
       "Overwrite existing file."},
   {"-noclobber", ARGV_CONSTANT, (char *) FALSE, (char *) &clobber,
       "Don't overwrite existing file (default)."},
   {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
       "Print out log messages (default)."},
   {"-quiet", ARGV_CONSTANT, (char *) FALSE, (char *) &verbose,
       "Do not print out log messages."},
   {"-copy_header", ARGV_CONSTANT, (char *) TRUE, (char *) &copy_all_header,
       "Copy all of the header from the first file."},
   {"-nocopy_header", ARGV_CONSTANT, (char *) FALSE, (char *) &copy_all_header,
       "Do not copy all of the header from the first file."},
   {"-filetype", ARGV_CONSTANT, (char *) NC_UNSPECIFIED, (char *) &datatype,
       "Use data type of first file (default)."},
   {"-byte", ARGV_CONSTANT, (char *) NC_BYTE, (char *) &datatype,
       "Write out byte data."},
   {"-short", ARGV_CONSTANT, (char *) NC_SHORT, (char *) &datatype,
       "Write out short integer data."},
   {"-long", ARGV_CONSTANT, (char *) NC_LONG, (char *) &datatype,
       "Write out long integer data."},
   {"-float", ARGV_CONSTANT, (char *) NC_FLOAT, (char *) &datatype,
       "Write out single-precision floating-point data."},
   {"-double", ARGV_CONSTANT, (char *) NC_DOUBLE, (char *) &datatype,
       "Write out double-precision floating-point data."},
   {"-signed", ARGV_CONSTANT, (char *) TRUE, (char *) &is_signed,
       "Write signed integer data."},
   {"-unsigned", ARGV_CONSTANT, (char *) FALSE, (char *) &is_signed,
       "Write unsigned integer data (default if type specified)."},
   {"-range", ARGV_FLOAT, (char *) 2, (char *) valid_range,
       "Valid range for output data."},
   {"-max_buffer_size_in_kb", ARGV_INT, (char *) 1, 
       (char *) &max_buffer_size_in_kb,
       "Specify the maximum size of the internal buffers (in kbytes)."},
   {"-check_dimensions", ARGV_CONSTANT, (char *) TRUE, 
       (char *) &check_dim_info,
       "Check that files have matching dimensions (default)."},
   {"-nocheck_dimensions", ARGV_CONSTANT, (char *) FALSE, 
       (char *) &check_dim_info,
       "Do not check that files have matching dimensions."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

/* Main program */

int main(int argc, char *argv[])
{
   char *prog, **input_files, **output_files, **infiles_pool;
   char *design_filename;
   char *path_in = NULL;
   char *path_end = "/";
   int num_input_files, num_output_files;
   int path_len;
   Loop_Options *loop_options;
   Matrix *design_matrix;
   Matrix *tmp, *tmp2, *singval, *buffer;
   Contrast_Mat_Array *contrast_matrices; 
   char *arg_string, *tmp_count;
   int i, j, pool_test, is_avg_fwhm;
   Contrast_Matrix *cur_matrix, *cur_matrix2;
   Avg_Info *cur_avg, *cur_avg2;
   Program_Data *program_data;
   Glm_Object *glm_obj;
   Initial_Data *initial_data;
   char **family_names, **link_names, **variance_function_names;
   char **corr_struct_names;
   Family family_input;
   Link link_input;
   Variance_Function variance_function_input;
   Corr_Struct corr_input;
   Lambda *lambda_buffer = NULL;
   VIO_Volume tmp_volume;
   volume_input_struct volume_input;
   int *sizes;
   int num_dim, itest, is_fwhm;
   int volume_sizes[VIO_MAX_DIMENSIONS];
   double volume_step[VIO_MAX_DIMENSIONS];
   Fwhm_Info *fwhm_general, *fwhm_gaussian;
   
   tmpfile_name = tempnam(NULL, "glim");

   tmpfile_list = fopen(tmpfile_name, "w");

   /* Save time stamp and args */
   arg_string = time_stamp(argc, argv); 

   /* Check arguments */
   prog = argv[0];
   if (ParseArgv(&argc, argv, argTable, 0) || (argc != 2)) {
      (void) fprintf(stderr, "\nUsage: %s [options] <des.file> \n", prog);
      (void) fprintf(stderr, 
        "       %s -help\n", prog);
      fprintf(stderr,"\nYour command line was (after parsing): \n\n      ");
      for(i=0; i<argc; i++) {
         fprintf(stderr,"%s ",argv[i]);
      }
      fprintf(stderr,"\n\n");
      fprintf(stderr,"Use '-' as des.file to read from stdin. Last line should end with 'END'.\n\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
   }
   
   create_family_names(&family_names, &link_names, &variance_function_names, 
                       &corr_struct_names); 

   family_input = get_family(family_str, family_names);
   link_input = get_link(link_str, link_names);
   variance_function_input = get_variance_function(variance_function_str,
                                                   variance_function_names);
   corr_input = get_corr_struct(corr_struct_str, corr_struct_names); 

   /* Check to see that path ends in '/' */

   if(path_tmp != NULL) {
      path_len = strlen(path_tmp);

      if(path_tmp[path_len-1] != path_end[0]) {
         path_in = GI_MALLOC(sizeof(char) * (path_len + 2));
         strcpy(path_in, path_tmp);
         strcat(path_in, path_end);
      }
   }

   /* Get design matrix from file design_filename */
   design_filename = argv[1];
   if( strcmp(design_filename,"-") == 0){
      design_filename = NULL;
      }

   if(get_design_matrix(design_filename, &input_files, &design_matrix,
                             path_in, num_response, max_num_row,
                             max_num_col) != TRUE) {
      fprintf(stderr,"Error: error in get_design_matrix_file. \n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
   }
   
   /* a bit of output */
   if(verbose){
      fprintf(stdout, "+++ design matrix from %s +++\n", design_filename);   
      for(i=0; i<design_matrix->num_rows; i++) {
         fprintf(stdout, " | [%s] - ", input_files[i]);   
         for(j=0; j<design_matrix->num_columns; j++) {
         fprintf(stdout, "%5g ", design_matrix->values[i][j]);
         }
         fprintf(stdout, "\n");
      }
   }
   
   contrast_matrices = GI_MALLOC(sizeof(Contrast_Mat_Array));
   contrast_matrices->num_columns = design_matrix->num_columns;
   contrast_matrices->num_contrasts = (contrast_raw_list.num_output - 
                                       contrast_raw_list.num_avg);
   contrast_matrices->num_avg = contrast_raw_list.num_avg;
   contrast_matrices->contrast_matrix_array =\
      GI_MALLOC(contrast_matrices->num_contrasts * sizeof(Matrix *));
   contrast_matrices->avg_array = GI_MALLOC(sizeof(Avg_Array));

   get_contrast_from_raw(&contrast_raw_list, &contrast_matrices);
   
   /* Check ParseArgv options */

   if ((fwhm_gaussian_file != NULL) || (fwhm_avg_file != NULL) || 
       (fwhm_simple_file != NULL) || (is_fwhm_flag == TRUE))
      is_fwhm = TRUE;

   if ((fwhm_gaussian_file != NULL) && (fwhm_simple_file != NULL)) {
      fprintf(stderr,"\nError: \"-fwhm_gauss\" and \"-fwhm_simple\" produce same output. \"-fwhm_gauss\" is only if data is Gaussian.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
   }

   if ((fwhm_gaussian_file != NULL) && (family_input != GAUSSIAN)) {
      fprintf(stderr,"\nError: \"-fwhm_gauss\" option called for when family is not gaussian, use \"-fwhm_simple\".\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
   }

   if ((contrast_matrices->avg_array == NULL) && (fwhm_avg_file != NULL)){
      fprintf(stderr,"\nError: \"-fwhm_avg\" option called for when no average test statistics being output.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
   }
   else if (contrast_matrices->avg_array != NULL) { 
      if ((contrast_matrices->avg_array->num_test <= 0)
          && (fwhm_avg_file != NULL)){
         fprintf(stderr,"\nError: \"-fwhm_avg\" option called for when no average test statistics being output.\n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }
   }

   if (family_input == GAUSSIAN) {
      if (contrast_matrices->avg_array == NULL) {
         if ((contrast_matrices->num_test <= 0) && (fwhm_gaussian_file != NULL)) {
            fprintf(stderr,"\nError: \"-fwhm_gauss\" option called for when no test statistics being output\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
      }
      else if (contrast_matrices->avg_array->num_test <= 0) {
         if ((contrast_matrices->num_test <= 0) && (fwhm_gaussian_file != NULL)) {
            fprintf(stderr,"\nError: \"-fwhm_gauss\" option called for when no test statistics being output\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
      }
   }
   else {
      if (contrast_matrices->avg_array == NULL) {
         if ((contrast_matrices->num_test <= 0) && (is_fwhm_flag == TRUE)) {
            fprintf(stderr,"\nError: \"-fwhm\" option called for when no test statistics being output\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
      }
      else if (contrast_matrices->avg_array->num_test <= 0) {
         if ((contrast_matrices->num_test <= 0) && (is_fwhm_flag == TRUE)) {
            fprintf(stderr,"\nError: \"-fwhm\" option called for when no test statistics being output\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
      }
   }
  
   /* Initialize lambda_buffer first because of its size */

   if (is_fwhm == TRUE) {

      if(start_volume_input(input_files[0], 3, NULL, NC_UNSPECIFIED, FALSE,
                            0.0, 0.0, TRUE, &tmp_volume,
                            (minc_input_options *) NULL, &volume_input) != VIO_OK){
         fprintf(stderr,"\nError opening %s to check for size[] and step[] for smoothness calculation.\n", input_files[0]);
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }

      num_dim = get_volume_n_dimensions(tmp_volume);
      get_volume_sizes(tmp_volume, volume_sizes);
      get_volume_separations(tmp_volume, volume_step);

      cancel_volume_input(tmp_volume, &volume_input);

      sizes = GI_MALLOC(sizeof(*sizes) * num_dim);

      for(i=0; i<num_dim; i++) {
         sizes[i] = volume_sizes[i];
      }

      if (contrast_matrices->avg_array != NULL) {
         if(contrast_matrices->avg_array->num_test > 0) 
            is_avg_fwhm = TRUE;
         else
            is_avg_fwhm = FALSE;
      }
      else
         is_avg_fwhm = FALSE;

      lambda_buffer = create_lambda_buffer(design_matrix->num_rows,
                                           contrast_matrices->num_test,
                                           num_dim, sizes, 
                                           (family_input == GAUSSIAN),
                                           is_avg_fwhm, fwhm_simple_file);

      if(lambda_buffer == NULL) {
         fprintf(stderr,"\nMemory error in create_lambda_buffer.\n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }

      lambda_buffer->constant = sqrt(4.0 * log(2.0));

      for(i=0; i<num_dim; i++) {
         lambda_buffer->constant *= pow(volume_step[i], 1.0/num_dim); 
         lambda_buffer->step[i] = volume_step[i];
      }

   }

   initial_data = GI_MALLOC(sizeof(*initial_data));

   initial_data->design_matrix = design_matrix;
   initial_data->contrast = contrast_matrices;
   initial_data->family = family_input;
   initial_data->link = link_input;
   initial_data->variance_function = variance_function_input;
   initial_data->corr_struct = corr_input;

   /* Get control parameters */

   initial_data->control = GI_MALLOC(sizeof(Glm_Control));

   if(pool_mask != NULL)
      initial_data->control->pool_id = TRUE;
   else
      initial_data->control->pool_id = FALSE;

   initial_data->control->status = CONT;

   initial_data->control->scylla = scylla;

   if(search_mask != NULL)
      initial_data->control->search_id = TRUE;
   else
      initial_data->control->search_id = FALSE;

   initial_data->control->fit_all = fit_all;
   if (deg_free_file != NULL) {
      initial_data->control->fit_all = TRUE;
      initial_data->control->deg_free_id = TRUE;
   }
   else
      initial_data->control->deg_free_id = -1;

   if(deviance_file!= NULL)
      initial_data->control->dev_id = TRUE;
   else 
      initial_data->control->dev_id = -1;

   if(avg_deviance_file!= NULL)
      initial_data->control->avg_dev_id = TRUE;
   else 
      initial_data->control->avg_dev_id = -1;

   if(scale_file!= NULL)
      initial_data->control->scale_id = TRUE;
   else 
      initial_data->control->scale_id = -1;

   if(error_file != NULL) {
      initial_data->control->err_id = TRUE;
   }
   else {
      initial_data->control->err_id = -1;
   }

   initial_data->control->maxVoxelSdBasedValue = maxVoxelSdBasedValue;

   initial_data->control->is_fwhm = is_fwhm;
   initial_data->control->use_mu_switch = use_mu_switch;
   initial_data->control->max_iter = max_iter;
   initial_data->control->tol = tol;
   if(scale_par >= 0.0) {
      initial_data->scale_est = scale_par;
      initial_data->control->is_scale = TRUE;
   }
   else 
      initial_data->control->is_scale = FALSE;
   initial_data->control->est_scale = est_scale;

   if(chi_sq_nu > 0) {
      initial_data->family = CHI_SQ;
      initial_data->scale_est = 2.0 / chi_sq_nu;
   }

   initial_data->control->eps_sing = eps_sing;
   initial_data->control->binomial_n = binomial_n;

   if ((initial_data->contrast->num_contrasts == 0) && 
       (initial_data->control->dev_id == -1) &&
       (initial_data->control->err_id == -1) &&
       (initial_data->control->scale_id == -1) &&
       (lambda_buffer == NULL))
      initial_data->control->do_fit = FALSE;
   else
      initial_data->control->do_fit = TRUE;

   glm_obj = GI_MALLOC(sizeof(Glm_Object));
   initialize_glm(glm_obj, initial_data);

   program_data = GI_MALLOC(sizeof(*program_data));
   
   if (lambda_buffer != NULL) {

      get_fwhm_filenames(lambda_buffer, glm_obj); 

      if (glm_obj->control->do_fit == TRUE) {
         if (glm_obj->family == GAUSSIAN)
            if (lambda_buffer->fwhm_gaussian != NULL)
               lambda_buffer->fwhm_gaussian->outfile = fwhm_gaussian_file;

         if(glm_obj->family != GAUSSIAN) {

            if(glm_obj->deg_free > 1) {
               for(itest=0; itest<glm_obj->contrast->num_test; itest++) {
                  fwhm_general = lambda_buffer->fwhm_general[itest];
                  fwhm_general->constant = 1.0;
               }
               if (lambda_buffer->fwhm_simple != NULL) {
                  lambda_buffer->fwhm_simple->constant = 1.0;
               }
            }
            else {
               fprintf(stderr,"\nError: to calculate FWHM degrees of freedom must be at least 2.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
         }
         else if (lambda_buffer->fwhm_gaussian != NULL){
            if(glm_obj->deg_free > 1) {
               lambda_buffer->fwhm_gaussian->constant = 1.0;
            }
            else {
               fprintf(stderr,"\nError: to calculate FWHM degrees of freedom must be at least 2.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
         }
      }

      if(lambda_buffer->is_avg == TRUE) {

         if(glm_obj->contrast->avg_array != NULL)
            if (lambda_buffer->fwhm_avg != NULL)
               lambda_buffer->fwhm_avg->outfile = fwhm_avg_file;

         if(glm_obj->response->num_rows > 2) {
            lambda_buffer->fwhm_avg->constant = 1.0;
         }
         else {
            fprintf(stderr,"\nError: to calculate FWHM degrees of freedom must be at least 2.\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
      }
   }

   glm_obj->family_names = family_names;
   glm_obj->link_names = link_names;
   glm_obj->variance_function_names = variance_function_names;
   glm_obj->corr_struct_names = corr_struct_names;
   glm_obj->deviance_file = deviance_file;
   if (glm_obj->contrast->avg_array != NULL)
      glm_obj->contrast->avg_array->deviance_file = avg_deviance_file;
   glm_obj->error_file = error_file;

   /*Calculate x_transpsose_x_inv to be used for Gaussian */

   if ((glm_obj->corr_struct == DIAGONAL_CORR) && 
       (glm_obj->control->do_fit == TRUE)) {
      switch(glm_obj->family) {
      case GAUSSIAN:
         tmp = create_matrix(design_matrix->num_columns, 
                                design_matrix->num_columns);

         multiply_matrix_sym(tmp, design_matrix, design_matrix, MATRIX_LEFT);

         if (invert_matrix_sym(glm_obj->info_inv, tmp) != TRUE) {
            fprintf(stderr,"\ndesign matrix:");
            print_matrix(stderr, design_matrix);
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
         
         for(i=0; i<glm_obj->contrast->num_contrasts; i++) {
            cur_matrix = glm_obj->contrast->contrast_matrix_array[i];
            if (cur_matrix->contrast_matrix != NULL) {
               multiply_matrix(glm_obj->tmp_con[i][1], glm_obj->info_inv,
                               cur_matrix->contrast_matrix, MATRIX_RIGHT);
               multiply_matrix_sym(glm_obj->tmp_con[i][2],
                                   cur_matrix->contrast_matrix, 
                                   glm_obj->tmp_con[i][1], MATRIX_NEITHER);
               invert_matrix_sym(glm_obj->tmp_con[i][3],
                                 glm_obj->tmp_con[i][2]);
            }
         }

         multiply_matrix(glm_obj->lm_matrix, glm_obj->info_inv,
                         glm_obj->design_matrix, MATRIX_RIGHT);
         delete_matrix(tmp);
         break;
         
      default:
         break;
      }
   }

   /* Checks to see if there is a need to create a temporary count file */

   glm_obj->control->count_id = -1;

   if(glm_obj->contrast->num_avg > 0)
      glm_obj->contrast->avg_array->num_pool = 0;

   for(i=0; i<glm_obj->contrast->num_avg; i++) {
      cur_avg = glm_obj->contrast->avg_array->avg_info[i];
      if (cur_avg->stdev_mode == POOLED_SD) {
         glm_obj->control->count_id = TRUE;
         glm_obj->contrast->avg_array->num_pool++;
      }
   }

   if (glm_obj->control->deg_free_id > 0)
      glm_obj->control->count_id = -1;

   if(glm_obj->control->count_id == TRUE) {
      tmp_count = tempnam(NULL, "glim");
      fprintf(tmpfile_list, "\n%s", tmp_count);
      num_tmpfiles++;
   }
   else
      tmp_count = NULL;

   /* Changing the value of num_test to make sure it's OK for gaussian */

   if ((glm_obj->family == GAUSSIAN) && (lambda_buffer != NULL)) {
      lambda_buffer->num_test = 1;
   }

   num_input_files = (design_matrix->num_rows + glm_obj->control->pool_id +
                      glm_obj->control->search_id);

   if(lambda_buffer != NULL) {
      num_output_files = (glm_obj->contrast->num_contrasts + 
                          (glm_obj->control->avg_dev_id > 0) +
                          (glm_obj->control->dev_id > 0) +
                          (glm_obj->control->scale_id > 0) +
                          (glm_obj->control->err_id > 0) +
                          (glm_obj->control->count_id > 0) +
                          (lambda_buffer->num_test * 
                           (!(lambda_buffer->is_gaussian))) + 
                          glm_obj->contrast->num_avg + 
                          ((lambda_buffer->fwhm_gaussian != NULL) &&
                           (fwhm_gaussian_file != NULL)) +
                          ((lambda_buffer->fwhm_avg != NULL) && 
                           (fwhm_avg_file != NULL)) +
                          ((lambda_buffer->fwhm_simple != NULL) &&
                           (fwhm_simple_file != NULL)) +
                          (glm_obj->control->deg_free_id > 0));
   }
   else {
      num_output_files = (glm_obj->contrast->num_contrasts + 
                          (glm_obj->control->avg_dev_id > 0) +
                          (glm_obj->control->dev_id > 0) +
                          (glm_obj->control->scale_id > 0) +
                          (glm_obj->control->err_id > 0) +
                          (glm_obj->control->count_id > 0)+
                          (glm_obj->control->deg_free_id > 0) +
                          glm_obj->contrast->num_avg);
      
   }

   /* Set pool_id to its index in input file list */

   i = design_matrix->num_rows;

   if(glm_obj->control->pool_id == TRUE) {
      glm_obj->control->pool_id = i;
      i++;
   }
   else
      glm_obj->control->pool_id = -1;

   if(glm_obj->control->search_id == TRUE) {
      glm_obj->control->search_id = i;
      i++;
   }
   else
      glm_obj->control->search_id = -1;

   
   /* Set up file lists  checking to see if to be removed */

   if (glm_obj->control->pool_id > 0)
      input_files[glm_obj->control->pool_id] = pool_mask;

   if (glm_obj->control->search_id > 0) {
      input_files[glm_obj->control->search_id] = search_mask;
      glm_obj->search_mask = search_mask;
   }
   else
      glm_obj->search_mask = NULL;
      
   if (num_output_files <= 0) {
      fprintf(stderr,"\nError: no output files specified.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
   }
         
   output_files = GI_MALLOC(num_output_files * sizeof(*output_files));

   for(i=0; i<glm_obj->contrast->num_contrasts; i++) {

      cur_matrix = glm_obj->contrast->contrast_matrix_array[i];
      if(cur_matrix->stdev_mode == VOXEL_SD)
         output_files[cur_matrix->out_id] = cur_matrix->outfile;
      else if(cur_matrix->stdev_mode == POOLED_SD) {
         if(cur_matrix->tmpfile != NULL) {
            output_files[cur_matrix->out_id] = cur_matrix->tmpfile;
         }
         else {
            fprintf(stderr,"\nError: no temporary file name for outfile %s.\n",
                    cur_matrix->outfile);
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
      }
   }
   
   for(i=0; i<glm_obj->contrast->num_avg; i++) {

      cur_avg = glm_obj->contrast->avg_array->avg_info[i];
      if(cur_avg->stdev_mode == VOXEL_SD)
         output_files[cur_avg->out_id] = cur_avg->outfile;
      else if(cur_avg->stdev_mode == POOLED_SD) {
         if(cur_avg->tmpfile != NULL) {
            output_files[cur_avg->out_id] = cur_avg->tmpfile;
         }
         else {
            fprintf(stderr,"\nError: no temporary file name for outfile %s.\n",
                    cur_avg->outfile);
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
      }
   }

   /* Check to see that all output files have different names and whether */
   /* pooling mask has been given unnecessarily.                          */

   pool_test = TRUE;

   for(i=0; i<glm_obj->contrast->num_contrasts; i++) {
      cur_matrix = glm_obj->contrast->contrast_matrix_array[i];
      if(cur_matrix->stdev_mode != VOXEL_SD)
         pool_test = FALSE;
      for(j=i+1; j<glm_obj->contrast->num_contrasts; j++) {
         cur_matrix2 = glm_obj->contrast->contrast_matrix_array[j];
         if(strcmp(cur_matrix->outfile, cur_matrix2->outfile) == 0) {
            fprintf(stderr,"\nError: output files %d and %d have the same filename.\n", cur_matrix->out_id+1, cur_matrix2->out_id+1);
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
      }
   }

   for(i=0; i<glm_obj->contrast->num_avg; i++) {
      cur_avg = glm_obj->contrast->avg_array->avg_info[i];
      if(cur_avg->stdev_mode != VOXEL_SD) {
         pool_test = FALSE;
      }
      if(glm_obj->contrast->num_avg > 1) {
         for(j=i+1; j<glm_obj->contrast->num_avg; j++) {
            cur_avg2 = glm_obj->contrast->avg_array->avg_info[j];
            if(strcmp(cur_avg->outfile, cur_avg2->outfile) == 0) {
               fprintf(stderr,"\nError: output files %d and %d have the same filename.\n", cur_avg->out_id+1, cur_avg2->out_id+1);
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
         }
      }
   }

   if((pool_test == TRUE) && (glm_obj->control->pool_id > 0)) {
      fprintf(stderr,"\nError: pooling mask file given when no output files are pooled stdev type.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
   }

   i = glm_obj->contrast->num_contrasts + glm_obj->contrast->num_avg;

   if (glm_obj->control->dev_id == TRUE) {
      output_files[i] = deviance_file;
      glm_obj->control->dev_id = i;
      i++;
   }

   if(glm_obj->control->deg_free_id == TRUE) {
      output_files[i] = deg_free_file;
      glm_obj->control->deg_free_id = i;
      i++;
   }
   else
      glm_obj->control->deg_free_id = -1;

   if (glm_obj->control->avg_dev_id == TRUE) {
      output_files[i] = avg_deviance_file;
      glm_obj->control->avg_dev_id = i;
      i++;
   }

   if (glm_obj->control->scale_id == TRUE) {
      output_files[i] = scale_file;
      glm_obj->control->scale_id = i;
      glm_obj->scale_file = scale_file;
      i++;
   }

   if (glm_obj->control->err_id == TRUE) {
      output_files[i] = error_file;
      glm_obj->control->err_id = i;
      i++;
   }

   if (glm_obj->control->count_id == TRUE) {
      output_files[i] = tmp_count;
      glm_obj->control->count_id = i;
      i++;
   }

   if ((glm_obj->control->is_fwhm == TRUE) && (glm_obj->family == GAUSSIAN)){
      if(lambda_buffer->fwhm_gaussian != NULL) {
         if(lambda_buffer->fwhm_gaussian->outfile != NULL) {
            output_files[i] = lambda_buffer->fwhm_gaussian->outfile;
            lambda_buffer->fwhm_gaussian->out_id = i;
            i++;
         }
      }
   }

   if (lambda_buffer != NULL) {
      if (glm_obj->contrast->avg_array != NULL) {
         if (lambda_buffer->fwhm_avg->outfile != NULL) {
            output_files[i] = lambda_buffer->fwhm_avg->outfile;
            lambda_buffer->fwhm_avg->out_id = i;
            i++;
         }
      }

      if (lambda_buffer->fwhm_simple != NULL) {
         if (lambda_buffer->fwhm_simple->outfile != NULL) {
            output_files[i] = lambda_buffer->fwhm_simple->outfile;
            lambda_buffer->fwhm_simple->out_id = i;
            i++;
         }
      }
   }
   

   if ((glm_obj->control->is_fwhm == TRUE) && (glm_obj->family != GAUSSIAN)){
      itest = 0;
      for(j=i; j<i+lambda_buffer->num_test; j++) {
         lambda_buffer->fwhm_general[itest]->out_id = j;
         output_files[j] = lambda_buffer->fwhm_general[itest]->outfile;
         itest++;
      }
   }

   /* Check to see if pooling wanted but no pool file */
   if (glm_obj->control->pool_id < 0) {
      for(j=0; j<glm_obj->contrast->num_contrasts; j++) {
         cur_matrix = glm_obj->contrast->contrast_matrix_array[j];
         if(cur_matrix->stdev_mode == POOLED_SD) {
            fprintf(stderr,"\nError: mask file needed in order to pool_sd.\n");
            fprintf(stderr,"\nEnter mask file as -pool_sd <mask.mnc> on command line.\n\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
      }
   }

   program_data = GI_MALLOC(sizeof(*program_data));
   
   program_data->glm_obj = glm_obj;
   program_data->lambda_buffer = lambda_buffer;

   /* Set up loop options */
   loop_options = create_loop_options();
   set_loop_verbose(loop_options, verbose);
   set_loop_clobber(loop_options, clobber);
   set_loop_datatype(loop_options, datatype, is_signed, 
                     valid_range[0], valid_range[1]);
   set_loop_copy_all_header(loop_options, copy_all_header);
   set_loop_buffer_size(loop_options, (long) 1024 * max_buffer_size_in_kb);
   set_loop_check_dim_info(loop_options, check_dim_info);

   /* Do voxel loop */

   voxel_loop(num_input_files, input_files, num_output_files, output_files, 
              arg_string, loop_options,
              voxel_function_glm, (void *) program_data);

   /* Free loop options for voxel_function_pool*/
   free_loop_options(loop_options);

   if (glm_obj->control->pool_id > 0) {

      if(glm_obj->num_pool > 0) 
         glm_obj->pooled_dev = glm_obj->pooled_dev / glm_obj->num_pool;
      if(glm_obj->avg_num_pool > 0)
         glm_obj->avg_pooled_dev = (glm_obj->avg_pooled_dev / \
                                    glm_obj->avg_num_pool);

      loop_options = create_loop_options();
      set_loop_verbose(loop_options, verbose);
      set_loop_clobber(loop_options, clobber);
      set_loop_datatype(loop_options, datatype, is_signed, 
                        valid_range[0], valid_range[1]);
      set_loop_copy_all_header(loop_options, copy_all_header);
      set_loop_buffer_size(loop_options, (long) 1024 * max_buffer_size_in_kb);
      set_loop_check_dim_info(loop_options, check_dim_info);

      infiles_pool = GI_MALLOC(sizeof(char*) *\
                            (1 + ((glm_obj->control->count_id > 0) ||
                                  (glm_obj->control->deg_free_id > 0))));

      for(i=0; i<glm_obj->contrast->num_contrasts; i++) {

         cur_matrix = glm_obj->contrast->contrast_matrix_array[i];

         cur_matrix->pooled_dev = glm_obj->pooled_dev;
         cur_matrix->deg_free = glm_obj->deg_free * 1.0;
         cur_matrix->maxVoxelSdBasedValue = glm_obj->control->maxVoxelSdBasedValue;

         cur_matrix->pooled_sd = sqrt(cur_matrix->pooled_dev);

         if(cur_matrix->stdev_mode == POOLED_SD) {

            infiles_pool[0] = cur_matrix->tmpfile;
            if (glm_obj->control->deg_free_id > 0)
               infiles_pool[1] = deg_free_file;

            voxel_loop(1 + (glm_obj->control->deg_free_id > 0), infiles_pool,
                       1, &cur_matrix->outfile, arg_string, loop_options,
                       voxel_function_pool_con, (void *) cur_matrix);

            if(remove(cur_matrix->tmpfile) != 0) {
               fprintf(stderr, "\nError removing file %s. Must be removed by user.\n", cur_matrix->tmpfile);  
            }
            else
               fprintf(stderr,"\nFile %s succesfully deleted.\n\n", 
                       cur_matrix->tmpfile); 
         }
      }

      for(i=0; i<glm_obj->contrast->num_avg; i++) {
         
         cur_avg = glm_obj->contrast->avg_array->avg_info[i];

         cur_avg->pooled_dev = glm_obj->avg_pooled_dev;
         cur_avg->deg_free = glm_obj->response->num_rows - 1.0;

         cur_avg->pooled_sd = sqrt(cur_avg->pooled_dev);
         cur_avg->maxVoxelSdBasedValue =glm_obj->control->maxVoxelSdBasedValue;

         if(cur_avg->stdev_mode == POOLED_SD) {

            infiles_pool[0] = cur_avg->tmpfile;
            if (glm_obj->control->count_id > 0) {
               infiles_pool[1] = tmp_count;
               cur_avg->df_offset = 0.0;
            }
            else { 
               infiles_pool[1] = deg_free_file;
               cur_avg->df_offset = cur_avg->deg_free - glm_obj->deg_free;
            }

            voxel_loop(2, infiles_pool, 1,
                       &cur_avg->outfile, arg_string, loop_options,
                       voxel_function_pool_avg, (void *) cur_avg);

            if(remove(cur_avg->tmpfile) != 0) {
               fprintf(stderr, "\nError removing file %s. Must be removed by user.\n", cur_matrix->tmpfile);  
            }
            else
               fprintf(stderr,"\nFile %s succesfully deleted.\n\n", 
                       cur_avg->tmpfile); 
         }
      }

      free_loop_options(loop_options);

      if (tmp_count != NULL)
         remove(tmp_count); 
   }

   if (glm_obj->control->is_fwhm == TRUE) {

      if ((glm_obj->family == GAUSSIAN) && (glm_obj->contrast->num_test > 0)){
         fwhm_gaussian = lambda_buffer->fwhm_gaussian;
      
         for(i=0; i<lambda_buffer->num_dim; i++) {
            for(j=0; j<=i; j++) {
               fwhm_gaussian->lambda_pool->values[i][j] =\
                  (fwhm_gaussian->lambda_pool->values[i][j] /
                   fwhm_gaussian->num_calc);
            }
         }

         for(i=0; i<lambda_buffer->num_dim; i++) {
            for(j=lambda_buffer->num_dim-1; j>i; j--) {
               fwhm_gaussian->lambda_pool->values[i][j] = fwhm_gaussian->lambda_pool->values[j][i];
            }
         }

         fwhm_gaussian->smoothness_pool = sqrt(determinant(fwhm_gaussian->lambda_pool));

         fwhm_gaussian->smoothness_pool *= fwhm_gaussian->constant;

         if(fwhm_gaussian->smoothness_pool > 0.0) {
            fwhm_gaussian->fwhm_pool = pow(fwhm_gaussian->smoothness_pool,
                                           1.0/lambda_buffer->num_dim);
         }
         else
            fwhm_gaussian->fwhm_pool = INVALID_DATA;
      
         if(fwhm_gaussian->fwhm_pool > 0.0) {
            fwhm_gaussian->fwhm_pool = (lambda_buffer->constant /
                                        fwhm_gaussian->fwhm_pool);
         }
         else
            fwhm_gaussian->fwhm_pool = INVALID_DATA;
         
         fwhm_gaussian->int_lambda = (fwhm_gaussian->int_lambda *
                                     volume_step[0] * volume_step[1] *
                                     volume_step[2]);
         fwhm_gaussian->search_vol =  (fwhm_gaussian->num_calc *
                                      volume_step[0] * volume_step[1] *
                                      volume_step[2]);
         fwhm_gaussian->smoothness_int = (fwhm_gaussian->int_lambda /
                                         fwhm_gaussian->search_vol);
            
         if(fwhm_gaussian->smoothness_int > 0.0) {
            fwhm_gaussian->fwhm_int = pow(fwhm_gaussian->smoothness_int,
                                          1.0/lambda_buffer->num_dim);
         }
         else
            fwhm_gaussian->fwhm_int = INVALID_DATA;
         
         if(fwhm_gaussian->fwhm_int > 0.0) {
            fwhm_gaussian->fwhm_int = (lambda_buffer->constant /
                                       fwhm_gaussian->fwhm_int);
         }

         fprintf(stderr,"fwhm %f %f %f", fwhm_gaussian->fwhm_pool, 
                 fwhm_gaussian->num_calc, fwhm_gaussian->int_lambda);
       /*  print_matrix(stderr, fwhm_gaussian->lambda_pool); */
      }
   }
       

   if ((glm_obj->control->is_fwhm == TRUE) && 
      (glm_obj->contrast->avg_array != NULL)) {
      if(glm_obj->contrast->avg_array->num_test > 0) {
      fwhm_gaussian = lambda_buffer->fwhm_avg;

      /* print_matrix(stderr, fwhm_gaussian->lambda_pool); */

      for(i=0; i<lambda_buffer->num_dim; i++) {
         for(j=0; j<=i; j++) {
            fwhm_gaussian->lambda_pool->values[i][j] =\
                (fwhm_gaussian->lambda_pool->values[i][j] /
                   fwhm_gaussian->num_calc);
         }
      }

      /* print_matrix(stderr, fwhm_gaussian->lambda_pool); */

      for(i=0; i<lambda_buffer->num_dim; i++) {
         for(j=lambda_buffer->num_dim-1; j>i; j--) {
            fwhm_gaussian->lambda_pool->values[i][j] = fwhm_gaussian->lambda_pool->values[j][i];
         }
      }

      fwhm_gaussian->smoothness_pool = sqrt(determinant(fwhm_gaussian->lambda_pool));

      if(fwhm_gaussian->smoothness_pool > 0.0) {
         fwhm_gaussian->fwhm_pool = pow(fwhm_gaussian->smoothness_pool,
                                        1.0/lambda_buffer->num_dim);
      }
      else
         fwhm_gaussian->fwhm_pool = INVALID_DATA;
      
      if(fwhm_gaussian->fwhm_pool > 0.0) {
         fwhm_gaussian->fwhm_pool = (lambda_buffer->constant /
                                     fwhm_gaussian->fwhm_pool);
      }
      else
         fwhm_gaussian->fwhm_pool = INVALID_DATA;

      fwhm_gaussian->int_lambda = (fwhm_gaussian->int_lambda *
                                  volume_step[0] * volume_step[1] *
                                  volume_step[2]);
      fwhm_gaussian->search_vol =  (fwhm_gaussian->num_calc *
                                   volume_step[0] * volume_step[1] *
                                   volume_step[2]);
      fwhm_gaussian->smoothness_int = (fwhm_gaussian->int_lambda /
                                      fwhm_gaussian->search_vol);
            
      if(fwhm_gaussian->smoothness_int > 0.0) {
         fwhm_gaussian->fwhm_int = pow(fwhm_gaussian->smoothness_int,
                                       1.0/lambda_buffer->num_dim);
      }
      else
         fwhm_gaussian->fwhm_int = INVALID_DATA;
      
      if(fwhm_gaussian->fwhm_int > 0.0) {
         fwhm_gaussian->fwhm_int = (lambda_buffer->constant /
                                    fwhm_gaussian->fwhm_int);
      }

      fprintf(stderr,"fwhm %f %f %f", fwhm_gaussian->fwhm_pool, 
              fwhm_gaussian->num_calc, fwhm_gaussian->int_lambda);
       /* print_matrix(stderr, fwhm_gaussian->lambda_pool); */
      }
   }

   if(glm_obj->control->is_fwhm == TRUE) {
      if ((glm_obj->family != GAUSSIAN) && (glm_obj->contrast->num_test > 0)){
         for(i=0; i<lambda_buffer->num_test; i++) {
            fwhm_general = lambda_buffer->fwhm_general[itest];
            fwhm_general->int_lambda = (fwhm_general->int_lambda *
                                        volume_step[0] * volume_step[1] *
                                        volume_step[2]);
            fwhm_general->search_vol =  (fwhm_general->num_calc *
                                        volume_step[0] * volume_step[1] *
                                        volume_step[2]);
            fwhm_general->smoothness_int = (fwhm_general->int_lambda /
                                            fwhm_general->search_vol);
            
            if(fwhm_general->smoothness_int > 0.0) {
               fwhm_general->fwhm_int = pow(fwhm_general->smoothness_int,
                                             1.0/lambda_buffer->num_dim);
            }
            else
               fwhm_general->fwhm_int = INVALID_DATA;
      
            if(fwhm_general->fwhm_int > 0.0) {
               fwhm_general->fwhm_int = (lambda_buffer->constant /
                                          fwhm_general->fwhm_int);
      }
      else
         fwhm_general->fwhm_int = INVALID_DATA;


         }
      }
   }

   delete_glm_matrices(program_data->glm_obj, FALSE);

   if (lambda_buffer != NULL)
      delete_lambda_buffer(program_data->lambda_buffer);

   create_stat_variable(program_data, input_files);

/*   delete_matrix(glm_obj->design_matrix); */

   delete_tmpfiles(&tmpfile_list); 
   return EXIT_SUCCESS;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : voxel_function_glm
@INPUT      : Standard for voxel loop
@OUTPUT     : Standard for voxel loop
@RETURNS    : (nothing)
@DESCRIPTION: Routine doing working on each input buffer
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 25, 1997 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void voxel_function_glm(void *caller_data, long num_voxels, 
                    int input_num_buffers, int input_vector_length, 
                    double *input_data[],
                    int output_num_buffers, int output_vector_length, 
                    double *output_data[],
                    Loop_Info *loop_info)
/* ARGSUSED */
{
   long ivox, num_values;
   long count;
   int i, l;
   int same_test;
   int is_same = FALSE;
   int in_search;
   int itest;
   double check;
   double tmp_sum;
   double tmp_sum_sq;
   double variance, norm_avg_resid;
   double mean;
   double tmp;
   Contrast_Matrix *cur_matrix;
   Avg_Info *cur_avg;
   Glm_Object *glm_obj;
   Program_Data *program_data;
   Lambda *lambda_buffer;
   Fwhm_Info *fwhm_general;

   /* Get pointer to glm_obj */

   num_values = num_voxels * input_vector_length;

   program_data = (Program_Data *) caller_data;
   glm_obj = program_data->glm_obj;
   lambda_buffer = program_data->lambda_buffer;

   for(ivox=0 ; ivox <  num_values ; ivox++) {

      if(glm_obj->control->search_id <= 0)
         in_search = TRUE;
      else if (glm_obj->control->fit_all == TRUE)
         in_search = TRUE;
      else if (input_data[glm_obj->control->search_id][ivox] == 1)
         in_search = TRUE;
      else if (input_data[glm_obj->control->search_id][ivox] != 0) {
         fprintf(stderr,"\nError: values in search_mask must be 0 or 1.\n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }
      else
         in_search = FALSE;

      tmp_sum = 0.0;
      tmp_sum_sq = 0.0;
      count = 0;
      glm_obj->num_missing = 0;

      if (in_search == TRUE) {

         if ((glm_obj->control->do_fit == TRUE) && 
             (glm_obj->response->num_rows > 0)) {
           is_same = TRUE;
         }

         for(i=0; i < glm_obj->response->num_rows; i++) {
         
            glm_obj->missing->values[i][0] = 1.0;
            if (input_data[i][ivox] == INVALID_DATA) { 
               glm_obj->num_missing++;
               glm_obj->missing->values[i][0] = 0.0;
            }
            else {
               tmp_sum += input_data[i][ivox];
               tmp_sum_sq += input_data[i][ivox] * input_data[i][ivox];
            }

            glm_obj->response->values[i][0] = input_data[0][ivox];

            if (glm_obj->control->do_fit == TRUE)  {
               
               if(i>0) {
                  if(input_data[i][ivox] == input_data[i-1][ivox])
                     same_test = TRUE;
                  else 
                     same_test = FALSE;
                  is_same = same_test * is_same;

                  if(glm_obj->control->binomial_n == 1)
                     glm_obj->response->values[i][0] = input_data[i][ivox]; 
                  else
                     glm_obj->response->values[i][0] = (input_data[i][ivox] /
                                                        glm_obj->control->binomial_n);
               }
            }
         }

         glm_obj->deg_free_tmp = 1.0 *(glm_obj->deg_free -
                                       glm_obj->num_missing);

         glm_obj->control->in_mask = input_data[input_num_buffers-1][ivox];
         
         if(glm_obj->contrast->num_avg > 0) {

            count = (glm_obj->response->num_rows - glm_obj->num_missing);

            if(count > 0)
               mean = tmp_sum / count;
            else 
               mean = INVALID_DATA;

            if(count >= 2) {
               variance = ((count * tmp_sum_sq - tmp_sum * tmp_sum) /
                           (count * (count-1)));
               norm_avg_resid = sqrt((count - 1) * variance);
               if(glm_obj->contrast->avg_array->num_test > 0) {
                  if(glm_obj->num_missing > 0) {
                     for(i=0; i<glm_obj->response->num_rows; i++) {
                        glm_obj->avg_resid->values[i][0] =\
                           INVALID_DATA;
                     }
                  }
                  else {
                     for(i=0; i<glm_obj->response->num_rows; i++) {
                        glm_obj->avg_resid->values[i][0] = \
                           (input_data[i][ivox] - mean) / norm_avg_resid;
                     }
                  }
               }
            }
            else {
               variance = INVALID_DATA;
               norm_avg_resid = INVALID_DATA;
            }
         }
      }

      if(is_same == TRUE) {
         glm_obj->control->status = IS_SAME;
      }

      if((glm_obj->num_missing > 0) && (glm_obj->control->deg_free_id < 0)) {
         glm_obj->control->status = MISSING_DATA;
      }

      if(in_search == FALSE)
         glm_obj->control->status = NOT_IN_SEARCH;

      if (glm_obj->deg_free_tmp <= 0)
         glm_obj->control->status = ZERO_DEG_FREE;

      if(((glm_obj->num_missing == 0) || (glm_obj->control->deg_free_id >= 0))
         && (is_same == FALSE) && (in_search == TRUE) && 
         (glm_obj->deg_free_tmp > 0)) {
         fit_glm(glm_obj);
      }

      num_fitted++;
      
      switch(glm_obj->control->status) {
      case CONVERGE:
         for (l=0; l<glm_obj->contrast->num_contrasts; l++) {
            cur_matrix = glm_obj->contrast->contrast_matrix_array[l];
            i = cur_matrix->out_id;
            if (glm_obj->output->values[l][0] != INVALID_DATA) {
               switch(cur_matrix->out_mode) {
               case T_STAT:
                  output_data[i][ivox] = INVALID_DATA;
                  if(cur_matrix->stdev_mode == VOXEL_SD) {
                     if (glm_obj->control->maxVoxelSdBasedValue >= 0) { 
                        if(fabs(glm_obj->output->values[l][0]) <= 
                           glm_obj->control->maxVoxelSdBasedValue) 
                           output_data[i][ivox]=glm_obj->output->values[l][0]; 
                        else
                           output_data[i][ivox]=(glm_obj->output->values[l][0]*
                                        glm_obj->control->maxVoxelSdBasedValue/
                                        fabs(glm_obj->output->values[l][0]));
                     }
                     else 
                        output_data[i][ivox] = glm_obj->output->values[l][0];
                  }
                  else
                     output_data[i][ivox] = glm_obj->output->values[l][0];
                  break;

               default:
                  output_data[i][ivox] = glm_obj->output->values[l][0];
                  break;
               }
            }
            else {
               output_data[i][ivox] = INVALID_DATA;
            }
         }
         break;

      default:
         for (l=0; l<glm_obj->contrast->num_contrasts; l++) {
            cur_matrix = glm_obj->contrast->contrast_matrix_array[l];
            i = cur_matrix->out_id;
            output_data[i][ivox]=INVALID_DATA;
         }
         break;
      }
      
      for(l=0; l<glm_obj->contrast->num_avg; l++) {

         cur_avg = glm_obj->contrast->avg_array->avg_info[l];
         i = cur_avg->out_id;

         switch(cur_avg->out_mode) {
         case T_STAT:
            output_data[i][ivox] = INVALID_DATA;
            if (cur_avg->stdev_mode == VOXEL_SD) {
               if (variance != INVALID_DATA) {
                  tmp = mean / sqrt(variance / count);
                  if (glm_obj->control->maxVoxelSdBasedValue > 0.0) {
                     if(fabs(tmp) <= glm_obj->control->maxVoxelSdBasedValue)
                        output_data[i][ivox] = tmp;
                     else
                        output_data[i][ivox] = (tmp *
                                glm_obj->control->maxVoxelSdBasedValue /
                                                fabs(tmp));
                  }
                  else 
                     output_data[i][ivox] = tmp;
               }
               else
                  output_data[i][ivox] = INVALID_DATA;
            }
            else if (cur_avg->stdev_mode == POOLED_SD) {
               if (variance != INVALID_DATA)
                  output_data[i][ivox] = mean * sqrt(count);
               else
                  output_data[i][ivox] = INVALID_DATA;
            }
            break;
                  
         case BETA_HAT:
            output_data[i][ivox] = mean;
            break;

         case STDEV_BETA:
            if (variance != INVALID_DATA)
               output_data[i][ivox] = sqrt(variance / count);
            else
               output_data[i][ivox] = INVALID_DATA;
            break;

         case P_CORR:
            output_data[i][ivox] = INVALID_DATA;
            if(cur_avg->stdev_mode == VOXEL_SD) {
               if(variance != INVALID_DATA) {
                  tmp = mean / sqrt(variance / count);
                  output_data[i][ivox] = tmp / sqrt(tmp*tmp + count - 1);
               }
               else
                  output_data[i][ivox] = INVALID_DATA;
            }
            else if (cur_avg->stdev_mode == POOLED_SD) {
               if(variance != INVALID_DATA)
                  output_data[i][ivox] = mean * sqrt(count);
               else
                  output_data[i][ivox] = INVALID_DATA;
            }
            break;

         default:
            fprintf(stderr,"\nError: improper out_mode for avg output.\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
            break;
         }
      }

      if (glm_obj->control->dev_id > 0) {
         if(glm_obj->control->status == CONVERGE)
            output_data[glm_obj->control->dev_id][ivox] = (glm_obj->deviance *
                                                           glm_obj->scale_est);
         else
            output_data[glm_obj->control->dev_id][ivox] = INVALID_DATA;
      }

      if (glm_obj->control->scale_id > 0) {
         if(glm_obj->control->status == CONVERGE)
            output_data[glm_obj->control->scale_id][ivox] = glm_obj->scale_est;
         else
            output_data[glm_obj->control->scale_id][ivox] = INVALID_DATA;
      }

      if (glm_obj->control->avg_dev_id > 0) {
            output_data[glm_obj->control->avg_dev_id][ivox] = norm_avg_resid;
      }

      if (glm_obj->control->err_id > 0) {
         output_data[glm_obj->control->err_id][ivox] =glm_obj->control->status;
      }

      if (glm_obj->control->count_id > 0) {
         output_data[glm_obj->control->count_id][ivox] = count-1;
      }

      if (glm_obj->control->deg_free_id > 0) { 
         if (glm_obj->deg_free_tmp > 0)
            output_data[glm_obj->control->deg_free_id][ivox] =\
               glm_obj->deg_free_tmp;
         else
            output_data[glm_obj->control->deg_free_id][ivox] = INVALID_DATA;
      }


      /* Sums deviance over voxels inside stdev mask */

      if (glm_obj->control->pool_id > 0 ) {
         if ((input_data[glm_obj->control->pool_id][ivox] != 0.0) &&
             (input_data[glm_obj->control->pool_id][ivox] != 1.0)) {
            fprintf(stderr,"\nError: pooling mask values must be 0 or 1.\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }
         if ((input_data[glm_obj->control->pool_id][ivox] == 1) &&
             (glm_obj->num_missing == 0)) {
            glm_obj->pooled_dev += glm_obj->pearson;
            glm_obj->num_pool++;
         }
         if(glm_obj->contrast->avg_array != NULL) {
            if (glm_obj->contrast->avg_array->num_pool > 0) {
               if((input_data[glm_obj->control->pool_id][ivox] == 1) &&
                   (glm_obj->num_missing == 0)){
                  if(variance != INVALID_DATA) {
                     glm_obj->avg_pooled_dev += variance;
                     glm_obj->avg_num_pool++;
                  }
               }
            }
         }
      }

      if(lambda_buffer != NULL) {

         if ((glm_obj->num_missing == 0) &&
             ((glm_obj->control->status == CONVERGE) || 
             (glm_obj->control->status == CONT))) {
            check = 1.0;
         }
         else
            check = 0.0;

         update_lambda_buffer(lambda_buffer, glm_obj, check);

         if ((glm_obj->control->status == CONVERGE) || 
             (glm_obj->control->status == CONT)) {

            if ((glm_obj->num_missing == 0) && 
                (lambda_buffer->total >= lambda_buffer->ring_size)) {
               calculate_smoothness(lambda_buffer);

               if ((glm_obj->family == GAUSSIAN) && 
                   (lambda_buffer->fwhm_gaussian != NULL)) {
                  fwhm_general = lambda_buffer->fwhm_gaussian;
                  if (fwhm_general->out_id > 0) {

                     if (fwhm_general->fwhm >= 0.0) {

                        if (fwhm_general->fwhm <= 100.0) {
                           output_data[fwhm_general->out_id][ivox] =\
                              fwhm_general->fwhm;
                        }
                        else {
                           output_data[fwhm_general->out_id][ivox] =\
                              100.0;
                        }
                     }
                     else {
                        output_data[fwhm_general->out_id][ivox] =\
                           INVALID_DATA;
                     }
                  }
               }

               if (lambda_buffer->fwhm_avg != NULL) {
                  fwhm_general = lambda_buffer->fwhm_avg;

                  if (fwhm_general->out_id > 0) {
                     if (fwhm_general->fwhm >= 0.0) {

                        if (fwhm_general->fwhm <= 100.0) {
                           output_data[fwhm_general->out_id][ivox] =\
                              fwhm_general->fwhm;
                        }
                        else {
                           output_data[fwhm_general->out_id][ivox] =\
                              100.0;
                        }
                     }
                     else {
                        output_data[fwhm_general->out_id][ivox] =\
                           INVALID_DATA;
                     }
                  }
               }

               if ((glm_obj->family != GAUSSIAN) && 
                   (lambda_buffer->fwhm_general != NULL)) {
                  for(itest=0; itest<lambda_buffer->num_test; itest++) {
                     fwhm_general = lambda_buffer->fwhm_general[itest];

                     if (fwhm_general->fwhm >= 0.0) {

                        if (fwhm_general->fwhm <= 100.0) {
                           output_data[fwhm_general->out_id][ivox] =\
                              fwhm_general->fwhm;
                        }
                        else {
                           output_data[fwhm_general->out_id][ivox] =\
                              100.0;
                        }
                     }
                     else {
                        output_data[fwhm_general->out_id][ivox] =\
                           INVALID_DATA;
                     }
                  }

                  if (lambda_buffer->fwhm_simple != NULL) {
                     fwhm_general = lambda_buffer->fwhm_simple;
                     if (fwhm_general->fwhm >= 0.0) {

                        if (fwhm_general->fwhm <= 100.0) {
                           output_data[fwhm_general->out_id][ivox] =\
                              fwhm_general->fwhm;
                        }
                        else {
                           output_data[fwhm_general->out_id][ivox] =\
                              100.0;
                        }
                     }
                     else {
                        output_data[fwhm_general->out_id][ivox] =\
                           INVALID_DATA;
                     }
                  }
               }
            }
            else {
               if ((glm_obj->family != GAUSSIAN) &&
                   (lambda_buffer->fwhm_general != NULL)){
                  for(itest=0; itest<lambda_buffer->num_test; itest++) {
                     fwhm_general = lambda_buffer->fwhm_general[itest];
                     output_data[fwhm_general->out_id][ivox] =\
                        INVALID_DATA;
                  }
               }

               if ((glm_obj->family == GAUSSIAN) && 
                   (lambda_buffer->fwhm_gaussian != NULL)) {
                  fwhm_general = lambda_buffer->fwhm_gaussian;
                  if (fwhm_general->out_id > 0) {
                     output_data[fwhm_general->out_id][ivox] =\
                        INVALID_DATA;
                  }
               }

               if (lambda_buffer->fwhm_avg != NULL) {
                  fwhm_general = lambda_buffer->fwhm_avg;
                  if (fwhm_general->out_id > 0) {
                     output_data[fwhm_general->out_id][ivox] =\
                        INVALID_DATA;
                  }
               }

            }
         }
         else {
            if ((glm_obj->family != GAUSSIAN) &&
                (lambda_buffer->fwhm_general != NULL)){
               for(itest=0; itest<lambda_buffer->num_test; itest++) {
                  fwhm_general = lambda_buffer->fwhm_general[itest];
                  if (fwhm_general->out_id > 0) 
                     output_data[fwhm_general->out_id][ivox] = INVALID_DATA;
               }
            }

            if (lambda_buffer->fwhm_gaussian != NULL) {
                  fwhm_general = lambda_buffer->fwhm_gaussian;
               if (fwhm_general->out_id > 0) {
                  output_data[fwhm_general->out_id][ivox] =\
                  INVALID_DATA;
               }
            }

            if (lambda_buffer->fwhm_avg != NULL) {
               fwhm_general = lambda_buffer->fwhm_avg;
               if (fwhm_general->out_id > 0) {
                  output_data[fwhm_general->out_id][ivox] =\
                  INVALID_DATA;
               }
            }
         }
      }

      glm_obj->control->last_status = glm_obj->control->status;

      if(glm_obj->family != GAUSSIAN) {
         glm_obj->control->status = CONT;
     }
      else
         glm_obj->control->status = CONVERGE;
   }

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : voxel_function_pool_con
@INPUT      : Standard for voxel loop
@OUTPUT     : Standard for voxel loop
@RETURNS    : (nothing)
@DESCRIPTION: Routine doing working on each input buffer
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 25, 1997 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void voxel_function_pool_con(void *caller_data, long num_voxels,
                    int input_num_buffers, int input_vector_length, 
                    double *input_data[],
                    int output_num_buffers, int output_vector_length, 
                    double *output_data[],
                    Loop_Info *loop_info)
/* ARGSUSED */
{
   long ivox, num_values;
   Contrast_Matrix *cur_matrix;
   double tmp;
   double pooled_dev, pooled_sd;
   int deg_free;

   /* Get pointer to cur_matrix */

   cur_matrix = (Contrast_Matrix *) caller_data;

   num_values = num_voxels * input_vector_length;

   pooled_dev = cur_matrix->pooled_dev;
   pooled_sd = cur_matrix->pooled_sd;

   switch(cur_matrix->out_mode) {
   case T_STAT:
      for(ivox=0; ivox < num_values; ivox++) {
         output_data[0][ivox] = INVALID_DATA;
         if (input_data[0][ivox] != INVALID_DATA) {
            tmp = input_data[0][ivox] / pooled_sd;
            if (cur_matrix->maxVoxelSdBasedValue > 0) {
               if(fabs(tmp) <= cur_matrix->maxVoxelSdBasedValue)
                  output_data[0][ivox] = tmp;
               else
                  output_data[0][ivox] = (cur_matrix->maxVoxelSdBasedValue *
                                          tmp / fabs(tmp));
            }
            else
               output_data[0][ivox] = tmp;
         }
         else
            output_data[0][ivox] = INVALID_DATA;
      }
      break;

   case F_STAT:
      for(ivox=0; ivox < num_values; ivox++) {
         if(input_data[0][ivox] != INVALID_DATA)
            output_data[0][ivox] = input_data[0][ivox] / pooled_dev;
         else
            output_data[0][ivox] = INVALID_DATA;
      }
      break;

   case P_CORR:
      for(ivox=0; ivox < num_values; ivox++) {
         if(input_data[0][ivox] != INVALID_DATA) {
            if (input_num_buffers == 2)
               deg_free = input_data[1][ivox];
            else
               deg_free = cur_matrix->deg_free * 1.0;
            tmp = input_data[0][ivox] / pooled_sd;
            output_data[0][ivox] = tmp / sqrt(tmp*tmp + deg_free);
         }
         else
            output_data[0][ivox] = INVALID_DATA;
      }
      break;

   default:
      break;
   }

   return;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : voxel_function_pool_avg
@INPUT      : Standard for voxel loop
@OUTPUT     : Standard for voxel loop
@RETURNS    : (nothing)
@DESCRIPTION: Routine doing working on each input buffer
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 25, 1997 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void voxel_function_pool_avg(void *caller_data, long num_voxels,
                    int input_num_buffers, int input_vector_length, 
                    double *input_data[],
                    int output_num_buffers, int output_vector_length, 
                    double *output_data[],
                    Loop_Info *loop_info)
/* ARGSUSED */
{
   long ivox, num_values;
   Avg_Info *cur_avg;
   double tmp;
   double pooled_dev, pooled_sd;
   int deg_free;

   /* Get pointer to cur_matrix */

   cur_avg = (Avg_Info *) caller_data;

   num_values = num_voxels * input_vector_length;

   pooled_dev = cur_avg->pooled_dev;
   pooled_sd = cur_avg->pooled_sd;

   switch(cur_avg->out_mode) {
   case T_STAT:
      for(ivox=0; ivox < num_values; ivox++) {
         output_data[0][ivox] = INVALID_DATA;
         if (input_data[0][ivox] != INVALID_DATA) {
            tmp = input_data[0][ivox] / pooled_sd;
            if (cur_avg->maxVoxelSdBasedValue > 0) {
               if(fabs(tmp) <= cur_avg->maxVoxelSdBasedValue)
                  output_data[0][ivox] = tmp;
               else
                  output_data[0][ivox] = (cur_avg->maxVoxelSdBasedValue *
                                          tmp / fabs(tmp));
            }
            else
               output_data[0][ivox] = tmp;
         }
         else
            output_data[0][ivox] = INVALID_DATA;
      }
      break;

   case F_STAT:
      for(ivox=0; ivox < num_values; ivox++) {
         if(input_data[0][ivox] != INVALID_DATA)
            output_data[0][ivox] = input_data[0][ivox] / pooled_dev;
         else
            output_data[0][ivox] = INVALID_DATA;
      }
      break;

   case P_CORR:
      for(ivox=0; ivox < num_values; ivox++) {
         if(input_data[0][ivox] != INVALID_DATA) {
            if (input_num_buffers == 2)
               deg_free = input_data[1][ivox] + cur_avg->df_offset;
            else
               deg_free = cur_avg->deg_free * 1.0;
            tmp = input_data[0][ivox] / pooled_sd;
            output_data[0][ivox] = tmp / sqrt(tmp*tmp + deg_free);
         }
         else
            output_data[0][ivox] = INVALID_DATA;
      }
      break;

   default:
      break;
   }

   return;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : start_function
@INPUT      : Standard for voxel loop
@OUTPUT     : Standard for voxel loop
@RETURNS    : (nothing)
@DESCRIPTION: Routine setting up output buffers
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 25, 1997 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void start_function(void *caller_data, long num_voxels, 
                    int output_num_buffers, int output_vector_length, 
                    double *output_data[], Loop_Info *loop_info)
/* ARGSUSED */
{
   long ivox, num_values;
   int ibuff;

   num_values = num_voxels * output_vector_length;

   for(ivox = 0; ivox<num_values; ivox++) {
      for(ibuff=0; ibuff<output_num_buffers; ibuff++) {
         output_data[ibuff][ivox] = INVALID_DATA;
      }
   }

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : end_function
@INPUT      : Standard for voxel loop
@OUTPUT     : Standard for voxel loop
@RETURNS    : (nothing)
@DESCRIPTION: Routine computing final values for output buffers
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 25, 1997 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void end_function(void *caller_data, long num_voxels, 
                  int output_num_buffers, int output_vector_length, 
                  double *output_data[], Loop_Info *loop_info)
/* ARGSUSED */
{
   long ivox, num_values;
   
   num_values = num_voxels * output_vector_length;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_f_stat_raw
@INPUT      : dst - client data passed by ParseArgv
              key - matching key in argv
              nextarg - argument following key in argv
@OUTPUT     : (none)
@RETURNS    : TRUE since nextarg is used.
@DESCRIPTION: Gets a contrast converts it to a string 
              and adds it to an array of contrasts to be parsed later
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 21, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int get_f_stat_raw(char *dst, char *key, int argc, char **argv)
{
   int num_output, i, length;
   Contrast_Raw *contrast;
   Contrast_Raw_Array *contrast_list;
   
   /* Check for the right amount of arguments */
   if (argc < 4) {
      (void) fprintf(stderr, 
                     "\"%s\" option requires three additional arguments\n",
                     key);
      exit(EXIT_FAILURE);
   }

   /* Set up pointers to filename, method and contrast string */
   contrast_list = (Contrast_Raw_Array *) dst;
   
   /* Get raw contrast */
   
   contrast = GI_MALLOC(sizeof(Contrast_Raw));

   length = strlen(argv[0]);
   contrast->outfile_name = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->outfile_name, argv[0]) == NULL) {
      fprintf(stderr,"\nError copying outfile name for outfile %s\n", argv[0]);
      return FALSE;
   }

   length = strlen(argv[1]);
   contrast->stdev_type = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->stdev_type, argv[1]) == NULL) {
      fprintf(stderr,"\nError copying in_type for outfile %s\n", argv[0]);
      return FALSE;
   }

   length = strlen(argv[2]);
   contrast->in_type = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->in_type, argv[2]) == NULL) {
      fprintf(stderr,"\nError copying in_type for outfile %s\n", argv[0]);
      return FALSE;
   }

   length = strlen(argv[3]);
   contrast->raw_contrast = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->raw_contrast, argv[3]) == NULL) {
      fprintf(stderr,"\nError copying raw contrast for outfile %s\n", argv[0]);
      return FALSE;
   }

   contrast->is_avg = FALSE;
   contrast->out_type = NULL; 
   
   /* Set up or reallocate memory for contrast_list */

   if(contrast_list->contrast_array == NULL) {
      contrast_list->contrast_array = GI_MALLOC((MAX_NUM_CONTRASTS 
                                    * sizeof(contrast)));
      contrast_list->num_output = 1;
      contrast_list->num_avg = 0;
   }
   else {
      contrast_list->num_output++;
   }

   num_output = contrast_list->num_output;
   contrast_list->contrast_array[num_output-1] = contrast;
     
   /* Collapse remaining arguments and free last three pointers of argv */

   for (i = 0; i < (argc - 4) ; i++) {
      argv[i] = argv[i+4];
   }

   return (argc - 4);
      
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_t_stat_raw
@INPUT      : dst - client data passed by ParseArgv
              key - matching key in argv
              nextarg - argument following key in argv
@OUTPUT     : (none)
@RETURNS    : TRUE since nextarg is used.
@DESCRIPTION: Gets a contrast converts it to a string 
              and adds it to an array of contrasts to be parsed later
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 21, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int get_t_stat_raw(char *dst, char *key, int argc, char **argv)
{
   int num_output, i, length;
   Contrast_Raw *contrast;
   Contrast_Raw_Array *contrast_list;
   
   /* Check for the right amount of arguments */
   if (argc < 5) {
      (void) fprintf(stderr, 
                     "\"%s\" option requires five additional arguments\n",
                     key);
      exit(EXIT_FAILURE);
   }

   /* Set up pointers to filename, method and contrast string */
   contrast_list = (Contrast_Raw_Array *) dst;
   
   /* Get raw contrast */
   
   contrast = GI_MALLOC(sizeof(Contrast_Raw));

   length = strlen(argv[0]);
   contrast->outfile_name = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->outfile_name, argv[0]) == NULL) {
      fprintf(stderr,"\nError copying outfile name for outfile %s\n", argv[0]);
      return FALSE;
   }

   length = strlen(argv[1]);
   contrast->out_type = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->out_type, argv[1]) == NULL) {
      fprintf(stderr,"\nError copying out_type for outfile %s\n", argv[0]);
      return FALSE;
   }

   length = strlen(argv[2]);
   contrast->stdev_type = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->stdev_type, argv[2]) == NULL) {
      fprintf(stderr,"\nError copying stdev_type for outfile %s\n", argv[0]);
      return FALSE;
   }

   length = strlen(argv[3]);
   contrast->in_type = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->in_type, argv[3]) == NULL) {
      fprintf(stderr,"\nError copying in_type for outfile %s\n", argv[0]);
      return FALSE;
   }

   length = strlen(argv[4]);
   contrast->raw_contrast = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->raw_contrast, argv[4]) == NULL) {
      fprintf(stderr,"\nError copying raw contrast for outfile %s\n", argv[0]);
      return FALSE;
   }

   contrast->is_avg = FALSE;

   /* Set up or reallocate memory for contrast_list */

   if(contrast_list->contrast_array == NULL) {
      contrast_list->contrast_array = GI_MALLOC((MAX_NUM_CONTRASTS 
                                    * sizeof(contrast)));
      contrast_list->num_output = 1;
      contrast_list->num_avg = 0;
   }
   else {
      contrast_list->num_output++;
   }

   num_output = contrast_list->num_output;
   contrast_list->contrast_array[num_output-1] = contrast;
     
   /* Collapse remaining arguments and free last three pointers of argv */

   for (i = 0; i < (argc - 5) ; i++) {
      argv[i] = argv[i+5];
   }

   return (argc - 5);
      
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_avg_raw
@INPUT      : dst - client data passed by ParseArgv
              key - matching key in argv
              nextarg - argument following key in argv
@OUTPUT     : (none)
@RETURNS    : TRUE since nextarg is used.
@DESCRIPTION: Gets a contrast converts it to a string 
              and adds it to an array of contrasts to be parsed later
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 21, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int get_avg_raw(char *dst, char *key, int argc, char **argv)
{
   int num_output, i, length;
   Contrast_Raw *contrast;
   Contrast_Raw_Array *contrast_list;
   
   /* Check for the right amount of arguments */
   if (argc < 3) {
      (void) fprintf(stderr, 
                     "\"%s\" option requires three additional arguments\n",
                     key);
      exit(EXIT_FAILURE);
   }

   /* Set up pointers to filename, method and contrast string */
   contrast_list = (Contrast_Raw_Array *) dst;
   
   /* Get raw contrast */
   
   contrast = GI_MALLOC(sizeof(Contrast_Raw));

   length = strlen(argv[0]);
   contrast->outfile_name = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->outfile_name, argv[0]) == NULL) {
      fprintf(stderr,"\nError copying outfile name for outfile %s\n", argv[0]);
      return FALSE;
   }

   length = strlen(argv[1]);
   contrast->out_type = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->out_type, argv[1]) == NULL) {
      fprintf(stderr,"\nError copying out_type for outfile %s\n", argv[0]);
      return FALSE;
   }

   length = strlen(argv[2]);
   contrast->stdev_type = GI_MALLOC(sizeof(char) * (length + 1));
   if(strcpy(contrast->stdev_type, argv[2]) == NULL) {
      fprintf(stderr,"\nError copying stdev_type for outfile %s\n", argv[0]);
      return FALSE;
   }

   contrast->in_type = NULL;
   contrast->is_avg = TRUE;

   /* Set up or reallocate memory for contrast_list */

   if(contrast_list->contrast_array == NULL) {
      contrast_list->contrast_array = GI_MALLOC((MAX_NUM_CONTRASTS 
                                    * sizeof(contrast)));
      contrast_list->num_output = 1;
      contrast_list->num_avg = 1;
   }
   else {
      contrast_list->num_output++;
      contrast_list->num_avg++;
   }

   num_output = contrast_list->num_output;
   contrast_list->contrast_array[num_output-1] = contrast;
     
   /* Collapse remaining arguments and free last three pointers of argv */

   for (i = 0; i < (argc - 3) ; i++) {
      argv[i] = argv[i+3];
   }

   return (argc - 3);
      
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_stat_variable
@INPUT      : list of files from glim_image 
@OUTPUT     : 
@RETURNS    : (nothing)
@DESCRIPTION: creates t_statistic or F_statistic variable for output files
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Oct. 14, 97 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void create_stat_variable(Program_Data *program_data, char **input_files)
{
   char *outfile = NULL;
   int mincid, varid;
   int ifile, itest;
   int icont;
   Contrast_Matrix *cur_matrix;
   Avg_Info *cur_avg;
   Glm_Object *glm_obj;
   Lambda *lambda_buffer;
   Fwhm_Info *fwhm_general;
   char *var_name= "glim_statistics";
   char *glim_par_str = "glim-parameters";
   char *t_stat_str = "MNI t-statistics variable";
   char *f_stat_str = "MNI F-statistics variable";
   char *p_corr_str = "MNI partial correlation variable";
   char *beta_hat_str = "MNI regression parameter variable";
   char *stdev_str = "MNI standard deviation variable";
   char *num_df_str = "numerator-df";
   char *den_df_str = "denominator-df";
   char *pooled_file_str = "pool-mask-file";
   char *pooled_sd_str = "pooled-stdev";
   char *pooled_dev_str = "pooled-deviance";
   char *contrast_str = "contrast-matrix";
   char *column_str = "regress-column-number";
   char *average_str = "comparison of means";
   char *glm_str = "generalized linear model (GLM)";
   char *model_str = "model-type";
   char *link_str = "link-function";
   char *variance_fun_str = "variance-function";
   char *family_str = "exponential-family";
   char *param_str = NULL, *contr_str = NULL;
   char *int_lambda_str = "integrated-det-lambda";
   char *search_vol_str = "search-volume";
   char *search_file_str = "search-region-file";
   char *pooled_fwhm_str = "pooled-FWHM-estimate";
   char *int_fwhm_str = "integrated-FWHM-estimate";
   char *fwhm_file_str = "local-FWHM-file";
   char *fwhm_matrix_str = "pooled-lambda-matrix";
   char *fwhm_matrix_gaussian_str = NULL;
   char *fwhm_matrix_avg_str = NULL;
   char *num_fwhm_str = "num-voxel-pool-FWHM";
   glm_obj = program_data->glm_obj;
   lambda_buffer = program_data->lambda_buffer;

   create_glim_parameters(glm_obj, input_files, &param_str);

   /* Get smoothness values and attach them to proper contrasts */

   if(lambda_buffer != NULL) {
      if (glm_obj->family != GAUSSIAN) {
         for(itest=0; itest<glm_obj->contrast->num_test; itest++) {
            fwhm_general = lambda_buffer->fwhm_general[itest];

            icont = glm_obj->contrast->test_map[itest];
            cur_matrix = glm_obj->contrast->contrast_matrix_array[icont];

            cur_matrix->int_lambda = fwhm_general->int_lambda;
            cur_matrix->num_fwhm = fwhm_general->num_calc;
            cur_matrix->search_vol = fwhm_general->search_vol;
            cur_matrix->fwhm_int = fwhm_general->fwhm_int;
            cur_matrix->fwhm_outfile = fwhm_general->outfile;
         }
      }
      else if ((glm_obj->control->do_fit == TRUE) && 
               (lambda_buffer->fwhm_gaussian != NULL)){
         sprint_matrix(&fwhm_matrix_gaussian_str, 
                       lambda_buffer->fwhm_gaussian->lambda_pool); 
         delete_matrix(lambda_buffer->fwhm_gaussian->lambda_pool);
      }

      if(lambda_buffer->fwhm_avg != NULL) {
         sprint_matrix(&fwhm_matrix_avg_str, 
                       lambda_buffer->fwhm_avg->lambda_pool);
         delete_matrix(lambda_buffer->fwhm_avg->lambda_pool);
      }
      
   }

   /* Loop over output files, defining variable glim_statistics */

   itest = 0;

   for(ifile = 0; ifile < glm_obj->contrast->num_contrasts; ifile++) {
      
      outfile = glm_obj->contrast->contrast_matrix_array[ifile]->outfile;
      fprintf(stderr,"\nUpdating %s's header.",outfile);
      mincid = miopen(outfile, NC_WRITE);

      (void) ncredef(mincid);

      varid = ncvardef(mincid, var_name, NC_LONG, 0, 0);

      cur_matrix = glm_obj->contrast->contrast_matrix_array[ifile];
      
      switch(cur_matrix->out_mode) {
      case T_STAT:
         (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
         (void) miattputstr(mincid, varid, MIvarid, t_stat_str);
         (void) miattputstr(mincid, varid, glim_par_str, param_str);
         (void) miattputint(mincid, varid, den_df_str, glm_obj->deg_free);
         if(cur_matrix->stdev_mode == POOLED_SD) {
            (void) miattputstr(mincid, varid, pooled_file_str, 
                               input_files[glm_obj->design_matrix->num_rows]);
            (void) miattputdbl(mincid, varid, pooled_sd_str,
                              cur_matrix->pooled_sd);
         }
         if(lambda_buffer != NULL) {
            if(glm_obj->family != GAUSSIAN) {
               (void) miattputdbl(mincid, varid, int_lambda_str,
                                  cur_matrix->int_lambda);
               (void) miattputdbl(mincid, varid, num_fwhm_str,
                                  cur_matrix->num_fwhm);
               if(glm_obj->search_mask != NULL) {
                  (void) miattputstr(mincid, varid, search_file_str,
                                     glm_obj->search_mask);
               }
               (void) miattputdbl(mincid, varid, search_vol_str,
                                  cur_matrix->search_vol);
               (void) miattputdbl(mincid, varid, int_fwhm_str,
                                  cur_matrix->fwhm_int);
               (void) miattputstr(mincid, varid, fwhm_file_str,
                                  cur_matrix->fwhm_outfile);
            }
            else {
               (void) miattputdbl(mincid, varid, int_lambda_str,
                                  lambda_buffer->fwhm_gaussian->int_lambda);
               (void) miattputdbl(mincid, varid, pooled_fwhm_str,
                                  lambda_buffer->fwhm_gaussian->fwhm_pool);
               (void) miattputdbl(mincid, varid, int_fwhm_str,
                                  lambda_buffer->fwhm_gaussian->fwhm_int);
               if(glm_obj->search_mask != NULL) {
                  (void) miattputstr(mincid, varid, search_file_str,
                                     glm_obj->search_mask);
               }
               (void) miattputdbl(mincid, varid, search_vol_str,
                                  lambda_buffer->fwhm_gaussian->search_vol);
               (void) miattputstr(mincid, varid, fwhm_matrix_str,
                               fwhm_matrix_gaussian_str);
               fprintf(stderr,"\n%s", fwhm_matrix_gaussian_str);
            }
         }
         (void) miattputstr(mincid, varid, model_str, glm_str);
         if (cur_matrix->contrast_matrix == NULL)
            (void) miattputint(mincid, varid, column_str, 
                               cur_matrix->column_num+1);
         else {
            sprint_matrix(&contr_str, cur_matrix->contrast_matrix);
            (void) miattputstr(mincid, varid, contrast_str, 
                               contr_str); 
            GI_FREE(contr_str);
         }
         (void) miattputstr(mincid, varid, family_str,
                            glm_obj->family_names[glm_obj->family]);
         (void) miattputstr(mincid, varid, link_str,
                            glm_obj->link_names[glm_obj->link]);
         (void) miattputstr(mincid, varid, variance_fun_str,
                            glm_obj->variance_function_names[glm_obj->variance_function]);
         break;

      case F_STAT:
         (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
         (void) miattputstr(mincid, varid, MIvarid, f_stat_str);
         (void) miattputstr(mincid, varid, glim_par_str, param_str);
         (void) miattputint(mincid, varid, den_df_str, glm_obj->deg_free);
         (void) miattputint(mincid, varid, num_df_str, 
                            cur_matrix->contrast_matrix->num_rows);

         if(cur_matrix->stdev_mode == POOLED_SD) {
            (void) miattputstr(mincid, varid, pooled_file_str, 
                               input_files[glm_obj->design_matrix->num_rows]);
            (void) miattputdbl(mincid, varid, pooled_dev_str,
                              cur_matrix->pooled_dev);
         }
         if (lambda_buffer != NULL) {
            (void) miattputdbl(mincid, varid, int_lambda_str,
                               lambda_buffer->fwhm_gaussian->int_lambda);
            (void) miattputdbl(mincid, varid, pooled_fwhm_str,
                               lambda_buffer->fwhm_gaussian->fwhm_pool);
            (void) miattputdbl(mincid, varid, int_fwhm_str,
                               lambda_buffer->fwhm_gaussian->fwhm_int);
            if(glm_obj->search_mask != NULL) {
               (void) miattputstr(mincid, varid, search_file_str,
                                  glm_obj->search_mask);
            }
            (void) miattputdbl(mincid, varid, search_vol_str,
                               lambda_buffer->fwhm_gaussian->search_vol);
            (void) miattputstr(mincid, varid, fwhm_matrix_str,
                               fwhm_matrix_gaussian_str);
            fprintf(stderr,"\n%s", fwhm_matrix_gaussian_str);
         }
         (void) miattputstr(mincid, varid, model_str, glm_str);
         sprint_matrix(&contr_str, cur_matrix->contrast_matrix);
         (void) miattputstr(mincid, varid, contrast_str, 
                            contr_str); 
         GI_FREE(contr_str);


         (void) miattputstr(mincid, varid, family_str,
                            glm_obj->family_names[glm_obj->family]);
         (void) miattputstr(mincid, varid, link_str,
                            glm_obj->link_names[glm_obj->link]);
         (void) miattputstr(mincid, varid, variance_fun_str,
                            glm_obj->variance_function_names[glm_obj->variance_function]);
         break;

      case P_CORR:
         (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
         (void) miattputstr(mincid, varid, MIvarid, p_corr_str);
         (void) miattputstr(mincid, varid, glim_par_str, param_str);
         (void) miattputint(mincid, varid, den_df_str, glm_obj->deg_free);
         if(cur_matrix->stdev_mode == POOLED_SD) {
            (void) miattputstr(mincid, varid, pooled_file_str, 
                               input_files[glm_obj->design_matrix->num_rows]);
            (void) miattputdbl(mincid, varid, pooled_sd_str,
                              cur_matrix->pooled_sd);
         }
         if (lambda_buffer != NULL) {
            if (glm_obj->family != GAUSSIAN) {
               (void) miattputdbl(mincid, varid, int_lambda_str,
                                  cur_matrix->int_lambda);
               (void) miattputdbl(mincid, varid, num_fwhm_str,
                                  cur_matrix->num_fwhm);
               if(glm_obj->search_mask != NULL) {
                  (void) miattputstr(mincid, varid, search_file_str,
                                     glm_obj->search_mask);
               }
               (void) miattputdbl(mincid, varid, search_vol_str,
                                  cur_matrix->search_vol);
               (void) miattputdbl(mincid, varid, int_fwhm_str,
                                  cur_matrix->fwhm_int);
               (void) miattputstr(mincid, varid, fwhm_file_str,
                                  cur_matrix->fwhm_outfile);
            }
            else {
               (void) miattputdbl(mincid, varid, int_lambda_str,
                                  lambda_buffer->fwhm_gaussian->int_lambda);
               (void) miattputdbl(mincid, varid, pooled_fwhm_str,
                                  lambda_buffer->fwhm_gaussian->fwhm_pool);
               (void) miattputdbl(mincid, varid, int_fwhm_str,
                                  lambda_buffer->fwhm_gaussian->fwhm_int);
               if(glm_obj->search_mask != NULL) {
                  (void) miattputstr(mincid, varid, search_file_str,
                                     glm_obj->search_mask);
               }
               (void) miattputdbl(mincid, varid, search_vol_str,
                                  lambda_buffer->fwhm_gaussian->search_vol);
               (void) miattputstr(mincid, varid, fwhm_matrix_str,
                                  fwhm_matrix_gaussian_str);
               (void) miattputstr(mincid, varid, fwhm_matrix_str,
                                  fwhm_matrix_gaussian_str);
            }
         }
         (void) miattputstr(mincid, varid, model_str, glm_str);
         if (cur_matrix->contrast_matrix == NULL)
            (void) miattputint(mincid, varid, column_str, 
                               cur_matrix->column_num+1);
         else {
            sprint_matrix(&contr_str, cur_matrix->contrast_matrix);
            (void) miattputstr(mincid, varid, contrast_str, 
                               contr_str); 
            GI_FREE(contr_str);
         }
         (void) miattputstr(mincid, varid, family_str,
                            glm_obj->family_names[glm_obj->family]);
         (void) miattputstr(mincid, varid, link_str,
                            glm_obj->link_names[glm_obj->link]);
         (void) miattputstr(mincid, varid, variance_fun_str,
                            glm_obj->variance_function_names[glm_obj->variance_function]);

         break;

      case BETA_HAT:
         (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
         (void) miattputstr(mincid, varid, MIvarid, beta_hat_str);
         (void) miattputstr(mincid, varid, glim_par_str, param_str);
         (void) miattputstr(mincid, varid, model_str, glm_str);
         if (cur_matrix->contrast_matrix == NULL)
            (void) miattputint(mincid, varid, column_str, 
                               cur_matrix->column_num+1);
         else {
            sprint_matrix(&contr_str, cur_matrix->contrast_matrix);
            (void) miattputstr(mincid, varid, contrast_str, 
                               contr_str); 
            GI_FREE(contr_str);
         }
         (void) miattputstr(mincid, varid, family_str,
                            glm_obj->family_names[glm_obj->family]);
         (void) miattputstr(mincid, varid, link_str,
                            glm_obj->link_names[glm_obj->link]);
         (void) miattputstr(mincid, varid, variance_fun_str,
                            glm_obj->variance_function_names[glm_obj->variance_function]);

         break;

      case STDEV_BETA:
         (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
         (void) miattputstr(mincid, varid, MIvarid, stdev_str);
         (void) miattputstr(mincid, varid, glim_par_str, param_str);
         (void) miattputint(mincid, varid, den_df_str, glm_obj->deg_free);
         (void) miattputstr(mincid, varid, model_str, glm_str);
         if (cur_matrix->contrast_matrix == NULL)
            (void) miattputint(mincid, varid, column_str, 
                               cur_matrix->column_num+1);
         else {
            sprint_matrix(&contr_str, cur_matrix->contrast_matrix);
            (void) miattputstr(mincid, varid, contrast_str, 
                               contr_str); 
            GI_FREE(contr_str);
         }
         (void) miattputstr(mincid, varid, family_str,
                            glm_obj->family_names[glm_obj->family]);
         (void) miattputstr(mincid, varid, link_str,
                            glm_obj->link_names[glm_obj->link]);
         (void) miattputstr(mincid, varid, variance_fun_str,
                            glm_obj->variance_function_names[glm_obj->variance_function]);

         break;

      default:
         break;
      }

      ncendef(mincid);

      miclose(mincid);

      fprintf(stderr," Done.\n");
   }

   for(ifile = 0; ifile < glm_obj->contrast->num_avg; ifile++) {
      
      outfile = glm_obj->contrast->avg_array->avg_info[ifile]->outfile;
      fprintf(stderr,"\nUpdating %s's header.",outfile);
      mincid = miopen(outfile, NC_WRITE);

      (void) ncredef(mincid);

      varid = ncvardef(mincid, var_name, NC_LONG, 0, 0);

      cur_avg = glm_obj->contrast->avg_array->avg_info[ifile];
      cur_avg->deg_free = glm_obj->design_matrix->num_rows - 1.0;

      switch(cur_avg->out_mode) {
      case T_STAT:
         (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
         (void) miattputstr(mincid, varid, MIvarid, t_stat_str);
         (void) miattputstr(mincid, varid, glim_par_str, param_str);
         (void) miattputint(mincid, varid, den_df_str, cur_avg->deg_free);
         if(cur_avg->stdev_mode == POOLED_SD) {
            (void) miattputstr(mincid, varid, pooled_file_str, 
                               input_files[glm_obj->control->pool_id]);
            (void) miattputdbl(mincid, varid, pooled_sd_str,
                              cur_avg->pooled_sd);
         }
         if (lambda_buffer != NULL) {
            (void) miattputdbl(mincid, varid, int_lambda_str,
                               lambda_buffer->fwhm_avg->int_lambda);
            (void) miattputdbl(mincid, varid, pooled_fwhm_str,
                                  lambda_buffer->fwhm_avg->fwhm_pool);
            (void) miattputdbl(mincid, varid, int_fwhm_str,
                               lambda_buffer->fwhm_avg->fwhm_int);
            if(glm_obj->search_mask != NULL) {
               (void) miattputstr(mincid, varid, search_file_str,
                                  glm_obj->search_mask);
            }
            (void) miattputdbl(mincid, varid, search_vol_str,
                               lambda_buffer->fwhm_avg->search_vol);
            (void) miattputstr(mincid, varid, fwhm_matrix_str,
                               fwhm_matrix_avg_str);
         }
         (void) miattputstr(mincid, varid, model_str, average_str);
         (void) miattputstr(mincid, varid, family_str,
                            glm_obj->family_names[glm_obj->family]);
         (void) miattputstr(mincid, varid, link_str,
                            glm_obj->link_names[glm_obj->link]);
         (void) miattputstr(mincid, varid, variance_fun_str,
                            glm_obj->variance_function_names[glm_obj->variance_function]);
         break;

      case F_STAT:
         (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
         (void) miattputstr(mincid, varid, MIvarid, f_stat_str);
         (void) miattputstr(mincid, varid, glim_par_str, param_str);
         (void) miattputint(mincid, varid, den_df_str, cur_avg->deg_free);
         (void) miattputint(mincid, varid, num_df_str, 1);

         if(cur_avg->stdev_mode == POOLED_SD) {
            (void) miattputstr(mincid, varid, pooled_file_str, 
                               input_files[glm_obj->control->pool_id]);
            (void) miattputdbl(mincid, varid, pooled_dev_str,
                              cur_avg->pooled_dev);
         }
         if (lambda_buffer != NULL) {
            (void) miattputdbl(mincid, varid, int_lambda_str,
                               lambda_buffer->fwhm_avg->int_lambda);
            (void) miattputdbl(mincid, varid, pooled_fwhm_str,
                               lambda_buffer->fwhm_avg->fwhm_pool);
            (void) miattputdbl(mincid, varid, int_fwhm_str,
                               lambda_buffer->fwhm_avg->fwhm_int);
            if(glm_obj->search_mask != NULL) {
               (void) miattputstr(mincid, varid, search_file_str,
                                  glm_obj->search_mask);
            }
            (void) miattputdbl(mincid, varid, search_vol_str,
                               lambda_buffer->fwhm_avg->search_vol);
            (void) miattputstr(mincid, varid, fwhm_matrix_str,
                               fwhm_matrix_avg_str);
            fprintf(stderr,"\n%s", fwhm_matrix_avg_str);
         }
         (void) miattputstr(mincid, varid, model_str, average_str);
         (void) miattputstr(mincid, varid, family_str,
                            glm_obj->family_names[glm_obj->family]);
         (void) miattputstr(mincid, varid, link_str,
                            glm_obj->link_names[glm_obj->link]);
         (void) miattputstr(mincid, varid, variance_fun_str,
                            glm_obj->variance_function_names[glm_obj->variance_function]);
         break;

      case P_CORR:
         (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
         (void) miattputstr(mincid, varid, MIvarid, p_corr_str);
         (void) miattputstr(mincid, varid, glim_par_str, param_str);
         (void) miattputint(mincid, varid, den_df_str, cur_avg->deg_free);
         if(cur_avg->stdev_mode == POOLED_SD) {
            (void) miattputstr(mincid, varid, pooled_file_str, 
                               input_files[glm_obj->control->pool_id]);
            (void) miattputdbl(mincid, varid, pooled_sd_str,
                              cur_avg->pooled_sd);
         }
         if (lambda_buffer != NULL) {
            (void) miattputdbl(mincid, varid, int_lambda_str,
                               lambda_buffer->fwhm_avg->int_lambda);
            (void) miattputdbl(mincid, varid, pooled_fwhm_str,
                               lambda_buffer->fwhm_avg->fwhm_pool);
            (void) miattputdbl(mincid, varid, int_fwhm_str,
                               lambda_buffer->fwhm_avg->fwhm_int);
            if(glm_obj->search_mask != NULL) {
               (void) miattputstr(mincid, varid, search_file_str,
                                  glm_obj->search_mask);
            }
            (void) miattputdbl(mincid, varid, search_vol_str,
                               lambda_buffer->fwhm_avg->search_vol);
            (void) miattputstr(mincid, varid, fwhm_matrix_str,
                               fwhm_matrix_avg_str);
            fprintf(stderr,"\n%s", fwhm_matrix_avg_str);
         }
         (void) miattputstr(mincid, varid, model_str, average_str);
         (void) miattputstr(mincid, varid, family_str,
                            glm_obj->family_names[glm_obj->family]);
         (void) miattputstr(mincid, varid, link_str,
                            glm_obj->link_names[glm_obj->link]);
         (void) miattputstr(mincid, varid, variance_fun_str,
                            glm_obj->variance_function_names[glm_obj->variance_function]);

         break;

      case BETA_HAT:
         (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
         (void) miattputstr(mincid, varid, MIvarid, beta_hat_str);
         (void) miattputstr(mincid, varid, glim_par_str, param_str);
         (void) miattputstr(mincid, varid, model_str, average_str);
         (void) miattputstr(mincid, varid, model_str, average_str);
         (void) miattputstr(mincid, varid, family_str,
                            glm_obj->family_names[glm_obj->family]);
         (void) miattputstr(mincid, varid, link_str,
                            glm_obj->link_names[glm_obj->link]);
         (void) miattputstr(mincid, varid, variance_fun_str,
                            glm_obj->variance_function_names[glm_obj->variance_function]);

         break;

      case STDEV_BETA:
         (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
         (void) miattputstr(mincid, varid, MIvarid, stdev_str);
         (void) miattputstr(mincid, varid, glim_par_str, param_str);
         (void) miattputint(mincid, varid, den_df_str, glm_obj->deg_free);
         (void) miattputstr(mincid, varid, model_str, average_str);
         (void) miattputstr(mincid, varid, family_str,
                            glm_obj->family_names[glm_obj->family]);
         (void) miattputstr(mincid, varid, link_str,
                            glm_obj->link_names[glm_obj->link]);
         (void) miattputstr(mincid, varid, variance_fun_str,
                            glm_obj->variance_function_names[glm_obj->variance_function]);

         break;

      default:
         break;
      }

      ncendef(mincid);

      miclose(mincid);

      fprintf(stderr," Done.\n");
   }

   if (glm_obj->control->dev_id > 0) {

      outfile = glm_obj->deviance_file;
      fprintf(stderr,"\nUpdating %s's header.",outfile);
      mincid = miopen(outfile, NC_WRITE);

      (void) ncredef(mincid);

      varid = ncvardef(mincid, var_name, NC_LONG, 0, 0);

      (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
      (void) miattputstr(mincid, varid, MIvarid,
                         "MNI deviance variable");
      (void) miattputstr(mincid, varid, glim_par_str, param_str);
      (void) miattputstr(mincid, varid, model_str, glm_str);
      (void) miattputstr(mincid, varid, family_str,
                         glm_obj->family_names[glm_obj->family]);
      (void) miattputstr(mincid, varid, link_str,
                         glm_obj->link_names[glm_obj->link]);
      (void) miattputstr(mincid, varid, variance_fun_str,
                         glm_obj->variance_function_names[glm_obj->variance_function]);

      ncendef(mincid);

      miclose(mincid);

      fprintf(stderr," Done.\n");
   }

   if (glm_obj->control->scale_id > 0) {

      outfile = glm_obj->scale_file;
      fprintf(stderr,"\nUpdating %s's header.",outfile);
      mincid = miopen(outfile, NC_WRITE);

      (void) ncredef(mincid);

      varid = ncvardef(mincid, var_name, NC_LONG, 0, 0);

      (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
      (void) miattputstr(mincid, varid, MIvarid,
                         "MNI scale parameter variable");
      (void) miattputstr(mincid, varid, glim_par_str, param_str);
      (void) miattputstr(mincid, varid, model_str, glm_str);
      (void) miattputstr(mincid, varid, family_str,
                         glm_obj->family_names[glm_obj->family]);
      (void) miattputstr(mincid, varid, link_str,
                         glm_obj->link_names[glm_obj->link]);
      (void) miattputstr(mincid, varid, variance_fun_str,
                         glm_obj->variance_function_names[glm_obj->variance_function]);

      ncendef(mincid);

      miclose(mincid);

      fprintf(stderr," Done.\n");
   }

   if ((glm_obj->control->avg_dev_id > 0) && 
       (glm_obj->contrast->avg_array != NULL)){

      outfile = glm_obj->contrast->avg_array->deviance_file;
      fprintf(stderr,"\nUpdating %s's header.",outfile);
      mincid = miopen(outfile, NC_WRITE);

      (void) ncredef(mincid);

      varid = ncvardef(mincid, var_name, NC_LONG, 0, 0);

      (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
      (void) miattputstr(mincid, varid, MIvarid,
                         "MNI deviance variable");
      (void) miattputstr(mincid, varid, glim_par_str, param_str);
      (void) miattputstr(mincid, varid, model_str, average_str);
      (void) miattputstr(mincid, varid, family_str,
                         glm_obj->family_names[glm_obj->family]);
      (void) miattputstr(mincid, varid, link_str,
                         glm_obj->link_names[glm_obj->link]);
      (void) miattputstr(mincid, varid, variance_fun_str,
                         glm_obj->variance_function_names[glm_obj->variance_function]);

      ncendef(mincid);

      miclose(mincid);

      fprintf(stderr," Done.\n");
   }

   if (lambda_buffer != NULL) {
      if(glm_obj->family != GAUSSIAN) {
         for(itest=0; itest<lambda_buffer->num_test; itest++) {
            fwhm_general = lambda_buffer->fwhm_general[itest];
            icont = glm_obj->contrast->test_map[itest];
            cur_matrix = glm_obj->contrast->contrast_matrix_array[itest];
            outfile = fwhm_general->outfile;
            fprintf(stderr,"\nUpdating %s's header.",outfile);
            mincid = miopen(outfile, NC_WRITE);

            (void) ncredef(mincid);

            varid = ncvardef(mincid, var_name, NC_LONG, 0, 0);

            (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
            (void) miattputstr(mincid, varid, MIvarid,
                         "MNI local FWHM variable");
            (void) miattputstr(mincid, varid, glim_par_str, param_str);
            (void) miattputstr(mincid, varid, model_str, glm_str);
            if (cur_matrix->contrast_matrix == NULL) {
               (void) miattputint(mincid, varid, column_str, 
                                  cur_matrix->column_num+1);
            }
            else {
               sprint_matrix(&contr_str, cur_matrix->contrast_matrix);
               (void) miattputstr(mincid, varid, contrast_str, 
                                  contr_str); 
               GI_FREE(contr_str);
            }

            (void) miattputstr(mincid, varid, family_str,
                               glm_obj->family_names[glm_obj->family]);
            (void) miattputstr(mincid, varid, link_str,
                               glm_obj->link_names[glm_obj->link]);
            (void) miattputstr(mincid, varid, variance_fun_str,
                               glm_obj->variance_function_names[glm_obj->variance_function]);

            ncendef(mincid);

            miclose(mincid);

            fprintf(stderr," Done.\n");
          }
       }   
             
       if (lambda_buffer->fwhm_gaussian != NULL) {
          if (lambda_buffer->fwhm_gaussian->out_id > 0) {
             fwhm_general = lambda_buffer->fwhm_gaussian;
             outfile = fwhm_general->outfile;
             fprintf(stderr,"\nUpdating %s's header.",outfile);
             mincid = miopen(outfile, NC_WRITE);

             (void) ncredef(mincid);

             varid = ncvardef(mincid, var_name, NC_LONG, 0, 0);

             (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
             (void) miattputstr(mincid, varid, MIvarid,
                          "MNI local FWHM variable");
             (void) miattputstr(mincid, varid, glim_par_str, param_str);
             (void) miattputstr(mincid, varid, model_str, glm_str);
             (void) miattputstr(mincid, varid, family_str,
                                glm_obj->family_names[glm_obj->family]);
             (void) miattputstr(mincid, varid, link_str,
                                glm_obj->link_names[glm_obj->link]);
             (void) miattputstr(mincid, varid, variance_fun_str,
                                glm_obj->variance_function_names[glm_obj->variance_function]);

             ncendef(mincid);

             miclose(mincid);

             fprintf(stderr," Done.\n");
          }
       }

       if (lambda_buffer->fwhm_avg != NULL) {
          if (lambda_buffer->fwhm_avg->out_id > 0) {
             fwhm_general = lambda_buffer->fwhm_avg;
             outfile = fwhm_general->outfile;
             fprintf(stderr,"\nUpdating %s's header.",outfile);
             mincid = miopen(outfile, NC_WRITE);

             (void) ncredef(mincid);

             varid = ncvardef(mincid, var_name, NC_LONG, 0, 0);

             (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
             (void) miattputstr(mincid, varid, MIvarid,
                          "MNI local FWHM variable");
             (void) miattputstr(mincid, varid, glim_par_str, param_str);
             (void) miattputstr(mincid, varid, model_str, average_str);
             (void) miattputstr(mincid, varid, family_str,
                                glm_obj->family_names[glm_obj->family]);
             (void) miattputstr(mincid, varid, link_str,
                                glm_obj->link_names[glm_obj->link]);
             (void) miattputstr(mincid, varid, variance_fun_str,
                                glm_obj->variance_function_names[glm_obj->variance_function]);

             ncendef(mincid);

             miclose(mincid);

             fprintf(stderr," Done.\n");
          }
       }

       if (lambda_buffer->fwhm_simple != NULL) {
          if (lambda_buffer->fwhm_simple->out_id > 0) {
             fwhm_general = lambda_buffer->fwhm_simple;
             outfile = fwhm_general->outfile;
             fprintf(stderr,"\nUpdating %s's header.",outfile);
             mincid = miopen(outfile, NC_WRITE);

             (void) ncredef(mincid);

             varid = ncvardef(mincid, var_name, NC_LONG, 0, 0);

             (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
             (void) miattputstr(mincid, varid, MIvarid,
                          "MNI local FWHM variable");
             (void) miattputstr(mincid, varid, glim_par_str, param_str);
             (void) miattputstr(mincid, varid, model_str, glm_str);
             (void) miattputstr(mincid, varid, family_str,
                                glm_obj->family_names[glm_obj->family]);
             (void) miattputstr(mincid, varid, link_str,
                                glm_obj->link_names[glm_obj->link]);
             (void) miattputstr(mincid, varid, variance_fun_str,
                                glm_obj->variance_function_names[glm_obj->variance_function]);

             ncendef(mincid);

             miclose(mincid);

             fprintf(stderr," Done.\n");
          }
       }
   
   }

   if (glm_obj->control->err_id > 0) {

      outfile = glm_obj->error_file;
      fprintf(stderr,"\nUpdating %s's header.",outfile);
      mincid = miopen(outfile, NC_WRITE);

      (void) ncredef(mincid);

      varid = ncvardef(mincid, var_name, NC_LONG, 0, 0);

      (void) miattputstr(mincid, varid, MIvartype, MI_GROUP);
      (void) miattputstr(mincid, varid, MIvarid,
                         "MNI convergence codes");
      (void) miattputstr(mincid, varid, glim_par_str, param_str);
      (void) miattputstr(mincid, varid, family_str,
                         glm_obj->family_names[glm_obj->family]);
      (void) miattputstr(mincid, varid, link_str,
                         glm_obj->link_names[glm_obj->link]);
      (void) miattputstr(mincid, varid, variance_fun_str,
                         glm_obj->variance_function_names[glm_obj->variance_function]);
      ncendef(mincid);

      miclose(mincid);

      fprintf(stderr," Done.\n");

   }

   GI_FREE(param_str);

   if (fwhm_matrix_avg_str != NULL)
      GI_FREE(fwhm_matrix_avg_str);

   if (fwhm_matrix_gaussian_str != NULL)
      GI_FREE(fwhm_matrix_gaussian_str);

   return;

}
 

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_glim_parameters
@INPUT      : glm_obj - Glm_Object with pointers to contrast matrices and
                        design matrix
              input_files - list of data files used to fit GLM
@OUTPUT     : param_str - array of strings containing info to be used
                          to modify header of each output file, contains
                          info about design matrix
              contr_str - array of strings containing info to be used
                          to modify header of each output file, contains
                          info about contrast_matrices
@RETURNS    : (nothing)
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Oct. 14, 97 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void create_glim_parameters(Glm_Object *glm_obj, char **input_files, 
                            char **param_str)
{
   int i, j;
   long length_tot;
   char *tmpstring;
   char **tmpstring2;

   tmpstring = GI_MALLOC(sizeof(*tmpstring) * 10000);
   tmpstring2 = GI_MALLOC(sizeof(char *) * glm_obj->design_matrix->num_rows);
   length_tot = 0;

   for(i=0; i<glm_obj->design_matrix->num_rows; i++) {
      tmpstring2[i] = GI_MALLOC(sizeof(char) * (strlen(input_files[i]) + 50
                                    + 15*glm_obj->design_matrix->num_columns +
                                    300));
      sprintf(tmpstring,"Using file: %s. Regressors: ", input_files[i]);
      strcpy(tmpstring2[i], tmpstring);
      for(j=0; j<glm_obj->design_matrix->num_columns; j++) {
         sprintf(tmpstring, "%5.2e  ",glm_obj->design_matrix->values[i][j]);
         strcat(tmpstring2[i], tmpstring);
      }
      strcat(tmpstring2[i], "\n");
      length_tot += strlen(tmpstring2[i]);
   }
      
   GI_FREE(tmpstring);

   *param_str = GI_MALLOC(sizeof(**param_str) * (length_tot + 1));
   strcpy(*param_str, "\n");

   for(i=0; i<glm_obj->design_matrix->num_rows; i++) {
      strcat(*param_str, tmpstring2[i]);
      GI_FREE(tmpstring2[i]);
   }

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_fwhm_filenames
@INPUT      : lambda_buffer - structure containing info for calculating
                              smoothness of output random fields
              glm_obj - structure with info about model
@OUTPUT     : list of names for the output of smoothness of random field
@RETURNS    : 
@DESCRIPTION: routine to get filenames from the contrast filenames
@CREATED    : Jan. 22, 1997, J. Taylor
@MODIFIED   : 
---------------------------------------------------------------------------- */
void get_fwhm_filenames(Lambda *lambda_buffer, Glm_Object *glm_obj)
{
   int num_test, icont, itest;
   Contrast_Matrix *cur_matrix;
   char *fwhm_str = "_fwhm.mnc\0";
   Fwhm_Info *fwhm_general;

   num_test = glm_obj->contrast->num_test;
   
   itest = 0;

   if(glm_obj->family != GAUSSIAN) {

      for(itest=0; itest<num_test; itest++) {

         fwhm_general = lambda_buffer->fwhm_general[itest];

         icont = glm_obj->contrast->test_map[itest];
         cur_matrix = glm_obj->contrast->contrast_matrix_array[icont];

         /* Note: just setting up for t-statistics now... later to be fixed
            for f-stats and p-corr (Feb. 1998)*/

         if (glm_obj->control->scylla == FALSE) {
            fwhm_general->outfile = GI_MALLOC(sizeof(char) * (strlen(cur_matrix->outfile + 20)));

            strncpy(fwhm_general->outfile, cur_matrix->outfile, 
                    strlen(cur_matrix->outfile) - 4);
            strcat(fwhm_general->outfile, fwhm_str);
         }
         else {
            fwhm_general->outfile = tempnam(NULL, "glim");
            fprintf(tmpfile_list, "\n%s", fwhm_general->outfile);
            num_tmpfiles++;
         }

         fprintf(stderr,"\n%s \n", fwhm_general->outfile); 
      }
   }

   return;

}
      
/* ----------------------------- MNI Header -----------------------------------
@NAME       : update_lambda_buffer
@INPUT      : lambda_buffer - ptr to ring buffers where values for calculation
                              of derivatives are stored
@OUTPUT     : nothing
@RETURNS    : 
@DESCRIPTION: routine to advance buffers for calculating smoothness of field
@CREATED    : Nov. 3, 1997, J. Taylor
@MODIFIED   : 
---------------------------------------------------------------------------- */
int update_lambda_buffer(Lambda *lambda_buffer, Glm_Object *glm_obj, 
                         double check)
{
   int ibuff, i, j, k, h;
   int itest, icont;
   double norm;
   double tmp_sum = 0.0;
   Contrast_Matrix *cur_contrast;
   Fwhm_Info *fwhm_general;

   /* Advance index on ring */

   lambda_buffer->current += 1;
   if (lambda_buffer->current == lambda_buffer->ring_size)
      lambda_buffer->current = 0;

   lambda_buffer->total++;

   /* Update coords */

   lambda_buffer->coord[lambda_buffer->num_dim-1] += 1;

   for(i=lambda_buffer->num_dim-1; i>0; i--) {
      if(lambda_buffer->coord[i] == lambda_buffer->sizes[i]) {
         lambda_buffer->coord[i] = 0;
         lambda_buffer->coord[i-1] += 1;
      }
   }
            
   /* Check to see if we have just finished last slice */

   if (lambda_buffer->coord[0] == lambda_buffer->sizes[0])
      return FALSE;

   /* Calculate the values of weight term to use in correction for
      changing weights */
               
   if ((glm_obj->control->status == CONVERGE) &&
       (glm_obj->family != GAUSSIAN)){ 

      for(itest = 0; itest<lambda_buffer->num_test; itest++) {
         icont = glm_obj->contrast->test_map[itest];
         cur_contrast = glm_obj->contrast->contrast_matrix_array[icont];
         
         if(cur_contrast->out_mode == T_STAT) {
            if(cur_contrast->contrast_matrix != NULL) {
               norm = sqrt(glm_obj->tmp_con[i][2]->values[0][0]);
               for(j=0; j<glm_obj->response->num_rows; j++) {
                  cur_contrast->deriv->values[j][0] = 0.0;
                  for(h=0; h<glm_obj->design_matrix->num_columns; h++) {
                     for(k=0; k<glm_obj->design_matrix->num_columns; k++){
                        tmp_sum += (glm_obj->design_matrix->values[j][k] *
                                    glm_obj->info_inv->values[h][k]);
                     }
                     cur_contrast->deriv->values[j][0] += (tmp_sum * cur_contrast->contrast_matrix->values[0][h]);
                     tmp_sum=0.0;
                  }
                  cur_contrast->deriv->values[j][0] *= \
                     (glm_obj->sqrt_weight->values[j][0] *
                      norm);
                  
                  /* Check with Keith to see if this is right */
                        
                  /*          if(cur_contrast->stdev_mode == VOXEL_SD)
                              cur_contrast->deriv->values[j][0]/=\
                              sqrt(glm_obj->scale_est);  */
               }
            }
            else {
               h = cur_contrast->column_num;
               norm = sqrt(glm_obj->info_inv->values[h][h]);
               for(j=0; j<glm_obj->response->num_rows; j++) {
                  cur_contrast->deriv->values[j][0] = 0.0;
                  for(k=0; k<glm_obj->design_matrix->num_columns; k++) {
                     cur_contrast->deriv->values[j][0] += (glm_obj->design_matrix->values[j][k] * glm_obj->info_inv->values[h][k]);
                  }
                  cur_contrast->deriv->values[j][0] *= \
                     (glm_obj->sqrt_weight->values[j][0] /
                      norm);
               }
            }
         }
      }
   }

   /* Update values in ring buffers, overwriting the "last" value in the ring 
      with the current value which we wish to use to calculate lambda */

   for (ibuff=0; ibuff<lambda_buffer->num_buffers; ibuff++) {
         
      if(lambda_buffer->is_avg == TRUE ) {
         lambda_buffer->fwhm_avg->data[ibuff][lambda_buffer->current] =\
            glm_obj->avg_resid->values[ibuff][0];
      }

      if ((glm_obj->contrast->num_test > 0) &&
          (glm_obj->family == GAUSSIAN)) {
         lambda_buffer->fwhm_gaussian->data[ibuff][lambda_buffer->current]=\
            glm_obj->resid->values[ibuff][0];
      }
         
      if(glm_obj->family != GAUSSIAN) {
         for(itest=0; itest<lambda_buffer->num_test; itest++) {
            fwhm_general = lambda_buffer->fwhm_general[itest];
            icont = glm_obj->contrast->test_map[itest];
            cur_contrast = glm_obj->contrast->contrast_matrix_array[icont];

            if ((glm_obj->control->status == CONVERGE) ||
                (glm_obj->control->status == CONT)) { 
               fwhm_general->data[ibuff][lambda_buffer->current] =\
                  (cur_contrast->deriv->values[ibuff][0] *
                   glm_obj->resid->values[ibuff][0]);
            }
            else {
               fwhm_general->data[ibuff][lambda_buffer->current] =\
                  INVALID_DATA;
            }
         }
         if (lambda_buffer->fwhm_simple != NULL) {
            if (lambda_buffer->fwhm_simple->outfile != NULL) {
               fwhm_general = lambda_buffer->fwhm_simple;
               fwhm_general->data[ibuff][lambda_buffer->current] =\
                  glm_obj->resid->values[ibuff][0];
            }
         }
      }
   }     

   lambda_buffer->test[lambda_buffer->current]\
      = check;

   return TRUE;

}
