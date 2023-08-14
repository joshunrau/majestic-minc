/* ----------------------------- MNI Header -----------------------------------
@NAME       : glim.c
@INPUT      : 
@OUTPUT     : (nothing)
@RETURNS    : 
@DESCRIPTION: Source file for generalized linear models
@CREATED    : Sept 11, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_CTYPE_H
#include <ctype.h>
#endif

#include "glim.h"
#include "volume_io.h"

extern FILE *tmpfile_list;
extern char *tmpfile_name;
extern int num_tmpfiles;
extern int verbose;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_tmpfiles
@INPUT      : 
@OUTPUT     : 
@RETURNS    : TRUE since nextarg is used.
@DESCRIPTION: deletes tmpfiles created by glim_image
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 18, 1998 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void delete_tmpfiles(FILE **tmpfile_list)
{
   char tmpstring[2056];
   FILE *tmp;
   int is_null;

   fclose(*tmpfile_list);

   if (num_tmpfiles > 0) {
      *tmpfile_list = fopen(tmpfile_name, "r");

      while(fscanf(*tmpfile_list, "%s", &tmpstring) == 1) {
         is_null = TRUE;
         tmp = fopen(tmpstring,"r");
         if (tmp != NULL) {
            is_null = FALSE;
            fclose(tmp);
         }
         if (is_null == FALSE)
            remove(tmpstring);
      }
      
      fclose(*tmpfile_list);
   }

   remove(tmpfile_name);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_contrast_from_raw
@INPUT      : Array of raw contrasts, pointer to final contrast matrices
@OUTPUT     : (none)
@RETURNS    : TRUE since nextarg is used.
@DESCRIPTION: Parses contrasts entered on command line to contrast matrices
@METHOD     : 
@GLOBALS    : 
@CALLS      : get_contrast_from_matrix
@CREATED    : June 21, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int get_contrast_from_raw( Contrast_Raw_Array *contrast_list,
                          Contrast_Mat_Array **contrast_matrices)
{
   Contrast_Raw *cur_raw;
   Contrast_Matrix *cur_matrix;
   Avg_Info *cur_avg;
   int ifile, icont, iavg, tmp_len;
   int itest;
   char *tmp_name;
   char *prefix = "glim";
   char *tmp_dir = NULL;

   iavg = 0;
   icont = 0;

   (*contrast_matrices)->num_avg = contrast_list->num_avg;
   (*contrast_matrices)->num_contrasts = (contrast_list->num_output - 
                                          contrast_list->num_avg);
   (*contrast_matrices)->num_test = 0;
   (*contrast_matrices)->avg_array->num_test = 0;

   if(contrast_list->num_avg > 0) {
      (*contrast_matrices)->avg_array->avg_info = GI_MALLOC(sizeof(Avg_Info*) *
                                                      contrast_list->num_avg);
   }
   else
      (*contrast_matrices)->avg_array = NULL;

   for(ifile=0 ; ifile<contrast_list->num_output ; ifile++) {

      /* Get pointer to contrast array */

      cur_raw = contrast_list->contrast_array[ifile];
   
      if(cur_raw->is_avg == FALSE) {
         cur_matrix = GI_MALLOC(sizeof(Contrast_Matrix));

         cur_matrix->outfile = cur_raw->outfile_name;

         cur_matrix->out_id = ifile;
         
         if(cur_raw->out_type == NULL) {
            cur_matrix->out_mode = F_STAT;
            (*contrast_matrices)->num_test++;
         }
         else if (strcmp(cur_raw->out_type,"stdev") == 0)
            cur_matrix->out_mode = STDEV_BETA;
         else if(strcmp(cur_raw->out_type,"beta") == 0)
            cur_matrix->out_mode = BETA_HAT;
         else if(strcmp(cur_raw->out_type,"t_stat") == 0) {
            cur_matrix->out_mode = T_STAT;
            (*contrast_matrices)->num_test++;
         }
         else if(strcmp(cur_raw->out_type,"corr") == 0) {
            cur_matrix->out_mode = P_CORR;
            (*contrast_matrices)->num_test++;
         }
         else { 
            fprintf(stderr,"\nError in get_contrast_from_raw.\n");
            fprintf(stderr,"\n%s is not an output type.\n",cur_raw->out_type);
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }

         GI_FREE(cur_raw->out_type);

         if(strcmp(cur_raw->stdev_type,"voxel") == 0)
            cur_matrix->stdev_mode = VOXEL_SD;
         else if(strcmp(cur_raw->stdev_type,"pool") == 0)
            cur_matrix->stdev_mode = POOLED_SD;
         else { 
            fprintf(stderr,"\nError in get_contrast_from_raw.\n");
            fprintf(stderr,"\n%s is not an stdev type.\n", 
                    cur_raw->stdev_type);
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }

         if ((cur_matrix->stdev_mode == POOLED_SD) && 
             (cur_matrix->out_mode == BETA_HAT)) {
            fprintf(stderr,"\nWarning: using pooled stdev for output of beta_hat for file %s.\nBeta_hat output option can't use pooled stdev. Pool option ignored.\n", cur_matrix->outfile);
            cur_matrix->stdev_mode = VOXEL_SD;
         }

         GI_FREE(cur_raw->stdev_type);

         if ((cur_matrix->stdev_mode == POOLED_SD) && 
             (cur_matrix->out_mode == STDEV_BETA)) {
            fprintf(stderr,"\nWarning: using pooled stdev for output of stdev for file %s.\nStdev output option can't use pooled stdev. Pool option ignored.\n", cur_matrix->outfile);
            cur_matrix->stdev_mode = VOXEL_SD;
         }

         if (cur_matrix->stdev_mode == VOXEL_SD)
            cur_matrix->tmpfile = NULL;
         else if (cur_matrix->stdev_mode == POOLED_SD) {
            tmp_name = tempnam(tmp_dir, prefix);
            fprintf(tmpfile_list, "%s\n", tmp_name);
            num_tmpfiles++;
            tmp_len = strlen(tmp_name);
            cur_matrix->tmpfile = GI_MALLOC(sizeof(char) * (tmp_len + 1));
            if (strcpy(cur_matrix->tmpfile, tmp_name) == NULL) {
               fprintf(stderr,"\nError in get_contrast_from_raw.\n");
               fprintf(stderr,"Can't copy tmp_name to outfile.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
         }

         if(strcmp(cur_raw->in_type,"matrix") == 0) {
            if(get_contrast_from_matrix(cur_raw->raw_contrast ,
                                        (*contrast_matrices)->num_columns,
                                        &cur_matrix) != TRUE) {
               fprintf(stderr,"\nError getting contrast matrix for outfile %s, type %s.\n", cur_matrix->outfile, cur_raw->in_type);
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
         }
         else if(strcmp(cur_raw->in_type,"file") == 0) {
            if(get_contrast_from_file(cur_raw->raw_contrast, 
                                      (*contrast_matrices)->num_columns, 
                                      &cur_matrix) != TRUE)  {
               fprintf(stderr,"\nError getting contrast matrix for outfile %s, type %s.\n", cur_matrix->outfile, cur_raw->in_type);
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
         }
         else if(strcmp(cur_raw->in_type,"column") == 0) {
            if(get_contrast_from_column(cur_raw->raw_contrast,
                                        (*contrast_matrices)->num_columns,
                                        &cur_matrix) != TRUE) {
               fprintf(stderr,"\nError getting contrast matrix for outfile %s, type %s.\n", cur_matrix->outfile, cur_raw->in_type);
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
         }  
         else {
            fprintf(stderr,"\nError in get_contrast_from_raw.");
            fprintf(stderr,"\n%s isn't a proper input method.\n",
                    cur_raw->in_type);
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }

         (*contrast_matrices)->contrast_matrix_array[icont] = cur_matrix;

         icont++;
         cur_matrix++; 

         GI_FREE(cur_raw->in_type);

      }
      else {
         cur_avg = GI_MALLOC(sizeof(Avg_Info));

         cur_avg->outfile = cur_raw->outfile_name;

         cur_avg->out_id = ifile;

         if(cur_raw->out_type == NULL) {
            cur_avg->out_mode = F_STAT;
            (*contrast_matrices)->avg_array->num_test++;
         }
         else if (strcmp(cur_raw->out_type,"stdev") == 0)
            cur_avg->out_mode = STDEV_BETA;
         else if(strcmp(cur_raw->out_type,"beta") == 0)
            cur_avg->out_mode = BETA_HAT;
         else if(strcmp(cur_raw->out_type,"t_stat") == 0) {
            cur_avg->out_mode = T_STAT;
            (*contrast_matrices)->avg_array->num_test++;
         }
         else if(strcmp(cur_raw->out_type,"corr") == 0) {
            cur_avg->out_mode = P_CORR;
            (*contrast_matrices)->avg_array->num_test++;
         }
         else { 
            fprintf(stderr,"\nError in get_contrast_from_raw.\n");
            fprintf(stderr,"\n%s is not an output type.\n",cur_raw->out_type);
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }

         GI_FREE(cur_raw->out_type);

         if(strcmp(cur_raw->stdev_type,"voxel") == 0)
            cur_avg->stdev_mode = VOXEL_SD;
         else if(strcmp(cur_raw->stdev_type,"pool") == 0)
            cur_avg->stdev_mode = POOLED_SD;
         else { 
            fprintf(stderr,"\nError in get_contrast_from_raw.\n");
            fprintf(stderr,"\n%s is not an stdev type.\n",
                    cur_raw->stdev_type);
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
         }

         if ((cur_avg->stdev_mode == POOLED_SD) && 
             (cur_avg->out_mode == BETA_HAT)) {
            fprintf(stderr,"\nWarning: using pooled stdev for output of beta_hat for file %s.\nBeta_hat output option can't use pooled stdev. Pool option ignored.\n", cur_avg->outfile);
            cur_avg->stdev_mode = VOXEL_SD;
         }

         GI_FREE(cur_raw->stdev_type);

         if ((cur_avg->stdev_mode == POOLED_SD) && 
             (cur_avg->out_mode == STDEV_BETA)) {
            fprintf(stderr,"\nWarning: using pooled stdev for output of stdev for file %s.\nStdev output option can't use pooled stdev. Pool option ignored.\n", cur_avg->outfile);
            cur_avg->stdev_mode = VOXEL_SD;
         }

         if (cur_avg->stdev_mode == VOXEL_SD)
            cur_avg->tmpfile = NULL;
         else if (cur_avg->stdev_mode == POOLED_SD) {
            tmp_name = tempnam(tmp_dir, prefix);
            fprintf(tmpfile_list, "%s\n", tmp_name);
            num_tmpfiles++;
            tmp_len = strlen(tmp_name);
            cur_avg->tmpfile = GI_MALLOC(sizeof(char) * (tmp_len + 1));
            if (strcpy(cur_avg->tmpfile, tmp_name) == NULL) {
               fprintf(stderr,"\nError in get_contrast_from_raw.\n");
               fprintf(stderr,"Can't copy tmp_name to outfile.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
         }
         (*contrast_matrices)->avg_array->avg_info[iavg] = cur_avg;
         
         iavg++;
         cur_avg++;

      }

   }

   /* Get the mapping of all the contrasts that are test-statistics */
      
   (*contrast_matrices)->test_map = GI_MALLOC((*contrast_matrices)->num_test
                                           * sizeof(int));
   itest = 0;

   for(icont=0 ; icont<(*contrast_matrices)->num_contrasts ; icont++) {
      cur_matrix = (*contrast_matrices)->contrast_matrix_array[icont];

      if ((cur_matrix->out_mode == T_STAT) || 
          (cur_matrix->out_mode == F_STAT)){
         (*contrast_matrices)->test_map[itest] = icont;
         itest++;
      }
   }

   return TRUE;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_contrast_from_matrix
@INPUT      : raw contrast string,
              num_columns 
              pointer to current contrast matrix in get_contrast_from_raw
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: gets a contrast matrix entered in the form "matrix"
              (i.e. '0 0 0 1 0 ; 0 1 0 0 0' etc.)
@METHOD     : 
@GLOBALS    : 
@CALLS      :  get_double_list, create_matrix
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int get_contrast_from_matrix(char *contrast_string, int num_columns,
                             Contrast_Matrix **cur_matrix)
{
   int num_rows;
   Double_Array *converted_values = NULL;
   double value;
   int i,j;
   char *cur, *end;

   num_rows = 1;

   /* Get number of rows of current contrast */

   cur = contrast_string;
   end = contrast_string + strlen(contrast_string);

   while(cur != end) {
      if(*cur == ROW_SEPARATOR)  
         num_rows++;
      cur++;
   }
   
   /* Initialize current matrix */

   (*cur_matrix)->contrast_matrix = create_matrix(num_rows, num_columns);

   get_double_list(contrast_string, &converted_values);
   if (converted_values->num_values != (num_columns * num_rows)) {
         fprintf(stderr,"\nMismatched number of elements for contrast %s \n", (*cur_matrix)->outfile);
         return FALSE;
      }

   for(i=0; i<num_rows; i++) {
      for(j=0 ; j<num_columns; j++) {
         value = converted_values->values[i*num_columns+j];
         (*cur_matrix)->contrast_matrix->values[i][j] = value;
      }
   }
         
   return TRUE;
}
      
/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_contrast_from_column
@INPUT      : raw contrast string,
              num_columns 
              pointer to current contrast matrix in get_contrast_from_raw
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: gets a contrast matrix entered in the form "column"
              where the contrast outputs a t-statistic testing the hypothesis
              glm_obj->beta_hat->values[column-1][0] == 0
@METHOD     : 
@GLOBALS    : 
@CALLS      :   create_matrix
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int get_contrast_from_column(char *contrast_column, int num_columns, 
                           Contrast_Matrix **cur_matrix)
{
   int column;
   char **ptr = NULL;

   column = strtol(contrast_column, ptr, 10);

   if(column > num_columns) {
      fprintf(stderr,"\nError in get_contrast_from_column: column number too large.\n");

printf( "column = %d  num_columns = %d\n", column, num_columns );
      return FALSE;
   }

   (*cur_matrix)->contrast_matrix = NULL;

   (*cur_matrix)->column_num = column - 1;

   return TRUE;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_contrast_from_file
@INPUT      : raw contrast string,
              num_columns 
              pointer to current contrast matrix in get_contrast_from_raw
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: gets a contrast matrix entered in the form "file"
              a file containing a contrast matrix with the appropriate number
              of rows
@METHOD     : 
@GLOBALS    : 
@CALLS      :   create_matrix
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int get_contrast_from_file(char *contrast_filename, int num_columns,
                           Contrast_Matrix **cur_matrix)
{
   FILE *contrast_file;
   int num_rows;
   int i,j,k;
   int initial_rows = 10;
   double **values;
   double value;

   contrast_file = fopen(contrast_filename, "r");
   if (contrast_file == NULL) {
      fprintf(stderr,"\nError: can't find file %s.\n", contrast_filename);
      return FALSE;
   }

   values = GI_MALLOC(sizeof(double *) * initial_rows);

   for(i=0; i<initial_rows; i++) 
      values[i] = GI_MALLOC(sizeof(double) * num_columns);
      
   i = 0;
   num_rows = 0;

   while(fscanf(contrast_file, "%lf", &(values[i][0])) == 1) {

      j=1;
      while(fscanf(contrast_file, "%lf", &(values[i][j])) == 1) {
         j++;
      }

      if(j != num_columns) {
         fprintf(stderr,"\nError: mismatched number of entries in line %d of %s.\n", i, contrast_filename);
         return FALSE;
      }
      i++;
      num_rows++;
      
      if(num_rows == (initial_rows - 1)) {
         initial_rows = initial_rows + 50;
         GI_REALLOC(values, sizeof(double *) * initial_rows);
         for(k=num_rows; k<initial_rows; k++) {
            values[k] = GI_MALLOC(sizeof(double) * num_columns); 
         }
      }
   }

   (*cur_matrix)->contrast_matrix = create_matrix(num_rows,num_columns);

   for(i=0; i<num_rows; i++) {
      for(j=0; j<num_columns; j++) {
         value = values[i][j];
         (*cur_matrix)->contrast_matrix->values[i][j] = value;
      }
   }

   return TRUE;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_double_list
@INPUT      : double_string - string to be converted
              double_array - destination of string
@OUTPUT     : (none)
@RETURNS    : TRUE since nextarg is used.
@DESCRIPTION: Gets a list (array) of double values.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : March 8, 1995 (Peter Neelin)
@MODIFIED   : June 22, 1997 (J. Taylor) changed function so it could take
              input other than from ParseArgv
---------------------------------------------------------------------------- */
int get_double_list(char *double_string, Double_Array **double_array)
{
   int num_elements;
   int num_alloc;
   double *double_list = NULL;
   double dvalue;
   char *cur, *end, *prev;

   /* Set up pointers to end of string and first non-space character */

   end = double_string + strlen(double_string);
   cur = double_string;
   while (isspace(*cur)) cur++;
   num_elements = 0;
   num_alloc = 0;
   double_list = NULL;

   /* Loop through string looking for doubles */
   while (cur!=end) {

      /* Get double */
      prev = cur;
      dvalue = strtod(prev, &cur);

      if(prev==cur) {
         (void) fprintf(stderr, 
            "Error: expected vector of doubles for but got \"%s\"\n", 
                        double_string);
         return FALSE;
      }

      /* Add the value to the list */
      num_elements++;
      if (num_elements > num_alloc) {
         num_alloc += 20;
         if (double_list == NULL) {
            double_list = 
               GI_MALLOC(num_alloc * sizeof(*double_list));
         }
         else {
            double_list = 
               GI_REALLOC(double_list, num_alloc * sizeof(*double_list));
         }
      }
      double_list[num_elements-1] = dvalue;

      /* Skip any spaces */
      while (isspace(*cur)) cur++;

      /* Skip row separator and an optional comma */
      if (*cur == ENTRY_SEPARATOR) cur++;
      if (*cur == ROW_SEPARATOR) cur++;
   }

   if((*double_array) == NULL) {
      (*double_array) = GI_MALLOC(sizeof(Double_Array *));
      (*double_array)->values = GI_MALLOC(sizeof(double) * num_elements);
   }

   /* Update the global variables */
   (*double_array)->num_values = num_elements;
   if ((*double_array)->values != NULL) {
      GI_FREE((*double_array)->values);
   }
   (*double_array)->values = double_list;

   return TRUE;
}
      
/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_glm_matrices
@INPUT      : pointer to a Glm_Object
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: deletes matrices of glm_obj
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void delete_glm_matrices(Glm_Object *glm_obj, int is_last) 
{
   int i,j;
   Contrast_Matrix *cur_contrast;

   if(glm_obj->control->do_fit == TRUE) {
      if(glm_obj->family == GAUSSIAN) {
         delete_matrix(glm_obj->info_inv);
         delete_matrix(glm_obj->lm_matrix);
         delete_matrix(glm_obj->response);
         delete_matrix(glm_obj->mu);
         delete_matrix(glm_obj->beta_hat);
         delete_matrix(glm_obj->output);
         delete_matrix(glm_obj->missing);
      }
      else if (glm_obj->family != GAUSSIAN) {
         delete_matrix(glm_obj->info_inv);
         delete_matrix(glm_obj->beta_hat);
         delete_matrix(glm_obj->lm_matrix);
         delete_matrix(glm_obj->mu);
         delete_matrix(glm_obj->response);
         delete_matrix(glm_obj->z_i); 
         delete_matrix(glm_obj->x_beta);
         delete_matrix(glm_obj->output);
         delete_matrix(glm_obj->g_prime);
         delete_matrix(glm_obj->info_mat);
         delete_matrix(glm_obj->sqrt_weight);
         delete_matrix(glm_obj->missing);

         if (glm_obj->corr_struct == DIAGONAL_CORR) {
            delete_matrix(glm_obj->weight); 
         }
         else {
            delete_matrix(glm_obj->sigma);
         }
      }

      if (is_last == TRUE) {
         for(j=0; j<glm_obj->contrast->num_contrasts; j++) {
            cur_contrast = glm_obj->contrast->contrast_matrix_array[j];
            if(cur_contrast->contrast_matrix != NULL) {
               for(i=0; i<6; i++) {
                  delete_matrix(glm_obj->tmp_con[j][i]);
               }
               delete_matrix(cur_contrast->contrast_matrix);
            }
         }
      }
   }

   return; 
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_glm
@INPUT      : pointer to Glm_Object
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: allocates memory for all matrices needed for fitting of glm and 
              contrasts
@METHOD     : 
@GLOBALS    : uses family_input, link_input, variance_function_input
@CALLS      : create_contrast
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void initialize_glm(Glm_Object *glm_obj, Initial_Data *initial_data)
{
   int num_rows, num_columns;
   
   glm_obj->pooled_dev = 0.0;
   glm_obj->avg_pooled_dev = 0.0;
   glm_obj->num_pool = 0;
   glm_obj->avg_num_pool = 0;

   glm_obj->design_matrix = initial_data->design_matrix;
   glm_obj->control = initial_data->control;

   glm_obj->scale_est = initial_data->scale_est;

   num_rows = glm_obj->design_matrix->num_rows;
   num_columns = glm_obj->design_matrix->num_columns;

   glm_obj->deg_free = 1.0 * (num_rows - num_columns); 
   glm_obj->num_var = num_columns;
   glm_obj->response = create_matrix(num_rows, 1);
   glm_obj->resid = create_matrix(num_rows, 1);
   glm_obj->missing = create_matrix(num_rows, 1);

   if(glm_obj->deg_free < 1) {
      fprintf(stderr,"Error: overspecified model, insufficient degrees of freedom\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
   }

   glm_obj->contrast = initial_data->contrast;
   if(glm_obj->contrast->avg_array != NULL) {
      if(glm_obj->contrast->avg_array->num_test > 0)
         glm_obj->avg_resid = create_matrix(num_rows,
                                            glm_obj->response->num_columns);
   }

   /* Get family information */

   glm_obj->family = initial_data->family;

   if ((glm_obj->family != QUASI) 
       && (initial_data->variance_function != VARIANCE_DEFAULT)) {
      fprintf(stderr, "\nError: If family type is not QUASI, variance_function is implied\n");
      glm_obj->family = QUASI;
   }

   /* Set up links and variances */
   
   /* Check to see if the link is the canonical link */

   if (initial_data->link == LINK_DEFAULT)
      glm_obj->control->is_canonical = TRUE;
   else
      glm_obj->control->is_canonical = FALSE;

   switch(glm_obj->family) {

   case BINOMIAL:
      if(glm_obj->control->est_scale == FALSE) {
         glm_obj->control->is_scale = TRUE;
         glm_obj->scale_est = 1.0 / glm_obj->control->binomial_n;
         glm_obj->pearson = glm_obj->scale_est / glm_obj->deg_free;
      }
      else
         glm_obj->control->is_scale = FALSE;

      glm_obj->variance_function = MU_ONE_MINUS_MU;
      switch(initial_data->link) {
      case LINK_DEFAULT:
         glm_obj->link = LOGIT;
         break;
      case LOGIT:
         glm_obj->link = LOGIT;
         glm_obj->control->is_canonical = TRUE;
         break;
      case PROBIT:
         glm_obj->link = PROBIT;
         break;
      case CLOGLOG:
         glm_obj->link = CLOGLOG;
         break;
      default:
         fprintf(stderr, "\nError: Improper choice of link function \n");
         fprintf(stderr, "\n For BINOMIAL, choices are: CLOGLOG, LOGIT, PROBIT\n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);  
         break;
      }
      if(initial_data->corr_struct == CORR_DEFAULT)
         glm_obj->corr_struct = DIAGONAL_CORR;

      break;

   case GAMMA:
      glm_obj->variance_function = MU_SQUARED;
         switch(initial_data->link) {
         case LINK_DEFAULT:
            glm_obj->link = INVERSE;
            break;
         case INVERSE:
            glm_obj->link = INVERSE;
            glm_obj->control->is_canonical = TRUE;
            break;
         case IDENTITY:
            glm_obj->link = IDENTITY;
            break;
         case LOG:
            glm_obj->link = LOG;
            break;
         default:
            fprintf(stderr, "\nError: Improper choice of link function \n");
            fprintf(stderr, "\n For GAMMA, choices are: IDENTITY, INVERSE, LOG\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
            break;
         }
      if(initial_data->corr_struct == CORR_DEFAULT)
         glm_obj->corr_struct = DIAGONAL_CORR;
      break;

   case EXPONENTIAL:
      glm_obj->variance_function = MU_SQUARED;
         switch(initial_data->link) {
         case LINK_DEFAULT:
            glm_obj->link = LOG;
            break;
         case LOG:
            glm_obj->link = LOG;
            break;
         default:
            fprintf(stderr, "\nError: Improper choice of link function \n");
            fprintf(stderr, "\n For EXPONENTIAL, link function is LOG.\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
            break;
         }
      glm_obj->scale_est = 1.0;
      glm_obj->pearson = glm_obj->scale_est / glm_obj->deg_free;
      if(initial_data->corr_struct == CORR_DEFAULT)
         glm_obj->corr_struct = DIAGONAL_CORR;
      break;

   case CHI_SQ:
      glm_obj->variance_function = MU_SQUARED;
         switch(initial_data->link) {
         case LINK_DEFAULT:
            glm_obj->link = INVERSE;
            break;
         case INVERSE:
            glm_obj->link = INVERSE;
            glm_obj->control->is_canonical = TRUE;
            break;
         default:
            fprintf(stderr, "\nError: Improper choice of link function \n");
            fprintf(stderr, "\n For CHI_SQ, link function is INVERSE.\n");
            delete_tmpfiles(&tmpfile_list);
            exit(EXIT_FAILURE);
            break;
         }
      if(initial_data->corr_struct == CORR_DEFAULT)
         glm_obj->corr_struct = DIAGONAL_CORR;
      break;

   case GAUSSIAN:
      glm_obj->variance_function = CONST;
      switch(initial_data->link) {
      case LINK_DEFAULT:
         glm_obj->link = IDENTITY;
         break;
      case IDENTITY:
         glm_obj->link = IDENTITY;
         glm_obj->control->is_canonical = TRUE;
         break;
      default:
         fprintf(stderr, "\nError: Improper choice of link function \n");
         fprintf(stderr, "\nError: For GAUSSIAN, choices are: IDENTITY\n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }
      if(initial_data->corr_struct == CORR_DEFAULT)
         glm_obj->corr_struct = DIAGONAL_CORR;
      break;
      
   case INV_GAUSS:
      glm_obj->variance_function = MU_CUBED;
      switch(initial_data->link) {
      case LINK_DEFAULT:
         glm_obj->link = INVERSE_SQ;
         break;
      case INVERSE_SQ:
         glm_obj->link = INVERSE_SQ;
         glm_obj->control->is_canonical = TRUE;
         break;
      default:
         fprintf(stderr, "\nError: Improper choice of link function \n");
         fprintf(stderr, "\nError: For INV_GAUSS, choices are: INVERSE_SQ\n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }
      if(initial_data->corr_struct == CORR_DEFAULT)
         glm_obj->corr_struct = DIAGONAL_CORR;
      break;
      
   case POISSON:
     if(glm_obj->control->est_scale == FALSE) {
           glm_obj->control->is_scale = TRUE;
           glm_obj->scale_est = 1.0;
           glm_obj->pearson = glm_obj->scale_est / glm_obj->deg_free;
        }
      else
         glm_obj->control->is_scale = FALSE;

      glm_obj->variance_function = MU;
      switch(initial_data->link) {
      case LINK_DEFAULT:
         glm_obj->link = LOG;
         break;
      case LOG:
         glm_obj->link = LOG;
         glm_obj->control->is_canonical = TRUE;
         break;
      case IDENTITY:
         glm_obj->link = IDENTITY;
         break;
      case SQRT:
         glm_obj->link = SQRT;
         break;
      default:
         fprintf(stderr, "\nError: Improper choice of link function \n");
         fprintf(stderr, "\nError: For POISSON, choices are: IDENTITY, LOG, SQRT\n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }
      if(initial_data->corr_struct == CORR_DEFAULT)
         glm_obj->corr_struct = DIAGONAL_CORR;
      break;

   case QUASI:
      glm_obj->variance_function = initial_data->variance_function;
      glm_obj->link = initial_data->link;
      if(initial_data->corr_struct == CORR_DEFAULT)
         glm_obj->corr_struct = DIAGONAL_CORR;
      break;
   }

   if(glm_obj->corr_struct != DIAGONAL_CORR)
      glm_obj->family = QUASI;

   /* Create necessary matrices */

   if(glm_obj->control->do_fit == TRUE) {

      if (glm_obj->family == GAUSSIAN) {
         glm_obj->beta_hat = create_matrix(glm_obj->design_matrix->num_columns,
                                           glm_obj->response->num_columns);
         glm_obj->lm_matrix = create_matrix(num_columns, num_rows);
         glm_obj->mu = create_matrix(num_rows, glm_obj->response->num_columns);
         glm_obj->info_inv = create_matrix(num_columns, num_columns);
         glm_obj->output = create_matrix(initial_data->contrast->num_contrasts,
                                         1);
         glm_obj->x_beta = NULL;
         glm_obj->z_i = NULL;
         glm_obj->g_prime = NULL;
         glm_obj->info_mat = NULL;
         glm_obj->sqrt_weight = NULL;

      }
      else {
         glm_obj->lm_matrix = create_matrix(num_columns, num_rows);
         glm_obj->beta_hat = create_matrix(glm_obj->design_matrix->num_columns,
                                           glm_obj->response->num_columns);
         glm_obj->mu = create_matrix(num_rows, glm_obj->response->num_columns);
         glm_obj->info_inv = create_matrix(num_columns, num_columns);
         glm_obj->x_beta = create_matrix(num_rows,
                                         glm_obj->response->num_columns);
         glm_obj->z_i = create_matrix(num_rows, 1);
         glm_obj->output = create_matrix(initial_data->contrast->num_contrasts,
                                         1);
         glm_obj->g_prime = create_matrix(num_rows,
                                       glm_obj->response->num_columns);
         glm_obj->info_mat = create_matrix(num_columns + 1, num_columns + 1);
         
         if(initial_data->corr_struct == CORR_DEFAULT)
            glm_obj->corr_struct = DIAGONAL_CORR;
         else
            glm_obj->corr_struct = initial_data->corr_struct;

         if(glm_obj->corr_struct != DIAGONAL_CORR) {
            glm_obj->sigma = create_matrix(num_rows, num_rows);
            glm_obj->weight= create_matrix(num_rows, glm_obj->response->num_columns);
         }
         else if(glm_obj->corr_struct == DIAGONAL_CORR) {
            glm_obj->sigma = create_matrix(num_rows, glm_obj->response->num_columns);
            glm_obj->weight= create_matrix(num_rows, glm_obj->response->num_columns);
         }
         glm_obj->sqrt_weight = create_matrix(num_rows,
                                             glm_obj->response->num_columns);
      }

      /* create tmp_matrices */

      create_contrast(glm_obj);
   }

   return; 
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_contrast
@INPUT      : pointer to Glm_Object
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: creates matrices needed for calculating contrasts
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void create_contrast(Glm_Object *glm_obj) 
{
   int icont, num_rows, need_tmp_con;
   Contrast_Matrix *cur_contrast;
   
   need_tmp_con = 0;

   glm_obj->tmp_con = GI_MALLOC(sizeof(Matrix **) *
                             glm_obj->contrast->num_contrasts);
  
   for(icont=0; icont < glm_obj->contrast->num_contrasts; icont++) { 

      cur_contrast = glm_obj->contrast->contrast_matrix_array[icont]; 

      if (glm_obj->family != GAUSSIAN)
         cur_contrast->deriv = create_matrix(glm_obj->response->num_rows, 1);

      if(cur_contrast->contrast_matrix != NULL) {
         need_tmp_con++;
         glm_obj->tmp_con[icont] = GI_MALLOC(sizeof(Matrix *) * 6);
         num_rows = cur_contrast->contrast_matrix->num_rows;

         glm_obj->tmp_con[icont][0]=create_matrix(num_rows,
                                               glm_obj->response->num_columns);
         glm_obj->tmp_con[icont][1]=create_matrix(glm_obj->design_matrix->num_columns, num_rows);
         glm_obj->tmp_con[icont][2]=create_matrix(num_rows, num_rows);
         glm_obj->tmp_con[icont][3]=create_matrix(num_rows, num_rows);
         glm_obj->tmp_con[icont][4]=create_matrix(num_rows,
                                             glm_obj->response->num_columns);
         glm_obj->tmp_con[icont][5]=create_matrix(1, 1);
         

      }
      else {
         glm_obj->tmp_con[icont] = NULL;
      }
   }

   if(need_tmp_con <= 0)  {    /* if no contrast matrices are being used */ 
      GI_FREE(glm_obj->tmp_con);  /* we don't need these temp matrices      */
      glm_obj->tmp_con = NULL;
   }

   return;                                 
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_design_matrix
@INPUT      : char *design_filename - name and/or path of file containing
                                      design filename
              char ***infiles - pointer to list of strings of input filenames
              Matrix **design_mat - pointer to design_matrix
              char *path - optional path for design file
@OUTPUT     : (none)
@RETURNS    : TRUE or FALSE
@DESCRIPTION: fixes pointer of type Matrix to the design matrix and
              a pointer to a list of input files
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int get_design_matrix(char *design_filename, char ***infiles,
                      Matrix ** design_matrix, char *path, int num_response,
                      int max_num_row, int max_num_col)
{
   FILE    *design_file;
   int      num_columns, num_rows;
   int      i, j;
   int      prev_j;
   int      path_len;
   double   tmp_values[max_num_row][max_num_col];

   VIO_STR  filename, fullname, tmp;
   double   value;

   /* open the design file */
   design_file = (design_filename == NULL) ? stdin : fopen(design_filename, "r");
   if(design_file == NULL) {
      fprintf(stderr, "\nError opening design filename %s.\n\n", design_filename);
      //delete_tmpfiles(&tmpfile_list); //VF: files will be deleted in main 
      return (FALSE);
   }

   /* get path length */
   path_len = (path == NULL) ? 0 : strlen(path);
   
   /* allocate memory */
   *infiles = GI_MALLOC(sizeof(char *) * (max_num_row + 2));
   filename = alloc_string(4096); 
   fullname = alloc_string(4096);

   if(verbose) {
      fprintf(stderr, "Reading Design file [%s]:", design_filename);
   }

   /* for each line of the input file... */

   char * tmp_string = alloc_string(10000);
   char * delims = " \t\n\0";

   /* for each line of the input file... 
      read the entire line with mni_input_string (will know about comment
      lines starting with #), then parse with strtok with spaces, tabs,
      new lines as delimiters. */

   i = 0;
   prev_j = -1;
   while(mni_input_string(design_file, &tmp_string, (char)0, (char)0) == VIO_OK) {

      filename = strtok( tmp_string, delims );

      /* allocate some space for the filename */
      (*infiles)[i] = GI_MALLOC(sizeof(char) * (strlen(filename) + path_len + 1));

      /* add in the path if required */
      if(path != NULL){
         
         fullname = concat_strings(path, filename);
         
         /* switch pointers about */
         tmp = filename;
         filename = fullname;
         fullname = tmp;
      }

      if(strcpy((*infiles)[i], filename) == NULL) {
         fprintf(stderr, "\nError getting name of infile %d in %s.",
                 i + 1, design_filename);
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }

      /* now get the values */
      j = 0;

      char * res = strtok( NULL, delims );
      while( res != NULL ) {
        tmp_values[i][j] = atof( res );
        j++;
        if(j > max_num_col) {
          fprintf(stderr,
                  "\nError: more columns in design_matrix [%d] than -max_num_col [%d]\n",
                  j, max_num_col);
          delete_tmpfiles(&tmpfile_list);
          exit(EXIT_FAILURE);
        }
        res = strtok( NULL, " \t\n\0" );
      }

      /* check number of values */
      if((prev_j != -1) && (j != prev_j)) {
         fprintf(stderr, "\nError: mismatched number of entries on line [%d]\n", i);
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }
      prev_j = j;

      i++;
      if(i > max_num_row) {
         fprintf(stderr,
                 "\nError: more rows in design_matrix [%d] than -max_num_row [%d]\n", i,
                 max_num_row);
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }
   }
   num_rows = i;
   num_columns = j;
   
   /* populate the design matrix */
   *design_matrix = create_matrix(num_rows, num_columns);
   for(i = 0; i < num_rows; i++) {
      for(j = 0; j < num_columns; j++) {
         (*design_matrix)->values[i][j] = tmp_values[i][j];
      }
   }

   /* clean up */
   fclose(design_file);
   delete_string(filename);
   delete_string(fullname);
   delete_string(tmp_string);
   
   return TRUE;
}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : fit_glm
@INPUT      : pointer to Glm_Object
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: fits a glm according to response from voxel_function, stores
              attributes in glm_obj
@METHOD     : 
@GLOBALS    : 
@CALLS      : general_least_sq, calculate_deviance, calculate_link,
              calculate_weight, calculate_contrasts
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void fit_glm(Glm_Object *glm_obj)
{
   double old_dev;
   int i;

   /* Reset distance and deviance of glm_obj */

   if(glm_obj->control->do_fit == TRUE) {

      glm_obj->distance = 100.0;
      glm_obj->deviance = INVALID_DATA;

      glm_obj->control->use_mu = ((glm_obj->control->use_mu_switch) && 
                                  (glm_obj->control->last_status == CONVERGE));

      if(glm_obj->family != GAUSSIAN) {

         glm_obj->control->num_iter = 1;

         initialize_mu(glm_obj);

         /* Check to see if old mu's were used or
            whether it was initialized with y*/

         if(glm_obj->control->use_mu == FALSE) {
            calculate_link(glm_obj->z_i, glm_obj, NORMAL);
            calculate_link(glm_obj->g_prime, glm_obj, DERIV);

            if(glm_obj->corr_struct != DIAGONAL_CORR) {
               calculate_sigma_inv(glm_obj->sigma, glm_obj);
            }
            else {
               calculate_weight(glm_obj->weight, glm_obj);
            }
            general_least_sq(glm_obj);

            multiply_matrix(glm_obj->x_beta, glm_obj->design_matrix,
                            glm_obj->beta_hat, MATRIX_NEITHER);

            calculate_link(glm_obj->mu, glm_obj, INVER);
            
            old_dev = calculate_deviance(glm_obj); 
         }

         while (glm_obj->control->status == CONT) {
         
            calculate_link(glm_obj->g_prime, glm_obj, DERIV); 
            
            for(i=0; i<glm_obj->response->num_rows; i++) {
               glm_obj->z_i->values[i][0] = ((glm_obj->g_prime->values[i][0] * 
                                              (glm_obj->response->values[i][0]-
                                               glm_obj->mu->values[i][0]))+
                                             glm_obj->x_beta->values[i][0]);
            }

            if(glm_obj->corr_struct != DIAGONAL_CORR) {
               calculate_sigma_inv(glm_obj->sigma, glm_obj);
            }
            else {
               calculate_weight(glm_obj->weight, glm_obj);
            }

            general_least_sq(glm_obj);
            
            multiply_matrix(glm_obj->x_beta, glm_obj->design_matrix,
                            glm_obj->beta_hat, MATRIX_NEITHER);
         
            calculate_link(glm_obj->mu, glm_obj, INVER);

            glm_obj->deviance = calculate_deviance(glm_obj);

            glm_obj->distance = ((fabs(glm_obj->deviance - old_dev)) /
                                 (old_dev + glm_obj->control->tol));
            
            if (glm_obj->distance <= glm_obj->control->tol) {
               glm_obj->control->status = CONVERGE;
               copy_matrix(glm_obj->info_inv,
                           glm_obj->info_mat, glm_obj->num_var,
                           glm_obj->num_var, (-1.0), MATRIX_NEITHER);
            }

            old_dev = glm_obj->deviance;

            glm_obj->control->num_iter++;

            if(glm_obj->control->num_iter > glm_obj->control->max_iter)
               glm_obj->control->status = MAX_ITER;
         }
      }
      
      else if(glm_obj->family == GAUSSIAN) {
         for(i=0; i<glm_obj->response->num_rows; i++) {
            glm_obj->response->values[i][0] *= glm_obj->missing->values[i][0];
         }
         multiply_matrix(glm_obj->beta_hat, glm_obj->lm_matrix,
                         glm_obj->response, MATRIX_NEITHER);
         multiply_matrix(glm_obj->mu, glm_obj->design_matrix,
                         glm_obj->beta_hat, MATRIX_NEITHER);
         glm_obj->deviance = calculate_deviance(glm_obj);
      }
   
      if (glm_obj->control->status == CONVERGE){
            calculate_contrasts(glm_obj);
      }      
   }

   return ;
   
}
/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_info_matrix
@INPUT      : pointer to Glm_Object
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: initializes a (p+1) X (p+1) matrix with t(x) %*% w %*% x in
              first p rows and p columns, and the last row and columns is
              that of the sums of squares and products of the response
               - to be used for Gaussian sweep operator by Beaton(1964)
               found in McCullagh and Nelder (1989) p. 82-3.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void initialize_info_matrix(Glm_Object *glm_obj) 
{
   int i,j,k;
   int n,p;
   double tmp_value;

   n = glm_obj->design_matrix->num_rows;
   p = glm_obj->design_matrix->num_columns;

   if(glm_obj->corr_struct == DIAGONAL_CORR) {
      for(i=0; i<p; i++) {
         for(j=0; j<=i; j++) {
            tmp_value = 0.0;
            for(k=0; k<n; k++) {
               tmp_value += (glm_obj->weight->values[k][0] *
                             glm_obj->design_matrix->values[k][i] *
                             glm_obj->design_matrix->values[k][j] *
                             glm_obj->missing->values[k][0]);
            }
            glm_obj->info_mat->values[i][j] = tmp_value;
            glm_obj->info_mat->values[j][i] = tmp_value;
         }
      }

      for(i=0; i<p; i++) {
         tmp_value = 0.0;
         for(k=0; k<n; k++) {
            tmp_value += (glm_obj->weight->values[k][0] *
                          glm_obj->z_i->values[k][0] *
                          glm_obj->design_matrix->values[k][i] *
                          glm_obj->missing->values[k][0]);
         }
         glm_obj->info_mat->values[i][p] = tmp_value;
         glm_obj->info_mat->values[p][i] = tmp_value;
      }
   }     

   tmp_value = 0.0;
   /* Mental note: shouldn't it be [k][0] below??? check to see if it affects
      results */

   for(k=0; k<n; k++) {
      tmp_value += (glm_obj->z_i->values[k][k] *
                    glm_obj->z_i->values[k][k] *
                    glm_obj->missing->values[k][0]);
   }
   glm_obj->info_mat->values[p][p] = tmp_value;
   
   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : general_least_sq
@INPUT      : pointer to Glm_Object
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: performs generalized least squares calculation, coefficients
              stored in glm_obj->beta_hat
@METHOD     : 
@GLOBALS    : 
@CALLS      : multiply_matrix, multiply_matrix_sym, invert_matrix_sym
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void general_least_sq(Glm_Object *glm_obj)
{
   int i, j, k, p;
   int i_tmp1, i_tmp2, j_tmp1, j_tmp2;
   double a_kk, abs_a_kk;
   double tol = 0.0001;

   p = glm_obj->design_matrix->num_columns;

   initialize_info_matrix(glm_obj);

   /* Perform sweeps described in McCullagh and Nelder (1989) pp. 82-3 
      Note: careful care must be taken to ensure that only values on one side
      of the diagonal are used at each step, (i.e. i_tmp1, i_tmp2, etc.) */

   /* Here we write over the values of beta_hat with the diagonal values of
      info_mat to be used for the test for collinearity instead of creating a
      matrix every time */

   for(k=0; k < p; k++) {
      if(fabs(glm_obj->beta_hat->values[k][0] /
              glm_obj->info_mat->values[k][k]) <= tol) {
      /*   fprintf(stderr,"\nWarning: Fit may be singular, column %d of design matrix may be aliased.\n",k+1); */
         glm_obj->control->status = SING_MATRIX;
      }
      glm_obj->beta_hat->values[k][0] = glm_obj->info_mat->values[k][k];
   }



   for(k = 0; k < p; k++) {

      a_kk = glm_obj->info_mat->values[k][k];
      abs_a_kk = fabs(a_kk);

      /* Test for collinearity as in  Clarke (1982) Appl. Stat. pp.166-7*/
      if (k >= 1) {
         for(i = 0; i< k; i++) {
            if (1.0 /(((glm_obj->info_mat->values[i][k] * 
                       glm_obj->info_mat->values[i][k]) / 
                       glm_obj->info_mat->values[k][k])  -
                      glm_obj->info_mat->values[i][i]) <= 
                glm_obj->beta_hat->values[i][0] * glm_obj->control->tol) {
/*               fprintf(stderr,"Warning: Fit may be singular, column %d of design matrix may be aliased. \n",i+1); */
               glm_obj->control->status = SING_MATRIX;
               return;
            }
         }
      }
               
      for(i=0; i < glm_obj->info_mat->num_rows; i++) {
         for(j=0; j <= i; j++) {
            if((i > k) && (j == k))
              glm_obj->info_mat->values[i][j]=(glm_obj->info_mat->values[j][i]/
                                               abs_a_kk);
            else if ((i == k) && (j != k))
              glm_obj->info_mat->values[i][j]=(glm_obj->info_mat->values[j][i]/
                                               abs_a_kk);
            else if ((i == k) && (j == k))
               glm_obj->info_mat->values[i][j] = -1.0/a_kk;
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
                            
              glm_obj->info_mat->values[i][j]=(glm_obj->info_mat->values[j][i]-
                                     glm_obj->info_mat->values[i_tmp1][j_tmp1]*
                                     glm_obj->info_mat->values[i_tmp2][j_tmp2]/
                                     a_kk);
            }
         }
      }
      for(i=0; i<glm_obj->info_mat->num_rows; i++) {
         for(j=0; j< i; j++) {
            glm_obj->info_mat->values[j][i] = glm_obj->info_mat->values[i][j];
         }
      }
   }

   for(i=0; i<p; i++)
      glm_obj->beta_hat->values[i][0] = glm_obj->info_mat->values[i][p];

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : calculate_sigma
@INPUT      : pointer to Glm_Object
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: calculates variance-covaraince matrix inverse
              to be used in generalized least squares procedure
@METHOD     : 
@GLOBALS    : 
@CALLS      : calculate_link
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void calculate_sigma_inv(Matrix *matrix_result, Glm_Object *glm_obj) 
{
   switch(glm_obj->corr_struct) {

   case AUTO_REG:
      fprintf(stderr,"\nError: Auto-Reg(1) correlation structure not ready yet.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
      calculate_weight(glm_obj->weight, glm_obj);
      break;

   default:
      fprintf(stderr,"\nError: Not a valid correlation structure.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
      break;
   }

   return;
}
/* ----------------------------- MNI Header -----------------------------------
@NAME       : calculate_weight
@INPUT      : pointer to Glm_Object
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: calculates weights to be used in weighted least
              squares procedure
@METHOD     : 
@GLOBALS    : 
@CALLS      : calculate_link
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void calculate_weight(Matrix *matrix_result, Glm_Object *glm_obj) 
{
   int i;
   double sign;

     /* Check if num_rows and num_columns are right */
   if((matrix_result->num_columns != glm_obj->response->num_columns) || 
      (matrix_result->num_rows != glm_obj->response->num_rows)) {
      fprintf(stderr,"\nError in calculate_weight:\n");
      fprintf(stderr,"Size of matrix_result doesn't match glm_obj->weight.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
   }
  
   /* If link function is the canonical link function, then the weights
      are equal to 1.0/g'(mu) */

   switch(glm_obj->control->is_canonical) {
   case TRUE:
      for(i=0; i<matrix_result->num_rows; i++) {
         matrix_result->values[i][0] =1.0/fabs(glm_obj->g_prime->values[i][0]);
         if (glm_obj->control->is_fwhm == TRUE) 
            glm_obj->sqrt_weight->values[i][0] =\
               sqrt(matrix_result->values[i][0]);
      }
      break;

   case FALSE:
      switch(glm_obj->variance_function) {
      case CONST:
         for(i=0; i<matrix_result->num_rows; i++)
            glm_obj->sqrt_weight->values[i][0] = 1.0;
         break;
   
      case MU:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(fabs(glm_obj->mu->values[i][0]) <= 0.01) {
               if(glm_obj->mu->values[i][0] <= 0)
                  sign = -1.0;
               else if(glm_obj->mu->values[i][0] > 0)
                  sign = 1.0;
               glm_obj->sqrt_weight->values[i][0] = 0.01 * sign;
            }
            else
               glm_obj->sqrt_weight->values[i][0] = glm_obj->mu->values[i][0];
         }
         break;
         
      case MU_SQUARED:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(fabs(glm_obj->mu->values[i][0]) <= 0.01) 
               glm_obj->sqrt_weight->values[i][0] = 0.0001;
            else
               glm_obj->sqrt_weight->values[i][0] = (glm_obj->mu->values[i][0]*
                                                    glm_obj->mu->values[i][0]);
         }
         break;
         
      case MU_CUBED:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(fabs(glm_obj->mu->values[i][0]) <= 0.01) {
               if(glm_obj->mu->values[i][0] <= 0)
                  sign = -1.0;
               if(glm_obj->mu->values[i][0] > 0)
                  sign = 1.0;
               glm_obj->sqrt_weight->values[i][0] = 0.000001 * sign;
            }
            else
               glm_obj->sqrt_weight->values[i][0] = (glm_obj->mu->values[i][0] *
                                                    glm_obj->mu->values[i][0] *
                                                    glm_obj->mu->values[i][0]);
         }
         break;
         
      case MU_ONE_MINUS_MU:
         for(i=0; i<matrix_result->num_rows; i++) {
            if((glm_obj->mu->values[i][0] <= 0.01) ||
               (glm_obj->mu->values[i][0] >= 0.99)) {
               glm_obj->sqrt_weight->values[i][0] = 0.0099;
            }
            else
               glm_obj->sqrt_weight->values[i][0] = (glm_obj->mu->values[i][0] *
                                              (1 - glm_obj->mu->values[i][0]));
         }
         break;
   
      default:
         fprintf(stderr,"\nError: improper variance function. \n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }   
   
      for(i=0; i<matrix_result->num_rows; i++) {
         matrix_result->values[i][0] = 1.0/(glm_obj->g_prime->values[i][0] *
                      glm_obj->g_prime->values[i][0] *
                      glm_obj->sqrt_weight->values[i][0]);
         if (glm_obj->control->is_fwhm == TRUE) 
            glm_obj->sqrt_weight->values[i][0] =\
               sqrt(matrix_result->values[i][0]);
      }
      break;

   default:
      fprintf(stderr,"\nError in calculate_weight.\n");
      fprintf(stderr,"\nImproper value for is_canonical.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
      break;
   }
   
   return; 
   
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_mu
@INPUT      : glm_obj - Pointer to glm_object
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: initializes glm_obj->mu with appropriate values
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void initialize_mu(Glm_Object *glm_obj) 
{
   int i;
   double eps;

   eps = glm_obj->control->eps_sing;

   switch(glm_obj->link) {
   case LOGIT:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if(glm_obj->response->values[i][0] >= 1.0 - eps) {
            glm_obj->response->values[i][0] = 1.0 - eps;
            if(glm_obj->control->use_mu == FALSE)
               glm_obj->mu->values[i][0] = 1.0 - eps;
         }
         else if(glm_obj->response->values[i][0] <= eps) {
            glm_obj->response->values[i][0] = eps;
            if(glm_obj->control->use_mu == FALSE)
               glm_obj->mu->values[i][0] = eps;
         }
         else if(glm_obj->control->use_mu == FALSE)
            glm_obj->mu->values[i][0] = glm_obj->response->values[i][0];
      }
      break;
      
   case CLOGLOG:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if(glm_obj->response->values[i][0] >= 1.0 - eps) {
            glm_obj->response->values[i][0] = 1.0 - eps;
            if(glm_obj->control->use_mu == FALSE)
               glm_obj->mu->values[i][0] = 1.0 - eps;
         }
         else if(glm_obj->response->values[i][0] <= eps) {
            glm_obj->response->values[i][0] = eps;
            if(glm_obj->control->use_mu == FALSE)
               glm_obj->mu->values[i][0] = eps;
         }
         else if(glm_obj->control->use_mu == FALSE)
            glm_obj->mu->values[i][0] = glm_obj->response->values[i][0];
      }
      break;

   case LOG:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if(glm_obj->response->values[i][0] <= eps) {
            glm_obj->response->values[i][0] = eps;
            if(glm_obj->control->use_mu == FALSE)
               glm_obj->mu->values[i][0] = eps;
         }
         else if(glm_obj->control->use_mu == FALSE)
               glm_obj->mu->values[i][0] = glm_obj->response->values[i][0];
      }
      break;

   case SQRT:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if(glm_obj->response->values[i][0] <= eps) {
            glm_obj->response->values[i][0] = eps;
            if(glm_obj->control->use_mu == FALSE)
               glm_obj->mu->values[i][0] = eps;
         }
         else if(glm_obj->control->use_mu == FALSE)
            glm_obj->mu->values[i][0] = glm_obj->response->values[i][0];
      }
      break;

   case INVERSE:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if(glm_obj->response->values[i][0] <= eps) {
            glm_obj->response->values[i][0] = eps;
            if(glm_obj->control->use_mu == FALSE)
               glm_obj->mu->values[i][0] = eps;
         }
         else if(glm_obj->control->use_mu == FALSE)
            glm_obj->mu->values[i][0] = glm_obj->response->values[i][0];
      }
      break;

   case INVERSE_SQ:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if(glm_obj->response->values[i][0] <= eps) {
            glm_obj->response->values[i][0] = eps;
            if(glm_obj->control->use_mu == FALSE)
               glm_obj->mu->values[i][0] = eps;
         }
         else if(glm_obj->control->use_mu == FALSE)
            glm_obj->mu->values[i][0] = glm_obj->response->values[i][0];
      }
      break;

   case IDENTITY:
      if(glm_obj->control->use_mu == FALSE) {
         copy_matrix(glm_obj->mu, glm_obj->response,
                     glm_obj->response->num_rows,
                     glm_obj->response->num_columns, 1.0, MATRIX_NEITHER);
      }
      break;

   default:
      fprintf(stderr,"\nError: other links aren't ready in initialize_mu.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
      break;
   }

   return;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : calculate_link
@INPUT      : pointer to Glm_Object
@OUTPUT     : (none)
@RETURNS    : void
@DESCRIPTION: calculates link function for use with iterative weighted
              least squares procedure
@METHOD     : 
@GLOBALS    : 
@CALLS      : math.h functions
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void calculate_link(Matrix *matrix_result, Glm_Object *glm_obj, Link_Mode mode) 
{
   double tmp_value;
   double eps;
   int i;

   /* Check if num_rows and num_columns are right */
   if((matrix_result->num_columns != glm_obj->mu->num_columns) || 
      (matrix_result->num_rows != glm_obj->mu->num_rows)) {
      fprintf(stderr,"\nError in calculate_link:\n");
      fprintf(stderr,"Size of matrix_result doesn't match glm_obj->mu.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
   }

   eps = (glm_obj->control->eps_sing) / 10.0;

   switch(glm_obj->link) {
   
   case LOGIT:
      switch(mode) {
      case NORMAL:
         for(i=0; i<glm_obj->mu->num_rows; i++) {
            if ((glm_obj->mu->values[i][0] < 0.0) ||
                (glm_obj->mu->values[i][0] > 1.0)) {
               fprintf(stderr,"\nError: values outside the range 0 and 1 in logit transform.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
            if (glm_obj->mu->values[i][0] <= eps) {
               matrix_result->values[i][0] = log(eps / (1.0 - eps));
            }
            else if (glm_obj->mu->values[i][0] >= 1.0 - eps) {
               matrix_result->values[i][0] = log((1.0 - eps) / eps);
            }
            else {
               matrix_result->values[i][0] = log(glm_obj->mu->values[i][0]/
                           (1.0-glm_obj->mu->values[i][0]));
            }
         }
         break;
         
      case DERIV:
         for(i=0; i<glm_obj->mu->num_rows; i++) {
            if ((glm_obj->mu->values[i][0] < 0.0) ||
                (glm_obj->mu->values[i][0] > 1.0)) {
               fprintf(stderr,"\nError: values outside the range 0 and 1 in logit transform.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
            if (glm_obj->mu->values[i][0] <= eps) {
               matrix_result->values[i][0] = 1.0/(eps * (1.0 - eps));
            }
            else if (glm_obj->mu->values[i][0] >= 1.0 - eps) {
               matrix_result->values[i][0] = 1.0/(0.0005*0.9995);
            }
            else {
               matrix_result->values[i][0] = 1.0/((glm_obj->mu->values[i][0]) *
                            (1.0 - glm_obj->mu->values[i][0]));
            }
         }
         break;

      case INVER:
         for(i=0; i<glm_obj->mu->num_rows; i++) {
            tmp_value = exp(glm_obj->x_beta->values[i][0]);
            matrix_result->values[i][0] = tmp_value / (1.0 + tmp_value);
        }
         break;
            
      default:
         fprintf(stderr,"\n Improper Link_Mode choice \n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }
      break;
      
   case CLOGLOG:
      switch(mode) {
      case NORMAL:
         for(i=0; i<glm_obj->mu->num_rows; i++) {
            if ((glm_obj->mu->values[i][0] < 0.0) ||
                (glm_obj->mu->values[i][0] > 1.0)) {
               fprintf(stderr,"\nError: values outside the range 0 and 1 in c-log-log transform.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
            if (glm_obj->mu->values[i][0] <= eps) {
               matrix_result->values[i][0] = log(-log(1.0 - eps));
            }
            else if (glm_obj->mu->values[i][0] >= 1.0 - eps) {
               matrix_result->values[i][0] = log(-log(eps));
            }
            else {
               matrix_result->values[i][0] = log(-log(1.0 - glm_obj->mu->values[i][0]));
            }
         }
         break;
 
      case DERIV:
         for(i=0; i<glm_obj->mu->num_rows; i++) {
            if ((glm_obj->mu->values[i][0] < 0.0) ||
                (glm_obj->mu->values[i][0] > 1.0)) {
               fprintf(stderr,"\nError: values outside the range 0 and 1 in c-log-log transform.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
            if (glm_obj->mu->values[i][0] <= eps) {
               matrix_result->values[i][0] = -1.0 / (log(1.0 - eps) * 
                                                     (1.0 - eps));
            }
            else if (glm_obj->mu->values[i][0] >= 1.0 - eps) {
               matrix_result->values[i][0] = -1.0 / (log(eps) * (eps));
            }
            else {
               tmp_value = 1.0 - glm_obj->mu->values[i][0];
               matrix_result->values[i][0] = -((1.0 / log(tmp_value)) *
                                               (1.0 / tmp_value));
            }
         }
         break;            
    
      case INVER:
         for(i=0; i<glm_obj->mu->num_rows; i++) {
            matrix_result->values[i][0] = 1 - exp(-exp(glm_obj->x_beta->values[i][0]));
         }
         break;
      
      default:
         fprintf(stderr,"\n Improper Link_Mode choice \n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }
      break;
      
   case IDENTITY:
      switch(mode) {
      case NORMAL:
         for(i=0; i<matrix_result->num_rows; i++)  {
            matrix_result->values[i][0] = glm_obj->mu->values[i][0];
         }
         break;
         
      case DERIV:
         for(i=0; i<matrix_result->num_rows; i++) {
            matrix_result->values[i][0] = 1.0;
         }
         break;
            
      case INVER:
         for(i=0; i<matrix_result->num_rows; i++) {
            matrix_result->values[i][0] = glm_obj->x_beta->values[i][0];
         }
         break;
            
      default:
         fprintf(stderr,"\n Improper Link_Mode choice \n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }
      break;
      
   case INVERSE:
      switch(mode) {
      case NORMAL:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(glm_obj->mu->values[i][0] < 0.0) {
               fprintf(stderr,"\nError: negative values in inverse link transform.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
            if(glm_obj->mu->values[i][0] < eps) {
               matrix_result->values[i][0] = 1.0 / eps;
            }
            else {
               matrix_result->values[i][0] = 1.0/(glm_obj->mu->values[i][0]);
            }
         }
         break;

      case DERIV:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(glm_obj->mu->values[i][0] < 0.0) {
               fprintf(stderr,"\nError: negative values in inverse link transform.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
            if(glm_obj->mu->values[i][0] < eps) {
               matrix_result->values[i][0] = -1.0 / (eps * eps);
            }
            else {
               matrix_result->values[i][0] = -1.0/((glm_obj->mu->values[i][0] *
                              glm_obj->mu->values[i][0]));
            }
         }
         break;

      case INVER:
         for(i=0; i<matrix_result->num_rows; i++) {
            matrix_result->values[i][0] = 1.0/(glm_obj->x_beta->values[i][0]);
         }
         break;
         
      default:
         fprintf(stderr,"\n Improper Link_Mode choice \n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }
      break;

   case INVERSE_SQ:
      switch(mode) {
      case NORMAL:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(glm_obj->mu->values[i][0] < 0.0) {
               fprintf(stderr,"\nError: negative values in inverse_sq link transform. \n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            } 
            if(fabs(glm_obj->mu->values[i][0]) < eps) {
               matrix_result->values[i][0] = 1.0 / (eps * eps);
            }
            else {
               matrix_result->values[i][0] = 1.0/((glm_obj->mu->values[i][0]) *
                            (glm_obj->mu->values[i][0]));
            }
         }
         break;
            
      case DERIV:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(glm_obj->mu->values[i][0] < 0.0) {
               fprintf(stderr,"\nError: negative values in inverse_sq link transform. \n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
            if(glm_obj->mu->values[i][0] < eps) {
               matrix_result->values[i][0] = -1.0 / (eps * eps * eps);
            }
            else {
               matrix_result->values[i][0] =-1.0/((glm_obj->mu->values[i][0]) *
                                                  (glm_obj->mu->values[i][0]) *
                                                  (glm_obj->mu->values[i][0]));
            }
         }
         break;
         
      case INVER:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(glm_obj->x_beta->values[i][0] < 0.0) {
               fprintf(stderr,"\nError: negative values in inverse_sq link transform. \n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }  
            matrix_result->values[i][0] = 1.0/(sqrt(glm_obj->x_beta->values[i][0]));
         }
         break;
         
      default:
         fprintf(stderr,"\n Improper Link_Mode choice \n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }
   
      break;

   case PROBIT:
      fprintf(stderr,"\nError: probit link not ready yet.\n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
      break;
      
   case LOG:
      switch(mode) {
      case NORMAL:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(glm_obj->mu->values[i][0] < 0.0) {
               fprintf(stderr,"\nError: negative values in log link transform.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
            if(glm_obj->mu->values[i][0] <= eps) {
               matrix_result->values[i][0] = log(eps);
            }
            else {
               matrix_result->values[i][0] = log(glm_obj->mu->values[i][0]);
            }
         }
         break;

      case DERIV:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(glm_obj->mu->values[i][0] < 0.0) {
               fprintf(stderr,"\nError: negative values in log link transform.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            } 
            if(glm_obj->mu->values[i][0] <= eps) {
               matrix_result->values[i][0] = 1.0 / eps;
            }
            else {
               matrix_result->values[i][0] = 1.0/(glm_obj->mu->values[i][0]);
            }
         }
         break;
         
      case INVER:
         for(i=0; i<matrix_result->num_rows; i++) {
            matrix_result->values[i][0] = exp(glm_obj->x_beta->values[i][0]);
         }
         break;
         
      default:
         fprintf(stderr,"\n Improper Link_Mode choice \n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }
      break;
      
   case SQRT:
      switch(mode) {
      case NORMAL:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(glm_obj->mu->values[i][0] < 0.0) {
               fprintf(stderr,"\nError: negative values in sqrt link transform.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
            if(glm_obj->mu->values[i][0] <= eps) {
               matrix_result->values[i][0] = sqrt(eps);
            }
            else {
               matrix_result->values[i][0] = sqrt(glm_obj->mu->values[i][0]);
            }
         }
         break;
         
      case DERIV:
         for(i=0; i<matrix_result->num_rows; i++) {
            if(glm_obj->mu->values[i][0] < 0.0) {
               fprintf(stderr,"\nError: negative values in sqrt link transform.\n");
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
            if(glm_obj->mu->values[i][0] <= eps) {
               matrix_result->values[i][0] = 1.0 / sqrt(eps);
            }
            else {
               matrix_result->values[i][0] = 1.0/sqrt(glm_obj->mu->values[i][0]);
            }
         }
         break;
         
      case INVER:
         for(i=0; i<matrix_result->num_rows; i++) {
            matrix_result->values[i][0] = ((glm_obj->x_beta->values[i][0]) * 
                     (glm_obj->x_beta->values[i][0]));
         }
         break;

      default:
         fprintf(stderr,"\n Improper Link_Mode choice \n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }
      break;
      
   default:
      fprintf(stderr,"\n Improper link function \n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
      break;
   }
      
   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : calculate_contrasts
@INPUT      : pointer to Glm_Object
@OUTPUT     : 
@RETURNS    : void
@DESCRIPTION: calculates values for output, based on contrast matrices entered
              on command line
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void calculate_contrasts(Glm_Object *glm_obj) 
{
   Contrast_Matrix *cur_contrast;
   double f_stat, t_stat;
   int i;

   for(i=0; i<glm_obj->contrast->num_contrasts; i++) {

      cur_contrast = glm_obj->contrast->contrast_matrix_array[i];

      switch(cur_contrast->out_mode) {

      case T_STAT:
         if(cur_contrast->contrast_matrix != NULL) {
            /* calculate c beta */
            multiply_matrix(glm_obj->tmp_con[i][0],
                            cur_contrast->contrast_matrix,
                            glm_obj->beta_hat, MATRIX_NEITHER);
            
            /* calculate (c (xt_w_x_inv) ct)_inv */

            if (glm_obj->family != GAUSSIAN) {
               multiply_matrix(glm_obj->tmp_con[i][1],
                               glm_obj->info_inv,
                               cur_contrast->contrast_matrix, MATRIX_RIGHT);
               multiply_matrix_sym(glm_obj->tmp_con[i][2],
                                   cur_contrast->contrast_matrix, 
                                   glm_obj->tmp_con[i][1], MATRIX_NEITHER);
            }

            t_stat = (glm_obj->tmp_con[i][0]->values[0][0] /
                      sqrt(glm_obj->tmp_con[i][2]->values[0][0]));

            if (cur_contrast->stdev_mode == VOXEL_SD){
               t_stat = t_stat / sqrt(glm_obj->pearson);
            }

            glm_obj->output->values[i][0] = t_stat; 
         }
         else {
            glm_obj->output->values[i][0] = \
               (glm_obj->beta_hat->values[cur_contrast->column_num][0] / \
                sqrt(glm_obj->info_inv->values[cur_contrast->column_num][cur_contrast->column_num]));
            if (cur_contrast->stdev_mode == VOXEL_SD){
               glm_obj->output->values[i][0] =
                  glm_obj->output->values[i][0] / sqrt(glm_obj->pearson);
            }
         }
         break;

      case P_CORR:
         if(cur_contrast->contrast_matrix != NULL) {
            /* calculate c beta */
            multiply_matrix(glm_obj->tmp_con[i][0],
                            cur_contrast->contrast_matrix,
                            glm_obj->beta_hat, MATRIX_NEITHER);
            
            /* calculate (c (xt_w_x_inv) ct)_inv */

            if (glm_obj->family != GAUSSIAN) {
               multiply_matrix(glm_obj->tmp_con[i][1],
                               glm_obj->info_inv,
                               cur_contrast->contrast_matrix, MATRIX_RIGHT);
               multiply_matrix_sym(glm_obj->tmp_con[i][2],
                                   cur_contrast->contrast_matrix, 
                                   glm_obj->tmp_con[i][1], MATRIX_NEITHER);
            }

            t_stat = (glm_obj->tmp_con[i][0]->values[0][0] /
                      sqrt(glm_obj->tmp_con[i][2]->values[0][0]));

            if (cur_contrast->stdev_mode == VOXEL_SD){
               t_stat = t_stat / sqrt(glm_obj->pearson);
            }

            if(cur_contrast->stdev_mode == VOXEL_SD)
               glm_obj->output->values[i][0] = t_stat / sqrt(t_stat * t_stat
                                                             + glm_obj->deg_free);
            else if (cur_contrast->stdev_mode == POOLED_SD)
               glm_obj->output->values[i][0] = t_stat;
         }
         else {
            t_stat = \
               (glm_obj->beta_hat->values[cur_contrast->column_num][0]/\
                sqrt(glm_obj->info_inv->values[cur_contrast->column_num][cur_contrast->column_num]));
            if (cur_contrast->stdev_mode == VOXEL_SD) {
               t_stat =  t_stat/sqrt(glm_obj->pearson);
               glm_obj->output->values[i][0] = t_stat / sqrt(t_stat * t_stat
                                                             + glm_obj->deg_free);
            }
            else if (cur_contrast->stdev_mode == POOLED_SD) {
               glm_obj->output->values[i][0] = t_stat;
            }
         }
         break;
         
      case F_STAT:
         /* calculate c beta */
         multiply_matrix(glm_obj->tmp_con[i][0],
                         cur_contrast->contrast_matrix, glm_obj->beta_hat, MATRIX_NEITHER);
         
         /* calculate (c (xt_w_x_inv) ct)_inv */

         if (glm_obj->family != GAUSSIAN) {
            multiply_matrix(glm_obj->tmp_con[i][1],glm_obj->info_inv,
                            cur_contrast->contrast_matrix, MATRIX_RIGHT);
            multiply_matrix_sym(glm_obj->tmp_con[i][2],
                                cur_contrast->contrast_matrix, 
                                glm_obj->tmp_con[i][1], MATRIX_NEITHER);
            invert_matrix_sym(glm_obj->tmp_con[i][3], glm_obj->tmp_con[i][2]);
         }

         multiply_matrix(glm_obj->tmp_con[i][4], glm_obj->tmp_con[i][3],
                         glm_obj->tmp_con[i][0], MATRIX_NEITHER);
         multiply_matrix(glm_obj->tmp_con[i][5], glm_obj->tmp_con[i][0],
                         glm_obj->tmp_con[i][4], MATRIX_LEFT);

         f_stat = (glm_obj->tmp_con[i][5]->values[0][0] /
                   cur_contrast->contrast_matrix->num_rows);
         
         if (cur_contrast->stdev_mode == VOXEL_SD){
            f_stat = f_stat / glm_obj->pearson;
         }
         if (f_stat < 0) {
         print_matrix(stderr, glm_obj->tmp_con[i][0]);
         print_matrix(stderr, glm_obj->tmp_con[i][4]);
         }

         glm_obj->output->values[i][0] = f_stat;

         break;

      case STDEV_BETA:
         if(cur_contrast->contrast_matrix != NULL) {
            if(cur_contrast->contrast_matrix->num_rows == 1) {
  
               /* calculate (c (xt_w_x_inv) ct)_inv */
               if (glm_obj->family != GAUSSIAN) {
                  multiply_matrix(glm_obj->tmp_con[i][1],glm_obj->info_inv,
                                  cur_contrast->contrast_matrix, MATRIX_RIGHT);
                  multiply_matrix_sym(glm_obj->tmp_con[i][2],
                                      cur_contrast->contrast_matrix, 
                                      glm_obj->tmp_con[i][1], MATRIX_NEITHER);
               }

               glm_obj->output->values[i][0]=sqrt(glm_obj->pearson * glm_obj->tmp_con[i][2]->values[0][0]);
               
            }
            else {
               fprintf(stderr,"\nError in calculate_contrasts.\n");
               fprintf(stderr,"Var(Contrast %d) is a matrix, it is scalar only if num_rows is 1.",i+1);
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
         }
         else {
            glm_obj->output->values[i][0] = sqrt(glm_obj->pearson * glm_obj->info_inv->values[cur_contrast->column_num][cur_contrast->column_num]); 
         }
         break;
         
      case BETA_HAT:
         if(cur_contrast->contrast_matrix != NULL) {
            if(cur_contrast->contrast_matrix->num_rows == 1) {
               /* calculate c beta */
               multiply_matrix(glm_obj->tmp_con[i][0],
                               cur_contrast->contrast_matrix,
                               glm_obj->beta_hat, MATRIX_NEITHER);
               glm_obj->output->values[i][0] = glm_obj->tmp_con[i][0]->values[0][0];
            }
            else {
               fprintf(stderr,"\nError in calculate_contrasts.\n");
               fprintf(stderr,"Beta_hat(Contrast %d) is a matrix, it is scalar only if num_rows is 1.",i+1);
               delete_tmpfiles(&tmpfile_list);
               exit(EXIT_FAILURE);
            }
         }
         else {
            glm_obj->output->values[i][0] = glm_obj->beta_hat->values[cur_contrast->column_num][0];
         }
         break;
         
      default:
         fprintf(stderr,"\nError in calculate_contrasts.\n");
         fprintf(stderr,"Improper choice of out_mode.\n");
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
         break;
      }

   }

   return;

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : calculate_deviance
@INPUT      : pointer to Glm_Object
@OUTPUT     : 
@RETURNS    : deviance of fitted model (-2 log (Likelihood))
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 22, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
double calculate_deviance(Glm_Object *glm_obj) 
{
   double deviance = 0.0;
   double y, mu, tmp;
   int i;

   glm_obj->pearson = 0.0;
   if(glm_obj->control->is_scale == FALSE)
      glm_obj->scale_est = 0.0;

   switch(glm_obj->family) {

   case GAUSSIAN:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         y = glm_obj->response->values[i][0];
         mu = glm_obj->mu->values[i][0];
         tmp = (y-mu)*(y-mu);
         if (glm_obj->missing->values[i][0] == 1.0) {
            deviance += tmp;
            glm_obj->resid->values[i][0] = (y-mu);
         }
         else {
            glm_obj->resid->values[i][0] = INVALID_DATA;
         }
      }

      glm_obj->scale_est = normalize_vector(glm_obj->resid);
      glm_obj->scale_est *= glm_obj->scale_est;
      glm_obj->pearson = glm_obj->scale_est / glm_obj->deg_free;
      break;

   case POISSON:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if (glm_obj->missing->values[i][0] == 1.0) {
            y = glm_obj->response->values[i][0];
            mu = glm_obj->mu->values[i][0];
            tmp = 2.0 * (y * log(y/mu) - (y - mu));
            deviance += tmp;
            glm_obj->resid->values[i][0] = ((y-mu) * glm_obj->g_prime->values[i][0] *
                                            sqrt(glm_obj->weight->values[i][0]));
         }
         else {
            glm_obj->resid->values[i][0] = INVALID_DATA;
         }
      }
      
      if(glm_obj->control->is_scale == FALSE) {
         glm_obj->scale_est = normalize_vector(glm_obj->resid);
         glm_obj->scale_est *= glm_obj->scale_est;
         glm_obj->pearson = glm_obj->scale_est / glm_obj->deg_free;
      }
      break;

   case BINOMIAL:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if (glm_obj->missing->values[i][0] == 1.0) {
            y = glm_obj->response->values[i][0];
            mu = glm_obj->mu->values[i][0];         
            tmp = 2.0 * (y * log(y/mu) + (1.0 - y) * log((1.0-y)/(1.0-mu)));
            deviance += tmp;
            glm_obj->resid->values[i][0] = ((y-mu) * glm_obj->g_prime->values[i][0] *
                                            sqrt(glm_obj->weight->values[i][0]));
         }
         else {
            glm_obj->resid->values[i][0] = INVALID_DATA;
         }
      }

      if(glm_obj->control->is_scale == FALSE) {
         glm_obj->scale_est = normalize_vector(glm_obj->resid);
         glm_obj->scale_est *= glm_obj->scale_est;
         glm_obj->pearson = glm_obj->scale_est / glm_obj->deg_free;
      }
      break;
      
   case GAMMA:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if (glm_obj->missing->values[i][0] == 1.0) {
            y = glm_obj->response->values[i][0];
            mu = glm_obj->mu->values[i][0];
            tmp = -2.0 * (log(y/mu) + (y-mu)/mu);
            deviance += tmp;
            glm_obj->resid->values[i][0] = ((y-mu) * glm_obj->g_prime->values[i][0] *
                                            sqrt(glm_obj->weight->values[i][0]));
         }
         else {
           glm_obj->resid->values[i][0] = INVALID_DATA;
         }
      }
      
      if(glm_obj->control->is_scale == FALSE) {
         glm_obj->scale_est = normalize_vector(glm_obj->resid);
         glm_obj->scale_est *= glm_obj->scale_est;
         glm_obj->pearson = glm_obj->scale_est / glm_obj->deg_free;
      }
      break;

   case INV_GAUSS:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if (glm_obj->missing->values[i][0] == 1.0) {
            y = glm_obj->response->values[i][0];
            mu = glm_obj->mu->values[i][0];
            tmp = (y - mu) * (y - mu) / (mu * mu * y);
            deviance += tmp;
            glm_obj->resid->values[i][0] = ((y-mu) * glm_obj->g_prime->values[i][0] *
                                            sqrt(glm_obj->weight->values[i][0]));
         }
         else {
            glm_obj->resid->values[i][0] = INVALID_DATA;
         }
      }
      
      if(glm_obj->control->is_scale == FALSE) {
         glm_obj->scale_est = normalize_vector(glm_obj->resid);
         glm_obj->scale_est *= glm_obj->scale_est;
         glm_obj->pearson = glm_obj->scale_est / glm_obj->deg_free;
      }
      break;

      /* Note: if family is Quasi, doesn't return deviance, rather it returns
         Pearson's X^2 Statistic */
      /* Sept. 26 1997, Pearson's not quite right yet */

   case QUASI:
      for(i=0; i<glm_obj->response->num_rows; i++) {
         if (glm_obj->missing->values[i][0] == 1.0) {
            y = glm_obj->response->values[i][0];
            mu = glm_obj->mu->values[i][0];
            glm_obj->resid->values[i][0] = ((y-mu) * glm_obj->g_prime->values[i][0] *
                                            sqrt(glm_obj->weight->values[i][0]));
         }
         else {
            glm_obj->resid->values[i][0] = INVALID_DATA;
         }
      }
      
      glm_obj->scale_est = normalize_vector(glm_obj->resid);
      glm_obj->scale_est *= glm_obj->scale_est;
      glm_obj->pearson = glm_obj->scale_est / glm_obj->deg_free;
      break;

   default:
      fprintf(stderr, "\n Error: improper exponential family in calculate_deviance. \n");
      delete_tmpfiles(&tmpfile_list);
      exit(EXIT_FAILURE);
      break;
   }

   if (glm_obj->family != QUASI) {
      if(deviance/glm_obj->pearson <= glm_obj->control->tol) {
         glm_obj->control->status = ZERO_DEV;
      }
   }

   if (glm_obj->family == QUASI)
      return glm_obj->scale_est;
   else
      return deviance / glm_obj->scale_est;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_family_names
@INPUT      : arrays to store family, link and variance function names
@OUTPUT     : 
@RETURNS    : (nothing)
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Oct. 14, 97 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void create_family_names(char ***family_names, char ***link_names, 
                         char ***variance_function_names,
                         char ***corr_struct_names)
{
   char *gauss = "gaussian";
   char *poisson = "poisson";
   char *binomial = "binomial";
   char *inv_gauss = "inverse_gaussian";
   char *gamma = "gamma";
   char *exp_rv = "exponential";
   char *chi_sq = "chi_squared";
   char *quasi = "quasi";
   char *inverse = "inverse";
   char *identity = "identity";
   char *logit = "logit";
   char *probit = "probit";
   char *log = "log";
   char *cloglog = "c-log-log";
   char *sqrt_link = "sqrt link";
   char *inverse_sq = "inv_squared";
   char *mu_one_minus_mu = "mu_one_minus_mu";
   char *mu_cubed = "mu_cubed";
   char *mu_squared = "mu_squared";
   char *constant_var = "constant";
   char *mu = "mu";
   char *diagonal = "diagonal";
   char *auto_reg = "auto_reg";
   int num_corr_struct = 2;
   int num_family = 8;
   int num_link = 8;
   int num_variance_function = 5;

   *family_names = GI_MALLOC(sizeof(char *) * num_family);
   *link_names = GI_MALLOC(sizeof(char *) * num_link);
   *variance_function_names = GI_MALLOC(sizeof(char *) * num_variance_function);
   *corr_struct_names = GI_MALLOC(sizeof(char *) * num_corr_struct);

   (*family_names)[GAUSSIAN] = gauss;
   (*family_names)[GAMMA] = gamma;
   (*family_names)[POISSON] = poisson;
   (*family_names)[BINOMIAL] = binomial;
   (*family_names)[INV_GAUSS] = inv_gauss;
   (*family_names)[QUASI] = quasi;
   (*family_names)[CHI_SQ] = chi_sq;
   (*family_names)[EXPONENTIAL] = exp_rv;

   (*link_names)[LOGIT] = logit;
   (*link_names)[PROBIT] = probit;
   (*link_names)[CLOGLOG]  = cloglog;
   (*link_names)[IDENTITY] = identity;
   (*link_names)[INVERSE] = inverse;
   (*link_names)[INVERSE_SQ] = inverse_sq;
   (*link_names)[LOG] = log;
   (*link_names)[SQRT] = sqrt_link;

   (*variance_function_names)[CONST] = constant_var;
   (*variance_function_names)[MU] = mu;
   (*variance_function_names)[MU_SQUARED] = mu_squared;
   (*variance_function_names)[MU_ONE_MINUS_MU] = mu_one_minus_mu;
   (*variance_function_names)[MU_CUBED] = mu_cubed;

   (*corr_struct_names)[DIAGONAL_CORR] = diagonal;
   (*corr_struct_names)[AUTO_REG] = auto_reg;

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_family
@INPUT      : family_str - string containing information on the exponential
                           family of the GLM
              family_names - arrays containing lists of valid family names
@OUTPUT     : 
@RETURNS    : family type for the GLM
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Nov. 12, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
Family get_family(char *family_str, char **family_names) 
{
   Family family;
   int num_family = 8;
   int i;
   int is_match = FALSE;

   if(family_str != NULL) {
      for(i=0; i<num_family; i++) {
         if (strcmp(family_str, family_names[i]) == 0) {
            family = (Family) i;
            is_match = TRUE;
         }
      }
      if (is_match == FALSE) {
         fprintf(stderr,"\nInvalid family choice, choices are:\n");
         for(i=0; i<num_family; i++) {
            fprintf(stderr,"%s\n", family_names[i]);
         }
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }
   }
   else
      family = GAUSSIAN;

   return family;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_link
@INPUT      : link_str - string containing information on the exponential
                           family of the GLM
              link_names - arrays containing lists of valid family names
@OUTPUT     : 
@RETURNS    : link type for the GLM
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Nov. 12, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
Link get_link(char *link_str, char **link_names) 
{
   Link link;
   int num_link = 8;
   int i;
   int is_match = FALSE;

   if(link_str != NULL) {
      for(i=0; i<num_link; i++) {
         if (strcmp(link_str, link_names[i]) == 0) {
            link = (Link) i;
            is_match = TRUE;
         }
      }
      if (is_match == FALSE) {
         fprintf(stderr,"\nInvalid link choice, choices are:\n");
         for(i=0; i<num_link; i++) {
            fprintf(stderr,"%s\n", link_names[i]);
         }
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }
   }
   else
      link = LINK_DEFAULT;

   return link;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_variance_function
@INPUT      : variance_function_str - string containing information 
                                      on the exponential family of the GLM
              variance_function_names - arrays containing lists of valid
                                        family names
@OUTPUT     : 
@RETURNS    : variance_function type for the GLM
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Nov. 12, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
Variance_Function get_variance_function(char *variance_function_str, 
                                        char **variance_function_names) 
{
   Variance_Function variance_function;
   int num_variance_function = 5;
   int i;
   int is_match = FALSE;

   if(variance_function_str != NULL) {
      for(i=0; i<num_variance_function; i++) {
         if (strcmp(variance_function_str, variance_function_names[i]) == 0) {
            variance_function = (Variance_Function) i;
            is_match = TRUE;
         }
      }
      if (is_match == FALSE) {
         fprintf(stderr,"\nInvalid variance_function choice, choices are:\n");
         for(i=0; i<num_variance_function; i++) {
            fprintf(stderr,"%s\n", variance_function_names[i]);
         }
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }
   }
   else
      variance_function = VARIANCE_DEFAULT;

   return variance_function;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_corr_struct
@INPUT      : corr_struct_str - string containing information on the 
                                exponential family of the GLM
              corr_struct_names - arrays containing lists of valid family names
@OUTPUT     : 
@RETURNS    : correlation structure type for the GLM
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Nov. 12, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
Corr_Struct get_corr_struct(char *corr_struct_str, char **corr_struct_names) 
{
   Corr_Struct corr_struct;
   int num_corr_struct = 2;
   int i;
   int is_match = FALSE;

   if(corr_struct_str != NULL) {
      for(i=0; i<num_corr_struct; i++) {
         if (strcmp(corr_struct_str, corr_struct_names[i]) == 0) {
            corr_struct = (Corr_Struct) i;
            is_match = TRUE;
         }
      }
      if (is_match == FALSE) {
         fprintf(stderr,"\nInvalid correlation structure choice, choices are:\n");
         for(i=0; i<num_corr_struct; i++) {
            fprintf(stderr,"%s\n", corr_struct_names[i]);
         }
         delete_tmpfiles(&tmpfile_list);
         exit(EXIT_FAILURE);
      }
   }
   else
      corr_struct = CORR_DEFAULT;

   return corr_struct;
}

