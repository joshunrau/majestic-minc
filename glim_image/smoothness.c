/* ----------------------------- MNI Header -----------------------------------
@NAME       : smoothness.c
@INPUT      :
@OUTPUT     : (nothing)
@RETURNS    : 
@DESCRIPTION: routines to estimate smoothness of random field for glim_image
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

extern *tmpfile_list;
extern *tmpfile_name;

#include "smoothness.h"
#include "glim.h"

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_lambda_buffer
@INPUT      : num_buffer - num_buffers to create, for glm_obj it is the length
                           of glm_obj->response, later perhaps to be changed to
                           twice the number to account for varying weights
              sizes - size of a volume being used                          
@OUTPUT     : buffer - pointer to buffers needed to calculate the numerical 
                       derivatives for calculating smoothness of field
@RETURNS    : 
@DESCRIPTION: routine to create buffers for calculating smoothness of field
@CREATED    : Nov. 3, 1997, J. Taylor
@MODIFIED   : 
---------------------------------------------------------------------------- */
Lambda *create_lambda_buffer(int num_buffers, int num_test, 
                             int num_dim, int *sizes, int is_gaussian, 
                             int is_avg_fwhm, char *fwhm_simple_file)
{
   int i, j;
   Lambda *lambda_buffer;
   int *coord, *coord_tmp;
   int *element_sizes;
   double *step;
   double *test;

   if (is_avg_fwhm + num_test <= 0) {
      fprintf(stderr,"\nError: local FWHM not needed, no test statistics being output.\n");
      return NULL;
   }

   /* The total size is the size necessary to ensure that all values for the 
      cube over which the derivative will be calculated will be stored for
      each point */

   lambda_buffer = GI_MALLOC(sizeof(*lambda_buffer));

   if(lambda_buffer == NULL) {
      fprintf(stderr,"\nError allocating buffers for smoothness calculation.\n");
      return NULL;
   }

   lambda_buffer->num_dim = num_dim;
   lambda_buffer->num_buffers = num_buffers;
   lambda_buffer->num_test = num_test;
   lambda_buffer->current = 0;
   lambda_buffer->total = 0;
   lambda_buffer->constant = 1.0;
   lambda_buffer->is_gaussian = is_gaussian;
   lambda_buffer->is_avg = is_avg_fwhm;
   lambda_buffer->num_check = 1;

   /* The array check contains the look up points for both checking whether
      all points are in the mask and the values themselves used to calculate */

   element_sizes = GI_MALLOC(sizeof(*element_sizes) * num_dim);
   if(element_sizes == NULL) {
      fprintf(stderr,"\nError allocating element_sizes for smoothness calculation.\n");
      return NULL;
   } 

   for(i=lambda_buffer->num_dim-1; i>=0; i--) {
      element_sizes[i] = 1.0;
      if (i<lambda_buffer->num_dim-1) {
         element_sizes[i] =\
            (sizes[i+1] * element_sizes[i+1]);
      }
   }
         
   lambda_buffer->element_sizes = element_sizes;
   lambda_buffer->ring_size = 0;

   coord = GI_MALLOC(sizeof(*coord) * num_dim);
   if(coord == NULL) {
      fprintf(stderr,"\nError allocating coord for smoothness calculation.\n");
      return NULL;
   } 
   lambda_buffer->coord = coord;

   coord_tmp = GI_MALLOC(sizeof(*coord_tmp) * num_dim);

   for(i=0; i<lambda_buffer->num_dim; i++) {
      coord_tmp[i] = 0;
      lambda_buffer->num_check *= 2;
      lambda_buffer->coord[i] = 1;
      lambda_buffer->ring_size += lambda_buffer->element_sizes[i];
   }

   lambda_buffer->ring_size += 1;
   lambda_buffer->current = lambda_buffer->ring_size - 1;

   lambda_buffer->check = GI_MALLOC(sizeof(*lambda_buffer->check) * 
                                 lambda_buffer->num_check);

   lambda_buffer->index = GI_MALLOC(sizeof(*lambda_buffer->index) * 
                                 lambda_buffer->num_check);

   for(i=0; i<lambda_buffer->num_check; i++) {
      lambda_buffer->check[i] = (lambda_buffer->current -
                                 convert_coord_to_index(coord_tmp,
                                                        lambda_buffer));
      coord_tmp[0]++;
      for(j=0; j<lambda_buffer->num_dim-1; j++) {
         if (coord_tmp[j] == 2) {
            coord_tmp[j+1]++;
            coord_tmp[j] = 0;
         }
      }
   }
   GI_FREE(coord_tmp);
   lambda_buffer->current = -1;

   step = GI_MALLOC(sizeof(*step) * num_dim);
   if(step == NULL) {
      fprintf(stderr,"\nError allocating step for smoothness calculation.\n");
      return NULL;
   } 

   test = GI_MALLOC(sizeof(*test) * lambda_buffer->ring_size);

   lambda_buffer->sizes = sizes;
   lambda_buffer->coord = coord;
   lambda_buffer->step = step;

   for(i=0; i<num_dim; i++) {
      lambda_buffer->sizes[i] = sizes[i];
      lambda_buffer->coord[i] = 0;
   }
   lambda_buffer->coord[2] = -1;

   if (is_gaussian == TRUE)
      lambda_buffer->fwhm_simple = NULL;

   if (num_test > 0) {
      if (is_gaussian == TRUE) 
         lambda_buffer->fwhm_gaussian = create_general_buffer(lambda_buffer,
                                                              TRUE);
      else {
         lambda_buffer->fwhm_gaussian = NULL;
         lambda_buffer->fwhm_general = GI_MALLOC(sizeof(*lambda_buffer->fwhm_general) * (lambda_buffer->num_test));
         for(i=0; i<lambda_buffer->num_test; i++) {
            lambda_buffer->fwhm_general[i]=create_general_buffer(lambda_buffer,
                                                                 FALSE);
         }
         if (fwhm_simple_file != NULL) {
            lambda_buffer->fwhm_simple = create_general_buffer(lambda_buffer,
                                                               TRUE);
            lambda_buffer->fwhm_simple->outfile = fwhm_simple_file;
         }
      }

      /* Above, fwhm_simple is for non-weight corrected FWHM estimates */
   }
   else {
      lambda_buffer->fwhm_gaussian = NULL;
      lambda_buffer->fwhm_general = NULL;
   }

   if(lambda_buffer->is_avg == TRUE) {
      lambda_buffer->fwhm_avg = create_general_buffer(lambda_buffer, TRUE);
   }
   else {
      lambda_buffer->fwhm_avg = NULL;
   }

   lambda_buffer->test = test;

   return lambda_buffer;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_coord_to_index
@INPUT      : coord - vector of coordinates of whose "address" we want to find
              lambda_buffer - pointer to buffers used to calculate smoothness
                              of output fields for glim_image
@OUTPUT     : converts coord values to index values for FWHM calc
@RETURNS    : 
@DESCRIPTION: 
@CREATED    : Apr. 29, 1998 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int convert_coord_to_index(int *coord, Lambda *lambda_buffer) 
{
   int idim;
   int dim2 = -1;
   int dim1 = -1;
   int diff = 0;
   int index;

   /* First check to make sure we've loaded this voxel into buffer 
      As long as dim1, the first dimension where coord[idim] is
      less than lambda_buffer->coord[idim] is less than the first
      one where it's bigger, than we will have "loaded" this
      data into lambda_buffer */

   for(idim = 0; idim < lambda_buffer->num_dim; idim++) {
      if ((coord[idim] < lambda_buffer->coord[idim]) && (dim1 < 0))
         dim1 = idim;
      else if ((coord[idim] > lambda_buffer->coord[idim]) && (dim2 < 0))
         dim2 = idim;
   }

   if ((dim1 > dim2) && (dim2 > 0)){
      fprintf(stderr,"\nError: this voxel hasn't been looped through yet. \ndim1: %d   dim2: %d\n", dim1, dim2);
   }

   if (dim1 == dim2) {
      return lambda_buffer->current;
   }

   /* Next check whether we've already written over the coordinate */

   for(idim = 0; idim < lambda_buffer->num_dim; idim++) {
      diff += ((lambda_buffer->coord[idim] - coord[idim]) *
               lambda_buffer->element_sizes[lambda_buffer->num_dim-(idim+1)]);
   }

   if(diff >= lambda_buffer->ring_size) {
      fprintf(stderr,"\nError: this voxel has been overwritten.\n");
   }

   index = lambda_buffer->current - diff;
   if (index < 0)
      index += lambda_buffer->ring_size;

   return index;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_index_to_coord
@INPUT      : index - vector of coordinates of whose coord we want to find
              lambda_buffer - pointer to buffers used to calculate smoothness
                              of output fields for glim_image
@OUTPUT     : converts coord values to index values for FWHM calc
@RETURNS    : 
@DESCRIPTION: 
@CREATED    : Apr. 29, 1998 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void convert_index_to_coord(int index, int *coord, Lambda *lambda_buffer) 
{
   int idim;
   int diff;

   diff = lambda_buffer->current - index;

   if (diff < 0)
      diff += lambda_buffer->ring_size;

   for(idim = 0; idim<lambda_buffer->num_dim; idim++) {
      coord[idim] = lambda_buffer->coord[idim];

      while(diff > 0) {
         diff -= lambda_buffer->element_sizes[idim];
         coord[idim]--;
      }
      
      diff += lambda_buffer->element_sizes[idim];
      coord[idim]++;
   }
      
   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_general_buffer
@INPUT      : lambda_buffer - pointer to buffers used to calculate smoothness
                            of output fields for glim_image
@OUTPUT     : nothing
@RETURNS    : 
@DESCRIPTION: routine to allocate space for Gaussian buffers for FWHM
@CREATED    : Nov. 3, 1997, J. Taylor
@MODIFIED   : 
---------------------------------------------------------------------------- */
Fwhm_Info *create_general_buffer(Lambda *lambda_buffer, int is_gaussian)
{
   int ibuff, i;
   Fwhm_Info *fwhm_info;

   fwhm_info = GI_MALLOC(sizeof(*fwhm_info));

   fwhm_info->lambda = create_matrix(lambda_buffer->num_dim,
                                         lambda_buffer->num_dim);

   fwhm_info->is_gaussian = is_gaussian;
   fwhm_info->num_calc = 0.0;
   fwhm_info->out_id = -1;
   
   if (is_gaussian == TRUE)
      fwhm_info->lambda_pool = create_matrix(lambda_buffer->num_dim,
                                             lambda_buffer->num_dim);
   else
      fwhm_info->lambda_pool = NULL;

   if(fwhm_info->lambda == NULL) {
      fprintf(stderr,"\nError creating lambda in create_general_buffer.\n");
      return NULL;
   }

   fwhm_info->deriv = GI_MALLOC(sizeof(*fwhm_info->deriv) * 
                                 lambda_buffer->num_dim);
   if(fwhm_info->deriv == NULL) {
      fprintf(stderr,"\nError allocating deriv in create_general_buffer.\n");
      return NULL;
   }

   fwhm_info->data = GI_MALLOC(sizeof(double) * 
                                 (lambda_buffer->num_buffers));

   fwhm_info->ring_size = lambda_buffer->ring_size;

   for(ibuff=0; ibuff<lambda_buffer->num_buffers; ibuff++) {
      fwhm_info->data[ibuff] = GI_MALLOC(sizeof(**(fwhm_info->data)) *
                                           fwhm_info->ring_size);
      if(fwhm_info->data[ibuff] == NULL) {
         fprintf(stderr,"\nError allocating data in create_general_buffer. ibuff: %d\n", ibuff);
         return NULL;
      }
   }

   return fwhm_info;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_lambda_buffer
@INPUT      : lambda_buffer - pointer to buffers used to calculate smoothness
                            of output fields for glim_image
@OUTPUT     : nothing
@RETURNS    : 
@DESCRIPTION: routine to delete buffers used for calculating
              smoothness of field
@CREATED    : Nov. 3, 1997, J. Taylor
@MODIFIED   : 
---------------------------------------------------------------------------- */
int delete_lambda_buffer(Lambda *lambda_buffer)
{
   int ibuff, itest;
   Fwhm_Info *fwhm_general;

   GI_FREE(lambda_buffer->sizes);
   GI_FREE(lambda_buffer->coord);
   GI_FREE(lambda_buffer->step);
   GI_FREE(lambda_buffer->test);

   if(lambda_buffer->is_gaussian == FALSE) {
      for(itest=0; itest<lambda_buffer->num_test; itest++) {
         fwhm_general = lambda_buffer->fwhm_general[itest];

         for(ibuff=0; ibuff<lambda_buffer->num_buffers; ibuff++) {
            GI_FREE(fwhm_general->data[ibuff]);
         }
         GI_FREE(fwhm_general->data);

         delete_matrix(fwhm_general->lambda);
         GI_FREE(fwhm_general->deriv);
      }
      if (lambda_buffer->fwhm_simple != NULL) {
         GI_FREE(lambda_buffer->fwhm_simple->deriv);

         for(ibuff=0; ibuff<lambda_buffer->num_buffers; ibuff++) {
            GI_FREE(lambda_buffer->fwhm_simple->data[ibuff]);
         }
         GI_FREE(lambda_buffer->fwhm_simple->data);
         delete_matrix(lambda_buffer->fwhm_simple->lambda);
      }
   }
   else if (lambda_buffer->fwhm_gaussian != NULL) {
      GI_FREE(lambda_buffer->fwhm_gaussian->deriv);

      for(ibuff=0; ibuff<lambda_buffer->num_buffers; ibuff++) {
         GI_FREE(lambda_buffer->fwhm_gaussian->data[ibuff]);
      }
      GI_FREE(lambda_buffer->fwhm_gaussian->data);
      delete_matrix(lambda_buffer->fwhm_gaussian->lambda);
   }

   if(lambda_buffer->is_avg == TRUE) {
      GI_FREE(lambda_buffer->fwhm_avg->deriv);

      for(ibuff=0; ibuff<lambda_buffer->num_buffers; ibuff++) {
         GI_FREE(lambda_buffer->fwhm_avg->data[ibuff]);
      }
      GI_FREE(lambda_buffer->fwhm_avg->data);
      delete_matrix(lambda_buffer->fwhm_avg->lambda);
   }

   return TRUE;

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : calculate_smoothness
@INPUT      : lambda_buffer - ptr to ring buffers where values for calculation
                              of derivatives are stored
@OUTPUT     : buffer - pointer to buffers needed to calculate the numerical 
                       derivatives for calculating smoothness of field
@RETURNS    : 
@DESCRIPTION: routine to create buffers for calculating smoothness of field
@CREATED    : Nov. 3, 1997, J. Taylor
@MODIFIED   : 
---------------------------------------------------------------------------- */
void calculate_smoothness(Lambda *lambda_buffer)
{
   int i, itest;
   Fwhm_Info *fwhm_general;

   /* find positions for values to be used for calculating derivatives 
      and check to see if all values of mask are VIO_OK                    */

   lambda_buffer->test_voxel = 1.0;

   for(i=0; i<lambda_buffer->num_check-1; i++) {
      lambda_buffer->index[i] = (lambda_buffer->current -
                                 lambda_buffer->check[i]);
      if(lambda_buffer->index[i] < 0)
         lambda_buffer->index[i] += lambda_buffer->ring_size;
      if (lambda_buffer->test[lambda_buffer->index[i]] != INVALID_DATA) {
         lambda_buffer->test_voxel *= lambda_buffer->test[lambda_buffer->index[i]];
      }
      else
         lambda_buffer->test_voxel = 0.0;
   }

   lambda_buffer->index[lambda_buffer->num_check-1] = lambda_buffer->current;

   if(!((lambda_buffer->coord[0]) > 0) && (lambda_buffer->coord[1] > 0) &&
      (lambda_buffer->coord[2]) > 0) {
      lambda_buffer->test_voxel = 0.0;
   }

   /* take spatial derivatives, centered around the cube                    */
   /* must be changed if higher dimensions are going to be used,            */
   /* i.e. searching scale space                                            */

   /* Loop through buffers */

   if(lambda_buffer->test_voxel == 1.0) {

      lambda_buffer->num_calc++;

      if ((lambda_buffer->is_gaussian == FALSE) &&
          (lambda_buffer->fwhm_general != NULL)){
         for(itest=0; itest<lambda_buffer->num_test; itest++) {
            fwhm_general = lambda_buffer->fwhm_general[itest];
            calculate_smoothness_ind(fwhm_general, lambda_buffer);
         }
         if (lambda_buffer->fwhm_simple != NULL) {
            if (lambda_buffer->fwhm_simple->outfile != NULL) {
               calculate_smoothness_ind(lambda_buffer->fwhm_simple,
                                        lambda_buffer);
            }
         }
      }
      else if (lambda_buffer->fwhm_gaussian != NULL){
         calculate_smoothness_ind(lambda_buffer->fwhm_gaussian,
                                  lambda_buffer);
      }

      if(lambda_buffer->is_avg == TRUE) {
         calculate_smoothness_ind(lambda_buffer->fwhm_avg,
                                  lambda_buffer);
      }
   }
   else {
      if ((lambda_buffer->is_gaussian == FALSE) && 
          (lambda_buffer->fwhm_general != NULL)){
         for(itest=0; itest<lambda_buffer->num_test; itest++) {
            fwhm_general = lambda_buffer->fwhm_general[itest];
            fwhm_general->fwhm = INVALID_DATA;
            fwhm_general->smoothness = INVALID_DATA;
         }
      }
      else if (lambda_buffer->fwhm_gaussian != NULL) {
         lambda_buffer->fwhm_gaussian->fwhm = INVALID_DATA;
         lambda_buffer->fwhm_gaussian->smoothness = INVALID_DATA;
      }

      if(lambda_buffer->is_avg == TRUE) {
         lambda_buffer->fwhm_avg->fwhm = INVALID_DATA;
         lambda_buffer->fwhm_avg->smoothness = INVALID_DATA;
      }
   }

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : calculate_smoothness_ind
@INPUT      : fwhm_info - pointer to data for calculating smoothness
                             in GLM model fit
              lambda_buffer - pointer to residuals and info about smoothness
@OUTPUT     : fwhm_gaussian->fwhm, smoothness, int_smoothness
@RETURNS    : nothing
@DESCRIPTION: Calculates the sample covariance matrix at a given voxel
@CREATED    : Apr. 20, 1997, J. Taylor
@MODIFIED   : 
---------------------------------------------------------------------------- */
void calculate_smoothness_ind(Fwhm_Info *fwhm_info, Lambda *lambda_buffer)
{
   int ibuff, i, j;
   double *data;
   double avg, avg_sq;

   avg_sq = 0.0;

   for(i=0; i<lambda_buffer->num_dim; i++) {
      for(j=0; j<=i; j++) {
         fwhm_info->lambda->values[i][j] = 0.0;
      }
   }

   for(ibuff=0; ibuff<lambda_buffer->num_buffers; ibuff++) {

      data = fwhm_info->data[ibuff];

      if(lambda_buffer->num_dim == 2) {
         fwhm_info->deriv[0] = 0.5* ((data[lambda_buffer->index[3]] +
                                 data[lambda_buffer->index[2]]) -
                                (data[lambda_buffer->index[1]] +
                                 data[lambda_buffer->index[0]]));

         fwhm_info->deriv[1] = 0.5 * ((data[lambda_buffer->index[3]] +
                                 data[lambda_buffer->index[1]]) -
                                (data[lambda_buffer->index[2]] +
                                 data[lambda_buffer->index[0]]));
      }
      else if(lambda_buffer->num_dim == 3) {

         fwhm_info->deriv[0] = 0.25 * ((data[lambda_buffer->index[7]] +
                                 data[lambda_buffer->index[6]] +
                                 data[lambda_buffer->index[5]] + 
                                 data[lambda_buffer->index[4]]) -
                                (data[lambda_buffer->index[3]] +
                                 data[lambda_buffer->index[2]] + 
                                 data[lambda_buffer->index[1]] +
                                 data[lambda_buffer->index[0]]));

         fwhm_info->deriv[1] = 0.25 * ((data[lambda_buffer->index[7]] +
                                 data[lambda_buffer->index[6]] +
                                 data[lambda_buffer->index[3]] + 
                                 data[lambda_buffer->index[2]]) -
                                (data[lambda_buffer->index[5]] +
                                 data[lambda_buffer->index[4]] + 
                                 data[lambda_buffer->index[1]] +
                                 data[lambda_buffer->index[0]]));
      
         fwhm_info->deriv[2] = 0.25 * ((data[lambda_buffer->index[7]] +
                                 data[lambda_buffer->index[5]] +
                                 data[lambda_buffer->index[3]] + 
                                 data[lambda_buffer->index[1]]) -
                                (data[lambda_buffer->index[6]] +
                                 data[lambda_buffer->index[4]] + 
                                 data[lambda_buffer->index[2]] +
                                 data[lambda_buffer->index[0]]));
      }
      else if (lambda_buffer->num_dim == 4) {
         fwhm_info->deriv[0] = 0.125 * ((data[lambda_buffer->index[15]] +
                                 data[lambda_buffer->index[14]] +
                                 data[lambda_buffer->index[13]] + 
                                 data[lambda_buffer->index[12]] +
                                 data[lambda_buffer->index[11]] +
                                 data[lambda_buffer->index[10]] + 
                                 data[lambda_buffer->index[9]] +
                                 data[lambda_buffer->index[8]]) -
                                (data[lambda_buffer->index[7]] +
                                 data[lambda_buffer->index[6]] +
                                 data[lambda_buffer->index[5]] + 
                                 data[lambda_buffer->index[4]] +
                                 data[lambda_buffer->index[3]] +
                                 data[lambda_buffer->index[2]] + 
                                 data[lambda_buffer->index[1]] +
                                 data[lambda_buffer->index[0]]));

         fwhm_info->deriv[1] = 0.125 * ((data[lambda_buffer->index[15]] +
                                 data[lambda_buffer->index[14]] +
                                 data[lambda_buffer->index[13]] + 
                                 data[lambda_buffer->index[12]] +
                                 data[lambda_buffer->index[7]] +
                                 data[lambda_buffer->index[6]] + 
                                 data[lambda_buffer->index[5]] +
                                 data[lambda_buffer->index[4]]) -
                                (data[lambda_buffer->index[11]] +
                                 data[lambda_buffer->index[10]] +
                                 data[lambda_buffer->index[9]] + 
                                 data[lambda_buffer->index[8]] +
                                 data[lambda_buffer->index[3]] +
                                 data[lambda_buffer->index[2]] + 
                                 data[lambda_buffer->index[1]] +
                                 data[lambda_buffer->index[0]]));

         fwhm_info->deriv[2] = 0.125 * ((data[lambda_buffer->index[15]] +
                                 data[lambda_buffer->index[14]] +
                                 data[lambda_buffer->index[11]] + 
                                 data[lambda_buffer->index[10]] +
                                 data[lambda_buffer->index[7]] +
                                 data[lambda_buffer->index[6]] + 
                                 data[lambda_buffer->index[3]] +
                                 data[lambda_buffer->index[2]]) -
                                (data[lambda_buffer->index[13]] +
                                 data[lambda_buffer->index[12]] +
                                 data[lambda_buffer->index[9]] + 
                                 data[lambda_buffer->index[8]] +
                                 data[lambda_buffer->index[5]] +
                                 data[lambda_buffer->index[4]] + 
                                 data[lambda_buffer->index[1]] +
                                 data[lambda_buffer->index[0]]));

         fwhm_info->deriv[3] = 0.125 * ((data[lambda_buffer->index[15]] +
                                 data[lambda_buffer->index[13]] +
                                 data[lambda_buffer->index[11]] + 
                                 data[lambda_buffer->index[9]] +
                                 data[lambda_buffer->index[7]] +
                                 data[lambda_buffer->index[5]] + 
                                 data[lambda_buffer->index[3]] +
                                 data[lambda_buffer->index[1]]) -
                                (data[lambda_buffer->index[14]] +
                                 data[lambda_buffer->index[12]] +
                                 data[lambda_buffer->index[10]] + 
                                 data[lambda_buffer->index[8]] +
                                 data[lambda_buffer->index[6]] +
                                 data[lambda_buffer->index[4]] + 
                                 data[lambda_buffer->index[2]] +
                                 data[lambda_buffer->index[0]]));

      }

      if (fwhm_info->is_gaussian == FALSE) {
         avg = 0.0;
         for(i=0; i<lambda_buffer->num_check; i++) {
            avg += data[lambda_buffer->index[i]];
         }
         avg /= lambda_buffer->num_check;            
      
         avg_sq += avg * avg;
      }

      for(i=0; i<lambda_buffer->num_dim; i++) {
         for(j=0; j<=i; j++) {

            fwhm_info->lambda->values[i][j] +=\
               fwhm_info->deriv[i] * fwhm_info->deriv[j];

            if (fwhm_info->lambda_pool != NULL) {
               fwhm_info->lambda_pool->values[i][j] +=\
                  fwhm_info->deriv[i] * fwhm_info->deriv[j];
            }
         }
      }

   }

   if (fwhm_info->is_gaussian == FALSE) {
      for(i=0; i<lambda_buffer->num_dim; i++) {
         for(j=0; j<=i; j++) {
            fwhm_info->lambda->values[i][j] /= avg_sq;
         }
      }
   }

   for(i=0; i<lambda_buffer->num_dim; i++) {
      for(j=lambda_buffer->num_dim-1; j>i; j--) {
         fwhm_info->lambda->values[i][j] = fwhm_info->lambda->values[j][i];
      }
   }

   fwhm_info->smoothness = sqrt(determinant(fwhm_info->lambda));

   if(fwhm_info->smoothness > 0.0)
      fwhm_info->int_lambda += fwhm_info->smoothness;

   fwhm_info->smoothness *= fwhm_info->constant;

   if(fwhm_info->smoothness > 0.0) {
      fwhm_info->fwhm = pow(fwhm_info->smoothness, 1.0/lambda_buffer->num_dim);
   }
   else
      fwhm_info->fwhm = INVALID_DATA;

   if(fwhm_info->fwhm != INVALID_DATA) {
      fwhm_info->fwhm = (lambda_buffer->constant / fwhm_info->fwhm);
   }
   else
      fwhm_info->fwhm = INVALID_DATA;

   if (fwhm_info->fwhm != INVALID_DATA)
      fwhm_info->num_calc++;

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : average_smoothness
@INPUT      : tmp_lambda_file - tmp file with unaveraged smoothness values
@OUTPUT     : lambda_file - file where final values are to be placed
@RETURNS    : nothing
@DESCRIPTION: routine to create buffers for calculating smoothness of field
@CREATED    : Nov. 3, 1997, J. Taylor
@MODIFIED   : 
---------------------------------------------------------------------------- */
int average_smoothness(char *tmp_lambda_file, char *lambda_file, double scalar)
{
   VIO_Volume volume;
   nc_type in_nc_type;
   int in_signed_flag = FALSE;
   volume_input_struct volume_input;
   int v1, v2, v3;
   int num_dim;
   int sizes[VIO_MAX_DIMENSIONS];
   double avg, values[8];
   int i;

   if(start_volume_input(tmp_lambda_file, 3, NULL, NC_UNSPECIFIED, FALSE,
                         0.0, 0.0, TRUE, &volume,
                         (minc_input_options *) NULL, &volume_input) != VIO_OK){
      fprintf(stderr,"\nError opening %s.\n", tmp_lambda_file); 
      exit(EXIT_FAILURE);
   }

   in_nc_type = get_volume_nc_data_type(volume, &in_signed_flag);
   get_volume_sizes(volume, sizes);
   num_dim = get_volume_n_dimensions(volume);

   cancel_volume_input(volume, &volume_input);

   if( input_volume(tmp_lambda_file, 3, NULL, NC_FLOAT, FALSE, 0.0, 0.0,
                    TRUE, &volume, (minc_input_options *) NULL) != VIO_OK) {
      (void) fprintf(stderr, "Error reading file \"%s\"\n", tmp_lambda_file);
      return FALSE;
   }

   for(v1=0; v1<sizes[0]; v1++) {
      for(v2=0; v2<sizes[1]; v2++) {
         for(v3=0; v3<sizes[2]; v3++) {

            avg = 0.0;

            if((v1<sizes[0]-1) && (v2<sizes[1]-1) && (v3<sizes[2])-1) {

               values[0] = get_volume_real_value(volume, v1+1, v2, v3, 0, 0); 
               values[1] = get_volume_real_value(volume, v1, v2+1, v3, 0, 0);
               values[2] = get_volume_real_value(volume, v1, v2, v3+1, 0, 0);
               values[3] = get_volume_real_value(volume,v1+1, v2+1, v3, 0, 0);
               values[4] = get_volume_real_value(volume, v1+1,v2, v3+1, 0, 0);
               values[5] = get_volume_real_value(volume, v1, v2+1,v3+1, 0, 0);
               values[6] = get_volume_real_value(volume, v1, v2, v3, 0, 0);
               values[7] = get_volume_real_value(volume, v1+1, v2+1,v3+1,0,0);
            
               for(i=0; i<8; i++) {
                  if(values[i] != INVALID_DATA)
                     avg += values[i];
                  else {
                     avg = INVALID_DATA;
                     i = 8;
                  }
               }
            }
            else
               avg = INVALID_DATA;

            if (avg == 0.0)
               avg = INVALID_DATA;
            if (avg != INVALID_DATA)
                avg = scalar * 8.0 / pow(avg, (1.0 / num_dim));

            set_volume_real_value(volume, v1, v2, v3, 0, 0, avg);

         }
      }
   }

   if( output_modified_volume(lambda_file, in_nc_type, in_signed_flag,
                     0.0, 0.0, volume, tmp_lambda_file, NULL,
                     (minc_output_options *) NULL) != VIO_OK) {
      (void) fprintf(stderr, "Error writing file \"%s\"\n", lambda_file);
      return FALSE;
   }

   return TRUE;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : pvalue_t
@INPUT      : t - value of t-statistic 
              int_fwhm - integrated fwhm, over search volume
              df - degrees of freedom of t-statistic
@OUTPUT     : nothing
@RETURNS    : prob - probability of exceedence of t
@DESCRIPTION: routine to calculate probability of exceedence for t field
@CREATED    : Dec. 17, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
double pvalue_t(double t, double int_fwhm, int df)
{
   double value;

   value = pow((1.0 + t*t/(df)),-0.5*(df-1)) * ((((df-1.0)*t*t)/df) - 1.0);
   value = value * int_fwhm/(2.0 * 3.141593 * 3.141593);

   return value;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : thresh_t
@INPUT      : pvalue - pvalue by which the threshold is determined
              int_fwhm - integrated fwhm, over search volume
              df - degrees of freedom of t-statistic
@OUTPUT     : nothing
@RETURNS    : prob - probability of exceedence of t
@DESCRIPTION: routine to calculate probability of exceedence for t field
@CREATED    : Dec. 17, 1997 (J. Taylor)
@MODIFIED   : 
---------------------------------------------------------------------------- */
double thresh_t(double pvalue, double int_fwhm, int df)
{
   double a,b,c;
   double p;

   a = 1.644854;
   b = 1000.0;

   while ((b-a) > 1.0e-4) {
      c = 0.5*(a+b);

      p = pvalue_t(c, int_fwhm, df);

      if (p < pvalue)
         b = c;

      if (p >= pvalue)
         a = c;

   }

   return c;

}
   
