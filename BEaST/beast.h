/*  beast.h
 *
 *  Copyright 2011  Simon Fristed Eskildsen, Vladimir Fonov,
 *   	      	    Pierrick Coupé, Jose V. Manjon
 *
 *  This file is part of mincbeast.
 *
 *  mincbeast is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  mincbeast is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with mincbeast.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  For questions and feedback, please contact:
 *  Simon Fristed Eskildsen <eskild@gmail.com> 
 */


#ifndef BEAST_H
#define BEAST_H

#include <stdio.h>

#ifdef HAVE_MINC
#ifndef DEF_VOLUME_IO
#include "mincio.h"
#endif
#endif //HAVE_MINC

#include "basic.h"

#define MAXLIBSIZE 1000
#define FILENAMELENGTH 255

#define VOXELSIZEMAX 4
#define VOXELSIZEMIN 1
#define PATCHSIZEMAX 10
#define PATCHSIZEMIN 0
#define SEARCHAREAMAX 10
#define SEARCHAREAMIN 0
#define ALPHAMAX 1.0
#define ALPHAMIN 0.0
#define BETAMAX 5.0
#define BETAMIN 0.0
#define THRESHOLDMAX 1.0
#define THRESHOLDMIN 0.0


extern const char LICENSE[];
extern const char REFERENCE[];

typedef struct {
  int index;
  float ssd;
} ssd_t;

typedef struct {
  int index;
  double ssd;
} ssd_td;


typedef struct { /*beast algorithm parameters*/
  int voxelsize;
  int patchsize;
  int searcharea;
  double alpha;
  double beta;
  double threshold;
  int selectionsize;
} beast_conf;



typedef struct  { /*micnbeast command line options*/
  char *input_file;
  char *output_file;
  char *libdir;
  char *history_label;
  
  double lambda1;
  double lambda2;
  int    sparse_mode;
  
  int sparse_stride;

  VIO_BOOL outputprob;
  VIO_BOOL flipimages;
  VIO_BOOL load_moments;
  VIO_BOOL fill_output;
  VIO_BOOL verbose;
  VIO_BOOL medianfilter;
  VIO_BOOL patchfilter;
  VIO_BOOL abspath;
  VIO_BOOL same_res;
  VIO_BOOL clobber   ;
  VIO_BOOL nomask    ;
  VIO_BOOL nopositive;
  VIO_BOOL use_sparse;
  VIO_BOOL v2;
  VIO_BOOL use_double;

  int voxelsize;
  int sizepatch;
  int searcharea;
  double alpha;
  double beta;
  double threshold;
  double kappa_limit;
  int selectionsize;
  
  char *positive_file;
  char *selection_file;
  char *count_file;
  char *conf_file;
  char *mask_file;
  char *library_prefix;
  char *compare_file;
  
} beast_options;


int fgetline(FILE *fp, char line[], int max);

int median_filter(float *volume, int *sizes, int filtersize);

int trilinear_interpolant(float *volume, int *sizes, point3D coord, float *result);
int resize_volume(float *input, int *sizes, int *sizes2, float *result);
int resize_trilinear(float *input, int *sizes, int *sizes2, float *result);

void cp_volume(float *data, float *copy, int *sizes);

int flip_data(float *data, float *result, int *sizes);
int combine_maps(float *data, float *map, float *mask, int *sizes);

int down_sample(float *subject, float *result, int factor, int *sizes);
int up_sample(float *subject, float *result, int factor, int *sizes, int *targetsizes);

int threshold_data(float *data,const int *sizes, float threshold);

int add_mask_data(float *data1, float *mask, int *sizes);
int wipe_data(float *data1, int *sizes, float value);
int update_mask(float *subject, float *mask, float *segmented, int *sizes, float min, float max);

int flood_fill_float(float *data, float *output, int *sizes, 
                     int sx, int sy, int sz, float fill_value, int connectivity);

int pre_selection(const float *subject, const float *mask, 
                  char **images,        const int *sizes, 
                  int librarysize,      int num_selected, 
                  int *selection, 
                  const char *outfile, 
                  VIO_BOOL verbose);

int pre_selection_double(const float *subject, const float *mask, 
                  char **images, const int *sizes, 
                  int librarysize, int num_selected, int *selection, 
                  const char *outfile, VIO_BOOL verbose);

int read_configuration(const char *filename, beast_conf *conf);
int read_list(const char *filename, char **list,const char *basedir);

image_metadata * read_volume(const char *filename, float **data,int *sizes);
int write_volume_generic(const char *filename,const float *data,const image_metadata *meta,VIO_BOOL binary_mask );

char* create_minc_timestamp(int argc,char *argv[]);

int get_arguments(int argc, char  *argv[] , beast_options * _options);
void cleanup_arguments(beast_options * _options);

int mincbeast_v1(beast_options * _options);
int mincbeast_v2(beast_options * _options);

double dice_kappa(const float *  I1,const float *  I2,const int *  sizes);

#endif

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
