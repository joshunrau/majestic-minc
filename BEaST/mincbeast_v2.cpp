/*  mincbeast_v1.c
 *
 *  Copyright 2011  Simon Fristed Eskildsen, Vladimir Fonov,
 *                  Pierrick Coupé, Jose V. Manjon
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif //HAVE_CONFIG_H

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#include "array_alloc.h"
#include "nlmseg.h"
#include "beast.h"
#include "label.h"


int mincbeast_v2(beast_options * _options)
{
  char imagelist[FILENAMELENGTH], masklist[FILENAMELENGTH],meanlist[FILENAMELENGTH],varlist[FILENAMELENGTH];
  char ***images, ***masks, ***means, ***vars;
  int num_images,i,sizes[3][5],tmpsizes[5],volumesize,*selection,steps=3,filled=0;
  float *imagedata,*maskdata,*meandata,*vardata,**subject,**mask,**positivemask=NULL,**segsubject,**patchcount,**filtered;
  float max,min;
  float **segmented;
  float *tempdata;
  
  int scale,scaledvolumesize,scales[3] = {1,2,4};
  int masksize=0,initialscale,targetscale,scalesteps;
  beast_conf input_conf[3],configuration[3];
  image_metadata **meta;
  image_metadata *mask_meta;
  image_metadata *temp_meta;
  int targetvoxelsize=1;
  
  char *selection_file=NULL;
  char *count_file=NULL;
  char *conf_file=NULL;

  time_t timer;
  
  /* Get the time, overwriting newline */
  timer = time(NULL);

  if (_options->mask_file==NULL) {
    _options->mask_file=(char*)malloc((strlen(_options->libdir)+20));
    sprintf(_options->mask_file,"%s/margin_mask.mnc",_options->libdir);
  }
  
  if ((!_options->nopositive) && (_options->positive_file==NULL)) {
    _options->positive_file=(char*)malloc((strlen(_options->libdir)+30));
    sprintf(_options->positive_file,"%s/intersection_mask.mnc",_options->libdir);
  }

  if(!_options->clobber) {
    if(!access(_options->output_file,F_OK)) {
      fprintf(stderr,"ERROR! File exists: %s , run with -clobber\n",_options->output_file);
      return STATUS_ERR;
    }
    if(_options->count_file && !access(_options->count_file,F_OK)) {
      fprintf(stderr,"ERROR! File exists: %s , run with -clobber\n",_options->count_file);
      return STATUS_ERR;
    }
  }

  if ((_options->voxelsize>4) || (_options->voxelsize<1) || (_options->voxelsize==3)) {
    fprintf(stderr,"ERROR! Initial voxel size must be either 4, 2, or 1\n");
    return STATUS_ERR;
  }

  meta = (image_metadata **)malloc(3*sizeof(image_metadata*));

  meta[0] = read_volume(_options->input_file, &tempdata, sizes[0]);
  if (meta[0] == NULL) {
    fprintf(stderr,"ERROR! Image not read: %s\n",_options->input_file);
    return STATUS_ERR;
  }
  volumesize=sizes[0][0]*sizes[0][1]*sizes[0][2];

  subject = alloc_2d_float(3, volumesize); /*VF:memory waste....*/
  cp_volume(tempdata, subject[0], sizes[0]);
  free(tempdata);

  if ((temp_meta=read_volume(_options->mask_file, &tempdata, tmpsizes)) == NULL) {
    fprintf(stderr,"ERROR! Image not read: %s\n",_options->mask_file);
    return STATUS_ERR;
  }
  free_meta(temp_meta);

  if ((tmpsizes[0]!=sizes[0][0]) || (tmpsizes[1]!=sizes[0][1]) || (tmpsizes[2]!=sizes[0][2])) {
    fprintf(stderr,"ERROR! Mask dimension does not match image dimension!\n");
    return STATUS_ERR;
  }
  mask = alloc_2d_float(3, volumesize);
  cp_volume(tempdata, mask[0], sizes[0]);
  free(tempdata);

  if (_options->nomask) {
    /* option for no segmentation mask - set the mask to all ones */
    wipe_data(mask[0],sizes[0],1.0);
  }

  if (_options->positive_file!=NULL) {
    image_metadata *positive_meta;
    if ((positive_meta=read_volume(_options->positive_file, &tempdata, tmpsizes)) == NULL) {
      fprintf(stderr,"ERROR! Image not read: %s\n", _options->positive_file);
      return STATUS_ERR;
    }
    if ((tmpsizes[0]!=sizes[0][0]) || (tmpsizes[1]!=sizes[0][1]) || (tmpsizes[2]!=sizes[0][2])) {
      fprintf(stderr, "ERROR! Positive mask dimension does not match image dimension!\n");
      return STATUS_ERR;
    }
    positivemask = alloc_2d_float(3,volumesize);
    cp_volume(tempdata, positivemask[0], sizes[0]);
    free(tempdata);
    free_meta(positive_meta);

    down_sample(positivemask[0], positivemask[1], 2, sizes[0]);
    down_sample(positivemask[0], positivemask[2], 4, sizes[0]);
  }

  segmented = alloc_2d_float(3,volumesize);

  /* downsample the subject and mask */
  down_sample(subject[0], subject[1], 2, sizes[0]);
  down_sample(subject[0], subject[2], 4, sizes[0]);
  down_sample(mask[0],    mask[1],    2, sizes[0]);
  down_sample(mask[0],    mask[2],    4, sizes[0]);
  for(i=1;i<3;i++) {
    int j;
    for(j=0;j<3;j++)
      sizes[i][j]=sizes[i-1][j]/2;
  }

  /* make the downsampled masks crisp */
  threshold_data( mask[1], sizes[1], 0.5);
  threshold_data( mask[2], sizes[2], 0.5);
  

  /* populate the entire configuration table for compatibility reasons */
  for (i=0; i<3; i++) {
    configuration[i].voxelsize = _options->voxelsize;
    configuration[i].patchsize = _options->sizepatch;
    configuration[i].searcharea = _options->searcharea;
    configuration[i].alpha      = _options->alpha;
    configuration[i].beta       = _options->beta;
    configuration[i].threshold  = _options->threshold;
    configuration[i].selectionsize = _options->selectionsize;
  }


  if (_options->conf_file != NULL) {
    steps=read_configuration(_options->conf_file, input_conf);
    if (steps==STATUS_ERR) {
      fprintf(stderr,"Error in configuration file. Values outside limits.\n");
      return STATUS_ERR;
    }
    _options->selectionsize=configuration[0].selectionsize;
    
    initialscale=-1;
    targetscale=4;
    for (i=0; i<steps; i++) {
      scale=(int)(input_conf[i].voxelsize/2);
      configuration[scale].voxelsize   = input_conf[i].voxelsize;
      configuration[scale].patchsize   = input_conf[i].patchsize;
      configuration[scale].searcharea  = input_conf[i].searcharea;
      configuration[scale].alpha       = input_conf[i].alpha;
      configuration[scale].beta        = input_conf[i].beta;
      configuration[scale].threshold   = input_conf[i].threshold;
      configuration[scale].selectionsize=input_conf[i].selectionsize;
      
      if (scale>initialscale)
        initialscale=scale;
      if (scale<targetscale)
        targetscale=scale;
    }
  } else {
    /* if no configuration file, apply user settings for single scale */
    targetscale=initialscale=(int)(_options->voxelsize/2);
  }

  scalesteps=initialscale-targetscale+1;

  fprintf(stderr,"%d scale steps:\n",scalesteps);

  for (i=initialscale; i>=targetscale; i--) {
    fprintf(stderr,"%d %d %d %4.2lf %4.2lf %4.2lf %d\n",
            configuration[i].voxelsize,                    configuration[i].patchsize, configuration[i].searcharea,
            configuration[i].alpha, configuration[i].beta, configuration[i].threshold, configuration[i].selectionsize);
  }

  images = alloc_3d_char(3, MAXLIBSIZE, FILENAMELENGTH);
  masks =  alloc_3d_char(3, MAXLIBSIZE, FILENAMELENGTH);
  means =  alloc_3d_char(3, MAXLIBSIZE, FILENAMELENGTH);
  vars =   alloc_3d_char(3, MAXLIBSIZE, FILENAMELENGTH);

  sprintf(imagelist,"%s/%s.stx.%dmm", _options->libdir, _options->library_prefix, 1);
  sprintf(masklist,"%s/%s.masks.%dmm",_options->libdir, _options->library_prefix, 1);
  sprintf(meanlist,"%s/%s.means.%dmm",_options->libdir, _options->library_prefix, 1);
  sprintf(varlist,"%s/%s.vars.%dmm",  _options->libdir, _options->library_prefix, 1);
  
  num_images=read_list( imagelist, images[0], _options->abspath?"":_options->libdir );
  if (read_list( masklist, masks[0], _options->abspath?"":_options->libdir )!=num_images) {
      fprintf(stderr,"ERROR! Number of images and masks does not match!\n");
      return STATUS_ERR;
  }
  if ( num_images< _options->selectionsize ) {
    fprintf(stderr,"ERROR! Cannot select more images than in the library!\n\tlibrary images: %d\n\tselection: %d\n",num_images,_options->selectionsize);
    return STATUS_ERR;
  }
  
  if ((mask_meta=read_volume(_options->mask_file, &tempdata, tmpsizes)) == NULL) {
    fprintf(stderr,"ERROR! Image not read: %s\n",_options->mask_file);
    return STATUS_ERR;
  }
  if ((tmpsizes[0]!=sizes[0][0]) || (tmpsizes[1]!=sizes[0][1]) || (tmpsizes[2]!=sizes[0][2])) {
    fprintf(stderr,"ERROR! Image dimension does not match library image dimension!\n");
    return STATUS_ERR;
  }
  free(tempdata);
  free_meta(mask_meta);

  /* make the downsampled masks crisp */
  segsubject     = alloc_2d_float(3, volumesize);
  patchcount     = alloc_2d_float(3, volumesize);
  filtered       = alloc_2d_float(3, volumesize);
  
  if (_options->verbose) fprintf(stderr,"Initial voxel size: %d\nTarget voxel size: %d\n", scales[initialscale], scales[targetscale]);
  
  /*perform selection before */
  selection = (int *)malloc(_options->selectionsize*sizeof(*selection));
  
  if(_options->use_double)
    pre_selection_double(subject[0], mask[0], images[0], sizes[0], num_images, _options->selectionsize, selection, selection_file, _options->verbose);
  else
    pre_selection(subject[0], mask[0], images[0], sizes[0], num_images, _options->selectionsize, selection, selection_file, _options->verbose);

  imagedata = (float *)malloc(_options->selectionsize*volumesize*sizeof(float));
  maskdata =  (float *)malloc(_options->selectionsize*volumesize*sizeof(float));
  
  /* read the library images, masks, and moments */
  for (i=0; i<_options->selectionsize; i++) {
    image_metadata *_meta;
    if (_options->verbose) fprintf(stderr,".");
    if ((_meta=read_volume(images[0][selection[i]], &tempdata, tmpsizes)) == NULL) {
      fprintf(stderr,"ERROR! Image not read: %s\n",images[0][selection[i]]);
      return STATUS_ERR;
    }
    cp_volume(tempdata, imagedata+i*volumesize, tmpsizes);
    free(tempdata);
    free_meta(_meta);
  }
  
  if (_options->verbose) fprintf(stderr,"*");
  for (i=0; i<_options->selectionsize; i++) {
    image_metadata *_meta;
    if (_options->verbose) fprintf(stderr,".");
    if ((_meta=read_volume(masks[0][selection[i]], &tempdata, tmpsizes)) == NULL) {
      fprintf(stderr,"ERROR! Image not read: %s\n",masks[0][selection[i]]);
      return STATUS_ERR;
    }
    cp_volume(tempdata, maskdata+i*volumesize, tmpsizes);
    free(tempdata);
    free_meta(_meta);
  }
  if (_options->verbose) fprintf(stderr,"*");
  
  
  for (scale=initialscale; scale>=targetscale; scale--) 
  {
    float * _imagedata;
    float * _maskdata;
    float * _meandata;
    float * _vardata;
    
    int scaledvolumesize = sizes[scale][0]*sizes[scale][1]*sizes[scale][2];
    
    printf("Scale:%d [%dx%dx%d]\n",scales[scale],sizes[scale][0],sizes[scale][1],sizes[scale][2]);
    
    if(scales[scale]>1) {
      _imagedata = (float *)malloc(_options->selectionsize*scaledvolumesize*sizeof(float));
      _maskdata =  (float *)malloc(_options->selectionsize*scaledvolumesize*sizeof(float));
      
      for (i=0; i<_options->selectionsize; i++) {
        down_sample(imagedata+i*volumesize, _imagedata+i*scaledvolumesize, scales[scale], sizes[0]);
        down_sample(maskdata+i*volumesize,  _maskdata+i*scaledvolumesize,  scales[scale], sizes[0]);
        /* make the downsampled masks crisp */
        threshold_data( _maskdata+i*scaledvolumesize, sizes[scale], 0.5);
      }
    } else {
      _imagedata = imagedata;
      _maskdata  = maskdata;
    }
    _meandata =  (float *)malloc(_options->selectionsize*scaledvolumesize*sizeof(float));
    _vardata =   (float *)malloc(_options->selectionsize*scaledvolumesize*sizeof(float));
    
    for (i=0; i<_options->selectionsize; i++) {
      if (_options->verbose) fprintf(stderr,"c");
      ComputeFirstMoment(_imagedata+i*scaledvolumesize, _meandata+i*scaledvolumesize, 
                         sizes[scale], configuration[scale].patchsize, &min, &max);
      
      ComputeSecondMoment(_imagedata+i*scaledvolumesize, _meandata+i*scaledvolumesize, 
                          _vardata+i*scaledvolumesize, sizes[scale], configuration[scale].patchsize, &min, &max);
    }
    
    /* remove any disconnected parts */
    masksize = getLargestObject_float(mask[scale], sizes[scale], 1, 0);

    if (_options->verbose) fprintf(stderr,"Mask size: %d\nAlpha: %f\n",masksize,configuration[scale].alpha);

    /* make sure we starting from a clean slate */
    wipe_data(segsubject[scale], sizes[scale] ,0.0);

    if(_options->use_sparse) {
#ifdef USE_SPAMS
      if(_options->use_double)
        max = nlmsegSparse4D_double(subject[scale], _imagedata, _maskdata, _meandata, _vardata, mask[scale], 
                        configuration[scale].patchsize, configuration[scale].searcharea, configuration[scale].beta, 
                        configuration[scale].threshold, sizes[scale], 
                        _options->selectionsize, segsubject[scale], patchcount[scale],
                        _options->lambda1, _options->lambda2, _options->sparse_mode, _options->sparse_stride
                        );
      else
        max = nlmsegSparse4D(subject[scale], _imagedata, _maskdata, _meandata, _vardata, mask[scale], 
                        configuration[scale].patchsize, configuration[scale].searcharea, configuration[scale].beta, 
                        configuration[scale].threshold, sizes[scale], 
                        _options->selectionsize, segsubject[scale], patchcount[scale],
                        _options->lambda1, _options->lambda2, _options->sparse_mode, _options->sparse_stride
                        );
#endif
    }  else {
      
      if(_options->use_double)
        max = nlmsegFuzzy4D_double(subject[scale], _imagedata, _maskdata, _meandata, _vardata, mask[scale], 
                        configuration[scale].patchsize, configuration[scale].searcharea, configuration[scale].beta, 
                        configuration[scale].threshold, sizes[scale], _options->selectionsize, segsubject[scale], patchcount[scale]
                           );
      else
        max = nlmsegFuzzy4D(subject[scale], _imagedata, _maskdata, _meandata, _vardata, mask[scale], 
                        configuration[scale].patchsize, configuration[scale].searcharea, configuration[scale].beta, 
                        configuration[scale].threshold, sizes[scale], _options->selectionsize, segsubject[scale], patchcount[scale]
                           );
    }
    if(scales[scale]>1) {
      free(_imagedata);
      free(_maskdata);
    }
    free(_meandata);
    free(_vardata);

    if (_options->positive_file!=NULL) {
      /* add the certain positive segmentation (inside mask) */
      add_mask_data(segsubject[scale], positivemask[scale], sizes[scale]);
    }

    /* add the certain segmentation from the previous scale */
    add_mask_data(segsubject[scale], segmented[scale], sizes[scale]);

    if (_options->medianfilter) {
      median_filter(segsubject[scale], sizes[scale], 3);
    }

    /* the patch filter is experimental */
    if (_options->patchfilter) {
      wipe_data(filtered[scale],sizes[scale],0.0);
      wipe_data(patchcount[scale],sizes[scale],0.0);
      max = nlmfilter(subject[scale], mask[scale], segsubject[scale], 2*configuration[scale].patchsize, 2*configuration[scale].searcharea, configuration[scale].beta, configuration[scale].threshold, sizes[scale], filtered[scale], patchcount[scale]);
      combine_maps(segsubject[scale], filtered[scale], mask[scale], sizes[scale]);
    }

    if (scale > targetscale) {
      /* if performing a higher resolution step, upsample the result and create new mask */
      resize_trilinear(segsubject[scale], sizes[scale], sizes[scale-1], segsubject[scale-1]);
      
      masksize=update_mask(segsubject[scale-1], mask[scale-1], segmented[scale-1], sizes[scale-1], 
                           configuration[scale].alpha, 1.0-configuration[scale].alpha);
    }

  } /* for each scale */
  free(selection);
  free(imagedata);
  free(maskdata);

  if (count_file!=NULL) {
    if(write_volume_generic(count_file, patchcount[targetscale], meta[targetscale],FALSE))
    {
      fprintf(stderr,"Can't save output to %s\n",count_file);
      return STATUS_ERR;
    }
  }

  if(targetscale!=0 && _options->same_res) { /* need to upsample final output */
    if (_options->verbose) fprintf(stderr,"Upsampling to input resolution, %dx%dx%d\n",sizes[0][0],sizes[0][1],sizes[0][2]);
    resize_trilinear(segsubject[targetscale], sizes[targetscale], sizes[0], segsubject[0]);
    masksize=update_mask(segsubject[0], mask[0], segmented[0], sizes[0], configuration[targetscale].alpha, 1.0-configuration[targetscale].alpha);
    targetscale=0;
    configuration[targetscale].alpha = _options->alpha;
  }

  if (!_options->outputprob) {
    if (_options->verbose) fprintf(stderr,"Thresholding estimator at %f\n",configuration[targetscale].alpha);
    threshold_data(segsubject[targetscale], sizes[targetscale], configuration[targetscale].alpha);
    getLargestObject_float(segsubject[targetscale], sizes[targetscale], 1, 0);

    if (_options->fill_output) {
      wipe_data(mask[targetscale], sizes[targetscale], 1.0);
      filled = flood_fill_float(segsubject[targetscale], mask[targetscale], sizes[targetscale], 0, 0, 0, 0, 6);
      //segsubject[targetscale]=mask[targetscale];
      cp_volume(mask[targetscale],segsubject[targetscale],sizes[targetscale]);
    }
  }

  meta[targetscale]->history=strdup(_options->history_label);
  if(write_volume_generic(_options->output_file, segsubject[targetscale], meta[targetscale],!_options->outputprob)) {
    fprintf(stderr,"Can't save output to %s\n",_options->output_file);
    return STATUS_ERR;
  }

  if(_options->compare_file)
  {
      image_metadata *ref_meta;
      float *ref_data;
      int ref_sizes[5];
      double kappa=0.0;
      
      if ((ref_meta=read_volume(_options->compare_file, &ref_data, ref_sizes)) == NULL) {
        fprintf(stderr,"ERROR! Image not read: %s\n",_options->compare_file);
        return STATUS_ERR;
      }
      
      /*TODO: check sizes!*/
      
      kappa=dice_kappa(ref_data,segsubject[targetscale],sizes[targetscale]);
      
      free(ref_data);
      free_meta(ref_meta);
      
      fprintf(stdout,"Comparision kappa:%f\n",kappa);
      
      if(kappa<_options->kappa_limit)
      {
        fprintf(stderr,"Kappa below threshold!\n");
        return STATUS_ERR;
      }
    
  }
  
  
  free_2d_float(mask);
  free_2d_float(subject);
  if (_options->positive_file!=NULL)
    free_2d_float(positivemask);
  
  free_2d_float(filtered);
  free_2d_float(segmented);
  free_2d_float(segsubject);
  free_2d_float(patchcount);

  free_3d_char(images);
  free_3d_char(masks);
  free_3d_char(means);
  free_3d_char(vars);
  
  free_meta(meta[0]);
  
  free(meta);
  

  return STATUS_OK;
}

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
