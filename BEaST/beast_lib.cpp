/*  beast_lib.c
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif //HAVE_CONFIG_H


#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "array_alloc.h"
#include "beast.h"

#ifdef HAVE_MINC
#include "mincio.h"
#endif 

#ifdef HAVE_NIFTI
#include "niftiio.h"
#endif 

#ifdef MT_USE_OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 1
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif


//#define DEBUG

#define LABEL_PUSH(stack,current,voxel) current++; \
                                        stack[current]=voxel;

#define LABEL_POP(stack,current,voxel) if(current>-1) { \
                                         voxel=stack[current]; \
                                         current--; } \
                                       else \
                                         current=-2;

int cmp_float(const void *vp, const void *vq){
  float diff = *(float *)vp - *(float *)vq;
  
  return ((diff>=0.0) ? ((diff>0.0) ? +1 : 0) : -1);  
}

/* Read one line from fp, */
/* copying it to line array (but no more than max chars). */
/* Does not place terminating \n in line array. */
/* Returns line length, or 0 for empty line, or EOF for end-of-file. */

int fgetline(FILE *fp, char line[], int max){
  int nch = 0;
  int c;
  max = max - 1;   /* leave room for '\0' */
  
  while((c = getc(fp)) != EOF)
    {
      if(c == '\n')
	break;

      if(nch < max)
	{
	  line[nch] = c;
	  nch = nch + 1;
	}
    }

  if(c == EOF && nch == 0)
    return EOF;
  
  line[nch] = '\0';
  return nch;
}


int median_filter(float *volume, int *sizes, int filtersize){
  int i,numelements,margin;
  float **_kernel;
  float *result;

  margin = filtersize / 2;
  
  numelements=filtersize*filtersize*filtersize;
  
  _kernel=(float**)malloc(omp_get_max_threads()*sizeof(float*));
  
  for(i=0;i<omp_get_max_threads();i++)
      _kernel[i] = (float*)malloc(numelements*sizeof(float));


  result = (float*)malloc(sizes[0]*sizes[1]*sizes[2]*sizeof(*result));

  #pragma omp parallel for shared(result) 
  for (i=margin;i<sizes[0]-margin;i++)
  {
    int j,k;
    float *kernel=_kernel[omp_get_thread_num()];
    
    for (j=margin;j<sizes[1]-margin;j++)
      for (k=margin;k<sizes[2]-margin;k++)
      {
        int ii,jj,kk,e,index,kindex;
        index = i*sizes[2]*sizes[1] + j*sizes[2] + k;
        e=0;
        for (ii=-margin;ii<=margin;ii++)
          for (jj=-margin;jj<=margin;jj++)
            for (kk=-margin;kk<=margin;kk++){
              kindex = index + (ii*sizes[2]*sizes[1] + jj*sizes[2] + kk);
              if(volume[kindex]>=0.0) /*DON't use negative values*/
                kernel[e++] = volume[kindex];
            }
        /* sort the values and assign the median */
        if(e>0) {
            qsort(kernel,e,sizeof(float),cmp_float);
            result[index]=kernel[numelements/2];
        }
      }
  }
  
  #pragma omp parallel for shared(result)
  for (i=margin;i<sizes[0]-margin;i++)
  {
    int j,k;
    for (j=margin;j<sizes[1]-margin;j++)
      for (k=margin;k<sizes[2]-margin;k++)
      {
        int index;
        index = i*sizes[2]*sizes[1] + j*sizes[2] + k;
        volume[index] = result[index];
      }
  }
  
  for(i=0;i<omp_get_max_threads();i++)
      free(_kernel[i]);
  
  free(_kernel);
  free(result);
  
  return STATUS_OK;
}

int trilinear_interpolant(float *volume, int *sizes, point3D coord, float *result)
{
   int slcind, rowind, colind, slcmax, rowmax, colmax;
   int slcnext, rownext, colnext;
   float f0, f1, f2, r0, r1, r2, r1r2, r1f2, f1r2, f1f2;
   float v000, v001, v010, v011, v100, v101, v110, v111;

   /* Check that the coordinate is inside the volume */
   slcmax = sizes[0] - 1;
   rowmax = sizes[1] - 1;
   colmax = sizes[2] - 1;
   if ((coord.z  < 0) || 
       (coord.z  > slcmax) ||
       (coord.y  < 0) || 
       (coord.y  > rowmax) ||
       (coord.x  < 0) || 
       (coord.x  > colmax)) {
      *result = 0;
      return FALSE;
   }

   /* Get the whole part of the coordinate */ 
   slcind = (int) coord.z;
   rowind = (int) coord.y;
   colind = (int) coord.x;
   if (slcind >= slcmax-1) slcind = slcmax-1;
   if (rowind >= rowmax-1) rowind = rowmax-1;
   if (colind >= colmax-1) colind = colmax-1;

   /* Get the next voxel up */
   slcnext = slcind+1;
   rownext = rowind+1;
   colnext = colind+1;

   /* Check for case of dimension of length one */
   if (slcmax == 0) {
      slcind = 0;
      slcnext = 0;
   }
   if (rowmax == 0) {
      rowind = 0;
      rownext = 0;
   }
   if (colmax == 0) {
      colind = 0;
      colnext = 0;
   }

   /* Get the relevant voxels */
   v000 = volume[slcind*sizes[1]*sizes[2] + rowind*sizes[2] + colind];
   v001 = volume[slcind*sizes[1]*sizes[2] + rowind*sizes[2] + colnext];
   v010 = volume[slcind*sizes[1]*sizes[2] + rownext*sizes[2] + colind];
   v011 = volume[slcind*sizes[1]*sizes[2] + rownext*sizes[2] + colnext];

   v100 = volume[slcnext*sizes[1]*sizes[2] + rowind*sizes[2] + colind];
   v101 = volume[slcnext*sizes[1]*sizes[2] + rowind*sizes[2] + colnext];
   v110 = volume[slcnext*sizes[1]*sizes[2] + rownext*sizes[2] + colind];
   v111 = volume[slcnext*sizes[1]*sizes[2] + rownext*sizes[2] + colnext];

   /* Get the fraction parts */
   f0 = coord.z  - slcind;
   f1 = coord.y  - rowind;
   f2 = coord.x - colind;
   r0 = 1.0 - f0;
   r1 = 1.0 - f1;
   r2 = 1.0 - f2;

   /* Do the interpolation */
   r1r2 = r1 * r2;
   r1f2 = r1 * f2;
   f1r2 = f1 * r2;
   f1f2 = f1 * f2;
   *result =
      r0 *  (r1r2 * v000 +
             r1f2 * v001 +
             f1r2 * v010 +
             f1f2 * v011);
   *result +=
      f0 *  (r1r2 * v100 +
             r1f2 * v101 +
             f1r2 * v110 +
             f1f2 * v111);
   
   return TRUE;

}

int resize_volume(float *input, int *sizes, int *sizes2, float *result){
  int zd,yd,xd,z2d,y2d,x2d;
  int i,outside=0;
  float x_ratio;
  float y_ratio;
  float z_ratio;

  /* original sizes */
  zd = sizes[0];
  yd = sizes[1];
  xd = sizes[2];

  /* new sizes */
  z2d = sizes2[0];
  y2d = sizes2[1];
  x2d = sizes2[2];

  /* the old to new ratio */
  x_ratio = (float)(xd)/(float)x2d ;
  y_ratio = (float)(yd)/(float)y2d ;
  z_ratio = (float)(zd)/(float)z2d ;

  /* for each voxel in the new resolution */
  #pragma omp parallel for reduction(+:outside)
  for (i=0;i<z2d;i++) {
    int j,k;
    for (j=0;j<y2d;j++) {
      for (k=0;k<x2d;k++) {
        int index;
        point3D p;
        float value;
        
        /* the position in the old resolution */
        p.x = x_ratio * (float)(k) ;
        p.y = y_ratio * (float)(j) ;
        p.z = z_ratio * (float)(i) ;
        
        if (!trilinear_interpolant(input, sizes, p, &value)){
          outside++;
        }

        index = i*y2d*x2d + j*x2d + k;
        
        result[index] = value;
      }
    }
  }

  //fprintf(stderr,"Outside: %d\n",outside);

  return STATUS_OK;
}

void cp_volume(float *data, float *copy, int *sizes){
  int i;
  
  #pragma omp parallel for
  for (i=0;i<sizes[0];i++)
  {
    int j,k;
    for (j=0;j<sizes[1];j++)
      for (k=0;k<sizes[2];k++)
        copy[i*sizes[1]*sizes[2]+j*sizes[2]+k] = data[i*sizes[1]*sizes[2]+j*sizes[2]+k];
  }
}


int cmp_ssd(const void *vp, const void *vq){
  const ssd_t *t1 = (ssd_t *)vp;
  const ssd_t *t2 = (ssd_t *)vq;
  const float p = t1->ssd;
  const float q = t2->ssd;
  float diff = p - q;
  
  return ((diff>=0.0) ? ((diff>0.0) ? +1 : 0) : -1);
}

int cmp_ssd_d(const void *vp, const void *vq){
  const ssd_td *t1 = (ssd_td *)vp;
  const ssd_td *t2 = (ssd_td *)vq;
  const double p = t1->ssd;
  const double q = t2->ssd;
  double diff = p - q;
  
  return ((diff>=0.0) ? ((diff>0.0) ? +1 : 0) : -1);
}


static inline float get_ssd(const float *  I1,const float *  I2,const float *  mask,const int *  sizes){
  int j;
  float count=0.0;
  float ssd=0.0;
  
  #pragma omp parallel for reduction(+:ssd) reduction(+:count)
  for (j=0;j<sizes[0];j++){
    int k,l;
    for (k=0;k<sizes[1];k++){
      
      #if _OPENMP>=201307
        #pragma omp simd
      #endif  
      for (l=0;l<sizes[2];l++){
        int index = j*sizes[1]*sizes[2] + k*sizes[2] + l;
        float msk=mask[index];
//         if (mask[index]>0.0){
//           ssd += SQR(I1[index] - I2[index]);
//           count++;
//         }
        ssd+=SQR(I1[index] - I2[index])*msk;
        count+=msk;
      }
    }
  }
  
  ssd /= count;
  return ssd;
}


double dice_kappa(const float *  I1,const float *  I2,const int *  sizes) {
  int j;
  
  double count1=0.0;
  double count2=0.0;
  double intersect=0.0;
  
  #pragma omp parallel for reduction(+:count1) reduction(+:count2) reduction(+:intersect)
  for (j=0;j<sizes[0];j++){
    int k,l;
    for (k=0;k<sizes[1];k++){
      
      #if _OPENMP>=201307
        #pragma omp simd
      #endif  
      for (l=0;l<sizes[2];l++){
        int index = j*sizes[1]*sizes[2] + k*sizes[2] + l;
        if(I1[index]>0.5) count1+=1.0;
        if(I2[index]>0.5) count2+=1.0;
        if(I1[index]>0.5 && I2[index]>0.5) intersect+=1.0;
      }
    }
  }
  
  count1+=count2;
  
  if(count1>0.0)
    return 2*intersect/count1;
  else
    return 0.0;
}


static inline float get_ssd_double(const float *  I1,const float *  I2,const float *  mask,const int *  sizes){
  int j;
  double count=0.0;
  double ssd=0;
  
  #pragma omp parallel for reduction(+:ssd) reduction(+:count)
  for (j=0;j<sizes[0];j++){
    int k,l;
    for (k=0;k<sizes[1];k++){
      
      #if _OPENMP>=201307
        #pragma omp simd
      #endif  
      for (l=0;l<sizes[2];l++){
        int index = j*sizes[1]*sizes[2] + k*sizes[2] + l;
        double msk=mask[index];
//         if (mask[index]>0.0){
//           ssd += SQR(I1[index] - I2[index]);
//           count++;
//         }
        ssd+=SQR(I1[index] - I2[index])*msk;
        count+=msk;
      }
    }
  }
  
  ssd /= count;
  return ssd;
}


int down_sample(float *subject, float *result, int factor, int *sizes){
  int i;
  int di,dj,dk;

  di = sizes[0]/factor;
  dj = sizes[1]/factor;
  dk = sizes[2]/factor;

  #pragma omp parallel for
  for (i=0;i<di;i++){
    int j,k;
    
    for (j=0;j<dj;j++){
      for (k=0;k<dk;k++){
        float av=0;
        int   c=0;
        int   l,n,m;
        
        for (l=0;l<factor;l++)
          for (m=0;m<factor;m++)
            for (n=0;n<factor;n++)
            {
              av+=subject[(i*factor+l)*sizes[2]*sizes[1] + (j*factor+m)*sizes[2] + k*factor+n];
              c++;
            }

        if(!c) continue;
        
        result[ i*dj*dk + j*dk + k] = av/c;
      }
    }
  }      
  
  return STATUS_OK;
}

int up_sample(float *subject, float *result, int factor, int *sizes, int *targetsizes){
  int i;
  int di,dj,dk;
  int adji,adjj,adjk;

  fprintf(stderr,"Upsampling...\nSource dimension: %d %d %d\nTarget dimensions: %d %d %d\n",sizes[0],sizes[1],sizes[2],targetsizes[0],targetsizes[1],targetsizes[2]);

  di = targetsizes[0];
  dj = targetsizes[1];
  dk = targetsizes[2];

  adji = sizes[0]*factor;
  adjj = sizes[1]*factor;
  adjk = sizes[2]*factor;

  #pragma omp parallel for 
  for (i=0;i<di;i++){
    int j,k;
    
    for (j=0;j<dj;j++){
      for (k=0;k<dk;k++){
        
        int index=i*dj*dk + j*dk + k;
        
        if ((i<adji) && (j<adjj) && (k<adjk)){
          int li = (i/factor)*sizes[2]*sizes[1] + (j/factor)*sizes[2] + k/factor;
          result[index] = subject[li];
        }else{
          result[index] = 0;
        }
      }
    }
  }

  return STATUS_OK;
}

int resize_trilinear(float *input, int *sizes, int *sizes2, float *result) {
  float x_ratio;
  float y_ratio;
  float z_ratio;
  int i;
  int zd,yd,xd,z2d,y2d,x2d;

  zd = sizes[0];
  yd = sizes[1];
  xd = sizes[2];

  z2d = sizes2[0];
  y2d = sizes2[1];
  x2d = sizes2[2];

  x_ratio = ((float)(xd))/x2d ;
  y_ratio = ((float)(yd))/y2d ;
  z_ratio = ((float)(zd))/z2d ;


  #pragma omp parallel for
  for (i=0;i<z2d;i++) {
    int j,k;
    
    for (j=0;j<y2d;j++) {
      
      for (k=0;k<x2d;k++) {
        
        float A, B, C, D, E, F, G, H;
        
        int x = (int)(x_ratio * (k - 0.5)) ;
        int y = (int)(y_ratio * (j - 0.5)) ;
        int z = (int)(z_ratio * (i - 0.5)) ;
        
        float x_diff = (x_ratio * (k - 0.5)) - x ;
        float y_diff = (y_ratio * (j - 0.5)) - y ;
        float z_diff = (z_ratio * (i - 0.5)) - z ;
        int  index = z*yd*xd + y*xd + x;
        int  offset = k+ j*x2d + i*x2d*y2d;
        
        // get the values
        A = input[index];
        B = input[index+1];
        C = input[index+xd];
        D = input[index+xd+1];

        E = input[index+yd*xd];
        F = input[index+yd*xd+1];
        G = input[index+yd*xd+xd];
        H = input[index+yd*xd+xd+1];
                
        result[offset] =
          A*(1-x_diff)*(1-y_diff)*(1-z_diff) +
          B*(x_diff)*(1-y_diff)*(1-z_diff) +
          C*(y_diff)*(1-x_diff)*(1-z_diff) +
          D*(x_diff*y_diff)*(1-z_diff) +
          E*(1-x_diff)*(1-y_diff)*z_diff +
          F*(x_diff)*(1-y_diff)*z_diff +
          G*(1-x_diff)*(y_diff)*z_diff +
          H*x_diff*y_diff*z_diff ;	
        }
      }
    }
    return STATUS_OK;
}


int threshold_data(float *data,const int *sizes, float threshold){
  int i;
  
  #pragma omp parallel for
  for (i=0;i<sizes[0];i++){
    int j,k;
    for (j=0;j<sizes[1];j++){
      for (k=0;k<sizes[2];k++){
        int index=i*sizes[2]*sizes[1] + j*sizes[2] + k;
        
        if (data[index]>threshold)
          data[index]=1;
        else
          data[index]=0;
      }
    }
  }
  return STATUS_OK;
}

int add_mask_data(float *data1, float *mask, int *sizes){
  int i;
  
  #pragma omp parallel for 
  for (i=0;i<sizes[0];i++){
    int j,k;
    
    for (j=0;j<sizes[1];j++){
      for (k=0;k<sizes[2];k++){
        int index=i*sizes[2]*sizes[1] + j*sizes[2] + k;
        if (mask[index])
          data1[index]=1;
      }
    }
  }

  return STATUS_OK;
}

int wipe_data(float *data1, int *sizes, float value){
  int i;
  
  #pragma omp parallel for 
  for (i=0;i<sizes[0];i++){
    int j,k;
    for (j=0;j<sizes[1];j++){
      for (k=0;k<sizes[2];k++){	  
        int index=i*sizes[2]*sizes[1] + j*sizes[2] + k;
        data1[index]=value;
      }
    }
  }

  return STATUS_OK;
}


int update_mask(float *subject, float *mask, float *segmented, int *sizes, float min, float max){
  int i;
  int count=0;
  
  #pragma omp parallel for reduction(+:count)
  for (i=0;i<sizes[0];i++){
    int j,k;
    for (j=0;j<sizes[1];j++){
      for (k=0;k<sizes[2];k++){
        int index=i*sizes[2]*sizes[1] + j*sizes[2] + k;
        
        if (mask[index])
        {
          mask[index] = 1;
          if (subject[index]<min){
            segmented[index] = 0;
            mask[index] = 0;
          }
          else if (subject[index]>max){
            segmented[index] = 1;
            mask[index] = 0;
          } else {
            count++;
          }
        }
        
      }
    }
  }
  return count;
}

int flip_data(float *data, float *result, int *sizes){
  int i;

  #pragma omp parallel for
  for (i=0;i<sizes[0];i++){
    int j,k;
    for (j=0;j<sizes[1];j++){
      for (k=0;k<sizes[2];k++){
        int index1=i*sizes[2]*sizes[1] + j*sizes[2] + k;
        int index2=i*sizes[2]*sizes[1] + j*sizes[2] + (sizes[2]-1-k);
        result[index1]=data[index2];
      }
    }
  }

  return STATUS_OK;
}

int combine_maps(float *data, float *map, float *mask, int *sizes){
  int i;

  #pragma omp parallel for
  for (i=0;i<sizes[0];i++){
    int j,k;
    for (j=0;j<sizes[1];j++){
      for (k=0;k<sizes[2];k++){	  
        int index=i*sizes[2]*sizes[1] + j*sizes[2] + k;	
        if (mask[index]){
          //data[index]=(data[index]+map[index])/2;
          data[index]=map[index];
        }
      }
    }
  }

  return STATUS_OK;
}

int flood_fill_float(float *data, float *output, int *sizes, int sx, int sy, int sz, float fill_value, int connectivity){
  int **mask,count,i,j,k,index;
  int marked, current, total_voxels;
  point3D *stack, current_voxel, test_voxel;
  float iso_value, original_value, output_value;

  /* Creating mask corresponding to connectivity */
  if(connectivity==26){
    mask = (int**)malloc(connectivity*sizeof(*mask));
    mask[0] = (int*)malloc(connectivity*3*sizeof(**mask));
    for(i=1;i<connectivity;i++)
      mask[i] = mask[0] + i*3;
    count=0;
    for(i=-1;i<2;i++)
      for(j=-1;j<2;j++)
        for(k=-1;k<2;k++){
          if(!((i==0)&&(j==0)&&(k==0))){
            mask[count][0] = i;
            mask[count][1] = j;
            mask[count++][2] = k;
          }
        }
    if(count!=connectivity) fprintf(stderr, "ERROR: error creating mask!\n");
  }
  else if(connectivity==18){
    mask = (int**)malloc(connectivity*sizeof(*mask));
    mask[0] = (int*)malloc(connectivity*3*sizeof(**mask));
    for(i=1;i<connectivity;i++)
      mask[i] = mask[0] + i*3;
    count=0;
    for(i=-1;i<2;i++)
      for(j=-1;j<2;j++)
        for(k=-1;k<2;k++){
          if(!((ABS(i)==1)&&(ABS(j)==1)&&(ABS(k)==1)) && !((i==0)&&(j==0)&&(k==0))) {
            mask[count][0] = i;
            mask[count][1] = j;
            mask[count++][2] = k;
          }
        }
    if(count!=connectivity) fprintf(stderr, "ERROR: error creating mask!\n");
  }
  else if(connectivity==6){
    mask = (int**)malloc(connectivity*sizeof(*mask));
    mask[0] = (int*)malloc(connectivity*3*sizeof(**mask));
    for(i=1;i<connectivity;i++)
      mask[i] = mask[0] + i*3;
    
    mask[0][0]=-1;mask[0][1]=0;mask[0][2]=0;
    mask[1][0]=0;mask[1][1]=-1;mask[1][2]=0;
    mask[2][0]=0;mask[2][1]=+1;mask[2][2]=0;
    mask[3][0]=0;mask[3][1]=0;mask[3][2]=-1;
    mask[4][0]=0;mask[4][1]=0;mask[4][2]=+1;
    mask[5][0]=1;mask[5][1]=0;mask[5][2]=0;
  }
  else {
    fprintf(stderr,"ERROR: the specified connectivity %d is not supported!\n",connectivity);
    return STATUS_ERR;
  }

  total_voxels = sizes[0]*sizes[1]*sizes[2];
  //fprintf(stderr,"Total voxels: %d\n",total_voxels);

  /* Allocate stack */
  stack = (point3D *)malloc(total_voxels*sizeof(*stack));

  marked=1;
  current = -1;
  SET_3DPOINT(current_voxel,sx,sy,sz);
  index=sz*sizes[2]*sizes[1] + sy*sizes[2] + sx;
  iso_value = data[index];
#ifdef DEBUG
  fprintf(stderr,"Iso value: %f\n", iso_value);
  fprintf(stderr,"Fill value: %f\n",fill_value);
#endif
  output[index] = fill_value;
  
  /* flood fill the connected component */
  do{
    for(i=0;i<connectivity;i++){
      SET_3DPOINT(test_voxel,current_voxel.x+mask[i][0],current_voxel.y+mask[i][1],current_voxel.z+mask[i][2]);
      if((test_voxel.x > -1) && (test_voxel.x < sizes[0]) && 
         (test_voxel.y > -1) && (test_voxel.y < sizes[1]) && 
         (test_voxel.z > -1) && (test_voxel.z < sizes[2])) {
        index=test_voxel.z*sizes[2]*sizes[1] + test_voxel.y*sizes[2] + test_voxel.x;
        original_value = data[index];
        output_value = output[index];
        /* if the voxel is connected to the component and the voxels has not already been filled */
        if ((original_value == iso_value) && (output_value != fill_value)) {
          output[index] = fill_value;
          LABEL_PUSH(stack,current,test_voxel);
          marked++;
        }
      }
    }

    if ((current > total_voxels - 1) || (marked > total_voxels )){
      fprintf(stderr,"ERROR! Stack overflow! (total_voxels=%d, current=%d, marked=%d)\n",total_voxels,current,marked);
      free(mask[0]);
      free(mask);
      free(stack);
      return marked;
    }

    LABEL_POP(stack,current,current_voxel);
    
  } while(current!=-2);

  free(mask[0]);
  free(mask);
  free(stack);
  
  return marked;
}


int read_configuration(const char *filename, beast_conf *conf){
  int i,size=0;
  FILE *fd;
  char line[FILENAMELENGTH];
  VIO_BOOL issane=TRUE;

  fd=fopen(filename,"r");

  if (fd==NULL){
    fprintf(stderr,"ERROR! Cannot open configuration file %s\n",filename);
    exit(-1);
  }

  while (fgetline(fd, line, FILENAMELENGTH)!=EOF){
    /* if not a comment */
    if (line[0] != 35){
      sscanf(line,"%d %d %d %lf %lf %lf %d",&conf[size].voxelsize,&conf[size].patchsize,&conf[size].searcharea,&conf[size].alpha,&conf[size].beta,&conf[size].threshold,&conf[size].selectionsize);
      size++;
    }

  }

  fclose(fd);

  /* configuration sanity check */
  for (i=0;i<size;i++){
    if ((conf[i].voxelsize > VOXELSIZEMAX) || (conf[i].voxelsize < VOXELSIZEMIN))
      issane = FALSE;
    if ((conf[i].patchsize > PATCHSIZEMAX) || (conf[i].patchsize < PATCHSIZEMIN))
      issane = FALSE;
    if ((conf[i].searcharea > SEARCHAREAMAX) || (conf[i].searcharea < SEARCHAREAMIN))
      issane = FALSE;
    if ((conf[i].alpha > ALPHAMAX) || (conf[i].alpha < ALPHAMIN))
      issane = FALSE;
    if ((conf[i].beta > BETAMAX) || (conf[i].beta < BETAMIN))
      issane = FALSE;
    if ((conf[i].threshold > THRESHOLDMAX) || (conf[i].threshold < THRESHOLDMIN))
      issane = FALSE;
  }

  if (issane)
    return size;
  else
    return STATUS_ERR;
}


int read_list(const char *filename, char **list,const char *basedir) {
  FILE *fd;
  char line[FILENAMELENGTH];
  int size=0;
  
  fd=fopen(filename,"r");

  while (fgetline(fd, line, FILENAMELENGTH)!=EOF){
    if(basedir && *basedir)
      sprintf(list[size],"%s/%s",basedir,line);
    else
      sprintf(list[size],"%s",line);
    size++;
  }

  fclose(fd);

  return size;
}


int pre_selection(const float *subject, const float *mask, 
                  char **images, const int *sizes, 
                  int librarysize, int num_selected, int *selection, 
                  const char *outfile, VIO_BOOL verbose)
{
  int i;
  int volumesize;
  ssd_t *ssd;
  FILE *fd=NULL;
  int _sizes[5];

  fprintf(stderr,"Performing pre-selection ");

  volumesize=sizes[0]*sizes[1]*sizes[2];
  ssd = (ssd_t *)malloc(librarysize*sizeof(ssd_t));

  fprintf(stderr,"(%ld MB/subject)",volumesize*sizeof(float)/(1024*1024));
  
  /*#pragma omp parallel for */
  for (i=0;i<librarysize;i++)
  {
    float *imagedata;
    image_metadata *_meta;
    fprintf(stderr,".");

    _meta=read_volume(images[i], &imagedata, _sizes);
    /*TODO:compare sizes and _sizes ?*/

    ssd[i].index=i;
    ssd[i].ssd=get_ssd(subject,imagedata,mask,sizes);
    
    free(imagedata);
    free_meta(_meta);
  }

  qsort(ssd, librarysize, sizeof(ssd_t), cmp_ssd);

  fprintf(stderr,"done\n");

  if (outfile!=NULL)
    fd = fopen(outfile,"a");
  else
    fd = stderr;

  for (i=0;i<num_selected;i++){   
    selection[i] = ssd[i].index;
    if ((verbose) || (outfile!=NULL)) fprintf(fd,"%s %f\n",images[selection[i]],ssd[i].ssd);
  }
  
  if (outfile!=NULL)
    fclose(fd);
  free(ssd);
  return STATUS_OK;
}

int pre_selection_double(const float *subject, const float *mask, 
                  char **images, const int *sizes, 
                  int librarysize, int num_selected, int *selection, 
                  const char *outfile, VIO_BOOL verbose)
{
  int i;
  int volumesize;
  ssd_td *ssd;
  FILE *fd=NULL;
  int _sizes[5];

  fprintf(stderr,"Performing pre-selection ");

  volumesize=sizes[0]*sizes[1]*sizes[2];
  ssd = (ssd_td *)calloc(librarysize,sizeof(ssd_td));

  fprintf(stderr,"(%ld MB/subject)",volumesize*sizeof(float)/(1024*1024));
  
  /*#pragma omp parallel for */
  for (i=0;i<librarysize;i++)
  {
    float *imagedata;
    image_metadata *_meta;
    fprintf(stderr,".");

    _meta=read_volume(images[i], &imagedata, _sizes);
    /*TODO:compare sizes and _sizes ?*/

    ssd[i].index=i;
    ssd[i].ssd=get_ssd_double(subject,imagedata,mask,sizes);
    
    free(imagedata);
    free_meta(_meta);
  }

  qsort(ssd, librarysize, sizeof(ssd_td), cmp_ssd_d);

  fprintf(stderr,"done\n");

  if (outfile!=NULL)
    fd = fopen(outfile,"a");
  else
    fd = stderr;

  for (i=0;i<num_selected;i++){   
    selection[i] = ssd[i].index;
    if ((verbose) || (outfile!=NULL)) fprintf(fd,"%s %f\n",images[selection[i]],ssd[i].ssd);
  }
  
  if (outfile!=NULL)
    fclose(fd);
  free(ssd);
  return STATUS_OK;
}


image_metadata * read_volume(const char *filename, float **data,int *sizes){
  image_metadata *meta=NULL;

  /* if minc format */
  if (!strcmp("mnc",filename + strlen(filename)-3) || !strcmp("mnc.gz",filename + strlen(filename)-6)){
  #ifdef HAVE_MINC
    meta = read_minc(filename, data, sizes);
  #else
    fprintf(stderr,"READ: unsupported file format (%s)\n", filename);
    return NULL;
  #endif

  }else{
    #ifdef HAVE_NIFTI
    if (!strcmp("nii",filename + strlen(filename)-3) || !strcmp("nii.gz",filename + strlen(filename)-6)){
      meta = read_nifti(filename, data, sizes);
    } else {
      fprintf(stderr,"READ: unsupported file format (%s)\n", filename);
      return NULL;
    }
   #else 
    fprintf(stderr,"READ: unsupported file format (%s)\n", filename);
    return NULL;
   #endif
  }

#ifdef DEBUG
  fprintf(stderr,"READ: Dimension sizes: %d, %d, %d\n",meta->length[0],meta->length[1],meta->length[2]);
  fprintf(stderr,"READ: Start coordinates: %f, %f, %f\n",meta->start[0],meta->start[1],meta->start[2]);
  fprintf(stderr,"READ: Step values: %f, %f, %f\n",meta->step[0],meta->step[1],meta->step[2]);
#endif
  
  if(meta)
    meta->history = NULL;

  return meta;
}

int write_volume_generic(const char *filename,const float *data,const image_metadata *meta,VIO_BOOL binary_mask){

#ifdef DEBUG
  fprintf(stderr,"WRITE: Dimension sizes: %d, %d, %d\n",meta->length[0],meta->length[1],meta->length[2]);
  fprintf(stderr,"WRITE: Start coordinates: %f, %f, %f\n",meta->start[0],meta->start[1],meta->start[2]);
  fprintf(stderr,"WRITE: Step values: %f, %f, %f\n",meta->step[0],meta->step[1],meta->step[2]);
#endif

  /* if minc format */
  if (!strcmp("mnc",filename + strlen(filename)-3)){
    #ifdef HAVE_MINC
    if(write_minc(filename, data, meta,binary_mask))
    {
      fprintf(stderr,"WRITE:Error writing file (%s)!\n", filename);
      return STATUS_ERR;
    }
    #else
    fprintf(stderr,"WRITE:Unsupported file format (%s)!\n", filename);
    #endif 
    
  }else{
    #ifdef HAVE_NIFTI
    /* assume nifti */
    write_nifti_generic(filename, data, meta);/*TODO: can nifti write binary masks?*/
    #else
    fprintf(stderr,"WRITE:Unsupported file format (%s)!\n", filename);
    #endif 
  }
  
  return STATUS_OK;
}




/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */

