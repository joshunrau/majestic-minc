/*  nlmsegFuzzy.c
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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <strings.h>
#include <string.h>
#include <float.h>
#include "nlmseg.h"

static const int MINCOUNT=100;
static const int MAXITER=5;

//#define INCREASE_RADIUS


#ifdef MT_USE_OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 1
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif


#ifdef USE_SPAMS

#include <spams.h>


float nlmsegSparse4D(const float *subject,const  float *imagedata, 
                     const float *maskdata,const  float *meandata,
                     const float *vardata, 
                     const float *mask, 
                     int sizepatch, int searcharea, float beta, float threshold, 
                     const int dims[3],  int librarysize, 
                     float *SegSubject, float *PatchCount,
                     float lambda1, float lambda2, int sparse_mode,int stride
                    )
{    
  float *MeansSubj, *VarsSubj, *localmask;
  int i,v,f,ndim;
  float min,max;
  int patch_volume;
  int patch_center_voxel;
  
  int sadims,volumesize,index;
  int mincount=MINCOUNT;
  int notfinished;
  int notconverged;
  int finished;
  double minidist;
  double epsi = 0.0001;
  time_t time1,time2;
  int   iterations=0;
  
  float constraint=1.0; /*TODO: find optimal value?*/
  
  /*per thread temp storage*/
  float  **_PatchImg  =(float **) malloc(omp_get_max_threads()*sizeof(float*));
  float  **_PatchMask =(float **) malloc(omp_get_max_threads()*sizeof(float*));
  float  **_PatchTemp =(float **) malloc(omp_get_max_threads()*sizeof(float*));
  
  
  float  **_SegAccum =(float **) malloc(omp_get_max_threads()*sizeof(float*));
  float  **_SegNorm  =(float **) malloc(omp_get_max_threads()*sizeof(float*));
  
  
  fprintf(stderr,"Running sparse segmentation\n");
  
  fprintf(stderr,"Patch size: %d\nSearch area: %d\nBeta: %f\nThreshold: %f\nSelection: %d\nLambda1:%f\nLambda2:%f\nSparse mode:%d\n",
        sizepatch,searcharea,beta,threshold,librarysize,lambda1,lambda2,sparse_mode);
  
  ndim = 3;
  volumesize=dims[0]*dims[1]*dims[2];
  
  /*Patch radius*/
  f = sizepatch;
  /*volume of patch*/
  patch_volume=(2*f+1)*(2*f+1)*(2*f+1);
  /*index of patch central voxel*/
  patch_center_voxel=f+f*(2*f+1)+f*(2*f+1)*(2*f+1);

  /*Search Area radius*/
  v = searcharea;
  
  sadims = pow(2*v+1,ndim)* librarysize;

  /* allocate memory for multithreading operation*/
  for(i=0;i<omp_get_max_threads();i++)
  {
    _PatchTemp[i]=(float*) malloc( patch_volume*sizeof(float) );
    
    _SegAccum[i]=(float*) malloc( volumesize*sizeof(float) );
    _SegNorm[i] =(float*) malloc( volumesize*sizeof(float) );
  }
  
  MeansSubj = (float *)calloc(volumesize,sizeof(*MeansSubj));
  VarsSubj =  (float *)calloc(volumesize,sizeof(*VarsSubj));
  localmask = (float *)calloc(volumesize,sizeof(*localmask));
  
  memmove(localmask,mask,volumesize*sizeof(*localmask));
  
  fprintf(stderr,"Dimensions: %d %d %d\n",dims[0],dims[1],dims[2]);
  
  ComputeFirstMoment(subject, MeansSubj, dims, f, &min, &max);
  ComputeSecondMoment(subject, MeansSubj, VarsSubj, dims, f, &min, &max);
  
  do {
    
    fprintf(stderr,"Segmenting              ");
    time2=time(NULL);
    notfinished=0;
    finished=0;
    notconverged=0;
    
    
    for(i=0;i<omp_get_max_threads();i++)
    {
      _PatchImg[i]  =(float*) malloc( patch_volume*sizeof(float)*sadims );
      _PatchMask[i] =(float*) malloc( patch_volume*sizeof(float)*sadims );
      
      /*reset weights*/
      memset(_SegAccum[i],0,sizeof(float)*volumesize);
      memset(_SegNorm[i],0,sizeof(float)*volumesize);
    }
    
    #pragma omp parallel for shared(_PatchImg,_PatchTemp,_PatchMask) reduction(+:notfinished) reduction(+:finished) reduction(+:notconverged) schedule(dynamic)
    for(i=0;i<dims[0];i+=stride) 
    { /*start parallel section*/
      int j,k;
      SpMatrix<float> _alpha(1, sadims, sadims);
      
      /*use thread-specific temp memory*/
      float * PatchImg     = _PatchImg[omp_get_thread_num()];
      float * PatchMask    =_PatchMask[omp_get_thread_num()];;
      float * PatchTemp    =_PatchTemp[omp_get_thread_num()];
      
      float * SegAccum     =_SegAccum[omp_get_thread_num()];
      float * SegNorm      =_SegNorm[omp_get_thread_num()];
      
       if( omp_get_thread_num()==0 )
         fprintf(stderr,"\b\b\b\b\b\b\b\b\b%3d / %3d", i+1, dims[0]);
      
      for(j=0;j<dims[1];j+=stride) 
      {
        for(k=0;k<dims[2];k+=stride) 
        {
          int index=i*(dims[2]*dims[1])+(j*dims[2])+k;
          
          /* mask check */
          if ( localmask[index]  > 0 )
          {
            int   ii,jj,kk;
            float proba=0.0;
            
            int   count = 0;
            
            float TMean,TVar;
            
            TMean = MeansSubj[index];
            TVar =  VarsSubj [index];
            
            
            //memset(PatchTemp,0,sizeof(float)*patch_volume);
            
            ExtractPatch_norm(subject, PatchTemp, i, j, k, f, dims[0], dims[1], dims[2], TMean);
            
            /*TODO:apply patch-normalization using Mean and var*/
            
            /* go through the search space  */
            for(ii=-v;ii<=v;ii++)
            {
              for(jj=-v;jj<=v;jj++)
              {
                for(kk=-v;kk<=v;kk++)
                {
                  int ni,nj,nk;
                  ni=i+ii;
                  nj=j+jj;
                  nk=k+kk;
                  
                  if( (ni>=0) && (nj>=0) && (nk>=0) && 
                      (ni<dims[0]) && (nj<dims[1]) && (nk<dims[2]) )
                  {
                    int t;
                    for(t=0;t<librarysize;t++)
                    {
                      int index2 = t*volumesize+ni*(dims[2]*dims[1])+(nj*dims[2])+nk;
                      
                      float Mean = meandata[index2];
                      float Var =  vardata [index2];
                      
                      /*Similar Luminance and contrast -> Cf Wang TIP 2004*/
                      float th = ((2 * Mean * TMean + epsi) / ( Mean*Mean + TMean*TMean + epsi))  *
                                 ((2 * sqrt(Var) * sqrt(TVar) + epsi) / (Var + TVar + epsi));

                      if(th > threshold)
                      {
                        //memset(&PatchImg[count*patch_volume],0,sizeof(float)*patch_volume);
                        
                        ExtractPatch4D_norm(imagedata, &PatchImg[count*patch_volume] ,ni,nj,nk, t, f, dims[0], dims[1], dims[2], Mean);
                        ExtractPatch4D(maskdata,  &PatchMask[count*patch_volume] ,ni,nj,nk, t, f, dims[0], dims[1], dims[2]);
                        
                        count++;
                        
                      }
                    }
                  }
                }
              }
            }
            
            /* require a minimum number of selected patches  */
            if (count >= mincount) {
                int p;
                
                
                Matrix<float> M(PatchImg,   patch_volume, count);
                Matrix<float> X(PatchTemp,  patch_volume, 1 );
                /*TEST*/
                lasso<float>(M, X, _alpha, 
                             MINCOUNT, lambda1, lambda2, (constraint_type)sparse_mode,  /*L1COEFFS*/
                             true, false, 1); /*TODO: maybe use MINCOUNT for max count?*/
                
                float minidist = FLT_MAX; /*FLT_MAX;*/
                
                if( _alpha.nnz()>0 )
                {
                    float totalweight=0.0;
                    
                    for(p=0;p<count;p++)
                    {
                        float w = _alpha[p];
                        if ( w>0.0 )
                        {
                            totalweight += w;
                        }
                    }
                    
                    for(p=0;p<count;p++)
                    {
                        float w = _alpha[p];
                        if ( w>0.0 )
                        {
                            AddWPatch(SegAccum,&PatchMask[p*patch_volume],w/totalweight, i, j, k, f, dims[0], dims[1], dims[2]);
                            AddW(SegNorm,w/totalweight,                                  i, j, k, f, dims[0], dims[1], dims[2]);
                        }
                    }                    
                    
                } else {
                    notconverged++;
                    notfinished++;
                }
            } else {
              /* Not enough similar patches */
              notfinished++;
            }
          
          }// mask check
          
        } // for k
      } // for j
      
      /*end of parallel section*/
    } // for i
    
    /*accumulate all segmentations across all threads*/
    notfinished=0;
    #pragma omp parallel for shared(SegSubject,PatchCount,localmask) reduction(+:notfinished) 
    for(i=0;i<dims[0];i++)
    {
      int j,k;
      for(j=0;j<dims[1];j++) 
      {
        for(k=0;k<dims[2];k++) 
        {
          int index=i*(dims[2]*dims[1])+(j*dims[2])+k;
          if ( localmask[index]  > 0 )
          {
            float avg_seg=0.0;
            float avg_weight=0.0;
            int t;
            for(t=0;t<omp_get_max_threads();t++)
            {
                avg_seg+=_SegAccum[t][index];
                avg_weight+=_SegNorm[t][index];
            }
            
            if(avg_weight>0.0)
            {
                SegSubject[index] = avg_seg/avg_weight;
                PatchCount[index] = avg_weight;
            } else {
                SegSubject[index] = -1;
                notfinished++;
            }
          }              
        }
      }
    }
    if ( notfinished>0 ) {
      /*relax preselection criteria*/
      int count=0;
      threshold=threshold*0.95;
      mincount=mincount*0.95;
#ifdef INCREASE_RADIUS      
      v=v+1;
#endif
      
      #pragma omp parallel for  shared(SegSubject,PatchCount,localmask) reduction(+:count) 
      for(i=0;i<dims[0];i++)
      {
        int j,k;
        for(j=0;j<dims[1];j++)
        {
          for(k=0;k<dims[2];k++)
          {
            int index=i*(dims[2]*dims[1])+(j*dims[2])+k;
            
            if ( SegSubject[index]<0 ){
              localmask[index] = 1;
              count++;
            } else {
              localmask[index] = 0;
            }
          }
        }
      }
      
      fprintf(stderr," (redoing %d voxels) t=%f, min=%d \n",count, threshold, mincount);
      notfinished=count;
      sadims = pow( 2*v+1,ndim) * librarysize;
    }
    
    for(i=0;i<omp_get_max_threads();i++)
    {
        free(_PatchImg[i]);
        free(_PatchMask[i]);
    }
    iterations++;
  } while (notfinished && mincount>1 && iterations<MAXITER);

  if(notfinished) /*stopped early*/
  {
      /*TODO: report?*/
  }
  
  for(i=0;i<omp_get_max_threads();i++)
  {
    free(_PatchTemp[i]);
    free(_SegAccum[i]);
    free(_SegNorm[i]);
  }
  
  free(_PatchImg);
  free(_PatchTemp);
  free(_PatchMask);
  free(_SegAccum);
  free(_SegNorm);
  
  
  fprintf(stderr," done (t=%f, min=%d)\n",threshold, mincount);
  
  free(MeansSubj);
  free(VarsSubj);
  free(localmask);
  
  return max;
}

float nlmsegSparse4D_double(const float *subject,const  float *imagedata, 
                    const float *maskdata, const float *meandata,const  float *vardata, 
                    const float *mask, 
                    int sizepatch, int searcharea, double beta, double threshold, 
                    const int dims[3],   int librarysize, float *SegSubject, float *PatchCount,
                    double lambda1, double lambda2,int sparse_mode, int stride)
{
  float *MeansSubj, *VarsSubj, *localmask;
  int i,v,f,ndim;
  float min,max;
  int patch_volume;
  int patch_center_voxel;
  
  int sadims,volumesize,index;
  int mincount=MINCOUNT;
  int notfinished;
  int notconverged;
  int finished;
  double minidist;
  double epsi = 0.0001;
  time_t time1,time2;
  int   iterations=0;
  
  float constraint=1.0; /*TODO: find optimal value?*/
  
  /*per thread temp storage*/
  double  **_PatchImg  =(double **) malloc(omp_get_max_threads()*sizeof(double*));
  double  **_PatchMask =(double **) malloc(omp_get_max_threads()*sizeof(double*));
  double  **_PatchTemp =(double **) malloc(omp_get_max_threads()*sizeof(double*));
  
  
  double  **_SegAccum =(double **) malloc(omp_get_max_threads()*sizeof(double*));
  double  **_SegNorm  =(double **) malloc(omp_get_max_threads()*sizeof(double*));
  
  
  fprintf(stderr,"Running sparse segmentation\n");
  
  fprintf(stderr,"Patch size: %d\nSearch area: %d\nBeta: %f\nThreshold: %f\nSelection: %d\nLambda1:%f\nLambda2:%f\nSparse mode:%d\n",
        sizepatch,searcharea,beta,threshold,librarysize,lambda1,lambda2,sparse_mode);
  
  ndim = 3;
  volumesize=dims[0]*dims[1]*dims[2];
  
  /*Patch radius*/
  f = sizepatch;
  /*volume of patch*/
  patch_volume=(2*f+1)*(2*f+1)*(2*f+1);
  /*index of patch central voxel*/
  patch_center_voxel=f+f*(2*f+1)+f*(2*f+1)*(2*f+1);

  /*Search Area radius*/
  v = searcharea;
  
  sadims = pow(2*v+1,ndim)* librarysize;

  /* allocate memory for multithreading operation*/
  for(i=0;i<omp_get_max_threads();i++)
  {
    _PatchTemp[i]= (double*) malloc( patch_volume*sizeof(double) );
    _SegAccum[i] = (double*) malloc( volumesize*sizeof(double) );
    _SegNorm[i]  = (double*) malloc( volumesize*sizeof(double) );
  }
  
  MeansSubj = (float *)calloc(volumesize,sizeof(*MeansSubj));
  VarsSubj =  (float *)calloc(volumesize,sizeof(*VarsSubj));
  localmask = (float *)calloc(volumesize,sizeof(*localmask));
  
  memmove(localmask,mask,volumesize*sizeof(*localmask));
  
  fprintf(stderr,"Dimensions: %d %d %d\n",dims[0],dims[1],dims[2]);
  
  ComputeFirstMoment(subject, MeansSubj, dims, f, &min, &max);
  ComputeSecondMoment(subject, MeansSubj, VarsSubj, dims, f, &min, &max);
  
  do {
    
    fprintf(stderr,"Segmenting              ");
    time2=time(NULL);
    notfinished=0;
    finished=0;
    notconverged=0;
    
    
    for(i=0;i<omp_get_max_threads();i++)
    {
      _PatchImg[i]  =(double*) malloc( patch_volume*sizeof(double)*sadims );
      _PatchMask[i] =(double*) malloc( patch_volume*sizeof(double)*sadims );
      
      /*reset weights*/
      memset(_SegAccum[i],0,sizeof(double)*volumesize);
      memset(_SegNorm[i],0,sizeof(double)*volumesize);
    }
    
    #pragma omp parallel for shared(_PatchImg,_PatchTemp,_PatchMask) reduction(+:notfinished) reduction(+:finished) reduction(+:notconverged) schedule(dynamic)
    for(i=0;i<dims[0];i+=stride) 
    { /*start parallel section*/
      int j,k;
      SpMatrix<double> _alpha(1, sadims, sadims);
      
      /*use thread-specific temp memory*/
      double * PatchImg     = _PatchImg[omp_get_thread_num()];
      double * PatchMask    =_PatchMask[omp_get_thread_num()];;
      double * PatchTemp    =_PatchTemp[omp_get_thread_num()];
      
      double * SegAccum     =_SegAccum[omp_get_thread_num()];
      double * SegNorm      =_SegNorm[omp_get_thread_num()];
      
       if( omp_get_thread_num()==0 )
         fprintf(stderr,"\b\b\b\b\b\b\b\b\b%3d / %3d", i+1, dims[0]);
      
      for(j=0;j<dims[1];j+=stride) 
      {
        for(k=0;k<dims[2];k+=stride) 
        {
          int index=i*(dims[2]*dims[1])+(j*dims[2])+k;
          
          /* mask check */
          if ( localmask[index]  > 0 )
          {
            int   ii,jj,kk;
            double proba=0.0;
            
            int   count = 0;
            
            double TMean,TVar;
            
            TMean = MeansSubj[index];
            TVar =  VarsSubj [index];
            
            
            //memset(PatchTemp,0,sizeof(float)*patch_volume);
            
            ExtractPatch_norm_double(subject, PatchTemp, i, j, k, f, dims[0], dims[1], dims[2], TMean);
            
            /*TODO:apply patch-normalization using Mean and var*/
            
            /* go through the search space  */
            for(ii=-v;ii<=v;ii++)
            {
              for(jj=-v;jj<=v;jj++)
              {
                for(kk=-v;kk<=v;kk++)
                {
                  int ni,nj,nk;
                  ni=i+ii;
                  nj=j+jj;
                  nk=k+kk;
                  
                  if( (ni>=0) && (nj>=0) && (nk>=0) && 
                      (ni<dims[0]) && (nj<dims[1]) && (nk<dims[2]) )
                  {
                    int t;
                    for(t=0;t<librarysize;t++)
                    {
                      int index2 = t*volumesize+ni*(dims[2]*dims[1])+(nj*dims[2])+nk;
                      
                      double Mean = meandata[index2];
                      double Var =  vardata [index2];
                      
                      /*Similar Luminance and contrast -> Cf Wang TIP 2004*/
                      double th = ((2 * Mean * TMean + epsi) / ( Mean*Mean + TMean*TMean + epsi))  *
                                 ((2 * sqrt(Var) * sqrt(TVar) + epsi) / (Var + TVar + epsi));

                      if(th > threshold)
                      {
                        //memset(&PatchImg[count*patch_volume],0,sizeof(float)*patch_volume);
                        
                        ExtractPatch4D_norm_double(imagedata, &PatchImg[count*patch_volume] ,ni,nj,nk, t, f, dims[0], dims[1], dims[2], Mean);
                        ExtractPatch4D_double(maskdata,  &PatchMask[count*patch_volume] ,ni,nj,nk, t, f, dims[0], dims[1], dims[2]);
                        
                        count++;
                        
                      }
                    }
                  }
                }
              }
            }
            
            /* require a minimum number of selected patches  */
            if (count >= mincount) {
                int p;
                
                
                Matrix<double> M(PatchImg,   patch_volume, count);
                Matrix<double> X(PatchTemp,  patch_volume, 1 );
                /*TEST*/
                lasso<double>(M, X, _alpha, 
                             MINCOUNT, lambda1, lambda2, (constraint_type)sparse_mode,  /*L1COEFFS*/
                             true, false, 1); /*TODO: maybe use MINCOUNT for max count?*/
                
                double minidist = DBL_MAX; /*FLT_MAX;*/
                
                if( _alpha.nnz()>0 )
                {
                    float totalweight=0.0;
                    
                    for(p=0;p<count;p++)
                    {
                        float w = _alpha[p];
                        if ( w>0.0 )
                        {
                            totalweight += w;
                        }
                    }
                    
                    for(p=0;p<count;p++)
                    {
                        float w = _alpha[p];
                        if ( w>0.0 )
                        {
                            AddWPatch_double(SegAccum,&PatchMask[p*patch_volume],w/totalweight, i, j, k, f, dims[0], dims[1], dims[2]);
                            AddW_double(SegNorm,w/totalweight,                                  i, j, k, f, dims[0], dims[1], dims[2]);
                        }
                    }                    
                    
                } else {
                    notconverged++;
                    notfinished++;
                }
            } else {
              /* Not enough similar patches */
              notfinished++;
            }
          
          }// mask check
          
        } // for k
      } // for j
      
      /*end of parallel section*/
    } // for i
    
    /*accumulate all segmentations across all threads*/
    notfinished=0;
    #pragma omp parallel for shared(SegSubject,PatchCount,localmask) reduction(+:notfinished) 
    for(i=0;i<dims[0];i++)
    {
      int j,k;
      for(j=0;j<dims[1];j++) 
      {
        for(k=0;k<dims[2];k++) 
        {
          int index=i*(dims[2]*dims[1])+(j*dims[2])+k;
          if ( localmask[index]  > 0 )
          {
            float avg_seg=0.0;
            float avg_weight=0.0;
            int t;
            for(t=0;t<omp_get_max_threads();t++)
            {
                avg_seg+=_SegAccum[t][index];
                avg_weight+=_SegNorm[t][index];
            }
            
            if(avg_weight>0.0)
            {
                SegSubject[index] = avg_seg/avg_weight;
                PatchCount[index] = avg_weight;
            } else {
                SegSubject[index] = -1;
                notfinished++;
            }
          }              
        }
      }
    }
    if ( notfinished>0 ) {
      /*relax preselection criteria*/
      int count=0;
      threshold=threshold*0.95;
      mincount=mincount*0.95;
#ifdef INCREASE_RADIUS      
      v=v+1;
#endif
      
      #pragma omp parallel for  shared(SegSubject,PatchCount,localmask) reduction(+:count) 
      for(i=0;i<dims[0];i++)
      {
        int j,k;
        for(j=0;j<dims[1];j++)
        {
          for(k=0;k<dims[2];k++)
          {
            int index=i*(dims[2]*dims[1])+(j*dims[2])+k;
            
            if ( SegSubject[index]<0 ){
              localmask[index] = 1;
              count++;
            } else {
              localmask[index] = 0;
            }
          }
        }
      }
      
      fprintf(stderr," (redoing %d voxels) t=%f, min=%d \n",count, threshold, mincount);
      notfinished=count;
      sadims = pow( 2*v+1,ndim) * librarysize;
    }
    
    for(i=0;i<omp_get_max_threads();i++)
    {
        free(_PatchImg[i]);
        free(_PatchMask[i]);
    }
    iterations++;
  } while (notfinished && mincount>1 && iterations<MAXITER);

  if(notfinished) /*stopped early*/
  {
      /*TODO: report?*/
  }
  
  for(i=0;i<omp_get_max_threads();i++)
  {
    free(_PatchTemp[i]);
    free(_SegAccum[i]);
    free(_SegNorm[i]);
  }
  
  free(_PatchImg);
  free(_PatchTemp);
  free(_PatchMask);
  free(_SegAccum);
  free(_SegNorm);
  
  
  fprintf(stderr," done (t=%f, min=%d)\n",threshold, mincount);
  
  free(MeansSubj);
  free(VarsSubj);
  free(localmask);
  
  return max;  
}

#else
float nlmsegSparse4D(float *subject, float *imagedata, 
                     float *maskdata, float *meandata, float *vardata, 
                     float *mask, 
                     int sizepatch, int searcharea, float beta, float threshold, 
                     int dims[3],  int librarysize, float *SegSubject, float *PatchCount)
{
    fprintf(stderr,"Compiled without SPAMS!\n");
    return 0.0;
}


float nlmsegSparse4D_double(const float *subject,const  float *imagedata, 
                    const float *maskdata, const float *meandata,const  float *vardata, 
                    const float *mask, 
                    int sizepatch, int searcharea, double beta, double threshold, 
                    const int dims[3],   int librarysize, float *SegSubject, float *PatchCount,
                    double lambda1, double lambda2,int sparse_mode, int stride)
{
    fprintf(stderr,"Compiled without SPAMS!\n");
    return 0.0;
}
#endif //USE_SPAMS



/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
