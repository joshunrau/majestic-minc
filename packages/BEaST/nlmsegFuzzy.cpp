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

#define MINCOUNT 100


#ifdef MT_USE_OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 1
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif

/*Can be used to store the distance and the location*/
/*Some compilation issues can occur here, use typedef strcut with window*/
typedef struct {
  double dist;
  int x;
  int y;
  int z;
  int t;
} data_t;


float nlmsegFuzzy4D(const float *subject,const  float *imagedata, 
                    const float *maskdata, const float *meandata,
                    const float *vardata, 
                    const float *mask, 
                    int sizepatch, int searcharea, float beta, float threshold, 
                    const int dims[3],  int librarysize, 
                    float *SegSubject, float *PatchCount)
{    
  float *MeansSubj;
  float *VarsSubj;
  float *localmask;
  int i,v,f,ndim;
  float min,max;
  
  int patch_volume;

  
  int sadims,volumesize,index;
  int mincount=MINCOUNT;
  int notfinished;
  double minidist;
  double epsi = 0.0001;
  time_t time1,time2;
  
  /*per thread temp storage*/
  data_t **_tab       =(data_t **)malloc(omp_get_max_threads()*sizeof(data_t*));
  float  **_PatchImg  =(float **) malloc(omp_get_max_threads()*sizeof(float*));
  float  **_PatchTemp =(float **) malloc(omp_get_max_threads()*sizeof(float*));
  

  fprintf(stderr,"Patch size: %d\nSearch area: %d\nBeta: %f\nThreshold: %f\nSelection: %d\n",sizepatch,searcharea,beta,threshold,librarysize);
  
  ndim = 3;
  volumesize=dims[0]*dims[1]*dims[2];
  
  /*Patch radius*/
  f = sizepatch;
  
  /*volume of patch*/
  patch_volume=(2*f+1)*(2*f+1)*(2*f+1);
  
  
  /*Search Area radius*/
  v = searcharea;
  
  sadims = pow(2*v+1,ndim);
  sadims = sadims * librarysize;

  /* allocate memory for multithreading operation*/
  for(i=0;i<omp_get_max_threads();i++)
  {
    _PatchImg[i] =(float*) malloc( (2*f+1)*(2*f+1)*(2*f+1)*sizeof(float) );
    _PatchTemp[i]=(float*) malloc( (2*f+1)*(2*f+1)*(2*f+1)*sizeof(float) );
  }
  
  
  MeansSubj = (float *)calloc(volumesize,sizeof(*MeansSubj));
  VarsSubj =  (float *)calloc(volumesize,sizeof(*VarsSubj));
  localmask = (float *)calloc(volumesize,sizeof(*localmask));
  
  memmove(localmask,mask,volumesize*sizeof(*localmask));
  
  fprintf(stderr,"Dimensions: %d %d %d\n",dims[0],dims[1],dims[2]);
  
  fprintf(stderr,"Computing first moment image...");
  time1=time(NULL);
  ComputeFirstMoment(subject, MeansSubj, dims, f, &min, &max);
  time2=time(NULL);
  fprintf(stderr,"done (%d sec)\nComputing second moment image...",(int)(time2-time1));
  ComputeSecondMoment(subject, MeansSubj, VarsSubj, dims, f, &min, &max);
  fprintf(stderr,"done");
  time1=time(NULL);
  
  do {
    fprintf(stderr," (%d sec)\nSegmenting      ",(int)(time1-time2));
    time2=time(NULL);
    notfinished=0;
    
    for(i=0;i<omp_get_max_threads();i++)
    {
      _tab[i]=(data_t*)malloc(sadims*sizeof(data_t));
    }
    
    #pragma omp parallel for shared(_tab,_PatchImg,_PatchTemp) reduction(+:notfinished)
    for(i=0;i<dims[0];i++)
    { /*start parallel section*/
      int j,k;
      
      float * PatchImg;
      float * PatchTemp;
      
      data_t  *tab;
      
      if( omp_get_thread_num()==0 )
        fprintf(stderr,"\b\b\b\b\b\b\b\b\b%3d / %3d", i*omp_get_num_threads()+1, dims[0]);
      /*use thread-specific temp memory*/
      
      tab          = _tab      [omp_get_thread_num()];
      PatchImg     = _PatchImg [omp_get_thread_num()];
      PatchTemp    = _PatchTemp[omp_get_thread_num()];
      
      for(j=0;j<dims[1];j++)
      {
        for(k=0;k<dims[2];k++)
        {
          int index=i*(dims[2]*dims[1])+(j*dims[2])+k;
          
          /* mask check */
          if ( localmask[index]  > 0 )
          {
            int ii,jj,kk;
            float proba=0.;
            float average=0;
            float totalweight=0;
            int   count = 0;
            float minidist = FLT_MAX; /*FLT_MAX;*/
            float TMean,TVar;
            
            //memset(PatchTemp,0,sizeof(float)*patch_volume);
            ExtractPatch(subject, PatchTemp, i, j, k, f, dims[0], dims[1], dims[2]);
            
            TMean = MeansSubj[index];
            TVar =  VarsSubj [index];
            
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
                  
                  if((ni>=0) && (nj>=0) && (nk>=0) && (ni<dims[0]) && (nj<dims[1]) && (nk<dims[2]))
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
                        float d;
                        data_t storage;
                        
                        //memset(PatchImg,0,sizeof(float)*patch_volume);
                        ExtractPatch4D(imagedata, PatchImg ,ni,nj,nk, t,f, dims[0], dims[1], dims[2]);
                        
                        d =  SSDPatch(PatchImg, PatchTemp, f);
                        
                        if (d < minidist) minidist = d;
                        
                        storage.dist  = d;
                        storage.z = ni;
                        storage.y = nj;
                        storage.x = nk;
                        storage.t = t;
                          
                        tab[count] = storage;
                        count ++;
                        
                      }
                    }
                  }
                }
              }
            }
            
            /* require a minimum number of selected patches  */
            if (count >= mincount) {
              int realcount=0;
              int p=0;
              /* Sort the distance*/
              /*This can be removed according to the chosen strategy*/
              /*qsort (tab, count, sizeof *tab, cmp);*/
              
              /*You can use the closest Patches (i.e. the 'lim' closest ones) or all the preselected patches (count)*/
              //lim = count; /*in this case, you take all the preselected patches into account*/
              
              if ( minidist<=epsi ) minidist = epsi; /*to avoid division by zero*/
                
                while (p < count)
                {
                  data_t storage = tab[p];
                  float w = exp(- ((storage.dist)/(beta*(minidist)) ) ); /*The smoothing parameter is the minimal distance*/
                  
                  if (w>0.0)  
                  {
                    average  = average + maskdata[(storage.t*volumesize)+(storage.z*(dims[2]*dims[1]))+(storage.y*dims[2])+storage.x]*w;			      
                    totalweight = totalweight + w;
                    realcount++;
                  }
                  
                  p++;
                } // while
                
                /* We compute the probability */
                proba = average / totalweight;
                
                SegSubject[index] = proba;
                PatchCount[index] = realcount;
                
            } else {
              /* Not enough similar patches */
              notfinished+=1;
              SegSubject[index] = -1;
            }
          }// mask check
        } // for k
      } // for j
      /*Freeing per-thread data*/
      /*end of parallel section*/
    } // for i
    
    time1=time(NULL);
    
    if ( notfinished>0 ) {
      /*relax preselection criteria*/
      int count;
      threshold=threshold*0.95;
      mincount=mincount*0.95;
      v=v+1;
      count=0;
      
      #pragma omp parallel for reduction(+:count)
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
      
      fprintf(stderr," (redoing %d voxels) t=%f, min=%d ",count, threshold, mincount);
      sadims = pow(2*v+1,ndim);
      sadims = sadims * librarysize;
    }
    
    for(i=0;i<omp_get_max_threads();i++)
      free(_tab[i]);

  } while (notfinished);

  for(i=0;i<omp_get_max_threads();i++)
  {
    free(_PatchImg[i]);
    free(_PatchTemp[i]);
  }
  
  free(_PatchImg);
  free(_PatchTemp);
  free(_tab);
  
  fprintf(stderr," done (%d sec, t=%f, min=%d)\n",(int)(time1-time2), threshold, mincount);
  
  free(MeansSubj);
  free(VarsSubj);
  free(localmask);
  
  return max;
}



float nlmsegFuzzy4D_double(const float *subject,const  float *imagedata, 
                    const float *maskdata, const float *meandata,
                    const float *vardata, 
                    const float *mask, 
                    int sizepatch, int searcharea, 
                    double beta, double threshold, 
                    const int dims[3], int librarysize, 
                    float *SegSubject, float *PatchCount)
{    
  float *MeansSubj;
  float *VarsSubj;
  float *localmask;
  int i,v,f,ndim;
  float min,max;
  
  int patch_volume;

  
  int sadims,volumesize,index;
  int mincount=MINCOUNT;
  int notfinished;
  double minidist;
  double epsi = 0.0001;
  time_t time1,time2;
  
  /*per thread temp storage*/
  data_t **_tab       =(data_t **)malloc(omp_get_max_threads()*sizeof(data_t*));
  float  **_PatchImg  =(float **) malloc(omp_get_max_threads()*sizeof(float*));
  float  **_PatchTemp =(float **) malloc(omp_get_max_threads()*sizeof(float*));
  

  fprintf(stderr,"Patch size: %d\nSearch area: %d\nBeta: %f\nThreshold: %f\nSelection: %d\n",sizepatch,searcharea,beta,threshold,librarysize);
  
  ndim = 3;
  volumesize=dims[0]*dims[1]*dims[2];
  
  /*Patch radius*/
  f = sizepatch;
  
  /*volume of patch*/
  patch_volume=(2*f+1)*(2*f+1)*(2*f+1);
  
  
  /*Search Area radius*/
  v = searcharea;
  
  sadims = pow(2*v+1,ndim);
  sadims = sadims * librarysize;

  /* allocate memory for multithreading operation*/
  for(i=0;i<omp_get_max_threads();i++)
  {
    _PatchImg[i] =(float*) malloc( (2*f+1)*(2*f+1)*(2*f+1)*sizeof(float) );
    _PatchTemp[i]=(float*) malloc( (2*f+1)*(2*f+1)*(2*f+1)*sizeof(float) );
  }
  
  
  MeansSubj = (float *)calloc(volumesize,sizeof(*MeansSubj));
  VarsSubj =  (float *)calloc(volumesize,sizeof(*VarsSubj));
  localmask = (float *)calloc(volumesize,sizeof(*localmask));
  
  memmove(localmask,mask,volumesize*sizeof(*localmask));
  
  fprintf(stderr,"Dimensions: %d %d %d\n",dims[0],dims[1],dims[2]);
  
  fprintf(stderr,"Computing first moment image...");
  time1=time(NULL);
  ComputeFirstMoment(subject, MeansSubj, dims, f, &min, &max);
  time2=time(NULL);
  fprintf(stderr,"done (%d sec)\nComputing second moment image...",(int)(time2-time1));
  ComputeSecondMoment(subject, MeansSubj, VarsSubj, dims, f, &min, &max);
  fprintf(stderr,"done");
  time1=time(NULL);
  
  do {
    fprintf(stderr," (%d sec)\nSegmenting      ",(int)(time1-time2));
    time2=time(NULL);
    notfinished=0;
    
    for(i=0;i<omp_get_max_threads();i++)
    {
      _tab[i]=(data_t*)malloc(sadims*sizeof(data_t));
    }
    
    #pragma omp parallel for shared(_tab,_PatchImg,_PatchTemp) reduction(+:notfinished)
    for(i=0;i<dims[0];i++)
    { /*start parallel section*/
      int j,k;
      
      float * PatchImg;
      float * PatchTemp;
      
      data_t  *tab;
      
      if( omp_get_thread_num()==0 )
        fprintf(stderr,"\b\b\b\b\b\b\b\b\b%3d / %3d", i*omp_get_num_threads()+1, dims[0]);
      /*use thread-specific temp memory*/
      
      tab          = _tab      [omp_get_thread_num()];
      PatchImg     = _PatchImg [omp_get_thread_num()];
      PatchTemp    = _PatchTemp[omp_get_thread_num()];
      
      for(j=0;j<dims[1];j++)
      {
        for(k=0;k<dims[2];k++)
        {
          int index=i*(dims[2]*dims[1])+(j*dims[2])+k;
          
          /* mask check */
          if ( localmask[index]  > 0 )
          {
            int ii,jj,kk;
            double proba=0.;
            double average=0;
            double totalweight=0;
            int   count = 0;
            double minidist = DBL_MAX; /*FLT_MAX;*/
            double TMean,TVar;
            
            //memset(PatchTemp,0,sizeof(float)*patch_volume);
            ExtractPatch(subject, PatchTemp, i, j, k, f, dims[0], dims[1], dims[2]);
            
            TMean = MeansSubj[index];
            TVar =  VarsSubj [index];
            
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
                  
                  if((ni>=0) && (nj>=0) && (nk>=0) && (ni<dims[0]) && (nj<dims[1]) && (nk<dims[2]))
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
                        double d;
                        data_t storage;
                        
                        //memset(PatchImg,0,sizeof(float)*patch_volume);
                        ExtractPatch4D(imagedata, PatchImg ,ni,nj,nk, t,f, dims[0], dims[1], dims[2]);
                        
                        d =  SSDPatch_double(PatchImg, PatchTemp, f);
                        
                        if (d < minidist) minidist = d;
                        
                        storage.dist  = d;
                        storage.z = ni;
                        storage.y = nj;
                        storage.x = nk;
                        storage.t = t;
                          
                        tab[count] = storage;
                        count ++;
                        
                      }
                    }
                  }
                }
              }
            }
            
            /* require a minimum number of selected patches  */
            if (count >= mincount) {
              int realcount=0;
              int p=0;
              /* Sort the distance*/
              /*This can be removed according to the chosen strategy*/
              /*qsort (tab, count, sizeof *tab, cmp);*/
              
              /*You can use the closest Patches (i.e. the 'lim' closest ones) or all the preselected patches (count)*/
              //lim = count; /*in this case, you take all the preselected patches into account*/
              
              if ( minidist<=epsi ) minidist = epsi; /*to avoid division by zero*/
                
                while (p < count)
                {
                  data_t storage = tab[p];
                  double w = exp(- ((storage.dist)/(beta*(minidist)) ) ); /*The smoothing parameter is the minimal distance*/
                  
                  if (w>0.0)  
                  {
                    average  = average + maskdata[(storage.t*volumesize)+(storage.z*(dims[2]*dims[1]))+(storage.y*dims[2])+storage.x]*w;			      
                    totalweight = totalweight + w;
                    realcount++;
                  }
                  
                  p++;
                } // while
                
                /* We compute the probability */
                proba = average / totalweight;
                
                SegSubject[index] = proba;
                PatchCount[index] = realcount;
                
            } else {
              /* Not enough similar patches */
              notfinished+=1;
              SegSubject[index] = -1;
            }
          }// mask check
        } // for k
      } // for j
      /*Freeing per-thread data*/
      /*end of parallel section*/
    } // for i
    
    time1=time(NULL);
    
    if ( notfinished>0 ) {
      /*relax preselection criteria*/
      int count;
      threshold=threshold*0.95;
      mincount=mincount*0.95;
      v=v+1;
      count=0;
      
      #pragma omp parallel for reduction(+:count)
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
      
      fprintf(stderr," (redoing %d voxels) t=%f, min=%d ",count, threshold, mincount);
      sadims = pow(2*v+1,ndim);
      sadims = sadims * librarysize;
    }
    
    for(i=0;i<omp_get_max_threads();i++)
      free(_tab[i]);

  } while (notfinished);

  for(i=0;i<omp_get_max_threads();i++)
  {
    free(_PatchImg[i]);
    free(_PatchTemp[i]);
  }
  
  free(_PatchImg);
  free(_PatchTemp);
  free(_tab);
  
  fprintf(stderr," done (%d sec, t=%f, min=%d)\n",(int)(time1-time2), threshold, mincount);
  
  free(MeansSubj);
  free(VarsSubj);
  free(localmask);
  
  return max;
}


