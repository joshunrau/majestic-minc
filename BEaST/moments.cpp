/*  moments.c
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


#include <float.h>
#include "basic.h"
#include "nlmseg.h"

void ComputeFirstMoment(const float* ima, float* means, const int* dims, int f, float *min, float *max)
{
  int i;
  float vmin=FLT_MAX;
  float vmax=-FLT_MAX;
  
#if _OPENMP >=  201107
  #pragma omp parallel for reduction(max:vmax) reduction(min:vmin)
#endif  
  for(i=0;i<dims[0];i++)
  {
    int j,k;
    
    for(j=0;j<dims[1];j++)
    {
      for(k=0;k<dims[2];k++)
      {
        float mean = 0.0;
        int indice=0;
        int ii,jj,kk;
        
        for(ii=-f;ii<=f;ii++)
        {
          for(jj=-f;jj<=f;jj++)
          {
            for(kk=-f;kk<=f;kk++)
            {
              int ni,nj,nk;
              
              ni=i+ii;
              nj=j+jj;
              nk=k+kk;
              
              if(ni<0) ni=-ni;
              if(nj<0) nj=-nj;
              if(nk<0) nk=-nk;
              if(ni>=dims[0]) ni=2*dims[0]-ni-1;
              if(nj>=dims[1]) nj=2*dims[1]-nj-1;
              if(nk>=dims[2]) nk=2*dims[2]-nk-1;
              
              
              mean += ima[ni*(dims[2]*dims[1])+(nj*dims[2])+nk];			    
              
              indice+=1;
              
            }
          }
        }
        
        mean=mean/indice;
        
        vmin=MIN(vmin,mean);
        vmax=MAX(vmax,mean);
        
        means[i*(dims[2]*dims[1])+(j*dims[2])+k]=mean;
        
      }
    }
  }
  
  *min=vmin;
  *max=vmax;
}


void ComputeSecondMoment(const float* ima,const  float* means, float* variance, const int* dims,int f, float *min, float *max)
{
  int i;
  float vmin=FLT_MAX,vmax=-FLT_MAX;

#if _OPENMP >=  201107
  #pragma omp parallel for reduction(max:vmax) reduction(min:vmin)
#endif  
  for(i=0;i<dims[0];i++)
  {
    int j,k;
    for(j=0;j<dims[1];j++)
    {
      for(k=0;k<dims[2];k++)
      {
        int indice =0;
        float var=0.0;
        int ii,jj,kk;
        int index1=i*(dims[2]*dims[1])+(j*dims[2])+k;
        
        for(ii=-f;ii<=f;ii++)
        {
          for(jj=-f;jj<=f;jj++)
          {
            for(kk=-f;kk<=f;kk++)
            {
              int ni,nj,nk;
              ni=i+ii;
              nj=j+jj;
              nk=k+kk;
              if(ni>=0 && nj>=0 && nk>0 && ni<dims[0] && nj<dims[1] && nk<dims[2])
              {
                int index2=ni*(dims[2]*dims[1])+(nj*dims[2])+nk;
                var +=  ((ima[index2]-means[index1])*(ima[index2]-means[index1]));
                indice+=1;
              }
            }
          }
        }
        var=var/(indice-1);
        
        vmin=MIN(vmin,var);
        vmax=MAX(vmax,var);
        
        variance[index1]=var;
        
      }
    }
  }

  *min=vmin;
  *max=vmax;
}


void ComputeFirstMoment4D(const float* ima,const float* atlas, float* means, float* Ameans,const int* dims, int f)
{
  int k,t;
  
  #pragma omp parallel for 
  for(k=0;k<dims[2];k++)
  {
    int i,j;
    for(i=0;i<dims[1];i++)
    {
      for(j=0;j<dims[0];j++)
      {
        int indice=0;
        float Amean = 0.0;
        int ii,jj,kk;
        
        for(ii=-f;ii<=f;ii++)
        {
          for(jj=-f;jj<=f;jj++)
          {
            for(kk=-f;kk<=f;kk++)
            {
              int ni,nj,nk;
              
              ni=i+ii;
              nj=j+jj;
              nk=k+kk;
              
              if(ni<0) ni=-ni;
              if(nj<0) nj=-nj;
              if(nk<0) nk=-nk;
              if(ni>=dims[1]) ni=2*dims[1]-ni-1;
              if(nj>=dims[0]) nj=2*dims[0]-nj-1;
              if(nk>=dims[2]) nk=2*dims[2]-nk-1;
              
              
              Amean += atlas[nk*(dims[0]*dims[1])+(ni*dims[0])+nj];
              indice+=1;
              
            }
          }
        }
        
        Amean=Amean/indice;
        Ameans[k*(dims[0]*dims[1])+(i*dims[0])+j]=Amean;
      }
    }
  }

  for(t=0;t<dims[3];t++)
  {
    #pragma omp parallel for 
    for(k=0;k<dims[2];k++)
    {
      int i,j;
      for(i=0;i<dims[1];i++)
      {
        for(j=0;j<dims[0];j++)
        {
          int indice=0;
          float mean = 0.0;
          int ii,jj,kk;
          
          for(ii=-f;ii<=f;ii++)
          {
            for(jj=-f;jj<=f;jj++)
            {
              for(kk=-f;kk<=f;kk++)
              {
                int ni,nj,nk;
                ni=i+ii;
                nj=j+jj;
                nk=k+kk;
                
                if(ni<0) ni=-ni;
                if(nj<0) nj=-nj;
                if(nk<0) nk=-nk;
                if(ni>=dims[1]) ni=2*dims[1]-ni-1;
                if(nj>=dims[0]) nj=2*dims[0]-nj-1;
                if(nk>=dims[2]) nk=2*dims[2]-nk-1;
                
                mean   +=  ima[t*(dims[0]*dims[1]*dims[2])+(nk*(dims[0]*dims[1])+(ni*dims[0])+nj)];
                indice +=1;
              }
            }
          }
          
          mean=mean/indice;
          means[t*(dims[0]*dims[1]*dims[2])+(k*(dims[0]*dims[1])+(i*dims[0])+j)]=mean;
          
        }
      }
    }
  }
}

void ComputeSecondMoment4D(const float* ima,const float* atlas,const float* means, const float* Ameans, float* variance, float* Avariance,const int* dims,int f)
{    
  int k,t;
  
  #pragma omp parallel for 
  for(k=0;k<dims[2];k++)
  {
    int i,j;
    for(i=0;i<dims[1];i++)
    {
      for(j=0;j<dims[0];j++)
      {
        int ii,jj,kk;
        float Avar=0.0;
        int indice =0;
        
        for(ii=-f;ii<=f;ii++)
        {
          for(jj=-f;jj<=f;jj++)
          {
            for(kk=-f;kk<=f;kk++)
            {
              int ni,nj,nk;
              
              ni=i+ii;
              nj=j+jj;
              nk=k+kk;
              if(ni>=0 && nj>=0 && nk>0 && ni<dims[1] && nj<dims[0] && nk<dims[2])
              {
                Avar += ((atlas[nk*(dims[0]*dims[1])+(ni*dims[0])+nj]-Ameans[k*(dims[0]*dims[1])+(i*dims[0])+j])*(atlas[nk*(dims[0]*dims[1])+(ni*dims[0])+nj]-Ameans[k*(dims[0]*dims[1])+(i*dims[0])+j]));
                indice+=1;
              }
            }
          }
        }
        
        Avar=Avar/(indice-1);
        Avariance[k*(dims[0]*dims[1])+(i*dims[0])+j]=Avar;
      }
    }
  }
  
  for(t=0;t<dims[3];t++)
  {
    #pragma omp parallel for 
    for(k=0;k<dims[2];k++)
    {
      int i,j;
      for(i=0;i<dims[1];i++)
      {
        for(j=0;j<dims[0];j++)
        {
          int ii,jj,kk;
          float var=0.0;
          int indice =0;
          
          for(ii=-f;ii<=f;ii++)
          {
            for(jj=-f;jj<=f;jj++)
            {
              for(kk=-f;kk<=f;kk++)
              {
                int ni,nj,nk;
                
                ni=i+ii;
                nj=j+jj;
                nk=k+kk;
                if(ni>=0 && nj>=0 && nk>0 && ni<dims[1] && nj<dims[0] && nk<dims[2])
                {
                  var += ((ima[t*(dims[0]*dims[1]*dims[2])+nk*(dims[0]*dims[1])+(ni*dims[0])+nj]-means[t*(dims[0]*dims[1]*dims[2])+k*(dims[0]*dims[1])+(i*dims[0])+j])*(ima[t*(dims[0]*dims[1]*dims[2])+nk*(dims[0]*dims[1])+(ni*dims[0])+nj]-means[t*(dims[0]*dims[1]*dims[2])+k*(dims[0]*dims[1])+(i*dims[0])+j]));
                  indice+=1;
                }
              }
            }
          }
          
          var=var/(indice-1);
          variance[t*(dims[0]*dims[1]*dims[2])+k*(dims[0]*dims[1])+(i*dims[0])+j]=var;
        }
      }
    }
  }
}
