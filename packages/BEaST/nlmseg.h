/*  nlmseg.h
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


#ifndef NLMSEG_H
#define NLMSEG_H

void ComputeFirstMoment(const float* ima, float* means, const int* dims, int f, float *min, float *max);
void ComputeSecondMoment(const float* ima,const  float* means, float* variance, const int* dims,int f, float *min, float *max);
void ComputeFirstMoment4D(const float* ima,const float* atlas, float* means, float* Ameans,const int* dims, int f);
void ComputeSecondMoment4D(const float* ima,const float* atlas,const  float* means,const  float* Ameans, float* variance, float* Avariance,const int* dims,int f);

float SSDPatch(const float*  PatchImg, const float*  PatchTemplate, int f);
double SSDPatch_double(const float*  PatchImg, const float*  PatchTemplate, int f);

void ExtractPatch(const float*  ima, float*  Patch, int x, int y, int z, int size, int sx, int sy, int sz);
void ExtractPatch4D(const float*  ima, float*  Patch, int x,int y, int z, int t,int size,int sx,int sy,int sz);

void ExtractPatch_norm(const float*  ima, float*  Patch, int x, int y, int z, int size, int sx, int sy, int sz,float mean);
void ExtractPatch4D_norm(const float*  ima, float*  Patch, int x,int y, int z, int t,int size,int sx,int sy,int sz,float mean);


void ExtractPatch_double(const float*  ima, double*  Patch, int x, int y, int z, int size, int sx, int sy, int sz);
void ExtractPatch4D_double(const float*  ima, double*  Patch, int x,int y, int z, int t,int size,int sx,int sy,int sz);

void ExtractPatch_norm_double(const float*  ima, double*  Patch, int x, int y, int z, int size, int sx, int sy, int sz,double mean);
void ExtractPatch4D_norm_double(const float*  ima, double*  Patch, int x,int y, int z, int t,int size,int sx,int sy,int sz,double mean);



void AddWPatch(float*  ima,const float*  Patch, float w, int x, int y, int z, int size, int sx, int sy, int sz);
void AddW(float*  ima,float w,                          int x, int y, int z, int size, int sx, int sy, int sz);

void AddWPatch_double(double*  ima,const double*  Patch, double w, int x, int y, int z, int size, int sx, int sy, int sz);
void AddW_double(double*  ima,double w,                          int x, int y, int z, int size, int sx, int sy, int sz);



float nlmsegFuzzy4D(const float *subject, const float *imagedata,const  float *maskdata, const float *meandata, 
                    const float *vardata, const float *mask, int sizepatch, int searcharea, 
                    float alpha, float threshold, const int sizes[3], 
                    int librarysize, float *SegSubject, float *PatchCount);

float nlmsegFuzzy4D_double(const float *subject,const  float *imagedata, 
                    const float *maskdata, const float *meandata,
                    const float *vardata, 
                    const float *mask, 
                    int sizepatch, int searcharea, 
                    double beta, double threshold, 
                    const int dims[3], int librarysize, 
                    float *SegSubject, float *PatchCount);

float nlmsegSparse4D(const float *subject,const  float *imagedata, 
                    const float *maskdata, const float *meandata,const  float *vardata, 
                    const float *mask, 
                    int sizepatch, int searcharea, float beta, float threshold, 
                    const int dims[3],   int librarysize, float *SegSubject, float *PatchCount,
                    float lambda1, float lambda2,int sparse_mode, int stride);

float nlmsegSparse4D_double(const float *subject,const  float *imagedata, 
                    const float *maskdata, const float *meandata,const  float *vardata, 
                    const float *mask, 
                    int sizepatch, int searcharea, double beta, double threshold, 
                    const int dims[3],   int librarysize, float *SegSubject, float *PatchCount,
                    double lambda1, double lambda2,int sparse_mode, int stride);

float nlmfilter(const float *subject,const float *mask,const float *maskdata, 
                int sizepatch, int searcharea, 
                float beta, float threshold,const int dims[3], 
                float *SegSubject, float *PatchCount);

#endif

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
