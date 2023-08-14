/*  mincbeast.c
 *
 *  Copyright 2011  Simon Fristed Eskildsen, Vladimir Fonov,
 *                Pierrick Coupé, Jose V. Manjon
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

#include "ParseArgv.h"
#include "array_alloc.h"
#include "nlmseg.h"
#include "beast.h"
#include "label.h"

const char LICENSE[]="Copyright (C) 2011\tSimon Fristed Eskildsen, Vladimir S. Fonov, \n\
\t\t\tPierrick Coupe, Jose V. Manjon\n\n\
This program comes with ABSOLUTELY NO WARRANTY; for details type 'cat COPYING'. \n\
This is free software, and you are welcome to redistribute it under certain\n\
conditions; type 'cat COPYING' for details.\n\
";

const char REFERENCE[]="Reference: \n\tEskildsen SF, Coupe P, Fonov V, Manjon JV, Leung KK,\n\tGuizard N, Wassef SN, Ostergaard LR, Collins DL;\n\tAlzheimer's Disease Neuroimaging Initiative.\n\
\tBEaST: brain extraction based on nonlocal segmentation technique.\n\
\tNeuroimage. 2012 Feb 1;59(3):2362-73.\n\
\thttp://dx.doi.org/10.1016/j.neuroimage.2011.09.012\n\n";

#ifdef MT_USE_OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 1
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif

int main(int argc, char  *argv[] )
{
  beast_options _options;
  int ret=STATUS_OK;
  
  fprintf(stderr,"\nmincbeast --\t\tan implementation of BEaST (Brain Extraction\n\t\t\tusing non-local Segmentation Technique) version %s\n\n",PACKAGE_VERSION);
  
#ifdef MT_USE_OPENMP
  fprintf(stderr,"Using OpenMP, max number of threads=%d\n",omp_get_max_threads());
#endif
  
  
  if( get_arguments(argc, argv, &_options)!=STATUS_OK )
    return STATUS_ERR;
  
  if(_options.v2)
    ret=mincbeast_v2( &_options );
  else
    ret=mincbeast_v1( &_options );
  
  cleanup_arguments(&_options);

  return ret;
}

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
