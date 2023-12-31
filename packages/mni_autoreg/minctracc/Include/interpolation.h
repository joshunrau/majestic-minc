#ifndef MINCTRACC_INTERPOLATION_H
#define MINCTRACC_INTERPOLATION_H

#include "minctracc_arg_data.h"
#include "minctracc_point_vector.h"


int trilinear_interpolant(VIO_Volume volume, 
                                 PointR *coord, double *result);

int tricubic_interpolant(VIO_Volume volume, 
                                PointR *coord, double *result);

void do_Ncubic_interpolation(VIO_Volume volume, 
                                    long index[], int cur_dim, 
                                    double frac[], double *result);

int nearest_neighbour_interpolant(VIO_Volume volume, 
                                         PointR *coord, double *result);

/* A point is not masked if it is a point we should consider.
   If the mask volume is NULL, we consider all points.
   Otherwise, consider a point if the mask volume value is > 0.
*/
int point_not_masked(VIO_Volume volume, 
                            VIO_Real wx, VIO_Real wy, VIO_Real wz);


int voxel_point_not_masked(VIO_Volume volume, 
                                  VIO_Real vx, VIO_Real vy, VIO_Real vz);


/*VF: not sure where to put this*/
void normalize_data_to_match_target(VIO_Volume d1, VIO_Volume m1, VIO_Real thresh1,
                                           VIO_Volume d2, VIO_Volume m2, VIO_Real thresh2,
                                           Arg_Data *globals);
#endif
