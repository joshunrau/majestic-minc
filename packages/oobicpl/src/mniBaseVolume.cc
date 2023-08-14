#include "mniBaseVolume.h"

VIO_Real* mniBaseVolume::convertVoxelToWorld(VIO_Real voxel[]) {
  VIO_Real *returnValue = new VIO_Real[3];
  convert_voxel_to_world(this->volume, 
                         voxel, 
                         &returnValue[0],
                         &returnValue[1],
                         &returnValue[2]);
  
  return returnValue;
}

VIO_Real* mniBaseVolume::convertWorldToVoxel(VIO_Real xWorld,
                                         VIO_Real yWorld,
                                         VIO_Real zWorld) 
{
  VIO_Real *voxel = new VIO_Real[3];
  convert_world_to_voxel(this->volume, xWorld, yWorld, zWorld, voxel);
  return voxel;
}

