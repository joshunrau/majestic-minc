#include "mniVolume.h"

// blank constructor - just initialises the sizes for now
mniVolume::mniVolume() {
  this->sizes = new int[VIO_MAX_DIMENSIONS];
}

// Constructor from file
mniVolume::mniVolume(VIO_STR filename, 
                     VIO_Real voxelMin, 
                     VIO_Real voxelMax,
                     int nDimensions,
                     VIO_STR dimensions[],
                     nc_type dataType,
                     VIO_BOOL volumeSigned,
                     VIO_BOOL createVolume,
                     minc_input_options *options ) {
  if ( input_volume(filename, nDimensions, dimensions, dataType, 
                    volumeSigned, voxelMin,
                    voxelMax, createVolume, &this->volume, options) != VIO_OK ) {
    throw loadException();
  }

  this->sizes = new int[VIO_MAX_DIMENSIONS];
  get_volume_sizes(this->volume, this->sizes);
  this->nDimensions = nDimensions;
  this->filename = filename;
  this->dimNames = dimensions;
  this->voxelMin = voxelMin;
  this->voxelMax = voxelMax;
  this->dataType = dataType;
  this->signedFlag = volumeSigned;

}

mniVolume::mniVolume(VIO_Volume volumeIO_volume) {

  this->sizes = new int[VIO_MAX_DIMENSIONS];
  this->volume = volumeIO_volume;
  this->voxelMin = get_volume_voxel_min(volumeIO_volume);
  this->voxelMax = get_volume_voxel_max(volumeIO_volume);
  this->dataType = get_volume_nc_data_type(volumeIO_volume, 
					   &this->signedFlag);
  get_volume_sizes(volumeIO_volume, this->sizes);
  this->dimNames = get_volume_dimension_names(volumeIO_volume);
}

mniVolume::mniVolume(mniBaseVolume *copyVolume, 
		     VIO_BOOL copyVolumeDefinitionOnly,
		     nc_type dataType,
		     VIO_BOOL signedFlag,
		     VIO_Real voxelMin,
		     VIO_Real voxelMax) {


  //initialise sizes
  this->sizes = new int[VIO_MAX_DIMENSIONS];
  
  // now copy relevant bits from other volume
  if (copyVolumeDefinitionOnly == TRUE) {
    this->volume = copy_volume_definition(copyVolume->getVolume(), dataType,
					  signedFlag, voxelMin, voxelMax);
    this->dataType = dataType;
    this->signedFlag = signedFlag;
    this->voxelMin = voxelMin;
    this->voxelMax = voxelMax;
  }
  else { // create an exact copy
    this->volume = copy_volume(copyVolume->getVolume());
    this->voxelMin = copyVolume->getVoxelMin();
    this->voxelMax = copyVolume->getVoxelMax();
    this->signedFlag = copyVolume->getSignedFlag();
    this->dataType = copyVolume->getDataType();
  }

    this->sizes = copyVolume->getSizes();
    this->dimNames = copyVolume->getDimNames();
}
  
  

mniVolume::~mniVolume() {
  delete_volume(this->volume);
  delete this->sizes;
}

VIO_Real mniVolume::getWorld(VIO_Real xWorld, VIO_Real yWorld, VIO_Real zWorld) {
  VIO_Real *voxelCoord;
  cout << "In get world" << endl;
  voxelCoord = this->convertWorldToVoxel(xWorld, yWorld, zWorld);
  cout << "Converted ... " << endl;
  return get_volume_real_value(this->volume, 
                               (int)rint(voxelCoord[0]), 
                               (int)rint(voxelCoord[1]),
                               (int)rint(voxelCoord[2]), 0, 0);
}

void mniVolume::output(VIO_STR file, int cropValue) {
  if (output_volume(file, this->dataType, this->signedFlag,
                    this->voxelMin, this->voxelMax, this->volume,
                    (char *) "mnipl-- test", NULL) != VIO_OK) {
    throw writeException();
  }
}



/*
  mniVolume& mniVolume::operator+(mniVolume *a, mniVolume *b) {
  // check to make sure that sizes are the same
  int *a_sizes = a->getSizes();
  int *b_sizes = b->getSizes();
  if (a_sizes[0] != b_sizes[0] ||
  a_sizes[1] != b_sizes[1] || 
  a_sizes[2] != b_sizes[2])
  throw differentSizesException();
  
  for (int v1 = 0; v1 < a_sizes[0]; v1++) {
  for (int v2 = 0; v2 < a_sizes[1]; v2++) {
  for (int v3 = 0; v3 < a_sizes[2]; v3++) {
*/
