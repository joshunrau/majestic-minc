#include "mniLabelVolume.h"

mniLabelVolume::mniLabelVolume(VIO_STR filename, 
			       VIO_Real voxelMin, 
			       VIO_Real voxelMax,
			       int nDimensions,
			       VIO_STR dimensions[],
			       nc_type dataType,
			       VIO_BOOL volumeSigned,
			       VIO_BOOL createVolume,
			       minc_input_options *options) {

  // initialise the sizes variable
  this->sizes = new int[VIO_MAX_DIMENSIONS];

  // load the volume headers first
  if ( input_volume_header_only(filename, nDimensions, dimensions,
				&this->volume, options) != VIO_OK ) {
    throw loadException();
  }
  // now create it as a label volume
  if ( create_label_volume_from_file(filename, this->volume, &this->volume)
       != VIO_OK ) {
    throw loadException();
  }

  get_volume_sizes(this->volume, this->sizes);
                               
}

mniLabelVolume::mniLabelVolume(VIO_Volume volumeIO_volume) {
  this->sizes = new int[VIO_MAX_DIMENSIONS];
  this->volume = volumeIO_volume;
  get_volume_sizes(this->volume, this->sizes);
}
                               
mniLabelVolume::mniLabelVolume(mniBaseVolume *copyVolume, 
			       nc_type dataType ) {

  // initialise sizes
  this->sizes = new int[VIO_MAX_DIMENSIONS];

  // now copy all the relevant bits from the other volume
  this->volume = create_label_volume(copyVolume->getVolume(), dataType);
  this->dimNames = copyVolume->getDimNames();
  *this->sizes = *copyVolume->getSizes(); // copy by value
  this->dataType = dataType;
}

/*
mniLabelVolume::mniLabelVolume(mniLabelVolume *copyVolume,
			       nc_type dataType = NC_SHORT) {

  // initialise sizes
  this->sizes = new int[VIO_MAX_DIMENSIONS];

  // now copy all the relevant bits from the other volume
  this->volume = create_label_volume(copyVolume->getVolume(), dataType);
  this->dimNames = copyVolume->getDimNames();
  *this->sizes = *copyVolume->getSizes(); // copy by value
  this->dataType = dataType;
}
*/
mniLabelVolume::mniLabelVolume() {
  // initialise sizes
  this->sizes = new int[VIO_MAX_DIMENSIONS];
}

mniLabelVolume::mniLabelVolume(VIO_STR filename,
			       int newVolume,
			       int nDimensions,
			       VIO_STR dimensions[],
			       nc_type dataType,
			       minc_input_options *options){

  // initialise sizes
  this->sizes = new int[VIO_MAX_DIMENSIONS];

  if (input_volume_header_only(filename, nDimensions, dimensions,
                               &this->volume, options) != VIO_OK) {
    throw loadException();
  }
  this->volume = create_label_volume(this->volume, dataType);

  get_volume_sizes(this->volume, this->sizes);
  this->nDimensions = nDimensions;
  this->dimNames = dimensions;
  this->filename = filename;
}

mniLabelVolume::~mniLabelVolume() {
  delete_volume(this->volume);
  delete this->sizes;
}  

void mniLabelVolume::output(VIO_STR file, int cropValue) {
  // should replace the constant with an argument option
    save_label_volume(file, this->filename, this->volume, cropValue);
}

  
  
