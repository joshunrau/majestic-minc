#ifndef __MNILABELVOLUME__
#define __MNILABELVOLUME__

extern "C" {
#include "bicpl.h"
#include "volume_io.h"
}

#include <iostream>
#include "mniBaseVolume.h"
#include "mniVolume.h"

using namespace std;

//! A class to deal with Label Volumes
/*!
  The mniLabelVolume class differs from the mniVolume class insofar as
  it uses the bicpl label_volume set of functions to deal with the
  volume as opposed to the volume_io functions. The general unit is also
  an int as opposed to a VIO_Real.
*/
class mniLabelVolume : public mniBaseVolume {
  
public:
  //! Empty Constructor
  /*!  
    An empty constructor which does nothing except initialise the
    sizes array. It does not do any volume handling, and therefore should
    only be used if the Volume variable will be assigned to it. And since
    it is a protected variable, this constructor should only be used in a
    subclass.
  */
  mniLabelVolume();
  //! Constructor from file, creating initialised volume
  mniLabelVolume(VIO_STR filename, 
		 VIO_Real voxelMin = 0.0, 
		 VIO_Real voxelMax = 0.0,
		 int nDimensions = 3,
		 VIO_STR dimensions[] = ZXYdimOrder,
		 nc_type dataType = NC_UNSPECIFIED, 
		 VIO_BOOL volumeSigned = FALSE,
		 VIO_BOOL createVolume = TRUE, 
		 minc_input_options *options = NULL
		 );
  //! Constructor from file, creating uninitialised volume
  /*!
    Creates a label volume with the same parameters as the filename
    passed in.
    \param newVolume Itself of no use, only here to differentiate it from
    the other constructors for function overloading.
  */
  mniLabelVolume(VIO_STR filename,
		 int newVolume,
		 int nDimensions = 3,
		 VIO_STR dimensions[] = ZXYdimOrder,
		 nc_type dataType = NC_UNSPECIFIED,
		 minc_input_options *options = NULL);
  //! Copy constructor from another volume
  mniLabelVolume(mniBaseVolume *copyVolume, nc_type dataType = NC_SHORT);
  //! Constructor from a volume_io volume struct.
  mniLabelVolume(VIO_Volume volumeIO_volume);
  //! Destructor to clean up memory
  virtual ~mniLabelVolume();
  //! Sets all voxels in the volume to one value
  /*!
    \bug Dubious behaviour that I have not quite figured out yet when the
    value is set to 0
  */
  void setAllVoxels(int value) { 
    set_all_volume_label_data(this->volume, value); };

  int getVoxel(int v1, int v2, int v3, int v4=0, int v5=0) {
    return get_volume_label_data_5d(this->volume, v1, v2, v3, v4, v5);
  };

  int getVoxel(int indices[3]) {
    return get_volume_label_data(this->volume, indices);
  };

  void setVoxel(int value, int v1, int v2, int v3, 
                int v4=0, int v5=0) {
    set_volume_label_data_5d(this->volume, v1, v2, v3, v4, v5, value);
  };

  void setVoxel(int value, int indices[3]) {
    set_volume_label_data(this->volume, indices, value);
  }

  virtual void output(VIO_STR file, int cropValue = 0);
};

#endif // __MNILABELVOLUME__

