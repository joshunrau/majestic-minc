#ifndef __MNIVOLUME__
#define __MNIVOLUME__

extern "C" {
#include "bicpl.h"
#include "volume_io.h"
}

#include <iostream>
#include <math.h>
#include "mniBaseVolume.h"

using namespace std;

//! A class for working with minc volumes
/*!

  This class uses the standard volume_io functions for handling
  volumes.  It should be used for most cases ... at least until
  classes which bypass volume_io come along.

  \todo Create a constructor for creating a new volume without reading
  header information from a file
*/
class mniVolume : public mniBaseVolume {

public:
  //! Empty Constructor
  /*!  
  An empty constructor - this constructor really should not be used by
  any calling code, as all that it does is initialiset the sizes
  variable. Moreover, the code to create a volume from scratch is not
  yet present in this class
  */

  mniVolume();
  //! Constructor from file
  /*!
  A constructor which initialises the class by reading in the volume
  information from the filename provided. The other arguments all have
  default values that can usually be left untouched.

  \exception loadException Thrown should there be an error reading the
  file
  */

  mniVolume(VIO_STR filename, 
            VIO_Real voxelMin = 0.0, 
            VIO_Real voxelMax = 0.0,
            int nDimensions = 3,
            VIO_STR dimensions[] = ZXYdimOrder,
            nc_type dataType = NC_UNSPECIFIED, 
            VIO_BOOL volumeSigned = FALSE,
            VIO_BOOL createVolume = TRUE, 
            minc_input_options *options = NULL
            );

  //! Copy constructor from mniVolume
  /*!
  Copies the volume definition, including the volume data itself, from
  the volume provided.
  
  \param copyVolumeDefinitionOnly If this parameter is set to TRUE,
  only the volume definition is copied and the subsequent arguments
  come into play. If it is false, an exact copy is created, and the
  susequent arguments are ignored

  */
  mniVolume(mniBaseVolume *copyVolume, 
	    VIO_BOOL copyVolumeDefinitionOnly=FALSE,
	    nc_type dataType = NC_UNSPECIFIED,
	    VIO_BOOL signedFlag = FALSE,
	    VIO_Real voxelMin = 0.0,
	    VIO_Real voxelMax = 0.0);

  //! Constructor from a volume_io volume struct.
  mniVolume(VIO_Volume volumeIO_volume);
  //! Destructor to free memory
  virtual ~mniVolume();

  //! Get voxel value
  VIO_Real getVoxel(int v1, int v2, int v3, int v4=0, int v5=0) {
    return get_volume_real_value(this->volume, v1, v2, v3, v4, v5);
  };
  //! Overloaded getVoxel, taking three dimensional array
  VIO_Real getVoxel(int indices[3]) {
    return get_volume_real_value(this->volume, indices[0], indices[1],
				 indices[2], 0, 0);
  };
  //! Get a point in world space
  /*!
    First converts the three world indices into voxel space, then
    uses the rint function to round them to integers before passing
    them onto the get_volume_real_value function
  */
  VIO_Real getWorld(VIO_Real xWorld, VIO_Real yWorld, VIO_Real zWorld);
  //! Set voxel value
  void setVoxel(VIO_Real value, int v1, int v2, int v3,
		int v4=0, int v5=0) {
    set_volume_real_value(this->volume, v1, v2, v3, v4, v5, value);
  };
  //! Overloaded setVoxel, taking three dimensional array
  void setVoxel(VIO_Real value, int indices[3]) {
    set_volume_real_value(this->volume, indices[0], indices[1],
			  indices[2], 0, 0, value);
  };
  //  virtual void output() { };
  virtual void output(VIO_STR file, int cropValue = 0);

};

#endif // __MNIVOLUME__
