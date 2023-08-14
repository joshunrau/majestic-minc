#ifndef __MNIBASEVOLUME__
#define __MNIBASEVOLUME__

extern "C" {
#include "bicpl.h"
#include "volume_io.h"
}

#include <iostream>

using namespace std;

// static value used as default for volume loading
static VIO_STR  ZXYdimOrder[] = {(char *) MIzspace, (char *) MIxspace, (char *) MIyspace};
static VIO_STR  ZYXdimOrder[] = {(char *) MIzspace, (char *) MIyspace, (char *) MIxspace};

static VIO_STR  XYZdimOrder[] = {(char *) MIxspace, (char *) MIyspace, (char *) MIzspace};
static VIO_STR  XZYdimOrder[] = {(char *) MIxspace, (char *) MIzspace, (char *) MIyspace};

static VIO_STR  YXZdimOrder[] = {(char *) MIyspace, (char *) MIxspace, (char *) MIzspace};
static VIO_STR  YZXdimOrder[] = {(char *) MIyspace, (char *) MIzspace, (char *) MIxspace};

//! An abstract baseclass for a minc volume
/*!

This is an abstract class, meaning that no members of this type can
ever be created. One of the derived classes, mniVolume or
mniLabelVolume has to be used instead.

\todo Incorporate the voxel to world transformation info

*/

class mniBaseVolume {
protected:
  //! Holds the volume_io volume
  VIO_Volume       volume;
  //! Holds the sizes - has to be instantiated before use
  int          *sizes; 
  //! Holds the number of dimensions
  int          nDimensions;
  //! Holds the dimension order
  VIO_STR*      dimNames;
  //! Holds the original filename - not yet used
  VIO_STR       filename;
  //! Holds the minc data type
  nc_type      dataType;
  //! Holds the minimum voxel value
  VIO_Real         voxelMin;
  //! Holds the maximum voxel value
  VIO_Real         voxelMax;
  //! Whether the data type is signed or not
  VIO_BOOL      signedFlag;

public:
  //! Load exception class
  class loadException { };
  //! Write exception class
  class writeException { };
  //! Set the filename
  void setFilename(VIO_STR file) { filename = file; };
  //! Return pointer to volume_io volume
  VIO_Volume getVolume() { return this->volume; };
  //! Get pointer to volume sizes
  int* getSizes() { return this->sizes; };
  //! Get one size from sizes array
  int getSize(int index) { return this->sizes[index]; };
  //! Get dimensions names
  VIO_STR *getDimNames() { return this->dimNames; };
  //! Return volume min
  VIO_Real getVoxelMin() { return this->voxelMin; };
  //! Retrun volume max
  VIO_Real getVoxelMax() { return this->voxelMax; };
  //! Set the volume real range
  void setRealRange(VIO_Real lower, VIO_Real upper) { set_volume_real_range(
                                          this->volume, lower, upper); }
  //! Return signed flag
  VIO_BOOL getSignedFlag() { return this->signedFlag; };
  //! Return data type
  nc_type getDataType() { return this->dataType; };

  // voxel to world coordinate stuff:

  //!converts a voxel to world space
  /*!
    Returns a voxel to from voxel space to world space
    \param voxel[] An array holding the voxel to be converted
    \return An array holding the world coordinates in X Y Z order
    \note You have to free the memory of the returned array yourself
  */
  VIO_Real* convertVoxelToWorld(VIO_Real voxel[]);
  //! Convert a world coordinate into a voxel
  /*!
    \return An array holding the voxel coordinates
  */
  VIO_Real* convertWorldToVoxel(VIO_Real xWorld, VIO_Real yWorld, VIO_Real zWorld);
  //! Gets interpolated value at indices
  /*!
    \bug Use with caution - the returned VIO_Real argument ought to be an array,
    since in some situations the underlying volume_io function is supposed
    to return more than one value. But I don't quite (yet) understand when
    and how this is supposed to happen.
  */
  VIO_Real getInterpolatedVoxel(VIO_Real indices[],
			    int degreesContinuity=2,
			    VIO_BOOL interpolatingDimensions[]=NULL,
			    int useLinearAtEdge=TRUE,
			    VIO_Real outsideValue=0,
			    VIO_Real **firstDerivative=NULL,
			    VIO_Real ***secondDerivative=NULL) {
    VIO_Real tmpReturnValue;
    evaluate_volume(this->volume,
		    indices,
		    interpolatingDimensions,
		    degreesContinuity,
		    useLinearAtEdge,
		    outsideValue,
		    &tmpReturnValue,
		    firstDerivative,
		    secondDerivative);
    return tmpReturnValue;
  };
  //! Overloaded version of getInterpolatedVoxel
  VIO_Real getInterpolatedVoxel(VIO_Real v1, VIO_Real v2, VIO_Real v3,
			    int degreesContinuity=2,
			    VIO_BOOL interpolatingDimensions[]=NULL,
			    int useLinearAtEdge=TRUE,
			    VIO_Real outsideValue=0,
			    VIO_Real **firstDerivative=NULL,
			    VIO_Real ***secondDerivative=NULL) {
    VIO_Real tmpReturnValue;
    VIO_Real indices[3] = {v1, v2, v3};
    evaluate_volume(this->volume,
		    indices,
		    interpolatingDimensions,
		    degreesContinuity,
		    useLinearAtEdge,
		    outsideValue,
		    &tmpReturnValue,
		    firstDerivative,
		    secondDerivative);
    return tmpReturnValue;
  };
  

  //! Output the volume
  virtual void output(VIO_STR file, int cropValue = 0) = 0;
};

#endif
