/*
 * create_validity_volume.cc:
 *
 * Goal: to check the validity of a white matter cortical extraction
 * by comparing the output of scan_object_to_volume with a manually
 * painted slice.
 *
 * Basic Procedure: the two volumes are loaded. A check is performed
 * at every voxel to test for overlap. Overlap is defined as any part
 * of a continuous series of voxels going in the dominant direction
 * overlapping with the other volume. This check is performed in both
 * cardinal directions of a two-dimensional slice.
 *
 * Output: a slice volume containing the areas where no overlap was
 * found. Possible values are:
 * 1: overlap found
 * 2: no overlap found in the X direction
 * 3: no overlap found in the Z direction
 * 4: no overlap found in either direction
 *
 * Bugs: direction of overlap search is still hardcoded
 *
 * Last modified: Sep 10, 2001
 * Jason Lerch <jason@bic.mni.mcgill.ca>
 */




#include <mniLabelVolume.h>
#include <iostream>
#include <iomanip>

int main ( int argc, char *argv[] ) {
  // load two volumes: painted slice, intersected slice
  static VIO_STR ZYXdimOrder[] = {(char *) MIzspace, (char *) MIyspace, (char *) MIxspace};

  int v0, v1, v2;

  if (argc != 5) {
    cout << "USAGE: " << argv[0] << " -horizontal|-sagital|-coronal painted_slice.mnc intersected_slice.mnc output.mnc" << endl << endl;
    cout << "This programme is desinged to check the validity of an extracted white matter cortex against a manually painted slice. The intersected_slice.mnc argument can be derived by using the scan_object_to_volume app with the cortex and the painted slice as an argument." << endl << endl;
    return 1;
  }
  
  // determine the slice ordering order
  VIO_STR *firstDimOrder = new VIO_STR[3];
  VIO_STR *secondDimOrder = new VIO_STR[3];

  if (strcmp(argv[1], "-sagital") == 0) {
    firstDimOrder = XYZdimOrder;
    secondDimOrder = XZYdimOrder;
  }
  else if (strcmp(argv[1], "-coronal") == 0) {
    firstDimOrder = YXZdimOrder;
    secondDimOrder = YZXdimOrder;
  }
  else if (strcmp(argv[1], "-horizontal") == 0) {
    firstDimOrder = ZXYdimOrder;
    secondDimOrder = ZYXdimOrder;
  }
  else {
    cout << "ERROR: " << argv[1] << " not a valid argument. Please use "
         << "either -sagital, -coronal, or -horizontal." << endl << endl;
    return 1;
  }

  // load the volumes
  mniLabelVolume *paintedSlice = 
    new mniLabelVolume( argv[2], 0.0, 0.0, 3, firstDimOrder );
  mniLabelVolume *intersectSlice = 
    new mniLabelVolume( argv[3], 0.0, 0.0, 3, firstDimOrder );
  mniLabelVolume *output = 
    new mniLabelVolume( argv[3], 0.0, 0.0, 3, firstDimOrder );

  output->setAllVoxels(1);

  for (v0 = 0; v0 < paintedSlice->getSize(0); v0++) {
    for (v1 = 0; v1 < paintedSlice->getSize(1); v1++) {
      for (v2 = 0; v2 < paintedSlice->getSize(2); v2++) {
        int startOfIntersect = v2;
        bool surfaceFound = false;
        bool paintFound = false;
        if (intersectSlice->getVoxel(v0, v1, v2) > 0.5) {
          // start of new surface intersection
          surfaceFound = true;
        }
        else {
          surfaceFound = false;
        }

        while (surfaceFound == true) {
          // keep increasing v2 until surface no longer found
          // or a painted label is found.
          if (paintedSlice->getVoxel(v0,v1,v2) > 0.5) {
            paintFound = true;
            break;
          }
          else if (intersectSlice->getVoxel(v0, v1, v2) == 0) {
            //            surfaceFound = false;
            break;
          }
          else {
            // still on surface but no painted voxel
            v2++;
          }
        }

        if (paintFound == false && surfaceFound == true) {
          for (int i=startOfIntersect; i <= v2; i++) {
            output->setVoxel(2,v0,v1,i);
          }
        }
      }
    }
  }

  // output the output volume - it will be reloaded with a different
  // dimension order
  output->output(argv[4]);

  // reload the volumes with different dimension order
  delete paintedSlice;
  delete intersectSlice;
  delete output;
  paintedSlice = 
    new mniLabelVolume( argv[2], 0.0, 0.0, 3, secondDimOrder );
  intersectSlice = 
    new mniLabelVolume( argv[3], 0.0, 0.0, 3, secondDimOrder );
  output = 
    new mniLabelVolume( argv[4], 0.0, 0.0, 3, secondDimOrder );


  for (v0 = 0; v0 < paintedSlice->getSize(0); v0++) {
    for (v1 = 0; v1 < paintedSlice->getSize(1); v1++) {
      for (v2 = 0; v2 < paintedSlice->getSize(2); v2++) {
        int startOfIntersect = v2;
        bool surfaceFound = false;
        bool paintFound = false;
        if (intersectSlice->getVoxel(v0, v1, v2) > 0.5) {
          // start of new surface intersection
          surfaceFound = true;
        }
        else {
          surfaceFound = false;
        }

        while (surfaceFound == true) {
          // keep increasing v2 until surface no longer found
          // or a painted label is found.
          if (paintedSlice->getVoxel(v0,v1,v2) > 0.5) {
            paintFound = true;
            break;
          }
          else if (intersectSlice->getVoxel(v0, v1, v2) == 0) {
            //            surfaceFound = false;
            break;
          }
          else {
            // still on surface but no painted voxel
            v2++;

          }
        }

        if (paintFound == false && surfaceFound == true) {
          for (int i=startOfIntersect; i <= v2; i++) {
            output->setVoxel(output->getVoxel(v0, v1, i)+2, v0,v1,i);
          }
        }
      }
    }
  }

  output->output(argv[4]);  
  return 0;
}
