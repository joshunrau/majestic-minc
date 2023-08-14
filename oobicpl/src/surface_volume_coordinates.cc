// takes an obj file and returns the vertices in volume (rather than
// world coordinates). Not particularly useful, but written to satisfy
// the folks in Minnesota.

// BIC library includes
extern "C" {
#include <bicpl.h>
#include <volume_io.h>
}
#include <mniVolume.h>

// STL includes
#include <iostream>
#include <string>

// argument parsing
#include <arguments.h>

using namespace std;

int main (int argc, char *argv[]) {
  // geometry definitions
  VIO_File_formats        format;
  int                 num_objects;
  object_struct**     object_list;
  polygons_struct*    polygons;
  object_struct*      object;
  VIO_Point*              points;
  int                 n_points;
  
  // filename definitions
  string              input_object;
  string              input_volume;
  string              output_object;
  
  // parse command line arguments
  Arguments cArg( "surface_volume_coordinates", "", "-");
  cArg.addOption("help", "display usage help");
  cArg.addArgument("input_object", "the input surface");
  cArg.addArgument("input_volume", "the volume to use for coordinate details");
  cArg.addArgument("output_object", "the surface with volume coordinates");

  // parse arguments
  if (!cArg.parse(argc, argv)) {
    return 1;
  }
  // print help message if -help
  if (cArg.getOption((char *) "help")) {
    cArg.usage();
    return 0;
  }

  // open the input object
  if (input_graphics_file( (char*) cArg[(char *) "input_object"].c_str(),
			   &format, &num_objects, &object_list) != VIO_OK) {
    cerr << "ERROR reading file " << cArg[(char *) "input_object"] << endl;
    return 1;
  }

  // complain if it is anything other than a polygon set
  if ( (object_list[0])->object_type != POLYGONS ) {
    cerr << "ERROR: can only read obj files containing polygons." << endl;
    return 1;
  }

  // open the volume which will be used to compute the world to voxel
  // coordinates.
  mniVolume *volume = new mniVolume( (char *) cArg[(char *) "input_volume"].c_str(),
				     0.0, 0.0, 3, XYZdimOrder);

  n_points = get_object_points(object_list[0], &points);
  VIO_Real x,y,z;
  for( int i=0; i < n_points; i++ ) {
    convert_3D_world_to_voxel(volume->getVolume(), 
			      (VIO_Real) Point_x(points[i]),
			      (VIO_Real) Point_y(points[i]),
			      (VIO_Real) Point_z(points[i]),
			      &x, &y, &z);
    fill_Point( points[i], x, y, z );
  }
  
  compute_polygon_normals( get_polygons_ptr(object_list[0]) );

  // output the revised obj file
  (void) output_graphics_file( (char *)cArg[(char *) "output_object"].c_str(), format,
			       num_objects, object_list );
  
  return 0;
}


