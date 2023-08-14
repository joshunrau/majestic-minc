/* A utility which prints out a texture map of a surface containing
   values of either 1 or 0 per vertex depending on whether that
   vertex' X coordinate is greater than or less than 0 */

extern "C" {
#include <bicpl.h>
}

#include <iostream>

using namespace std;

int main( int argc, char **argv ) {
  
  VIO_File_formats    format;
  int             num_objects;
  object_struct** object_list;
  char*           filename;
  polygons_struct p;

  if (argc == 2) {
    filename = argv[1];
  }
  else {
    cerr << "USAGE: " << argv[0] << " filename.obj" << endl << endl
	 << "Prints out a value of either 1 or 0 per vertex depending " << endl
	 << "on whether that vertex is across X=0 or not." << endl << endl;
    return 0;
  }

  // open the file
  if ( input_graphics_file( filename, &format, &num_objects, &object_list )
       != VIO_OK ) {
    cerr << "ERROR reading file " << filename << endl;
    return 0;
  }

  // can't deal with anything other than polygons for now
  if ((object_list[0])->object_type != POLYGONS ) {
    cerr << "ERROR: can only read obj files containing polygons." << endl;
    return 0;
  }
  
  p = (object_list[0])->specific.polygons;
  
  // check each polygons coordinate against X=0
  for ( int i = 0; i < p.n_points; ++i ) {
    if (p.points[i].coords[0] < 0)
      cout << 0 << endl;
    else
      cout << 1 << endl;
  }
  return 1;
}

  
    
  
