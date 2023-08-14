#include <iostream>
#include <arguments.h>

extern "C" {
#include <bicpl.h>
}

using namespace std;

int main (int argc, char *argv[]) {
  Arguments cArg("create_lines", "(c) jharlap@bic", "-");
  cArg.addOption("help", "display usage help");
  cArg.addArgument("white_surface_mesh_file", "white matter surface mesh");
  cArg.addArgument("gray_surface_mesh_file", "gray matter surface mesh");
  cArg.addArgument("output_object_file", "output object");

  if(!cArg.parse(argc, argv))
    return 1;

  if(cArg.getOption((char *)"help")) {
    cArg.usage();
    return 0;
  }

  // initialize variables for surfaces
  VIO_File_formats    format;
  int             num_objects;
  object_struct** object_list_white;
  object_struct** object_list_gray;
  polygons_struct *polygons_white;
  polygons_struct *polygons_gray;

  object_struct* object_list_lines;
  object_list_lines = create_object(LINES);
  lines_struct *lines;
  lines = get_lines_ptr(object_list_lines);
  initialize_lines(lines, 1);

  // read in the white surface
  if ( input_graphics_file( (char*) cArg[(char *) "white_surface_mesh_file"].c_str(), &format, &num_objects, &object_list_white )
       != VIO_OK ) {
    cerr << "ERROR reading file " << cArg[(char *) "white_surface_mesh_file"] << endl;
    return 1;
  }

  // can't deal with anything other than polygons for now
  if ((object_list_white[0])->object_type != POLYGONS ) {
    cerr << "ERROR: can only read obj files containing polygons." << endl;
    return 1;
  }
  
  // read in the gray surface
  if ( input_graphics_file( (char*) cArg[ (char *) "gray_surface_mesh_file"].c_str(), &format, &num_objects, &object_list_gray )
       != VIO_OK ) {
    cerr << "ERROR reading file " << cArg[(char *)"gray_surface_mesh_file"] << endl;
    return 1;
  }

  // can't deal with anything other than polygons for now
  if ((object_list_gray[0])->object_type != POLYGONS ) {
    cerr << "ERROR: can only read obj files containing polygons." << endl;
    return 1;
  }
  
  // initialize points
  VIO_Point *points_white;
  VIO_Point *points_gray;

  polygons_white = get_polygons_ptr(object_list_white[0]);
  polygons_gray = get_polygons_ptr(object_list_gray[0]);

  points_white = polygons_white->points;
  points_gray = polygons_gray->points;

  for(int pidx=0; pidx<polygons_white->n_points; ++pidx) {
    // add a line to the lines object
    start_new_line(lines);
    add_point_to_line(lines, &points_white[pidx]);
    add_point_to_line(lines, &points_gray[pidx]);
  }

  if( output_graphics_file( (char*) cArg[(char *) "output_object_file"].c_str(), format, 1, &object_list_lines) != VIO_OK) {
    cerr << "ERROR writing file " << cArg[(char *) "output_object_file"] << endl;
    return 1;
  }
  cout << "Done" << endl;
  return 0;
}
