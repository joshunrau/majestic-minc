// takes a coloured obj file and turns the colours into a texture map.
// for the moment works only with per_item colours (which is what Display produces by
// default when painting surfaces) in a really awkward sort of way which does somehow 
// appear to work.

extern "C" {
#include <bicpl.h>
#include <ParseArgv.h>
}

#include <iostream>
#include <vector>
#include <string>
#include <mniVertstatsFile.h>

// variables for argument parsing
float ignoredValue = 6;

ArgvInfo argTable[] = 
  { {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, "VIO_Colour handling" },
    {"-ignore_value", ARGV_FLOAT, (char *) 1, (char *) &ignoredValue,
     "Set this value in the coloured obj file to 0 in the output" },
    {NULL, ARGV_END, NULL, NULL, NULL}
  };

using namespace std;

int main (int argc, char *argv[]) {
  VIO_File_formats format;
  int num_objects;
  object_struct **object_list;
  polygons_struct p;
  string outfile;
  mniVertstatsFile out;

  out.addToHistory(argc, argv);

  if ( ParseArgv( &argc, argv, argTable, 0 ) || argc != 3) {
    cerr << "USAGE: " << argv[0] << " input.obj " << " output.vertstats" << endl;
    return 1;
  }

  // load the file
  if ( input_graphics_file( argv[1] , &format, &num_objects, &object_list) != VIO_OK ) {
    return 1;
  }

  outfile = argv[2];

  if (num_objects > 1) {
    cerr << "WARNING: only dealing with first object in obj file" << endl;
  }

  if ( (object_list[0])->object_type != POLYGONS ) {
    cerr << "ERROR: can only deal with POLYGONS!" << endl;
    return 1;
  }

  p = (object_list[0])->specific.polygons;

  // make sure that it actually is a coloured obj file
  
  if (p.colour_flag == PER_VERTEX_COLOURS) {
    cerr << "PER_VERTEX_COLOURS not yet supported" << endl;
    return 1;
  }
  else if (p.colour_flag == ONE_COLOUR) {
    cerr << "obj file has no colours in it" << endl;
    return 1;
  }

  std::vector<float> v(p.n_points);

  // now actually get to the colours - there has to be a smarter way of doing this, but
  // it appears to work.
  for ( int i = 0; i < p.n_items; ++i ) {
    int f_size = GET_OBJECT_SIZE( p, i );
    for (int j = 0; j < f_size; ++j) {
      float r = get_Colour_r_0_1(p.colours[i]) * 1;
      float g = get_Colour_g_0_1(p.colours[i]) * 2;
      float b = get_Colour_b_0_1(p.colours[i]) * 3;
      
      // add the three colour components together for the final texture.
      float value = r + g + b;
      if (value == ignoredValue) {
        v[p.indices[POINT_INDEX(p.end_indices, i, j)]] = 0;
      }
      else {
        v[p.indices[POINT_INDEX(p.end_indices, i, j)]] = value;
      }
    }
  }
  
  out.putDataColumn(v, "labels");
  out.writeFile(outfile);
  // print the per-vertex colours to stdout
  //  for (int i = 0; i < p.n_points; ++i) {
  //    cout << v[i] << endl;
  // }
  return 0;
}
