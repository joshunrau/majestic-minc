/* computes the surface area of a region of interest on a 
polygonal surface, as defined by a vertstat file.

$Id: surface_area_roi.cc,v 1.3 2013-02-07 19:52:48 claude Exp $
*/

#include <mniVertstatsFile.h>

extern "C" {
#include <bicpl.h>
}

#include <iostream>
#include <list>
#include <arguments.h>
#include <algorithm>

using namespace std;

int main (int argc, char *argv[]) {
  vertexColumn  regions;
  list<int>     included_vertices;
  int           included_polygons = 0;
  int           roi;

  VIO_File_formats    format;
  int             num_objects;
  object_struct** object_list;
  polygons_struct *polygons;

  VIO_Point points[10];
  int n_points;
  int poly;

  VIO_Real area;
  area = 0.0;

  // parse command line arguments
  Arguments cArg( "surface_area_roi", "(c) jharlap@bic", "-" );
  cArg.addOption("help", "display usage help");
  Arguments::Option cOptRegion("region", "define the region of interest (defaults to all regions > 0)");
  cOptRegion.addArgument("region", "region of interest");
  cArg.addOption(cOptRegion);

  Arguments::Option cOptColumn("column", "define the column in the vertstats file containing the ROI segmentation");
  cOptColumn.addArgument("column", "name of the column");
  cArg.addOption(cOptColumn);

  Arguments::Option cOptInclusion("include", "inclusion rule - either one or all vertices of a polygon must be in the ROI to include the polygon in the surface area");
  cOptInclusion.addArgument("include", "'one' or 'all' (default is 'one')");
  cArg.addOption(cOptInclusion);

  cArg.addArgument("surface_file", "surface object file");
  cArg.addArgument("vertstats_file", "vertstats file containing ROI segmentation");

  if(!cArg.parse(argc, argv)) 
    return 1;

  if(cArg.getOption( "help")) {
    cArg.usage();
    return 0;
  }

  if(cArg.getOption( "region")) {
    cout << "Using ROI " << cArg.getOption( "region")["region"] << endl;
  }

  if(cArg.getOption( "column")) {
    cout << "Using column " << cArg.getOption( "column")[ "column"] << endl;
  }

  // open the surface file
  cout << "Loading: " << cArg[ "surface_file"] << endl;
  if ( input_graphics_file( (char*) cArg[ "surface_file"].c_str(), &format, &num_objects, &object_list )
       != VIO_OK ) {
    cerr << "ERROR reading file " << cArg[ "surface_file"] << endl;
    return 0;
  }

  // can't deal with anything other than polygons for now
  if ((object_list[0])->object_type != POLYGONS ) {
    cerr << "ERROR: can only read obj files containing polygons." << endl;
    return 0;
  }
  
  polygons = get_polygons_ptr(object_list[0]);
  
  // open the verstat file
  mniVertstatsFile stats(cArg[ "vertstats_file"].c_str());
  cout << "Loading: " << cArg[ "vertstats_file"] << endl;

  // read the user-defined column (or the first column if not defined)
  if(cArg.getOption( "column"))
    regions = stats.getDataColumn(cArg.getOption( "column")[ "column"].c_str());
  else
    regions = stats.getDataColumn(0);

  // get the region of interest
  if(cArg.getOption( "region")) {
    roi = atoi(cArg.getOption( "region")[ "region"].c_str());
  } else {
    roi = -1;
  }

  // find vertices that match the desired ROI
  vertexColumn::iterator regions_iterator;
  int vert_index;
  vert_index = 0;
  cout << "Starting loop through region list" << endl;
  for (regions_iterator = regions.begin() ; regions_iterator != regions.end() ; regions_iterator++) {
    // if this vertex is part of the ROI
    if ((roi == -1 && *regions_iterator > 0) || *regions_iterator == roi) {
      // add the vertex to the list of vertices to include in the surface area tally
      included_vertices.push_back(vert_index);
    }
    vert_index++;
  }
  cout << "Found " << included_vertices.size() << " vertices for ROI " << roi << endl;

  list<int>::iterator point_idx;
  int included_vertex_count = 0;
  int required_vertex_count = 1;

  if(cArg.getOption( "include") && cArg.getOption( "include")[ "include"] == "all")
    required_vertex_count = n_points;
    

  // loop through all the polygons
  cout << "Starting loop through polygons" << endl;
  for (poly = 0 ; poly < polygons->n_items ; ++poly) {
    // find the points for the polygon
    n_points = get_polygon_points(polygons, poly, points);
    if(n_points > 0) {
      // reset the included vertex counter and set the number of required vertices
      included_vertex_count = 0;
      if(cArg.getOption( "include") && cArg.getOption( "include")[ "include"] == "all")
        required_vertex_count = n_points;

      // loop over all the points in the polygon
      for (int i = 0; i < n_points; ++i) {
        // include the point in the included vertex count if it's in
        // the included vertex list
        point_idx = std::find(included_vertices.begin(), included_vertices.end(), POINT_INDEX(polygons->end_indices,poly,i));
        if(point_idx != included_vertices.end()) {
          included_vertex_count++;
        }

        // if we found enough vertices to include the poly
        if(included_vertex_count >= required_vertex_count) {
          // add the polygon's surface area to the tally
          area += get_polygon_surface_area (n_points, points);
          included_polygons++;
          
          // stop looping over points if we've included enough points
          break;
        }
      } // end looping over points
    }
  } // end looping through polygons

  cout << "Surface area: " << area << " using " << included_polygons << " polygons" << endl;
  
  return 0;
}
