//! \file examples/Arrangement_2/ex_postscript.cpp
// Using the arrangement Postscript output operator.

#include <CGAL/basic.h>

#ifndef CGAL_USE_LEDA
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs LEDA ..." << std::endl; 
  return 0;
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_postscript_file_stream.h>

#include "arr_rational_nt.h"
#include "point_location_utils.h"

typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;

int main ()
{
  // Construct the arrangement:
  Arrangement_2 arr;
  construct_segments_arr(arr);

  // Extract dimensions:
  double x_min = DBL_MAX, y_min = DBL_MAX, x_max = DBL_MIN, y_max = DBL_MIN;
  Arrangement_2::Vertex_const_iterator vi;
  for (vi = arr.vertices_begin(); vi != arr.vertices_end(); ++vi) {
    double x = CGAL::to_double(vi->point().x());
    double y = CGAL::to_double(vi->point().y());
    if (x < x_min) x_min = x;
    if (x_max < x) x_max = x;
    if (y < y_min) y_min = y;
    if (y_max < y) y_max = y;
  }
  // Create margin:
  float width = x_max - x_min;
  float height = y_max - y_min;
  float margin = ((width < height) ? height : width) * 0.1f;

  // Write the arrangement in postscript format to a file:
  CGAL::Arr_postscript_file_stream ps_file;
  ps_file.init(x_min - margin, x_max + margin, y_min - margin);
  ps_file << arr;
  ps_file.close();
  return 0;
}

#endif
