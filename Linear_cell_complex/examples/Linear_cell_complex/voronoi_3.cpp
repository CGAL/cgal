// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <iostream>
#include <fstream>
// #include "cgal_map_viewer_qt.h"
// #include "cgal_map_viewer_vtk.h"


typedef CGAL::Linear_cell_complex<3, 3> Map_3;
typedef Map_3::Dart_handle              Dart_handle;
typedef Map_3::Point                    Point;

typedef CGAL::Delaunay_triangulation_3<Map_3::Traits> Triangulation;

int main(int argc, char** argv)
{
   if (argc!=2)
   {
      std::cout<<"Usage : voronoi_3 filename"<<std::endl   
               <<"   filename being a fine containing 3D points used to "
               <<" compute the Delaunay_triangulation_3."<<std::endl;
      exit(0);         
   }   

  // 1) Compute the Delaunay_triangulation_3 
  Triangulation T;

  std::ifstream iFile(argv[1], std::ios::in);
  if (!iFile)
  {
     std::cout << "Problem reading file " << argv[1] << std::endl;
     exit(0);
  }
  else
  {
     std::istream_iterator<Point> begin(iFile), end;
     T.insert(begin, end);
  }
  assert(T.is_valid(false));
 
  // 2) Convert the triangulation into a 3D combinatorial map
  Map_3 map;
  CGAL::import_from_triangulation_3<Map_3, Triangulation>(map, T);

  std::cout<<"Delaunay triangulation :"<<std::endl<<"  ";
  map.display_characteristics(std::cout) << ", valid=" 
					 << map.is_valid() << std::endl;

  //  display_map(map); // Possible if you use cgal_map_viewer_qt or _vtk

  // 3) Compute the dual map.
  Map_3 dual_map;
  CGAL::dual<Map_3>(map,dual_map);
  
  // 4) Display the map characteristics.
  std::cout<<"Voronoi subdvision :"<<std::endl<<"  ";
  dual_map.display_characteristics(std::cout) << ", valid=" 
					      << dual_map.is_valid()
					      << std::endl;

  return 1;
}

