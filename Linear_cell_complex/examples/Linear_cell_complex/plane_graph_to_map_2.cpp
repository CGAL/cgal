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
#include <CGAL/Dart.h>
#include <iostream>
#include <fstream>
#include <cstring>

typedef CGAL::Linear_cell_complex<2,2> Map_2;
typedef Map_2::Dart_handle             Dart_handle;
typedef Map_2::Point                   Point;

int main(int narg, char** argv)
{
  if ( narg!=2 || !strcmp(argv[1],"-h") || !strcmp(argv[1],"-?") )
    {
      std::cout<<"Usage: a.out filename"<<std::endl
               <<"  with filename the name of the file containing a 2D "
		   "plane graph."<<std::endl<<std::endl
	       <<"File must be in text mode, respecting the following format:"
	       <<std::endl
	       <<"***********************************************"<<std::endl
	       <<"OFF2D"<<std::endl
	       <<"nbvertices nbedges"<<std::endl
	       <<"x y //coordinates, one line for each vertex"<<std::endl
	       <<"..."<<std::endl
	       <<"i j //edge betwen vertices number i and j,"
	       <<" one line for each edge"<<std::endl
	       <<"..."<<std::endl
	       <<"***********************************************"<<std::endl
	       <<std::endl;
      exit(EXIT_FAILURE);
    }
  
  Map_2 m;

  std::ifstream is(argv[1]);
  std::cout<<"Import plane graph from "<<argv[1]<<std::endl;
  CGAL::import_from_plane_graph(m, is);
  
  // Display the m characteristics.
  std::cout<<"Map characteristics:"<<std::endl<<"  ";
  m.display_characteristics(std::cout) 
    << ", valid=" << m.is_valid() << std::endl;
  
  return 1;
}

