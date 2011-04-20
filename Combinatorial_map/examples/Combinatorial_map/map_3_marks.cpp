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
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<3> CMap_3;
typedef CMap_3::Dart_handle Dart_handle;

int main()
{
  CMap_3 cm;
  
  // Reserve a mark
  int mark = cm.get_new_mark();
  if ( mark==-1 )
    {
      std::cerr<<"No more free mark, exit."<<std::endl;
      exit(-1);
    }
  
  // Create two tetrahedra.
  Dart_handle dh1 = make_combinatorial_tetrahedron(cm);  
  Dart_handle dh2 = make_combinatorial_tetrahedron(cm);

  // 3-sew them.
  cm.sew<3>(dh1, dh2);
  
  // Mark the darts belonging to the first tetrahedron.
  for  (CMap_3::Dart_of_cell_range<3>::iterator 
          it(cm.darts_of_cell<3>(dh1).begin()),
          itend(cm.darts_of_cell<3>(dh1).end()); it!=itend; ++it)
    cm.mark(it, mark);

  // Remove the common 2-cell between the two cubes:
  // the two tetrahedra are merged.
  CGAL::remove_cell<CMap_3, 2>(cm, dh1);

  // Thanks to the mark, we know which darts come from the first tetrahedron.
  unsigned int res=0;
  for (CMap_3::Dart_range::iterator it(cm.darts().begin()),
	 itend(cm.darts().end()); it!=itend; ++it)
    {
      if ( cm.is_marked(it, mark) )
	++res;
    }
  
  std::cout<<"Number of darts from the first tetrahedron: "<<res<<std::endl;
  cm.free_mark(mark);

  return EXIT_SUCCESS;
}

