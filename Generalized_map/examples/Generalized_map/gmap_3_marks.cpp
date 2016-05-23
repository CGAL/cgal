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
#include <CGAL/Generalized_map.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Generalized_map<3> GMap_3;
typedef GMap_3::Dart_handle Dart_handle;

int main()
{
  GMap_3 gm;

  // Reserve a mark
  int mark = gm.get_new_mark();
  if ( mark==-1 )
    {
      std::cerr<<"No more free mark, exit."<<std::endl;
      exit(-1);
    }

  // Create two tetrahedra.
  Dart_handle dh1 = gm.make_combinatorial_tetrahedron();
  Dart_handle dh2 = gm.make_combinatorial_tetrahedron();

  CGAL_assertion( gm.is_valid() );
  CGAL_assertion( gm.is_volume_combinatorial_tetrahedron(dh1) );
  CGAL_assertion( gm.is_volume_combinatorial_tetrahedron(dh2) );

  // 3-sew them.
  gm.sew<3>(dh1, dh2);

  // Mark the darts belonging to the first tetrahedron.
  for  (GMap_3::Dart_of_cell_range<3>::iterator
          it(gm.darts_of_cell<3>(dh1).begin()),
          itend(gm.darts_of_cell<3>(dh1).end()); it!=itend; ++it)
    gm.mark(it, mark);

  // Remove the common 2-cell between the two cubes:
  // the two tetrahedra are merged.
  //TODO  gm.remove_cell<2>(dh1);

  // Thanks to the mark, we know which darts come from the first tetrahedron.
  unsigned int res=0;
  for (GMap_3::Dart_range::iterator it(gm.darts().begin()),
         itend(gm.darts().end()); it!=itend; ++it)
  {
    if ( gm.is_marked(it, mark) )
      ++res;
  }

  std::cout<<"Number of darts from the first tetrahedron: "<<res<<std::endl;
  gm.free_mark(mark);

  return EXIT_SUCCESS;
}
