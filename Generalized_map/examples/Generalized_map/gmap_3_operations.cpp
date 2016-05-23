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
typedef GMap_3::Dart_handle      Dart_handle;

int main()
{
  GMap_3 gm;

  // Create one hexahedron.
  Dart_handle d1 = gm.make_combinatorial_hexahedron();
  CGAL_assertion( gm.is_valid() );
  CGAL_assertion( gm.is_volume_combinatorial_hexahedron(d1) );

  // Add two edges along two opposite facets.
  //  CGAL::insert_cell_1_in_cell_2(cm,d1->alpha(1),d1->alpha(0));
  CGAL_assertion( gm.is_valid() );

  Dart_handle d2=d1->alpha(2)->alpha(1)->alpha(1)->alpha(2);
  //  CGAL::insert_cell_1_in_cell_2(cm,d2,d2->alpha(1)->alpha(1));
  CGAL_assertion( gm.is_valid() );

  // Insert a facet along these two new edges plus two initial edges of the cube.
  std::vector<Dart_handle> path;
  path.push_back(d1->alpha(1));
  path.push_back(d1->alpha(0)->alpha(2)->alpha(1));
  path.push_back(d2->alpha(0));
  path.push_back(d2->alpha(2)->alpha(1));

  //Dart_handle d3=CGAL::insert_cell_2_in_cell_3(cm,path.begin(),path.end());
  CGAL_assertion( gm.is_valid() );

  // Display the m characteristics.
  gm.display_characteristics(std::cout) << ", valid=" <<
    gm.is_valid() << std::endl;

  // We use the removal operations to get back to the initial cube.
  //  gm.remove_cell<2>(d3);
  CGAL_assertion( gm.is_valid() );

  //  gm.remove_cell<1>(d1->alpha(1));
  CGAL_assertion( gm.is_valid() );

  //  gm.remove_cell<1>(d2->alpha(0));
  CGAL_assertion( gm.is_valid() );

  // Display the m characteristics.
  gm.display_characteristics(std::cout) << ", valid="
					<< gm.is_valid() << std::endl;

  return EXIT_SUCCESS;
}
