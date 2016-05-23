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

typedef CGAL::Generalized_map<4> GMap_4;
typedef GMap_4::Dart_handle Dart_handle;

int main()
{
  GMap_4 gm;
  Dart_handle d1 = gm.make_combinatorial_tetrahedron();
  Dart_handle d2 = gm.make_combinatorial_tetrahedron();

  CGAL_assertion(gm.is_valid());

  gm.sew<4>(d1,d2);

  gm.display_characteristics(std::cout);
  std::cout<<", valid="<<gm.is_valid()<<std::endl;

  return EXIT_SUCCESS;
}
