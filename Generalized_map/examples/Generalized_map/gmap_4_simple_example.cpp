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
#include<CGAL/Generalized_map_constructors.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Generalized_map<4> CMap_4;
typedef CMap_4::Dart_handle Dart_handle;

int main()
{
  CMap_4 cm;
  Dart_handle d1 = CGAL::make_combinatorial_tetrahedron(cm);
  Dart_handle d2 = CGAL::make_combinatorial_tetrahedron(cm);

  CGAL_assertion(cm.is_valid());

  cm.sew<4>(d1,d2);

  cm.display_characteristics(std::cout);
  std::cout<<", valid="<<cm.is_valid()<<std::endl;

  return EXIT_SUCCESS;
}
