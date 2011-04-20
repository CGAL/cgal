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
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<4> CMap_4;
typedef CMap_4::Dart_handle Dart_handle;

Dart_handle make_triangle_2(CMap_4& amap)
{
 Dart_handle d1 = amap.create_dart();
 Dart_handle d2 = amap.create_dart();
 Dart_handle d3 = amap.create_dart();
 amap.link_beta<1>(d1,d2);
 amap.link_beta<1>(d2,d3);
 amap.link_beta<1>(d3,d1);
 return d1;
}
 
Dart_handle make_tetrahedral_3(CMap_4& amap)
{
  Dart_handle d1 = make_triangle_2(amap);
  Dart_handle d2 = make_triangle_2(amap);
  Dart_handle d3 = make_triangle_2(amap);
  Dart_handle d4 = make_triangle_2(amap);
  amap.link_beta<2>(d1, d2);
  amap.link_beta<2>(d3, d2->beta(0));
  amap.link_beta<2>(d1->beta(1), d3->beta(0));
  amap.link_beta<2>(d4, d2->beta(1));
  amap.link_beta<2>(d4->beta(0), d3->beta(1));
  amap.link_beta<2>(d4->beta(1), d1->beta(0));
  return d1;
}

int main()
{
  CMap_4 cm;
  Dart_handle d1 = make_tetrahedral_3(cm);
  Dart_handle d2 = make_tetrahedral_3(cm);
  
  cm.sew<4>(d1,d2);
  
  cm.display_characteristics(std::cout);
  std::cout<<", valid="<<cm.is_valid()<<std::endl;

  return EXIT_SUCCESS;
}

