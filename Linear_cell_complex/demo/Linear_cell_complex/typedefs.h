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
// $URL: svn+ssh://gdamiand@scm.gforge.inria.fr/svn/cgal/branches/features/Linear_cell_complex-gdamiand/Linear_cell_complex/demo/Linear_cell_complex/typedefs.h $
// $Id: typedefs.h 65446 2011-09-20 16:55:42Z gdamiand $
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/Color.h>
#include <CGAL/Timer.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#define COLOR_VOLUME 1 // Pour activer la couleur des volumes

#ifdef COLOR_VOLUME
template<class Cell>
struct Average_functor : public std::binary_function<Cell,Cell,void>
  {
    void operator()(Cell& acell1,Cell& acell2)
    { 
      acell1.attribute()=
	CGAL::Color((acell1.attribute().r()+acell2.attribute().r())/2,
		    (acell1.attribute().g()+acell2.attribute().g())/2,
		    (acell1.attribute().b()+acell2.attribute().b())/2);
    }
  };
class Myitems
{
public:
  // typedef CGAL::Exact_predicates_inexact_constructions_kernel Traits;  
  
  template < class Refs >
  struct Dart_wrapper 
  {
    typedef CGAL::Dart<3, Refs > Dart;
    
    typedef CGAL::Cell_attribute_with_point< Refs > Vertex_attrib;
    typedef CGAL::Cell_attribute< Refs, CGAL::Color > Volume_attrib;
    
    typedef CGAL::cpp0x::tuple<Vertex_attrib,CGAL::Disabled,CGAL::Disabled,Volume_attrib>
    Attributes;
  };
};
#else // COLOR_VOLUME
typedef CGAL::Combinatorial_map_with_points_min_items<3,3> Myitems;
#endif

typedef CGAL::Linear_cell_complex_traits<3,CGAL::Exact_predicates_inexact_constructions_kernel> Mytraits;

typedef CGAL::Combinatorial_map_with_points<3,3,Mytraits,Myitems> Map;
typedef Map::Dart_handle      Dart_handle;
typedef Map::Vertex_attribute Vertex;

typedef Map::Point    Point_3;
typedef Map::Vector   Vector_3;
typedef Map::Traits::Iso_cuboid_3 Iso_cuboid_3;

typedef CGAL::Timer Timer;

struct Scene {
  Map* map;
};





#endif
