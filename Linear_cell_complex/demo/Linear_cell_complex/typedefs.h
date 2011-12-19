// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
  template < class Refs >
  struct Dart_wrapper 
  {
    typedef CGAL::Dart<3, Refs > Dart;
    
    typedef CGAL::Cell_attribute_with_point< Refs > Vertex_attrib;
    typedef CGAL::Cell_attribute< Refs, CGAL::Color > Volume_attrib;
    
    typedef CGAL::cpp0x::tuple<Vertex_attrib,void,void,
                               Volume_attrib> Attributes;
  };
};

typedef CGAL::Linear_cell_complex_traits
<3,CGAL::Exact_predicates_inexact_constructions_kernel> Mytraits;

typedef CGAL::Linear_cell_complex<3,3,Mytraits,Myitems> LCC;
typedef LCC::Dart_handle      Dart_handle;
typedef LCC::Vertex_attribute Vertex;

typedef LCC::Point    Point_3;
typedef LCC::Vector   Vector_3;

typedef CGAL::Timer Timer;

struct Scene {
  LCC* lcc;
};

#endif
