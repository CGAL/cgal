// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_INTERSECTION_LINE_2_CIRCLE_2_MAP_H
#define CGAL_INTERSECTION_LINE_2_CIRCLE_2_MAP_H

#include <CGAL/license/Circular_kernel_2.h>


#include <map>
#include <vector>

namespace CGAL {
namespace internal {

class Intersection_line_2_circle_2_map {

typedef struct inter_map_pair {
  int x, y;
  inter_map_pair(int xx=0, int yy=0) : x(xx), y(yy) {}
  inter_map_pair(const inter_map_pair &i) : x(i.x), y(i.y) {}
  bool operator<(const inter_map_pair &i) const {
    if(x < i.x) return true;
    if(x > i.x) return false; 
    if(y < i.y) return true;
    return false;
  }
} inter_map_pair;
typedef std::map< inter_map_pair , CGAL::Object > Table;

private:
  Table intersection_map;
  unsigned int id_gen;

public:
  Intersection_line_2_circle_2_map() : id_gen(0) { intersection_map.clear(); }
  ~Intersection_line_2_circle_2_map() { intersection_map.clear(); }
  
  unsigned int get_new_id() {
    return ++id_gen;
  } 

  template < class T >
  bool find(int id1, int id2, T& res) const {
    Table::const_iterator p = intersection_map.find(
      inter_map_pair(id1,id2));
    if(p == intersection_map.end()) return false;
    assign(res, p->second);
    return true;
  }

  template < class T >
  void put(const int id1, const int id2, const T& res) {
    intersection_map[inter_map_pair(id1,id2)] = CGAL::make_object(res);
  }
};

} // endof internal cgal namespace
} //endof cgal namespace

#endif // CGAL_INTERSECTION_LINE_2_CIRCLE_2_MAP_H
