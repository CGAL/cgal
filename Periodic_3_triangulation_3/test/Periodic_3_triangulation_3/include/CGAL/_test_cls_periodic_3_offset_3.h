// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
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
// 
//
// Author(s)     : Manuel Caroli
//                 Francois Rebufat

#include <sstream>
#include <cassert>

template <class Offset>
void _test_cls_periodic_3_offset_3(const Offset &)
{
  // Creation
  Offset o;

  Offset o0(0,0,0);
  Offset o1(1,2,3);
  Offset o2(1,2,3);
  Offset o3(-3,2,-1);
  Offset o5(1,-2,3);

  assert(o0.is_null());
  assert(!o1.is_null());
  assert(!o2.is_null());
  assert(!o3.is_null());
  assert(o1 == o2);
  assert(o2 != o3);
  
  // Operations
  o = o1+o3;
  assert(o == Offset(-2,4,2));
  
  o = o1-o2;
  assert(o.is_null());
  o = o1-o3;
  assert(o == Offset(4,0,4));
  
  assert(-o == Offset(-4,0,-4));

  o += o2;
  assert(o == Offset(5,2,7));

  o -= o2;
  assert(o == Offset(4,0,4));

  assert(o1 == o2);

  assert(o0 != o1);

  assert(o3 < o0);
  assert(o0 < o5);
  assert(o5 < o1);

  // Access
  assert(o3[0] == -3);
  assert(o3[1] ==  2);
  assert(o3[2] == -1);

  assert(o3.x() == -3);
  assert(o3.y() ==  2);
  assert(o3.z() == -1);

  assert(o0.is_null());
  assert(!o1.is_null());

  std::stringstream ss;
  ss << o1;
  assert(ss.str() == "1 2 3");

  ss >> o;
  assert(o == o1);
}
