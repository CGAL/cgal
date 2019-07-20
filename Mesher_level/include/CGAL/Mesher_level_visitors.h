// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESHER_LEVEL_VISITORS_H
#define CGAL_MESHER_LEVEL_VISITORS_H

#include <CGAL/disable_warnings.h>

namespace CGAL {

template <typename Previous_level>
class Null_mesh_visitor_level
{
  Previous_level& previous;
public:
  typedef Previous_level Previous_visitor;

  Null_mesh_visitor_level(Previous_level& p) : previous(p)
  {
  }

  Previous_level& previous_level()
  {
    return previous;
  }

  template <typename E, typename P>
  void before_conflicts(E, P) const {}

  template <typename E, typename P, typename Z>
  void before_insertion(E, P, Z) const {}

  template <typename V>
  void after_insertion(V) const {}

  template <typename E, typename P, typename Z>
  void after_no_insertion(E, P, Z) const {}
}; // end class Null_mesh_visitor_level

class Null_mesh_visitor
{
public:
  typedef Null_mesh_visitor Previous_visitor;

  const Null_mesh_visitor& previous_level() const 
  {
    return *this;
  }

  template <typename E, typename P>
  void before_conflicts(E, P) const {}

  template <typename E, typename P, typename Z>
  void before_insertion(E, P, Z) const {}

  template <typename V>
  void after_insertion(V) const {}

  template <typename E, typename P, typename Z>
  void after_no_insertion(E, P, Z) const {}
}; // end class Null_mesh_visitor

template <typename V1, typename V2>
struct Combine_mesh_visitor {
  typedef Combine_mesh_visitor<typename V1::Previous_visitor,
			       typename V2::Previous_visitor>
                                                       Previous_visitor;

  V1 v1;
  V2 v2;

  Combine_mesh_visitor(V1 v1_, V2 v2_)
    : v1(v1_), v2(v2_)
  {
  }

  Previous_visitor previous_level()
  { 
    return Previous_visitor(v1.previous_level(),
			    v2.previous_level());
  }

  template <typename E, typename P>
  void before_conflicts(E e, P p) 
  {
    v1.before_conflicts(e, p);
    v2.before_conflicts(e, p);    
  }

  template <typename E, typename P, typename Z>
  void before_insertion(E e, P p, Z z)
  {
    v1.before_insertion(e, p, z);
    v2.before_insertion(e, p, z);
  }

  template <typename V>
  void after_insertion(V v)
  {
    v1.after_insertion(v);
    v2.after_insertion(v);
  } 

  template <typename E, typename P, typename Z>
  void after_no_insertion(E e, P p, Z z)
  {
    v1.after_no_insertion(e, p, z);
    v2.after_no_insertion(e, p, z);
  }
}; // end class Combine_mesh_visitor

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESHER_LEVEL_VISITORS_H
