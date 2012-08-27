// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_TRIANGLE_ACCESSOR_3_H
#define CGAL_TRIANGLE_ACCESSOR_3_H

#include <CGAL/Polyhedron_3.h>

namespace CGAL {

template <typename Polyhedron, typename K>
class Triangle_accessor_3 {
  // This class template should never be instantiated, but only its
  // specializations.
  typedef typename Polyhedron::Error_bad_match Error_bad_match;
};


template < class K,class Items,
           template < class T, class I, class A>
           class T_HDS,
           class Alloc>
class Triangle_accessor_3<Polyhedron_3<K,Items,T_HDS,Alloc>, K >
{
  typedef Polyhedron_3<K,Items,T_HDS,Alloc> Polyhedron;
public:
  typedef typename K::Triangle_3                    Triangle_3;
  typedef typename Polyhedron::Facet_const_iterator Triangle_iterator;
  typedef typename Polyhedron::Facet_const_handle   Triangle_handle;

  Triangle_accessor_3() { }

  Triangle_iterator triangles_begin(const Polyhedron& p) const
  {
    return p.facets_begin();
  }

  Triangle_iterator triangles_end(const Polyhedron& p) const
  {
    return p.facets_end();
  }

  Triangle_3 triangle(const Triangle_handle& handle) const
  {
    typedef typename K::Point_3 Point;
    const Point& a = handle->halfedge()->vertex()->point();
    const Point& b = handle->halfedge()->next()->vertex()->point();
    const Point& c = handle->halfedge()->next()->next()->vertex()->point();
    return Triangle_3(a,b,c);
  }
};


} // end namespace CGAL


#endif // POLYHEDRON_TRIANGLE_ACCESSOR_H
