// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>

#ifndef CGAL_POLYGON_2_EDGE_CIRCULATOR_H
#define CGAL_POLYGON_2_EDGE_CIRCULATOR_H

#include <iterator>
#include <CGAL/circulator.h>
#include <CGAL/Polygon_2/Polygon_2_vertex_circulator.h>
#include <CGAL/Polygon_2/polygon_assertions.h>

namespace CGAL {
#ifndef DOXYGEN_RUNNING //to avoid conflicts
template <class _Traits, class _Container> class Polygon_2;
#endif
template <class _Traits, class _Container>
class Polygon_2_const_edge_circulator {
  public:
    typedef _Traits Traits;
    typedef typename _Traits::Segment_2 Segment_2;
    typedef _Container Container;
//    typedef Bidirectional_const_circulator_from_container<_Container>
    typedef Polygon_circulator<_Container>
            Vertex_const_circulator;

    typedef Segment_2                            value_type;
    typedef typename _Container::difference_type difference_type;
    typedef typename _Container::size_type       size_type;
    typedef Segment_2*                           pointer;
    typedef const Segment_2*                     const_pointer;
    typedef Segment_2&                           reference;
    typedef const Segment_2&                     const_reference;
    typedef Bidirectional_circulator_tag    iterator_category;

  private:
    Vertex_const_circulator first_vertex;
    // carry around the segment to be a proper iterator
    Segment_2               segment;

  public:
    Polygon_2_const_edge_circulator() {}

    Polygon_2_const_edge_circulator(Vertex_const_circulator f)
      : first_vertex(f) {}

  bool operator==( Nullptr_t CGAL_assertion_code(p) ) const {
      CGAL_polygon_assertion( p == 0);
      return (first_vertex == 0);
    }

    bool operator!=( Nullptr_t p ) const
    {
      return !(*this == p);
    }

    bool
    operator==(
      const Polygon_2_const_edge_circulator<_Traits, _Container>& x) const
    {
      return first_vertex == x.first_vertex;
    }

    bool
    operator!=(
      const Polygon_2_const_edge_circulator<_Traits, _Container>& x) const
    {
      return !(first_vertex == x.first_vertex);
    }

    const Segment_2& operator*() const
    {
      Vertex_const_circulator second_vertex = first_vertex;
      ++second_vertex;
      typename Traits::Construct_segment_2 construct_segment_2 = 
            Traits().construct_segment_2_object();
      const_cast<Polygon_2_const_edge_circulator*>(this)->segment = 
          construct_segment_2(*first_vertex, *second_vertex);
      return segment;
    }

    const Segment_2* operator->() const 
    {
      return &(**this);
    }

    Polygon_2_const_edge_circulator<_Traits, _Container>& operator++()
    {
      ++first_vertex;
      return *this;
    }

    Polygon_2_const_edge_circulator<_Traits, _Container> operator++(int)
    {
      Polygon_2_const_edge_circulator<_Traits, _Container> tmp = *this;
      ++*this;
      return tmp;
    }

    Polygon_2_const_edge_circulator<_Traits, _Container>& operator--()
    {
      --first_vertex;
      return *this;
    }

    Polygon_2_const_edge_circulator<_Traits, _Container> operator--(int)
    {
      Polygon_2_const_edge_circulator<_Traits, _Container> tmp = *this;
      --*this;
      return tmp;
    }

// random access iterator requirements
    Polygon_2_const_edge_circulator<_Traits, _Container>&
    operator+=(difference_type n)
    {
      first_vertex += n;
      return *this;
    }

    Polygon_2_const_edge_circulator<_Traits, _Container>
    operator+(difference_type n) const
    {
      return Polygon_2_const_edge_circulator<_Traits, _Container>(
        first_vertex + n);
    }

    Polygon_2_const_edge_circulator<_Traits, _Container>&
    operator-=(difference_type n)
    {
      return (*this) -= n;
    }

    Polygon_2_const_edge_circulator<_Traits, _Container>
    operator-(difference_type n) const
    {
      return Polygon_2_const_edge_circulator<_Traits, _Container>(
        first_vertex - n);
    }

    difference_type
    operator-(
      const Polygon_2_const_edge_circulator<_Traits, _Container>& a) const
    {
      return first_vertex - a.first_vertex;
    }

    Segment_2 operator[](int n) const
    {
      return *Polygon_2_const_edge_circulator<_Traits, _Container>(
        first_vertex+n);
    }

    bool operator<(
      const Polygon_2_const_edge_circulator<_Traits, _Container>& a) const
    {
      return first_vertex < a.first_vertex;
    }

    bool operator>(
      const Polygon_2_const_edge_circulator<_Traits, _Container>& a) const
    {
      return first_vertex > a.first_vertex;
    }

    bool operator<=(
      const Polygon_2_const_edge_circulator<_Traits, _Container>& a) const
    {
      return first_vertex <= a.first_vertex;
    }

    bool operator>=(
      const Polygon_2_const_edge_circulator<_Traits, _Container>& a) const
    {
      return first_vertex >= a.first_vertex;
    }

};

/*
template <class _Traits, class _Container>
typename _Container::difference_type
distance_type(
  const Polygon_2_const_edge_circulator<_Traits,_Container>&)
{
  return typename _Container::difference_type();
}

template <class _Traits, class _Container>
typename _Traits::Segment_2*
value_type(const Polygon_2_const_edge_circulator<_Traits,_Container>&)
{
  return (typename _Traits::Segment_2*)(0);
}
*/

//-----------------------------------------------------------------------//
//                          implementation
//-----------------------------------------------------------------------//

//--------------------------------------------------------------------//
// I don't know how to implement the following function:
//
// template <class _Traits, class _Container>
// inline
// Polygon_2_const_edge_circulator<_Traits, _Container>
// operator+(_Container::difference_type n,
//           Polygon_2_const_edge_circulator<_Traits, _Container>& a)
// { return a+n; }
//--------------------------------------------------------------------//

} //namespace CGAL

#endif
