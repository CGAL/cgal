// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_DUMMY_ITERATOR_H
#define CGAL_VORONOI_DIAGRAM_2_DUMMY_ITERATOR_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/iterator.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

template<class Value_t>
class Dummy_iterator : public Emptyset_iterator
{
 private:
  typedef Dummy_iterator<Value_t>  Self;

 public:
  typedef Value_t                         value_type;
  typedef value_type&                     reference;
  typedef value_type*                     pointer;
  typedef const value_type&               const_reference;
  typedef const value_type*               const_pointer;
  typedef std::size_t                     size_type;
  typedef std::ptrdiff_t                  difference_type;
  typedef std::bidirectional_iterator_tag iterator_category;

  Dummy_iterator() {}
  Dummy_iterator(const Dummy_iterator&) {}

  template< class T >
  Self& operator=(const T&) { return *this; }

  Self& operator++()        { return *this; }
  Self& operator++(int)     { return *this; }

  Self& operator--()        { return *this; }
  Self& operator--(int)     { return *this; }

  reference operator*()              { return *dummy_pointer(); }
  pointer   operator->()             { return dummy_pointer(); }

  const_reference operator*()  const { return *dummy_pointer(); }
  const_pointer   operator->() const { return dummy_pointer(); }

  bool operator==(const Self&) const { return true; }
  bool operator!=(const Self&) const { return false; }

  bool operator<(const Self& other) const {
    return this < &other;
  }

  static const_reference dummy_reference() {
    static value_type dummy_reference_static;
    return dummy_reference_static;
  }
 private:
  static pointer dummy_pointer() {
    static value_type dummy_pointer_static;
    return &dummy_pointer_static;
  }
};

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL


#endif // CGAL_VORONOI_DIAGRAM_2_DUMMY_ITERATOR_H
