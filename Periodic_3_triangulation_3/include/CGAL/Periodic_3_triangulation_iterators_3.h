// Copyright (c) 1999,2008-2009   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_TRIANGULATION_ITERATORS_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_ITERATORS_3_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/array.h>
#include <CGAL/iterator.h>

namespace CGAL {

template < class T >
class Periodic_3_triangulation_tetrahedron_iterator_3 {
// Iterates over the primitives in a periodic triangulation.
// Options:
// - STORED: output each primitive from the TDS exactly once
// - UNIQUE: output exactly one periodic copy of each primitive, no matter
//     whether the current tds stores a n-sheeted covering for n!=1.
// - STORED_COVER_DOMAIN: output each primitive whose intersection with the
//     actually used periodic domain is non-zero.
// - UNIQUE_COVER_DOMAIN: output each primitive whose intersection
//     with the original domain that the user has given is non-zero
//
// Comments:
// When computing in 1-sheeted covering, there will be no difference in the
// result of STORED and UNIQUE as well as STORED_COVER_DOMAIN and
// UNIQUE_COVER_DOMAIN.

public:

  typedef typename T::Periodic_tetrahedron                value_type;
  typedef const typename T::Periodic_tetrahedron *        pointer;
  typedef const typename T::Periodic_tetrahedron &        reference;
  typedef std::size_t                                     size_type;
  typedef std::ptrdiff_t                                  difference_type;
  typedef std::bidirectional_iterator_tag                 iterator_category;

  typedef typename T::Periodic_tetrahedron                Periodic_tetrahedron;
  typedef Periodic_3_triangulation_tetrahedron_iterator_3<T>
                                                 Periodic_tetrahedron_iterator;
  typedef typename T::Cell                                Cell;
  typedef typename T::Cell_iterator                       Cell_iterator;

  typedef typename T::Offset                              Offset;
  typedef typename T::Iterator_type                       Iterator_type;

  Periodic_3_triangulation_tetrahedron_iterator_3(Iterator_type it = T::STORED)
    : _t(NULL), _it(it), _off(0) {}

  Periodic_3_triangulation_tetrahedron_iterator_3(const T * t,
      Iterator_type it = T::STORED)
    : _t(t), pos(_t->cells_begin()), _it(it), _off(0) {
    if (_it == T::UNIQUE || _it == T::UNIQUE_COVER_DOMAIN) {
      while (pos != _t->cells_end() && !is_canonical() )
	++pos;
    }
  }

  // used to initialize the past-the-end iterator
  Periodic_3_triangulation_tetrahedron_iterator_3(const T* t, int,
					      Iterator_type it = T::STORED)
    : _t(t), pos(_t->cells_end()), _it(it), _off(0) {}

  Periodic_tetrahedron_iterator& operator++() {
    switch (_it) {
    case T::STORED:
      ++pos;
      break;
    case T::UNIQUE:
      do { ++pos; } while (pos != _t->cells_end() && !is_canonical());
      break;
    case T::STORED_COVER_DOMAIN:
    case T::UNIQUE_COVER_DOMAIN:
      increment_domain();
      break;
    default:
      CGAL_triangulation_assertion(false);
    };
    return *this;
  }

  Periodic_tetrahedron_iterator& operator--() {
    switch (_it) {
    case T::STORED:
      --pos;
      break;
    case T::UNIQUE:
      do { --pos; } while (pos != _t->cells_begin() && !is_canonical());
      break;
    case T::STORED_COVER_DOMAIN:
    case T::UNIQUE_COVER_DOMAIN:
      decrement_domain();
    };
    return *this;
  }

  Periodic_tetrahedron_iterator operator++(int)
    {
      Periodic_tetrahedron_iterator tmp(*this);
      ++(*this);
      return tmp;
    }

  Periodic_tetrahedron_iterator operator--(int)
    {
      Periodic_tetrahedron_iterator tmp(*this);
      --(*this);
      return tmp;
    }

  bool operator==(const Periodic_tetrahedron_iterator& ti) const
    {
      CGAL_triangulation_assertion(_it == ti._it);
      return _t == ti._t && pos == ti.pos && _off == ti._off;
    }

  bool operator!=(const Periodic_tetrahedron_iterator& ti) const
    {
      return !(*this == ti);
    }

  reference operator*() const
    {
      periodic_tetrahedron = construct_periodic_tetrahedron();
      return periodic_tetrahedron;
    }

  pointer operator->() const
    {
      periodic_tetrahedron = construct_periodic_tetrahedron();
      return &periodic_tetrahedron;
    }

  Cell_iterator get_cell() const
  {
    return pos;
  }

private:
  const T*  _t;
  Cell_iterator pos; // current cell.
  Iterator_type _it;
  int _off; // current offset
  mutable Periodic_tetrahedron periodic_tetrahedron; // current tetrahedron.

private:
  // check whether pos points onto a unique edge or not.
  // If we are computing in 1-sheeted covering this should
  // always be true.
  bool is_canonical() {
    // fetch all offsets
    Offset off0, off1, off2, off3;
    get_edge_offsets(off0, off1, off2, off3);
    
    if (_t->number_of_sheets() != make_array(1,1,1)) {
      // If there is one offset with entries larger than 1 then we are
      // talking about a vertex that is too far away from the original
      // domain to belong to a canonical triangle.
      if (off0.x() > 1) return false;
      if (off0.y() > 1) return false;
      if (off0.z() > 1) return false;
      if (off1.x() > 1) return false;
      if (off1.y() > 1) return false;
      if (off1.z() > 1) return false;
      if (off2.x() > 1) return false;
      if (off2.y() > 1) return false;
      if (off2.z() > 1) return false;
      if (off3.x() > 1) return false;
      if (off3.y() > 1) return false;
      if (off3.z() > 1) return false;
    }

    // If there is one direction of space for which all offsets are
    // non-zero then the edge is not canonical because we can
    // take the copy closer towards the origin in that direction.
    int offx = off0.x() & off1.x() & off2.x() & off3.x();
    int offy = off0.y() & off1.y() & off2.y() & off3.y();
    int offz = off0.z() & off1.z() & off2.z() & off3.z();

    return (offx == 0 && offy == 0 && offz == 0);
  }

  // Artificial incrementation function that takes periodic
  // copies into account.
  void increment_domain() {
    int off = get_drawing_offsets();
    CGAL_triangulation_assertion(_off <= off);
    if (_off == off) {
      _off = 0;
      do { ++pos; } while (_it == T::UNIQUE_COVER_DOMAIN
			   && pos != _t->cells_end() && !is_canonical());
    } else {
      do {
	++_off;
      } while ((((~_off)|off)&7)!=7); // Increment until a valid
				      // offset has been found
    }
  }

  // Artificial decrementation function that takes periodic
  // copies into account.
  void decrement_domain() {
    if (_off == 0) {
      if (pos == _t->cells_begin()) return;
      do { --pos; } while (_it == T::UNIQUE_COVER_DOMAIN && !is_canonical());
      _off = get_drawing_offsets();
    } else {
      int off = get_drawing_offsets();
      do {
	--_off;
      } while ((((~_off)|off)&7)!=7); // Decrement until a valid
				      // offset has been found
    }
  }

  // Get the canonicalized offsets of an edge.
  // This works in any cover that is encoded in _t->combine_offsets
  void get_edge_offsets(Offset &off0, Offset &off1,
			Offset &off2, Offset &off3) const {
    Offset cell_off0 = _t->int_to_off(pos->offset(0));
    Offset cell_off1 = _t->int_to_off(pos->offset(1));
    Offset cell_off2 = _t->int_to_off(pos->offset(2));
    Offset cell_off3 = _t->int_to_off(pos->offset(3));
    Offset diff_off((cell_off0.x() == 1 
		     && cell_off1.x() == 1 
		     && cell_off2.x() == 1
		     && cell_off3.x() == 1)?-1:0,
		    (cell_off0.y() == 1 
		     && cell_off1.y() == 1
		     && cell_off2.y() == 1
		     && cell_off3.y() == 1)?-1:0,
		    (cell_off0.z() == 1 
		     && cell_off1.z() == 1
		     && cell_off2.z() == 1
		     && cell_off3.z() == 1)?-1:0);
    off0 = _t->combine_offsets(_t->get_offset(pos,0), diff_off);
    off1 = _t->combine_offsets(_t->get_offset(pos,1), diff_off);
    off2 = _t->combine_offsets(_t->get_offset(pos,2), diff_off);
    off3 = _t->combine_offsets(_t->get_offset(pos,3), diff_off);
  }

  // return an integer that encodes the translations which have to be
  // applied to the edge *pos
  int get_drawing_offsets() {
    Offset off0, off1, off2, off3;
    // Choose edges that are to be duplicated. These are edges that
    // intersect the boundary of the periodic domain. In UNIQUE mode
    // this means that the offset with respect to drawing should
    // differ in some entries. Otherwise we consider the offsets
    // internally stored inside the cell telling us that this cell
    // wraps around the domain.
    if (_it == T::UNIQUE_COVER_DOMAIN)
      get_edge_offsets(off0,off1,off2,off3);
    else {
      CGAL_triangulation_assertion(_it == T::STORED_COVER_DOMAIN);
      off0 = _t->int_to_off(pos->offset(0));
      off1 = _t->int_to_off(pos->offset(1));
      off2 = _t->int_to_off(pos->offset(2));
      off3 = _t->int_to_off(pos->offset(3));
    }

    CGAL_triangulation_assertion(off0.x() == 0 || off0.x() == 1);
    CGAL_triangulation_assertion(off0.y() == 0 || off0.y() == 1);
    CGAL_triangulation_assertion(off0.z() == 0 || off0.z() == 1);
    CGAL_triangulation_assertion(off1.x() == 0 || off1.x() == 1);
    CGAL_triangulation_assertion(off1.y() == 0 || off1.y() == 1);
    CGAL_triangulation_assertion(off1.z() == 0 || off1.z() == 1);
    CGAL_triangulation_assertion(off2.x() == 0 || off2.x() == 1);
    CGAL_triangulation_assertion(off2.y() == 0 || off2.y() == 1);
    CGAL_triangulation_assertion(off2.z() == 0 || off2.z() == 1);
    CGAL_triangulation_assertion(off3.x() == 0 || off3.x() == 1);
    CGAL_triangulation_assertion(off3.y() == 0 || off3.y() == 1);
    CGAL_triangulation_assertion(off3.z() == 0 || off3.z() == 1);
    
    int offx = ( ((off0.x() == 0 && off1.x() == 0 
		   && off2.x() == 0 && off3.x() == 0)
	       || (off0.x() == 1 && off1.x() == 1 
		   && off2.x() == 1 && off3.x() == 1)) ? 0 : 1);
    int offy = ( ((off0.y() == 0 && off1.y() == 0 
		   && off2.y() == 0 && off3.y() == 0)
	       || (off0.y() == 1 && off1.y() == 1 
		   && off2.y() == 1 && off3.y() == 1)) ? 0 : 1);
    int offz = ( ((off0.z() == 0 && off1.z() == 0 
		   && off2.z() == 0 && off3.z() == 0)
	       || (off0.z() == 1 && off1.z() == 1 
		   && off2.z() == 1 && off3.z() == 1)) ? 0 : 1);
    
    return( 4*offx + 2*offy + offz );
  }

  Periodic_tetrahedron construct_periodic_tetrahedron() const {
    CGAL_triangulation_assertion(pos != typename T::Cell_handle());
    Offset off0, off1, off2, off3;
    get_edge_offsets(off0, off1, off2, off3);
    Offset transl_off = Offset((((_off>>2)&1)==1 ? -1:0),
			       (((_off>>1)&1)==1 ? -1:0),
    			       (( _off    &1)==1 ? -1:0));
    if (_it == T::STORED_COVER_DOMAIN) {
      off0 = _t->combine_offsets(off0,transl_off);
      off1 = _t->combine_offsets(off1,transl_off);
      off2 = _t->combine_offsets(off2,transl_off);
      off3 = _t->combine_offsets(off3,transl_off);
    }
    if (_it == T::UNIQUE_COVER_DOMAIN) {
      off0 += transl_off;
      off1 += transl_off;
      off2 += transl_off;
      off3 += transl_off;
    }
    return make_array(
	std::make_pair(pos->vertex(0)->point(),off0),
	std::make_pair(pos->vertex(1)->point(),off1),
	std::make_pair(pos->vertex(2)->point(),off2),
	std::make_pair(pos->vertex(3)->point(),off3));
  }
};

template < class T >
class Periodic_3_triangulation_triangle_iterator_3 {
// Iterates over the primitives in a periodic triangulation.
// Options:
// - STORED: output each primitive from the TDS exactly once
// - UNIQUE: output exactly one periodic copy of each primitive, no matter
//     whether the current tds stores a n-sheeted covering for n!=1.
// - STORED_COVER_DOMAIN: output each primitive whose intersection with the
//     actually used periodic domain is non-zero.
// - UNIQUE_COVER_DOMAIN: output each primitive whose intersection
//     with the original domain that the user has given is non-zero
//
// Comments:
// When computing in 1-sheeted covering, there will be no difference in the
// result of STORED and UNIQUE as well as STORED_COVER_DOMAIN and
// UNIQUE_COVER_DOMAIN.

public:

  typedef typename T::Periodic_triangle                   value_type;
  typedef const typename T::Periodic_triangle *           pointer;
  typedef const typename T::Periodic_triangle &           reference;
  typedef std::size_t                                     size_type;
  typedef std::ptrdiff_t                                  difference_type;
  typedef std::bidirectional_iterator_tag                 iterator_category;

  typedef typename T::Periodic_triangle                   Periodic_triangle;
  typedef Periodic_3_triangulation_triangle_iterator_3<T>
                                                    Periodic_triangle_iterator;
  typedef typename T::Facet                               Facet;
  typedef typename T::Facet_iterator                      Facet_iterator;

  typedef typename T::Offset                              Offset;
  typedef typename T::Iterator_type                       Iterator_type;

  Periodic_3_triangulation_triangle_iterator_3(Iterator_type it = T::STORED)
    : _t(NULL), _it(it), _off(0) {}

  Periodic_3_triangulation_triangle_iterator_3(const T * t,
					      Iterator_type it = T::STORED)
    : _t(t), pos(_t->facets_begin()), _it(it), _off(0) {
    if (_it == T::UNIQUE || _it == T::UNIQUE_COVER_DOMAIN) {
      while (pos != _t->facets_end() && !is_canonical() )
	++pos;
    }
  }

  // used to initialize the past-the-end iterator
  Periodic_3_triangulation_triangle_iterator_3(const T* t, int,
					      Iterator_type it = T::STORED)
    : _t(t), pos(_t->facets_end()), _it(it), _off(0) {}

  Periodic_triangle_iterator& operator++() {
    switch (_it) {
    case T::STORED:
      ++pos;
      break;
    case T::UNIQUE:
      do { ++pos; } while (pos != _t->facets_end() && !is_canonical());
      break;
    case T::STORED_COVER_DOMAIN:
    case T::UNIQUE_COVER_DOMAIN:
      increment_domain();
      break;
    default:
      CGAL_triangulation_assertion(false);
    };
    return *this;
  }

  Periodic_triangle_iterator& operator--() {
    switch (_it) {
    case T::STORED:
      --pos;
      break;
    case T::UNIQUE:
      do { --pos; } while (pos != _t->facets_begin() && !is_canonical());
      break;
    case T::STORED_COVER_DOMAIN:
    case T::UNIQUE_COVER_DOMAIN:
      decrement_domain();
    };
    return *this;
  }

  Periodic_triangle_iterator operator++(int)
    {
      Periodic_triangle_iterator tmp(*this);
      ++(*this);
      return tmp;
    }

  Periodic_triangle_iterator operator--(int)
    {
      Periodic_triangle_iterator tmp(*this);
      --(*this);
      return tmp;
    }

  bool operator==(const Periodic_triangle_iterator& ti) const
    {
      CGAL_triangulation_assertion(_it == ti._it);
      return _t == ti._t && pos == ti.pos && _off == ti._off;
    }

  bool operator!=(const Periodic_triangle_iterator& ti) const
    {
      return !(*this == ti);
    }

  reference operator*() const
    {
      periodic_triangle = construct_periodic_triangle();
      return periodic_triangle;
    }

  pointer operator->() const
    {
      periodic_triangle = construct_periodic_triangle();
      return &periodic_triangle;
    }

  Facet_iterator get_facet() const
  {
    return pos;
  }

private:
  const T*  _t;
  Facet_iterator pos; // current facet.
  Iterator_type _it;
  int _off; // current offset
  mutable Periodic_triangle periodic_triangle; // current segment.

private:
  // check whether pos points onto a unique edge or not.
  // If we are computing in 1-sheeted covering this should
  // always be true.
  bool is_canonical() {
    // fetch all offsets
    Offset off0, off1, off2;
    get_edge_offsets(off0, off1, off2);
    
    if (_t->number_of_sheets() != make_array(1,1,1)) {
      // If there is one offset with entries larger than 1 then we are
      // talking about a vertex that is too far away from the original
      // domain to belong to a canonical triangle.
      if (off0.x() > 1) return false;
      if (off0.y() > 1) return false;
      if (off0.z() > 1) return false;
      if (off1.x() > 1) return false;
      if (off1.y() > 1) return false;
      if (off1.z() > 1) return false;
      if (off2.x() > 1) return false;
      if (off2.y() > 1) return false;
      if (off2.z() > 1) return false;
    }

    // If there is one direction of space for which all offsets are
    // non-zero then the edge is not canonical because we can
    // take the copy closer towards the origin in that direction.
    int offx = off0.x() & off1.x() & off2.x();
    int offy = off0.y() & off1.y() & off2.y();
    int offz = off0.z() & off1.z() & off2.z();

    return (offx == 0 && offy == 0 && offz == 0);
  }

  // Artificial incrementation function that takes periodic
  // copies into account.
  void increment_domain() {
    int off = get_drawing_offsets();
    CGAL_triangulation_assertion(_off <= off);
    if (_off == off) {
      _off = 0;
      do { ++pos; } while (_it == T::UNIQUE_COVER_DOMAIN
			   && pos != _t->facets_end() && !is_canonical());
    } else {
      do {
	++_off;
      } while ((((~_off)|off)&7)!=7); // Increment until a valid
				      // offset has been found
    }
  }

  // Artificial decrementation function that takes periodic
  // copies into account.
  void decrement_domain() {
    if (_off == 0) {
      if (pos == _t->facets_begin()) return;
      do { --pos; } while (_it == T::UNIQUE_COVER_DOMAIN && !is_canonical());
      _off = get_drawing_offsets();
    } else {
      int off = get_drawing_offsets();
      do {
	--_off;
      } while ((((~_off)|off)&7)!=7); // Decrement until a valid
				      // offset has been found
    }
  }

  // Get the canonicalized offsets of an edge.
  // This works in any cover that is encoded in _t->combine_offsets
  void get_edge_offsets(Offset &off0, Offset &off1, Offset &off2) const {
    Offset cell_off0 = _t->int_to_off(pos->first->offset((pos->second+1)&3));
    Offset cell_off1 = _t->int_to_off(pos->first->offset((pos->second+2)&3));
    Offset cell_off2 = _t->int_to_off(pos->first->offset((pos->second+3)&3));
    Offset diff_off((cell_off0.x() == 1 
		     && cell_off1.x() == 1 
		     && cell_off2.x() == 1)?-1:0,
		    (cell_off0.y() == 1 
		     && cell_off1.y() == 1
		     && cell_off2.y() == 1)?-1:0,
		    (cell_off0.z() == 1 
		     && cell_off1.z() == 1
		     && cell_off2.z() == 1)?-1:0);
    off0 = _t->combine_offsets(_t->get_offset(pos->first,
						     (pos->second+1)&3),
			       diff_off);
    off1 = _t->combine_offsets(_t->get_offset(pos->first,
						     (pos->second+2)&3),
			       diff_off);
    off2 = _t->combine_offsets(_t->get_offset(pos->first,
						     (pos->second+3)&3),
			       diff_off);
  }

  // return an integer that encodes the translations which have to be
  // applied to the edge *pos
  int get_drawing_offsets() {
    Offset off0, off1, off2;
    // Choose edges that are to be duplicated. These are edges that
    // intersect the boundary of the periodic domain. In UNIQUE mode
    // this means that the offset with respect to drawing should
    // differ in some entries. Otherwise we consider the offsets
    // internally stored inside the cell telling us that this cell
    // wraps around the domain.
    if (_it == T::UNIQUE_COVER_DOMAIN)
      get_edge_offsets(off0,off1,off2);
    else {
      CGAL_triangulation_assertion(_it == T::STORED_COVER_DOMAIN);
      off0 = _t->int_to_off(pos->first->offset((pos->second+1)&3));
      off1 = _t->int_to_off(pos->first->offset((pos->second+2)&3));
      off2 = _t->int_to_off(pos->first->offset((pos->second+3)&3));
    }

    CGAL_triangulation_assertion(off0.x() == 0 || off0.x() == 1);
    CGAL_triangulation_assertion(off0.y() == 0 || off0.y() == 1);
    CGAL_triangulation_assertion(off0.z() == 0 || off0.z() == 1);
    CGAL_triangulation_assertion(off1.x() == 0 || off1.x() == 1);
    CGAL_triangulation_assertion(off1.y() == 0 || off1.y() == 1);
    CGAL_triangulation_assertion(off1.z() == 0 || off1.z() == 1);
    CGAL_triangulation_assertion(off2.x() == 0 || off2.x() == 1);
    CGAL_triangulation_assertion(off2.y() == 0 || off2.y() == 1);
    CGAL_triangulation_assertion(off2.z() == 0 || off2.z() == 1);
    
    int offx = ( ((off0.x() == 0 && off1.x() == 0 && off2.x() == 0)
	       || (off0.x() == 1 && off1.x() == 1 && off2.x() == 1)) ? 0 : 1);
    int offy = ( ((off0.y() == 0 && off1.y() == 0 && off2.y() == 0)
	       || (off0.y() == 1 && off1.y() == 1 && off2.y() == 1)) ? 0 : 1);
    int offz = ( ((off0.z() == 0 && off1.z() == 0 && off2.z() == 0)
	       || (off0.z() == 1 && off1.z() == 1 && off2.z() == 1)) ? 0 : 1);
    
    return( 4*offx + 2*offy + offz );
  }

  Periodic_triangle construct_periodic_triangle() const {
    CGAL_triangulation_assertion(pos->first != typename T::Cell_handle());
    Offset off0, off1, off2;
    get_edge_offsets(off0, off1, off2);
    Offset transl_off = Offset((((_off>>2)&1)==1 ? -1:0),
			       (((_off>>1)&1)==1 ? -1:0),
    			       (( _off    &1)==1 ? -1:0));
    if (_it == T::STORED_COVER_DOMAIN) {
      off0 = _t->combine_offsets(off0,transl_off);
      off1 = _t->combine_offsets(off1,transl_off);
      off2 = _t->combine_offsets(off2,transl_off);
    }
    if (_it == T::UNIQUE_COVER_DOMAIN) {
      off0 += transl_off;
      off1 += transl_off;
      off2 += transl_off;
    }
    return make_array(
	std::make_pair(pos->first->vertex((pos->second+1)&3)->point(),off0),
	std::make_pair(pos->first->vertex((pos->second+2)&3)->point(),off1),
	std::make_pair(pos->first->vertex((pos->second+3)&3)->point(),off2));
  }
};

template < class T >
class Periodic_3_triangulation_segment_iterator_3 {
// Iterates over the primitives in a periodic triangulation.
// Options:
// - STORED: output each primitive from the TDS exactly once
// - UNIQUE: output exactly one periodic copy of each primitive, no matter
//     whether the current tds stores a n-sheeted covering for n!=1.
// - STORED_COVER_DOMAIN: output each primitive whose intersection with the
//     actually used periodic domain is non-zero.
// - UNIQUE_COVER_DOMAIN: output each primitive whose intersection
//     with the original domain that the user has given is non-zero
//
// Comments:
// When computing in 1-sheeted covering, there will be no difference in the
// result of STORED and UNIQUE as well as STORED_COVER_DOMAIN and
// UNIQUE_COVER_DOMAIN.

public:

  typedef typename T::Periodic_segment                    value_type;
  typedef const typename T::Periodic_segment *            pointer;
  typedef const typename T::Periodic_segment &            reference;
  typedef std::size_t                                     size_type;
  typedef std::ptrdiff_t                                  difference_type;
  typedef std::bidirectional_iterator_tag                 iterator_category;

  typedef typename T::Periodic_segment                    Periodic_segment;
  typedef Periodic_3_triangulation_segment_iterator_3<T>
                                                     Periodic_segment_iterator;
  typedef typename T::Edge                                Edge;
  typedef typename T::Edge_iterator                       Edge_iterator;

  typedef typename T::Offset                              Offset;
  typedef typename T::Iterator_type                       Iterator_type;

  Periodic_3_triangulation_segment_iterator_3(Iterator_type it = T::STORED)
    : _t(NULL), _it(it), _off(0) {}

  Periodic_3_triangulation_segment_iterator_3(const T * t,
					      Iterator_type it = T::STORED)
    : _t(t), pos(_t->edges_begin()), _it(it), _off(0) {
    if (_it == T::UNIQUE || _it == T::UNIQUE_COVER_DOMAIN) {
      while (pos != _t->edges_end() && !is_canonical() )
	++pos;
    }
  }

  // used to initialize the past-the-end iterator
  Periodic_3_triangulation_segment_iterator_3(const T* t, int,
					      Iterator_type it = T::STORED)
    : _t(t), pos(_t->edges_end()), _it(it), _off(0) {}

  Periodic_segment_iterator& operator++() {
    switch (_it) {
    case T::STORED:
      ++pos;
      break;
    case T::UNIQUE:
      do { ++pos; } while (pos != _t->edges_end() && !is_canonical());
      break;
    case T::STORED_COVER_DOMAIN:
    case T::UNIQUE_COVER_DOMAIN:
      increment_domain();
      break;
    default:
      CGAL_triangulation_assertion(false);
    };
    return *this;
  }

  Periodic_segment_iterator& operator--() {
    switch (_it) {
    case T::STORED:
      --pos;
      break;
    case T::UNIQUE:
      do { --pos; } while (pos != _t->edges_begin() && !is_canonical());
      break;
    case T::STORED_COVER_DOMAIN:
    case T::UNIQUE_COVER_DOMAIN:
      decrement_domain();
    };
    return *this;
  }

  Periodic_segment_iterator operator++(int)
    {
      Periodic_segment_iterator tmp(*this);
      ++(*this);
      return tmp;
    }

  Periodic_segment_iterator operator--(int)
    {
      Periodic_segment_iterator tmp(*this);
      --(*this);
      return tmp;
    }

  bool operator==(const Periodic_segment_iterator& ti) const
    {
      CGAL_triangulation_assertion(_it == ti._it);
      return _t == ti._t && pos == ti.pos && _off == ti._off;
    }

  bool operator!=(const Periodic_segment_iterator& ti) const
    {
      return !(*this == ti);
    }

  reference operator*() const
    {
      periodic_segment = construct_periodic_segment();
      return periodic_segment;
    }

  pointer operator->() const
    {
      periodic_segment = construct_periodic_segment();
      return &periodic_segment;
    }

  Edge_iterator get_edge() const
  {
    return pos;
  }
private:
  const T*  _t;
  Edge_iterator pos; // current edge.
  Iterator_type _it;
  int _off; // current offset
  mutable Periodic_segment periodic_segment; // current segment.

private:
  // check whether pos points onto a unique edge or not.
  // If we are computing in 1-sheeted covering this should
  // always be true.
  bool is_canonical() {
    // fetch all offsets
    Offset off0, off1;
    get_edge_offsets(off0, off1);
    
    if (_t->number_of_sheets() != make_array(1,1,1)) {
      // If there is one offset with entries larger than 1 then we are
      // talking about a vertex that is too far away from the original
      // domain to belong to a canonical triangle.
      if (off0.x() > 1) return false;
      if (off0.y() > 1) return false;
      if (off0.z() > 1) return false;
      if (off1.x() > 1) return false;
      if (off1.y() > 1) return false;
      if (off1.z() > 1) return false;
    }

    // If there is one direction of space for which all offsets are
    // non-zero then the edge is not canonical because we can
    // take the copy closer towards the origin in that direction.
    int offx = off0.x() & off1.x();
    int offy = off0.y() & off1.y();
    int offz = off0.z() & off1.z();

    return (offx == 0 && offy == 0 && offz == 0);
  }

  // Artificial incrementation function that takes periodic
  // copies into account.
  void increment_domain() {
    int off = get_drawing_offsets();
    CGAL_triangulation_assertion(_off <= off);
    if (_off == off) {
      _off = 0;
      do { ++pos; } while (_it == T::UNIQUE_COVER_DOMAIN
			   && pos != _t->edges_end() && !is_canonical());
    } else {
      do {
	++_off;
      } while ((((~_off)|off)&7)!=7); // Increment until a valid
				      // offset has been found
    }
  }

  // Artificial decrementation function that takes periodic
  // copies into account.
  void decrement_domain() {
    if (_off == 0) {
      if (pos == _t->edges_begin()) return;
      do { --pos; } while (_it == T::UNIQUE_COVER_DOMAIN && !is_canonical());
      _off = get_drawing_offsets();
    } else {
      int off = get_drawing_offsets();
      do {
	--_off;
      } while ((((~_off)|off)&7)!=7); // Decrement until a valid
				      // offset has been found
    }
  }

  // Get the canonicalized offsets of an edge.
  // This works in any cover that is encoded in _t->combine_offsets
  void get_edge_offsets(Offset &off0, Offset &off1) const {
    Offset cell_off0 = _t->int_to_off(pos->first->offset(pos->second));
    Offset cell_off1 = _t->int_to_off(pos->first->offset(pos->third));
    Offset diff_off((cell_off0.x()==1 && cell_off1.x()==1)?-1:0,
		    (cell_off0.y()==1 && cell_off1.y()==1)?-1:0,
		    (cell_off0.z()==1 && cell_off1.z()==1)?-1:0);
    off0 = _t->combine_offsets(_t->get_offset(pos->first,pos->second),
			       diff_off);
    off1 = _t->combine_offsets(_t->get_offset(pos->first,pos->third),
			       diff_off);
  }

  // return an integer that encodes the translations which have to be
  // applied to the edge *pos
  int get_drawing_offsets() {
    Offset off0, off1;
    // Choose edges that are to be duplicated. These are edges that
    // intersect the boundary of the periodic domain. In UNIQUE mode
    // this means that the offset with respect to drawing should
    // differ in some entries. Otherwise we consider the offsets
    // internally stored inside the cell telling us that this cell
    // wraps around the domain.
    if (_it == T::UNIQUE_COVER_DOMAIN)
      get_edge_offsets(off0,off1);
    else {
      CGAL_triangulation_assertion(_it == T::STORED_COVER_DOMAIN);
      off0 = _t->int_to_off(pos->first->offset(pos->second));
      off1 = _t->int_to_off(pos->first->offset(pos->third));
    }
    Offset diff_off = off0 - off1;
    
    CGAL_triangulation_assertion(diff_off.x() >= -1 || diff_off.x() <= 1);
    CGAL_triangulation_assertion(diff_off.y() >= -1 || diff_off.y() <= 1);
    CGAL_triangulation_assertion(diff_off.z() >= -1 || diff_off.z() <= 1);

    return( 4*(diff_off.x() == 0 ? 0:1)
	    + 2*(diff_off.y() == 0 ? 0:1)
	    + (diff_off.z() == 0 ? 0:1) );
  }

  Periodic_segment construct_periodic_segment() const {
    CGAL_triangulation_assertion(pos->first != typename T::Cell_handle());
    Offset off0, off1;
    get_edge_offsets(off0, off1);
    Offset transl_off = Offset((((_off>>2)&1)==1 ? -1:0),
			       (((_off>>1)&1)==1 ? -1:0),
    			       (( _off    &1)==1 ? -1:0));
    if (_it == T::STORED_COVER_DOMAIN) {
      off0 = _t->combine_offsets(off0,transl_off);
      off1 = _t->combine_offsets(off1,transl_off);
    }
    if (_it == T::UNIQUE_COVER_DOMAIN) {
      off0 += transl_off;
      off1 += transl_off;
    }
    return make_array(
	std::make_pair(pos->first->vertex(pos->second)->point(),off0),
	std::make_pair(pos->first->vertex(pos->third)->point(),off1));
  }
};

template < class T >
class Periodic_3_triangulation_point_iterator_3 {
// Iterates over the primitives in a periodic triangulation.
// Options:
// - STORED: output each primitive from the TDS exactly once
// - UNIQUE: output exactly one periodic copy of each primitive, no matter
//     whether the current tds stores a n-sheeted covering for n!=1.
// - STORED_COVER_DOMAIN: output each primitive whose intersection with the
//     actually used periodic domain is non-zero.
// - UNIQUE_COVER_DOMAIN: output each primitive whose intersection
//     with the original domain that the user has given is non-zero
//
// Comments:
// When computing in 1-sheeted covering, there will be no difference in the
// result of STORED and UNIQUE as well as STORED_COVER_DOMAIN and
// UNIQUE_COVER_DOMAIN.

public:
  typedef typename T::Periodic_point                      value_type;
  typedef const typename T::Periodic_point *              pointer;
  typedef const typename T::Periodic_point &              reference;
  typedef std::size_t                                     size_type;
  typedef std::ptrdiff_t                                  difference_type;
  typedef std::bidirectional_iterator_tag                 iterator_category;

  typedef typename T::Periodic_point                      Periodic_point;
  typedef Periodic_3_triangulation_point_iterator_3<T>  Periodic_point_iterator;

  typedef typename T::Vertex                              Vertex;
  typedef typename T::Vertex_iterator                     Vertex_iterator;

  typedef typename T::Offset                              Offset;
  typedef typename T::Iterator_type                       Iterator_type;

  Periodic_3_triangulation_point_iterator_3(Iterator_type it = T::STORED)
    : _t(NULL), _it(it) {}

  Periodic_3_triangulation_point_iterator_3(const T * t,
					    Iterator_type it = T::STORED)
    : _t(t), pos(_t->vertices_begin()), _it(it) {
    if (_it == T::UNIQUE || _it == T::UNIQUE_COVER_DOMAIN) {
      while (pos != _t->vertices_end() && !is_canonical() )
	++pos;
    }
  }

  // used to initialize the past-the-end iterator
  Periodic_3_triangulation_point_iterator_3(const T* t, int,
					    Iterator_type it = T::STORED)
    : _t(t), pos(_t->vertices_end()), _it(it) {}

  Periodic_point_iterator& operator++() {
    switch (_it) {
    case T::STORED:
    case T::STORED_COVER_DOMAIN:
      ++pos;
      break;
    case T::UNIQUE:
    case T::UNIQUE_COVER_DOMAIN:
      do { ++pos; } while (pos != _t->vertices_end() && !is_canonical());
      break;
    default:
      CGAL_triangulation_assertion(false);
    };
    return *this;
  }

  Periodic_point_iterator& operator--() {
    switch (_it) {
    case T::STORED:
    case T::STORED_COVER_DOMAIN:
      --pos;
      break;
    case T::UNIQUE:
    case T::UNIQUE_COVER_DOMAIN:
      do { --pos; } while (pos != _t->vertices_begin() && !is_canonical());
      break;
    default:
      CGAL_triangulation_assertion(false);
    };
    return *this;
  }

  Periodic_point_iterator operator++(int)
  {
    Periodic_point_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
  
  Periodic_point_iterator operator--(int)
  {
    Periodic_point_iterator tmp(*this);
    --(*this);
    return tmp;
  }
  
  bool operator==(const Periodic_point_iterator& pi) const
  {
    CGAL_triangulation_assertion(_it == pi._it);
    return _t == pi._t && pos == pi.pos;
  }
  
  bool operator!=(const Periodic_point_iterator& pi) const
  {
    return !(*this == pi);
  }
  
  reference operator*() const
  {
    periodic_point = construct_periodic_point();
    return periodic_point;
  }
  
  pointer operator->() const
  {
    periodic_point = construct_periodic_point();
    return &periodic_point;
  }
  
  Vertex_iterator get_vertex() const
  {
    return pos;
  }
private:
  const T*  _t;
  Vertex_iterator pos; // current vertex.
  Iterator_type _it;
  int _off; // current offset
  mutable Periodic_point periodic_point; // current point.

private:
  // check whether pos points onto a vertex inside the original
  // domain. If we are computing in 1-sheeted covering this should
  // always be true.
  bool is_canonical() {
    return (_t->get_offset(pos) == Offset(0,0,0));
  }

  Periodic_point construct_periodic_point() const {
    CGAL_triangulation_assertion(pos != typename T::Vertex_handle());
    Offset off = _t->get_offset(pos);
    return std::make_pair(pos->point(),off);
  }
};

template <class T>
class Domain_tester {  
  const T *t;

public:
  Domain_tester() {}
  Domain_tester(const T *tr) : t(tr) {}

  bool operator()(const typename T::Vertex_iterator & v) const {
    return (t->get_offset(v) != typename T::Offset(0,0,0));
  }
};

// Iterates over the vertices in a periodic triangulation that are
// located inside the original cube.
// Derives from Filter_iterator in order to add a conversion to handle
//
// Comments:
// When computing in 1-sheeted covering, there will be no difference
// between a normal Vertex_iterator and this iterator
template <class T>
class Periodic_3_triangulation_unique_vertex_iterator_3
  : public Filter_iterator<typename T::Vertex_iterator, Domain_tester<T> > {

  typedef typename T::Vertex_handle Vertex_handle;
  typedef typename T::Vertex_iterator Vertex_iterator;

  typedef Filter_iterator<Vertex_iterator, Domain_tester<T> > Base;
  typedef Periodic_3_triangulation_unique_vertex_iterator_3 Self;
public:

  Periodic_3_triangulation_unique_vertex_iterator_3() : Base() {}
  Periodic_3_triangulation_unique_vertex_iterator_3(const Base &b) : Base(b) {}

  Self & operator++() { Base::operator++(); return *this; }
  Self & operator--() { Base::operator--(); return *this; }
  Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
  Self operator--(int) { Self tmp(*this); --(*this); return tmp; }

  operator Vertex_handle() const { return Base::base(); }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_ITERATORS_3_H
