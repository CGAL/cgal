// Copyright (c) 1999-2016   INRIA Nancy - Grand Est (France).
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
//                 Iordan Iordanov  <Iordan.Iordanov@loria.fr>

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_ITERATORS_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_ITERATORS_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/array.h>
#include <CGAL/iterator.h>

namespace CGAL {

template < class T >
class Periodic_4_hyperbolic_triangulation_triangle_iterator_2 {

public:

  typedef typename T::Periodic_triangle                   				value_type;
  typedef const typename T::Periodic_triangle *           				pointer;
  typedef const typename T::Periodic_triangle &           				reference;
  typedef std::size_t                                     				size_type;
  typedef std::ptrdiff_t                                  				difference_type;
  typedef std::bidirectional_iterator_tag                 				iterator_category;

  typedef typename T::Periodic_triangle                   				Periodic_triangle;
  typedef Periodic_4_hyperbolic_triangulation_triangle_iterator_2<T>	Periodic_triangle_iterator;
  typedef typename T::Facet                               				Facet;
  typedef typename T::Facet_iterator                      				Facet_iterator;

  typedef typename T::Offset                              				Offset;
  typedef typename T::Iterator_type                       				Iterator_type;

  Periodic_4_hyperbolic_triangulation_triangle_iterator_2(Iterator_type it = T::STORED)
    : _t(NULL), _it(it), _off() {}

  Periodic_4_hyperbolic_triangulation_triangle_iterator_2(const T * t, Iterator_type it = T::STORED)
    : _t(t), pos(_t->facets_begin()), _it(it), _off() {
    if (_it == T::UNIQUE || _it == T::UNIQUE_COVER_DOMAIN) {
      while (pos != _t->facets_end() && !is_canonical() )
		++pos;
    }
  }

  // used to initialize the past-the-end iterator
  Periodic_4_hyperbolic_triangulation_triangle_iterator_2(const T* t, int, Iterator_type it = T::STORED)
    : _t(t), pos(_t->facets_end()), _it(it), _off() {}

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
  Offset  _off; // current offset
  mutable Periodic_triangle periodic_triangle; // current segment.

private:
  // check whether pos points onto a unique edge or not.
  // If we are computing in 1-sheeted covering this should
  // always be true.
  bool is_canonical() {
    // fetch all offsets
    Offset off0, off1, off2;
    get_edge_offsets(off0, off1, off2);

    return (off0.is_identity() && off1.is_identity() && off2.is_identity());
  }

  // Get the canonicalized offsets of an edge.
  // This works in any cover that is encoded in _t->combine_offsets
  void get_edge_offsets(Offset &off0, Offset &off1, Offset &off2) const {
    off0 = (pos+0)->second;
    off1 = (pos+1)->second;
    off2 = (pos+2)->second;
  }


  Periodic_triangle construct_periodic_triangle() const {
    CGAL_triangulation_assertion(pos->first != typename T::Cell_handle());
    Offset off0, off1, off2;
    get_edge_offsets(off0, off1, off2);
    return make_array(
	std::make_pair(pos->first->vertex(0)->point(),off0),
	std::make_pair(pos->first->vertex(1)->point(),off1),
	std::make_pair(pos->first->vertex(2)->point(),off2));
  }
};

template < class T >
class Periodic_4_hyperbolic_triangulation_segment_iterator_2 {

public:

  typedef typename T::Periodic_segment                    value_type;
  typedef const typename T::Periodic_segment *            pointer;
  typedef const typename T::Periodic_segment &            reference;
  typedef std::size_t                                     size_type;
  typedef std::ptrdiff_t                                  difference_type;
  typedef std::bidirectional_iterator_tag                 iterator_category;

  typedef typename T::Periodic_segment                    Periodic_segment;
  typedef Periodic_4_hyperbolic_triangulation_segment_iterator_2<T>
                                                     	  Periodic_segment_iterator;
  typedef typename T::Edge                                Edge;
  typedef typename T::Edge_iterator                       Edge_iterator;

  typedef typename T::Offset                              Offset;
  typedef typename T::Iterator_type                       Iterator_type;

  Periodic_4_hyperbolic_triangulation_segment_iterator_2(Iterator_type it = T::STORED)
    : _t(NULL), _it(it), _off() {}

  Periodic_4_hyperbolic_triangulation_segment_iterator_2(const T * t,
					      Iterator_type it = T::STORED)
    : _t(t), pos(_t->edges_begin()), _it(it), _off() {
    if (_it == T::UNIQUE || _it == T::UNIQUE_COVER_DOMAIN) {
      while (pos != _t->edges_end() && !is_canonical() )
	++pos;
    }
  }

  // used to initialize the past-the-end iterator
  Periodic_4_hyperbolic_triangulation_segment_iterator_2(const T* t, int,
					      Iterator_type it = T::STORED)
    : _t(t), pos(_t->edges_end()), _it(it), _off() {}

  Periodic_segment_iterator& operator++() {
    switch (_it) {
    case T::STORED:
      ++pos;
      break;
    case T::UNIQUE:
      do { ++pos; } while (pos != _t->edges_end() && !is_canonical());
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
  Offset _off; // current offset
  mutable Periodic_segment periodic_segment; // current segment.

private:
  // check whether pos points onto a unique edge or not.
  // If we are computing in 1-sheeted covering this should
  // always be true.
  bool is_canonical() {
    // fetch all offsets
    Offset off0, off1;
    get_edge_offsets(off0, off1);

    return (off0.is_identity() && off1.is_identity());
  }

  // Get the canonicalized offsets of an edge.
  // This works in any cover that is encoded in _t->combine_offsets
  void get_edge_offsets(Offset &off0, Offset &off1) const {
    off0 = (pos+0)->second;
    off1 = (pos+1)->second;
  }


  Periodic_segment construct_periodic_segment() const {
    CGAL_triangulation_assertion(pos->first != typename T::Cell_handle());
    Offset off0, off1;
    get_edge_offsets(off0, off1);
    
    return make_array(
	std::make_pair(pos->first->vertex(0)->point(),off0),
	std::make_pair(pos->first->vertex(1)->point(),off1));
  }
};

template < class T >
class Periodic_4_hyperbolic_triangulation_point_iterator_2 {

public:
  typedef typename T::Periodic_point                      value_type;
  typedef const typename T::Periodic_point *              pointer;
  typedef const typename T::Periodic_point &              reference;
  typedef std::size_t                                     size_type;
  typedef std::ptrdiff_t                                  difference_type;
  typedef std::bidirectional_iterator_tag                 iterator_category;

  typedef typename T::Periodic_point                      Periodic_point;
  typedef Periodic_4_hyperbolic_triangulation_point_iterator_2<T>    
  														  Periodic_point_iterator;

  typedef typename T::Vertex                              Vertex;
  typedef typename T::Vertex_iterator                     Vertex_iterator;

  typedef typename T::Offset                              Offset;
  typedef typename T::Iterator_type                       Iterator_type;

  Periodic_4_hyperbolic_triangulation_point_iterator_2(Iterator_type it = T::STORED)
    : _t(NULL), _it(it) {}

  Periodic_4_hyperbolic_triangulation_point_iterator_2(const T * t, Iterator_type it = T::STORED)
    : _t(t), pos(_t->vertices_begin()), _it(it) {
    if (_it == T::UNIQUE || _it == T::UNIQUE_COVER_DOMAIN) {
      while (pos != _t->vertices_end() && !is_canonical() )
		++pos;
    }
  }

  // used to initialize the past-the-end iterator
  Periodic_4_hyperbolic_triangulation_point_iterator_2(const T* t, int, Iterator_type it = T::STORED)
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
  Offset _off; // current offset
  mutable Periodic_point periodic_point; // current point.

private:
  // check whether pos points onto a vertex inside the original
  // domain. If we are computing in 1-sheeted covering this should
  // always be true.
  bool is_canonical() {
    return (pos.second.is_identity());
  }

  Periodic_point construct_periodic_point() const {
    CGAL_triangulation_assertion(pos != typename T::Vertex_handle());
    Offset off = pos->second;
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
    return (!t->v.offset().is_identity());
  }
};

// Iterates over the vertices in a periodic triangulation that are
// located inside the original domain.
// Derives from Filter_iterator in order to add a conversion to handle
//
template <class T>
class Periodic_4_hyperbolic_triangulation_unique_vertex_iterator_2
  : public Filter_iterator<typename T::Vertex_iterator, Domain_tester<T> > {

  typedef typename T::Vertex_handle 	Vertex_handle;
  typedef typename T::Vertex_iterator 	Vertex_iterator;

  typedef Filter_iterator<Vertex_iterator, Domain_tester<T> > 			Base;
  typedef Periodic_4_hyperbolic_triangulation_unique_vertex_iterator_2 	Self;
public:

  Periodic_4_hyperbolic_triangulation_unique_vertex_iterator_2() : Base() {}
  Periodic_4_hyperbolic_triangulation_unique_vertex_iterator_2(const Base &b) : Base(b) {}

  Self & operator++() { Base::operator++(); return *this; }
  Self & operator--() { Base::operator--(); return *this; }
  Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
  Self operator--(int) { Self tmp(*this); --(*this); return tmp; }

  operator Vertex_handle() const { return Base::base(); }
};

}  // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_ITERATORS_2_H