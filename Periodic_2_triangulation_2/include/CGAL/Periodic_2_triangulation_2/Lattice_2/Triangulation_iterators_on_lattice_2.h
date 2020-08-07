// Copyright (c) 2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//                 Georg Osang

#ifndef CGAL_P2T2_TRIANGULATION_ITERATORS_ON_LATTICE_2_H
#define CGAL_P2T2_TRIANGULATION_ITERATORS_ON_LATTICE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>

#include <CGAL/iterator.h>

namespace CGAL {

// All these iterators are basically TDS iterators but f->neighbor(i) is replaced by tr->neighbor(f, i)
// because while in multiple cover, a DT2 is used and the neighbor is not trivial to compute.
// See also Tr::neighbor(Face_handle, int)

template <class Tr>
class Triangulation_ds_vertex_iterator_on_lattice_2
{
  typedef typename Tr::Triangulation_data_structure             Tds;
  typedef Triangulation_ds_vertex_iterator_on_lattice_2<Tr>     Self;

public:
  typedef typename Tds::Vertex                                  Vertex;
  typedef typename Tds::Vertex_handle                           Vertex_handle;
  typedef typename Tds::Vertex_iterator                         Base_vertex_iterator;

  typedef Vertex                                                value_type;
  typedef Vertex*                                               pointer;
  typedef Vertex&                                               reference;
  typedef std::size_t                                           size_type;
  typedef std::ptrdiff_t                                        difference_type;
  typedef std::bidirectional_iterator_tag                       iterator_category;

private:
  Base_vertex_iterator _pos;

  const Tr* _tr;

public:
  Triangulation_ds_vertex_iterator_on_lattice_2() : _pos(nullptr), _tr(nullptr) { }

  Triangulation_ds_vertex_iterator_on_lattice_2(const Tr* tr)
    : _tr(tr)
  {
    if(_tr->dimension() == -2)
    {
      _pos = _tr->tds().vertices().end(); // there is no vertex
      return;
    }

    _pos = _tr->tds().vertices().begin();
    while(_pos != _tr->tds().vertices().end())
      increment();
  }

  Triangulation_ds_vertex_iterator_on_lattice_2(const Tr* tr, int) // used to initialize an iterator at 'end()'
    : _tr(tr)
  {
    _pos = _tr->tds().vertices().end();
  }

  bool operator==(const Self& fi) const
  {
    return _tr == fi._tr && _pos == fi._pos;
  }

  bool operator!=(const Self& fi) const
  {
    return !(*this == fi);
  }

  Self& operator++()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Vertex_handle() && _pos != _tr->tds().vertices().end());

    do
      increment();
    while(_pos != _tr->tds().vertices().end());

    return *this;
  }

  Self& operator--()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Vertex_handle() && *this != Self(_tr));

    do
      decrement();
    while(*this != Self(_tr));

    return *this;
  }

  Self operator++(int)
  {
    Self tmp(*this);
    ++(*this);
    return tmp;
  }

  Self operator--(int)
  {
    Self tmp(*this);
    --(*this);
    return tmp;
  }

  Vertex* operator->() const
  {
    return _pos;
  }

  Vertex& operator*() const
  {
    return *_pos;
  }

private:
  void increment()
  {
    CGAL_triangulation_precondition(_tr->tds().dimension() == 2);

    do {
      ++_pos;
    } while(_pos != _tr->tds().vertices().end() && !_tr->is_canonical(_pos));
  }

  void decrement()
  {
    CGAL_triangulation_precondition(_tr->tds().dimension() == 2);

    do {
      --_pos;
    } while(*this != Self(_tr) && !_tr->is_canonical(_pos));
  }
};

template <class Tr>
class Triangulation_ds_face_iterator_on_lattice_2
{
  typedef typename Tr::Triangulation_data_structure             Tds;
  typedef Triangulation_ds_face_iterator_on_lattice_2<Tr>       Self;

public:
  typedef typename Tds::Face                                    Face;
  typedef typename Tds::Face_handle                             Face_handle;
  typedef typename Tds::Face_iterator                           Face_iterator;

  typedef Face                                                  value_type;
  typedef Face*                                                 pointer;
  typedef Face&                                                 reference;
  typedef std::size_t                                           size_type;
  typedef std::ptrdiff_t                                        difference_type;
  typedef std::bidirectional_iterator_tag                       iterator_category;

private:
  Face_iterator _pos;

  const Tr* _tr;

public:
  Triangulation_ds_face_iterator_on_lattice_2() : _pos(nullptr), _tr(nullptr) { }

  Triangulation_ds_face_iterator_on_lattice_2(const Tr* tr)
    : _tr(tr)
  {
    if(_tr->dimension() == -2)
    {
      _pos = _tr->tds().faces().end(); // there is no face
      return;
    }

    _pos = _tr->tds().faces().begin();
    while(_pos != _tr->tds().faces().end())
      increment();
  }

  Triangulation_ds_face_iterator_on_lattice_2(const Tr* tr, int) // used to initialize an iterator at 'end()'
    : _tr(tr)
  {
    _pos = _tr->tds().faces().end();
  }

  bool operator==(const Self& fi) const
  {
    return _tr == fi._tr && _pos == fi._pos;
  }

  bool operator!=(const Self& fi) const
  {
    return !(*this == fi);
  }

  Self& operator++()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _pos != _tr->tds().faces().end());

    do
      increment();
    while(_pos != _tr->tds().faces().end());

    return *this;
  }

  Self& operator--()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && *this != Self(_tr));

    do
      decrement();
    while(*this != Self(_tr));

    return *this;
  }

  Self operator++(int)
  {
    Self tmp(*this);
    ++(*this);
    return tmp;
  }

  Self operator--(int)
  {
    Self tmp(*this);
    --(*this);
    return tmp;
  }

  Face* operator->() const
  {
    return _pos;
  }

  Face& operator*() const
  {
    return *_pos;
  }

private:
  void increment()
  {
    CGAL_triangulation_precondition(_tr->tds().dimension() == 2);

    do {
      ++_pos;
    } while(_pos != _tr->tds().faces().end() && !_tr->is_canonical(_pos));
  }

  void decrement()
  {
    CGAL_triangulation_precondition(_tr->tds().dimension() == 2);

    do {
      --_pos;
    } while(*this != Self(_tr) && !_tr->is_canonical(_pos));
  }
};

template <class Tr>
class Triangulation_ds_edge_iterator_on_lattice_2
{
  typedef typename Tr::Triangulation_data_structure             Tds;
  typedef Triangulation_ds_edge_iterator_on_lattice_2<Tr>       Self;

public:
  typedef typename Tds::Edge                                    Edge;
  typedef typename Tds::Face_iterator                           Face_iterator;
  typedef typename Tds::Face_handle                             Face_handle;

  typedef Edge                                                  value_type;
  typedef Edge*                                                 pointer;
  typedef Edge&                                                 reference;
  typedef std::size_t                                           size_type;
  typedef std::ptrdiff_t                                        difference_type;
  typedef std::bidirectional_iterator_tag                       iterator_category;

private:
  Face_iterator _pos;
  mutable Edge _edge;

  const Tr* _tr;

public:
  Triangulation_ds_edge_iterator_on_lattice_2() : _pos(nullptr), _edge(), _tr(nullptr) { }

  Triangulation_ds_edge_iterator_on_lattice_2(const Tr* tr)
    : _tr(tr)
  {
    _edge.second = 0;
    if(_tr->dimension() <= 0)
    {
      _pos = _tr->tds().faces().end(); // there is no edge
      return;
    }

    _pos = _tr->tds().faces().begin();
    while(_pos != _tr->tds().faces().end() && !_tr->is_canonical(_pos) && !associated_edge())
      increment();
  }

  Triangulation_ds_edge_iterator_on_lattice_2(const Tr* tr, int) // used to initialize an iterator at 'end()'
    : _tr(tr)
  {
    _pos = _tr->tds().faces().end();
    _edge.second = 0;
  }

  bool operator==(const Self& fi) const
  {
    return _tr == fi._tr && _pos == fi._pos && _edge.second == fi._edge.second;
  }

  bool operator!=(const Self& fi) const
  {
    return !(*this == fi);
  }

  Self& operator++()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _pos != _tr->tds().faces().end());

    do
      increment();
    while(!associated_edge() && _pos != _tr->tds().faces().end());

    return *this;
  }

  Self& operator--()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && *this != Self(_tr));

    do
      decrement();
    while(!associated_edge() && *this != Self(_tr));

    return *this;
  }

  Self operator++(int)
  {
    Self tmp(*this);
    ++(*this);
    return tmp;
  }

  Self operator--(int)
  {
    Self tmp(*this);
    --(*this);
    return tmp;
  }

  Edge* operator->() const
  {
    _edge.first = _pos;
    return &_edge;
  }

  Edge& operator*() const
  {
    _edge.first = _pos;
    return _edge;
  }

private:
  void increment()
  {
    CGAL_triangulation_precondition(_tr->tds().dimension() == 2);

    if(_edge.second == 2)
    {
      do {
        ++_pos;
      } while(_pos != _tr->tds().faces().end() && !_tr->is_canonical(_pos));

      _edge.second = 0;
    }
    else
    {
      _edge.second += 1;
    }
  }

  void decrement()
  {
    CGAL_triangulation_precondition(_tr->tds().dimension() == 2);

    if(_edge.second == 0)
    {
      do {
        --_pos;
      } while(!_tr->is_canonical(_pos) && *this != Self(_tr));

      _edge.second = 2;
    }
    else
    {
      _edge.second -= 1;
    }
  }

  // Because (f, i) is really a halfedge, not an edge, so ignore half of them
  bool associated_edge()
  {
    return Face_handle(_pos) < _pos->neighbor(_edge.second);
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Tr>
struct Construct_periodic_point
{
  typedef typename Tr::Vertex                                   argument_type;
  typedef typename Tr::Periodic_point                           result_type;

  result_type operator()(const argument_type& v) const { return result_type(v.point(), v.offset()); }
};

template <class Tr>
struct Construct_periodic_segment
{
  typedef typename Tr::Vertex_handle                            Vertex_handle;
  typedef typename Tr::Periodic_point                           Periodic_point;

  typedef typename Tr::Edge                                     argument_type;
  typedef typename Tr::Periodic_segment                         result_type;

  result_type operator()(const argument_type& e) const
  {
    const Vertex_handle v1 = e.first->vertex(Triangulation_cw_ccw_2::cw(e.second));
    const Vertex_handle v2 = e.first->vertex(Triangulation_cw_ccw_2::ccw(e.second));

    return result_type(Periodic_point(v1->point(), v1->offset()),
                       Periodic_point(v2->point(), v2->offset()));
  }
};

template <class Tr>
struct Construct_periodic_triangle
{
  typedef typename Tr::Vertex_handle                            Vertex_handle;
  typedef typename Tr::Periodic_point                           Periodic_point;

  typedef typename Tr::Face                                     argument_type;
  typedef typename Tr::Periodic_triangle                        result_type;

  result_type operator()(const argument_type& f) const
  {
    const Vertex_handle v0 = f->vertex(0);
    const Vertex_handle v1 = f->vertex(1);
    const Vertex_handle v2 = f->vertex(2);

    return result_type(Periodic_point(v0->point(), v0->offset()),
                       Periodic_point(v1->point(), v1->offset()),
                       Periodic_point(v2->point(), v2->offset()));
  }
};

} // namespace CGAL

#endif // CGAL_P2T2_TRIANGULATION_ITERATORS_ON_LATTICE_2_H
