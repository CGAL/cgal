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

#ifndef CGAL_P2T2_TRIANGULATION_CIRCULATORS_ON_LATTICE_2_H
#define CGAL_P2T2_TRIANGULATION_CIRCULATORS_ON_LATTICE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>

#include <CGAL/iterator.h>

#include <cstddef>

namespace CGAL {

// All these circulators are basically TDS circulators but f->neighbor(i) is replaced by tr->neighbor(f, i)
// because while in multiple cover, a DT2 is used and the neighbor is not trivial to compute.
// See also Tr::neighbor(Face_handle, int)

template <class Tr>
class Triangulation_ds_face_circulator_on_lattice_2
  : public Bidirectional_circulator_base<typename Tr::Triangulation_data_structure::Face,
                                         std::ptrdiff_t,
                                         std::size_t>,
    public Triangulation_cw_ccw_2
{
private:
  typedef
  Bidirectional_circulator_base<typename Tr::Triangulation_data_structure::Face,
                                std::ptrdiff_t,
                                std::size_t>                  Base_circulator;
  typedef typename Tr::Triangulation_data_structure           Tds;
  typedef Triangulation_ds_face_circulator_on_lattice_2<Tds>  Self;

public:
  typedef typename Tds::Face                                  Face;
  typedef typename Tds::Vertex                                Vertex;
  typedef typename Tds::Face_handle                           Face_handle;
  typedef typename Tds::Vertex_handle                         Vertex_handle;

private:
  Vertex_handle _v;
  Face_handle _pos;

  const Tr* _tr;

public:
  Triangulation_ds_face_circulator_on_lattice_2() : _v(), _pos(), _tr(nullptr) { }

  Triangulation_ds_face_circulator_on_lattice_2(Vertex_handle v,
                                                Face_handle f,
                                                const Tr* tr)
    : _v(v), _pos(f), _tr(tr)
  {
    if(_v == Vertex_handle())
      _pos = Face_handle();
    else if(_pos == Face_handle())
      _pos = v->face();

    CGAL_assertion(_tr->is_canonical(_pos));

    if(_pos == Face_handle())
    {
      _v = Vertex_handle();
      _pos = Face_handle();
    }
    else
    {
      CGAL_triangulation_precondition(_pos->has_vertex(v));
    }
  }

  Self& operator++()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());

    const int i = _pos->index(_v);
    _pos = _tr->neighbor(_pos, ccw(i));

    return *this;
  }

  Self operator++(int)
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());
    Self tmp(*this);
    ++(*this);
    return tmp;
  }

  Self& operator--()
  {
     CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());
     int i = _pos->index(_v);
     _pos = _tr->neighbor(_pos, cw(i));
     return *this;
  }

  Self operator--(int)
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());
    Self tmp(*this);
    --(*this);
    return tmp;
  }

  bool operator==(const Self &fc) const
  {
    return (_v == fc._v) && (_pos == fc._pos);
  }

  bool operator!=(const Self &fc) const
  {
  return !(*this == fc);
  }

  bool operator==(const Face_handle f) const { return _pos == f; }
  bool operator!=(const Face_handle f) const { return _pos != f; }

  bool is_empty() const { return (_tr == nullptr || _v == Vertex_handle() || _pos == Face_handle()); }
  bool operator==(std::nullptr_t CGAL_triangulation_assertion_code(n)) const
  {
    CGAL_triangulation_assertion(n == nullptr);
    return (_tr == nullptr || _v == Vertex_handle() || _pos == Face_handle());
  }

  bool operator!=(std::nullptr_t CGAL_triangulation_assertion_code(n)) const
  {
    CGAL_triangulation_assertion(n == nullptr);
    return !(*this == nullptr);
  }

  Face& operator*() const
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());
    return *_pos;
  }

  Face* operator->() const
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());
    return &*_pos;
  }

  Face_handle base() const { return _pos; }
  operator Face_handle() const { return _pos; }
};

template <class Tr_>
bool operator==(typename Tr_::Triangulation_data_structure::Face_handle f,
                Triangulation_ds_face_circulator_on_lattice_2<Tr_> fc)
{
  return (fc == f);
}

template <class Tr_>
bool operator!=(typename Tr_::Triangulation_data_structure::Face_handle f,
                Triangulation_ds_face_circulator_on_lattice_2<Tr_> fc)
{
  return (fc != f);
}

template <class Tr>
class Triangulation_ds_vertex_circulator_on_lattice_2
  : public Bidirectional_circulator_base<typename Tr::Triangulation_data_structure::Vertex,
                                         std::ptrdiff_t,
                                         std::size_t>,
    public Triangulation_cw_ccw_2
{
  typedef Triangulation_ds_vertex_circulator_on_lattice_2<Tr> Self;
  typedef typename Tr::Triangulation_data_structure           Tds;

public:
  typedef typename Tds::Face                                  Face;
  typedef typename Tds::Vertex                                Vertex;
  typedef typename Tds::Face_handle                           Face_handle;
  typedef typename Tds::Vertex_handle                         Vertex_handle;

private:
  Vertex_handle _v;
  Face_handle _pos;
  int _ri;

  const Tr* _tr;

public:
  Triangulation_ds_vertex_circulator_on_lattice_2() : _v(), _pos(), _tr(nullptr) { }

  Triangulation_ds_vertex_circulator_on_lattice_2(Vertex_handle v,
                                                  Face_handle f,
                                                  const Tr* tr)
    : _v( v ), _pos(f), _tr(tr)
  {
    if(_v == Vertex_handle())
     _pos = Face_handle();
    else if(_pos == Face_handle())
      _pos = v->face();

    if(_pos == Face_handle())
    {
      _v = Vertex_handle();
      _pos = Face_handle();
      return;
    }

    CGAL_assertion(_tr->is_canonical(_pos));

    _ri = ccw(_pos->index(_v));
  }

  Self& operator++()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());

    _pos = _tr->neighbor(_pos, ccw(_pos->index(_v)));
    _ri = ccw(_pos->index(_v));

    return *this;
  }

  Self operator++(int)
  {
    Self tmp(*this);
    ++(*this);
    return tmp;
  }

  Self& operator--()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());

    _pos = _pos->neighbor(cw(_pos->index(_v)));
    _ri = ccw(_pos->index(_v));

    return *this;
  }

  Self operator--(int)
  {
    Self tmp(*this);
    --(*this);
    return tmp;
  }

  bool operator==(const Self& vc) const { (_tr == vc._tr) && (_v == vc._v) && (_ri == vc._ri) && (_pos == vc._pos); }
  bool operator!=(const Self& vc) const { !(*this == vc); }

  bool operator==(const Vertex_handle vh) const { return _pos->vertex(_ri) == vh; }
  bool operator!=(const Vertex_handle vh) const { return _pos->vertex(_ri) != vh; }

  bool is_empty() const { return (_tr == nullptr || _v == Vertex_handle() || _pos == Face_handle()); }

  bool operator==(std::nullptr_t CGAL_triangulation_assertion_code(n)) const
  {
    CGAL_triangulation_assertion(n == nullptr);
    return (_tr == nullptr || _v == Vertex_handle() || _pos == Face_handle());
  }

  bool operator!=(std::nullptr_t CGAL_triangulation_assertion_code(n)) const
  {
    CGAL_triangulation_assertion(n == nullptr);
    return !(*this == nullptr);
  }

  Vertex& operator*() const
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());
    return *(_pos->vertex(_ri));
  }

  Vertex* operator->() const
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());
    return &*(_pos->vertex(_ri));
  }

  Vertex_handle base() const { return _pos->vertex(_ri) ;}
  operator Vertex_handle() const { return _pos->vertex(_ri) ;}
};

template <class Tr_>
inline bool operator==(typename Tr_::Triangulation_data_structure::Vertex_handle vh,
                       Triangulation_ds_vertex_circulator_on_lattice_2<Tr_> vc)
{
  return (vc == vh);
}

template <class Tr_>
inline bool operator!=(typename Tr_::Triangulation_data_structure::Vertex_handle vh,
                       Triangulation_ds_vertex_circulator_on_lattice_2<Tr_> vc)
{
  return !(vc == vh);
}

template <class Tr>
class Triangulation_ds_edge_circulator_on_lattice_2
  : public Bidirectional_circulator_base<typename Tr::Triangulation_data_structure::Edge,
                                         std::ptrdiff_t,
                                         std::size_t>,
    public Triangulation_cw_ccw_2
{
  typedef Triangulation_ds_edge_circulator_on_lattice_2<Tr>   Self;
  typedef typename Tr::Triangulation_data_structure           Tds;

public:
  typedef typename Tds::Face                                  Face;
  typedef typename Tds::Vertex                                Vertex;
  typedef typename Tds::Edge                                  Edge;
  typedef typename Tds::Face_handle                           Face_handle;
  typedef typename Tds::Vertex_handle                         Vertex_handle;

private:
  int _ri;
  Vertex_handle _v;
  Face_handle _pos;
  mutable Edge _edge;

  const Tr* _tr;

public:
  Triangulation_ds_edge_circulator_on_lattice_2() : _ri(0), _v(), _pos(), _tr(nullptr) { }

  Triangulation_ds_edge_circulator_on_lattice_2(Vertex_handle v,
                                                Face_handle f,
                                                const Tr* tr)
    : _v(v), _pos(f), _tr(tr)
  {
    if(_v == Vertex_handle())
      _pos = Face_handle();
    else if(_pos == Face_handle())
      _pos = v->face();

    if(_pos == Face_handle())
    {
      _ri = 0;
      _v = Vertex_handle();
      _pos = Face_handle();
      return;
    }

    CGAL_assertion(_tr->is_canonical(_pos));

    _ri = ccw(_pos->index(_v));
  }

  Self& operator++()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());

    _pos = _pos->neighbor(ccw(_pos->index(_v)));
    _ri = ccw(_pos->index(_v));

    return *this;
  }

  Self operator++(int)
  {
    Self tmp(*this);
    ++(*this);
    return tmp;
  }

  Self& operator--()
  {
    CGAL_triangulation_precondition(_tr != nullptr && _pos != Face_handle() && _v != Vertex_handle());
    _pos = _pos->neighbor(cw(_pos->index(_v)));
    _ri = ccw(_pos->index(_v));

    return *this;
  }

  Self operator--(int)
  {
    Self tmp(*this);
    --(*this);
    return tmp;
  }

  bool operator==(const Self& vc) const { return (_v == vc._v) &&  (_ri == vc._ri) && (_pos == vc._pos); }
  bool operator!=(const Self& vc) const { return ! (*this == vc); }
  bool is_empty() const { return (_tr == nullptr || _v == Vertex_handle() || _pos == Face_handle()); }
  bool operator==(std::nullptr_t CGAL_triangulation_assertion_code(n)) const
  {
    CGAL_triangulation_assertion(n == nullptr);
    return (_tr == nullptr || _v == Vertex_handle() || _pos == Face_handle());
  }

  bool operator!=(std::nullptr_t CGAL_triangulation_assertion_code(n)) const
  {
    CGAL_triangulation_assertion(n == nullptr);
    return !(*this == nullptr);
  }

  Edge* operator->() const
  {
    _edge.first = _pos;
    _edge.second = _ri;
    return &_edge;
  }

  Edge& operator*() const
  {
    _edge.first = _pos;
    _edge.second = _ri;
    return _edge;
  }
};

} // namespace CGAL

#endif // CGAL_P2T2_TRIANGULATION_CIRCULATORS_ON_LATTICE_2_H
