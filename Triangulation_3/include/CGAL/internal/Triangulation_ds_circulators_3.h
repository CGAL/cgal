// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

#ifndef CGAL_INTERNAL_TRIANGULATION_DS_CIRCULATORS_3_H
#define CGAL_INTERNAL_TRIANGULATION_DS_CIRCULATORS_3_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/circulator_bases.h>

namespace CGAL { namespace internal {

template < class Tds_ >
class Triangulation_ds_cell_circulator_3
 : public Bidirectional_circulator_base<typename Tds_::Cell,
                                        std::ptrdiff_t, std::size_t>
{
  // circulates on cells around a given edge

  typedef Tds_                         Tds;
  typedef typename Tds::Cell           Cell;
  typedef typename Tds::Edge           Edge;
  typedef typename Tds::Vertex_handle  Vertex_handle;
  typedef typename Tds::Cell_handle    Cell_handle;

  typedef Triangulation_ds_cell_circulator_3<Tds> Cell_circulator;

public:

  Triangulation_ds_cell_circulator_3()
    : _s(), _t(), pos()
  {}

  Triangulation_ds_cell_circulator_3(Cell_handle c, int s, int t)
    : _s(c->vertex(s)), _t(c->vertex(t)), pos(c)
  {
    CGAL_triangulation_precondition( c != Cell_handle() &&
				     s >= 0 && s < 4 &&
				     t >= 0 && t < 4 );
  }

  Triangulation_ds_cell_circulator_3(const Edge & e)
    : _s(e.first->vertex(e.second)), _t(e.first->vertex(e.third)), pos(e.first)
  {
    CGAL_triangulation_precondition( e.first != Cell_handle() &&
				     e.second >=0 && e.second < 4 &&
				     e.third  >=0 && e.third  < 4);
  }

  Triangulation_ds_cell_circulator_3(Cell_handle c, int s, int t,
	                             Cell_handle start)
    : _s(c->vertex(s)), _t(c->vertex(t)), pos(start)
  {
    CGAL_triangulation_precondition( c != Cell_handle() &&
				     s >= 0 && s < 4 &&
				     t >= 0 && t < 4 &&
                                     start->has_vertex( _s ) &&
	                             start->has_vertex( _t ) );
  }

  Triangulation_ds_cell_circulator_3(const Edge & e, Cell_handle start)
    : _s(e.first->vertex(e.second)), _t(e.first->vertex(e.third)), pos(start)
  {
    CGAL_triangulation_precondition( e.first != Cell_handle() &&
				     e.second >=0 && e.second < 4 &&
				     e.third  >=0 && e.third  < 4 &&
                                     start->has_vertex( _s ) &&
	                             start->has_vertex( _t ) );
  }

  Cell_circulator & operator++()
  {
    CGAL_triangulation_precondition( pos != Cell_handle() );
    //then dimension() cannot be < 3

    pos = pos->neighbor(next_around_edge(pos->index(_s), pos->index(_t)));
    return *this;
  }

  Cell_circulator operator++(int)
  {
    Cell_circulator tmp(*this);
    ++(*this);
    return tmp;
  }

  Cell_circulator & operator--()
  {
    CGAL_triangulation_precondition( pos != Cell_handle() );

    pos = pos->neighbor(next_around_edge(pos->index(_t), pos->index(_s)));
    return *this;
  }

  Cell_circulator operator--(int)
  {
    Cell_circulator tmp(*this);
    --(*this);
    return tmp;
  }

  Cell& operator*() const
  {
    return *pos;
  }

  Cell* operator->() const
  {
    return &*pos;
  }

  bool operator==(const Cell_circulator & ccir) const
  {
    return pos == ccir.pos && _s == ccir._s && _t == ccir._t;
  }

  bool operator!=(const Cell_circulator & ccir) const
  {
    return ! (*this == ccir);
  }

  bool operator==(Cell_handle ch) const
  {
    return ch == pos;
  }

  bool operator!=(Cell_handle ch) const
  {
    return ch != pos;
  }

  bool operator==(Nullptr_t CGAL_triangulation_assertion_code(n)) const
  {
    CGAL_triangulation_assertion( n == NULL);
    return pos == Cell_handle();
  }

  bool operator!=(Nullptr_t n) const
  {
    return ! (*this == n);
  }

  // For TDS's private use only.
  Cell_handle base() const { return pos; }
  operator Cell_handle() const { return pos; }

private:
  Vertex_handle _s;    // source vertex of the edge
  Vertex_handle _t;    // target vertex of the edge
  Cell_handle pos;     // current cell

  static int next_around_edge(const int i, const int j)
  {
      return Triangulation_utils_3::next_around_edge(i,j);
  }
};

template < class Tds_ >
inline
bool
operator==(typename Tds_::Cell_handle ch, Triangulation_ds_cell_circulator_3<Tds_> cc)
{
  return (cc==ch);
}

template < class Tds_ >
inline
bool
operator!=(typename Tds_::Cell_handle ch, Triangulation_ds_cell_circulator_3<Tds_> cc)
{
  return !(cc==ch);
}

template < class Tds_ >
class Triangulation_ds_facet_circulator_3
  : public Bidirectional_circulator_base<typename Tds_::Facet,
                                         std::ptrdiff_t, std::size_t>
{
  // circulates on facets around a given edge

  typedef Tds_                         Tds;
  typedef typename Tds::Cell           Cell;
  typedef typename Tds::Cell_handle    Cell_handle;
  typedef typename Tds::Facet          Facet;
  typedef typename Tds::Edge           Edge;
  typedef typename Tds::Vertex_handle  Vertex_handle;

  typedef Triangulation_ds_facet_circulator_3<Tds> Facet_circulator;

public:

  Triangulation_ds_facet_circulator_3()
    : _s(), _t(), pos()
  {}

  Triangulation_ds_facet_circulator_3(Cell_handle c, int s, int t)
    : _s(c->vertex(s)), _t(c->vertex(t)), pos(c)
  {
    CGAL_triangulation_precondition( c != Cell_handle() &&
				     s >= 0 && s < 4 &&
				     t >= 0 && t < 4 );
  }

  Triangulation_ds_facet_circulator_3(const Edge & e)
    : _s(e.first->vertex(e.second)), _t(e.first->vertex(e.third)), pos(e.first)
  {
    CGAL_triangulation_precondition( e.first != Cell_handle() &&
				     e.second >= 0 && e.second < 4 &&
				     e.third  >= 0 && e.third  < 4);
  }

  Triangulation_ds_facet_circulator_3(Cell_handle c, int s, int t,
	                              Cell_handle start, int f)
    : _s(c->vertex(s)), _t(c->vertex(t))
  {
    CGAL_triangulation_precondition( c != Cell_handle() &&
				     s >= 0 && s < 4 &&
				     t >= 0 && t < 4 &&
				     f >= 0 && f < 4 &&
                                     start->has_vertex( _s ) &&
	                             start->has_vertex( _t ) );

    int i = start->index( _s );
    int j = start->index( _t );

    CGAL_triangulation_precondition( f!=i && f!=j );

    if ( f == next_around_edge(i,j) )
	pos = start;
    else
	pos = start->neighbor(f); // other cell with same facet
  }

  Triangulation_ds_facet_circulator_3(Cell_handle c, int s, int t,
	                              const Facet & start)
    : _s(c->vertex(s)), _t(c->vertex(t))
  {
    CGAL_triangulation_precondition( c != Cell_handle() &&
				     s >= 0 && s < 4 &&
				     t >= 0 && t < 4 &&
                                     start.first->has_vertex( _s ) &&
	                             start.first->has_vertex( _t ) );

    int i = start.first->index( _s );
    int j = start.first->index( _t );

    CGAL_triangulation_precondition( start.second !=i && start.second !=j );

    if ( start.second == next_around_edge(i,j) )
	pos = start.first;
    else
      pos = start.first->neighbor(start.second); // other cell with same facet
  }

  Triangulation_ds_facet_circulator_3(const Edge & e, Cell_handle start, int f)
    : _s(e.first->vertex(e.second)), _t(e.first->vertex(e.third))
  {
    CGAL_triangulation_precondition( e.first != Cell_handle() &&
				     e.second >= 0 && e.second < 4 &&
				     e.third  >= 0 && e.third  < 4 &&
				     f >= 0 && f < 4 &&
                                     start->has_vertex( _s ) &&
	                             start->has_vertex( _t ) );

    int i = start->index( _s );
    int j = start->index( _t );

    CGAL_triangulation_precondition( f!=i && f!=j );

    if ( f == next_around_edge(i,j) )
	pos = start;
    else
	pos = start->neighbor(f); // other cell with same facet
  }

  Triangulation_ds_facet_circulator_3(const Edge & e, const Facet & start)
    : _s(e.first->vertex(e.second)), _t(e.first->vertex(e.third))
  {
    CGAL_triangulation_precondition( e.first != Cell_handle() &&
				     e.second >= 0 && e.second < 4 &&
				     e.third  >= 0 && e.third  < 4 &&
                                     start.first->has_vertex( _s ) &&
	                             start.first->has_vertex( _t ) );

    int i = start.first->index( _s );
    int j = start.first->index( _t );

    if ( start.second == next_around_edge(i,j) )
	pos = start.first;
    else
	pos = start.first->neighbor(start.second);
  }

  Facet_circulator & operator++()
  {
    CGAL_triangulation_precondition( pos != Cell_handle() );
    //then dimension() cannot be < 3

    pos = pos->neighbor( next_around_edge( pos->index(_s), pos->index(_t) ) );
    return *this;
  }

  Facet_circulator operator++(int)
  {
    Facet_circulator tmp(*this);
    ++(*this);
    return tmp;
  }

  Facet_circulator & operator--()
  {
    CGAL_triangulation_precondition( pos != Cell_handle() );

    pos = pos->neighbor( next_around_edge( pos->index(_t), pos->index(_s) ) );
    return *this;
  }

  Facet_circulator operator--(int)
  {
    Facet_circulator tmp(*this);
    --(*this);
    return tmp;
  }

  Facet operator*() const
  {
    return Facet(pos, next_around_edge( pos->index(_s), pos->index(_t) ) );
  }

  struct Proxy_Facet {
    Proxy_Facet(const Facet & ff) : f(ff) {}
    Facet f;
    const Facet* operator->() const { return &f; }
  };

  Proxy_Facet operator->() const
  {
    return Proxy_Facet(* *this);
  }

  bool operator==(const Facet_circulator & ccir) const
  {
    return pos == ccir.pos && _s == ccir._s && _t == ccir._t;
  }

  bool operator!=(const Facet_circulator & ccir) const
  {
    return ! (*this == ccir);
  }

  bool operator==(Nullptr_t CGAL_triangulation_assertion_code(c)) const
  {
    CGAL_triangulation_assertion(c == NULL);
    return pos == Cell_handle();
  }

  bool operator!=(Nullptr_t c) const
  {
    return ! (*this == c);
  }

private:
  Vertex_handle _s; // source vertex of the edge
  Vertex_handle _t; // target vertex of the edge
  Cell_handle pos;  // current cell
  // the current facet is the facet of pos numbered
  // next_around_edge( pos->index(_c->vertex(_s)),
  //                   pos->index(_c->vertex(_t)) )

  static int next_around_edge(const int i, const int j)
  {
      return Triangulation_utils_3::next_around_edge(i,j);
  }
};

template < class Tds_ >
class Triangulation_ds_face_circulator_3
 : public Bidirectional_circulator_base<typename Tds_::Cell,
                                        std::ptrdiff_t, std::size_t>
{
  // circulates on faces (Cell) around a given vertex,
  // valid in dimension 2 only.

  typedef Tds_                         Tds;
  typedef typename Tds::Cell           Cell;
  typedef typename Tds::Vertex         Vertex;
  typedef typename Tds::Vertex_handle  Vertex_handle;
  typedef typename Tds::Cell_handle    Cell_handle;

  typedef Triangulation_ds_face_circulator_3<Tds> Face_circulator;

public:

  Triangulation_ds_face_circulator_3()
    : _s(), pos() {}

  Triangulation_ds_face_circulator_3(Vertex_handle v, Cell_handle c)
    : _s(v), pos(c) {}

  Face_circulator & operator++()
  {
    CGAL_triangulation_precondition( pos != Cell_handle() );
    //then dimension() cannot be < 3

    pos = pos->neighbor(ccw(pos->index(_s)));
    return *this;
  }

  Face_circulator operator++(int)
  {
    Face_circulator tmp(*this);
    ++(*this);
    return tmp;
  }

  Face_circulator & operator--()
  {
    CGAL_triangulation_precondition( pos != Cell_handle() );

    pos = pos->neighbor(cw(pos->index(_s)));
    return *this;
  }

  Face_circulator operator--(int)
  {
    Face_circulator tmp(*this);
    --(*this);
    return tmp;
  }

  Cell& operator*() const
  {
    return *pos;
  }

  Cell* operator->() const
  {
    return &*pos;
  }

  bool operator==(const Face_circulator & ccir) const
  {
    return pos == ccir.pos && _s == ccir._s;
  }

  bool operator!=(const Face_circulator & ccir) const
  {
    return ! (*this == ccir);
  }

  bool operator==(Nullptr_t CGAL_triangulation_assertion_code(c)) const
  {
    CGAL_triangulation_assertion(c == NULL);
    return pos == Cell_handle();
  }

  bool operator!=(Nullptr_t c) const
  {
    return ! (*this == c);
  }

  // For TDS's private use only.
  Cell_handle base() const { return pos; }
  operator Cell_handle() const { return pos; }

private:
  Vertex_handle _s;    // source vertex
  Cell_handle pos;     // current cell

  static int cw(int i)
  {
      return Triangulation_utils_3::cw(i);
  }
  static int ccw(int i)
  {
      return Triangulation_utils_3::ccw(i);
  }
};

}} // namespace CGAL::internal

#endif // CGAL_INTERNAL_TRIANGULATION_DS_CIRCULATORS_3_H
