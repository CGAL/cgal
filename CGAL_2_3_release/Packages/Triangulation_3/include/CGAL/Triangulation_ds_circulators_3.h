// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_ds_circulators_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_CIRCULATORS_3_H
#define CGAL_TRIANGULATION_DS_CIRCULATORS_3_H

#include <utility>

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/circulator.h>
#include <CGAL/Triangulation_utils_3.h>

CGAL_BEGIN_NAMESPACE

template < class Tds_ >
class Triangulation_ds_cell_circulator_3
  : public Bidirectional_circulator_base<typename Tds_::Cell,
    ptrdiff_t, size_t>,
    public Triangulation_utils_3
{
  // circulates on cells around a given edge
public:

  typedef Tds_                  Tds;
  typedef typename Tds::Cell    Cell;
  typedef typename Tds::Edge    Edge;
  typedef typename Tds::Vertex  Vertex;

  typedef Triangulation_ds_cell_circulator_3<Tds> Cell_circulator;

  Triangulation_ds_cell_circulator_3()
    : _tds(NULL), _c(NULL), _s(0), _t(0), pos(NULL)
  {}

  Triangulation_ds_cell_circulator_3(const Tds * tds, Cell* c, int s, int t)
    : _tds(const_cast<Tds *>(tds)), _c(c), _s(s), _t(t)
  {
    CGAL_triangulation_precondition( c != NULL &&
				     s >= 0 && s < 4 &&
				     t >= 0 && t < 4 );
    //     if ( _tds->dimension() <3 ) useless since precondition in tds
    //     incident_cells
    pos = c;
  }

  Triangulation_ds_cell_circulator_3(const Tds * tds, const Edge & e)
    : _tds(const_cast<Tds *>(tds)), _c(e.first), _s(e.second), _t(e.third)
  {
    CGAL_triangulation_precondition( e.first != NULL &&
				     e.second >=0 && e.second < 4 &&
				     e.third  >=0 && e.third  < 4);
    //     if ( _tds->dimension() <3 ) useless since precondition in tds
    //     incident_cells
      pos = e.first;
  }

  Triangulation_ds_cell_circulator_3(const Tds * tds, Cell* c, int s, int t,
				     Cell* start)
    : _tds(const_cast<Tds *>(tds)), _c(c), _s(s), _t(t)
  {
    CGAL_triangulation_precondition( c != NULL &&
				     s >= 0 && s < 4 &&
				     t >= 0 && t < 4 &&
                                     start->has_vertex( c->vertex(s) ) &&
	                             start->has_vertex( c->vertex(t) ) );
    pos = start;
  }

  Triangulation_ds_cell_circulator_3(const Tds * tds, const Edge & e,
	                             Cell* start)
    : _tds(const_cast<Tds *>(tds)), _c(e.first), _s(e.second), _t(e.third)
  {
    CGAL_triangulation_precondition( e.first != NULL &&
				     e.second >=0 && e.second < 4 &&
				     e.third  >=0 && e.third  < 4 &&
                            start->has_vertex( e.first->vertex(e.second) ) &&
	                    start->has_vertex( e.first->vertex(e.third) ) );
    pos = start;
  }

  Triangulation_ds_cell_circulator_3(const Cell_circulator & ccir)
    : _tds(ccir._tds), _c(ccir._c), _s(ccir._s), _t(ccir._t), pos(ccir.pos)
  {}

  Cell_circulator & operator++()
  {
    CGAL_triangulation_precondition( pos != NULL );
    //then dimension() cannot be < 3

    pos = pos->neighbor( next_around_edge( pos->index(_c->vertex(_s)),
					   pos->index(_c->vertex(_t)) ) );
    return *this;
  }

  Cell_circulator operator++(int)
  {
    CGAL_triangulation_precondition( pos != NULL );
    Cell_circulator tmp(*this);
    ++(*this);
    return tmp;
  }

  Cell_circulator & operator--()
  {
    CGAL_triangulation_precondition( pos != NULL );

    pos = pos->neighbor( next_around_edge( pos->index(_c->vertex(_t)),
					   pos->index(_c->vertex(_s)) ) );
    return *this;
  }

  Cell_circulator operator--(int)
  {
    CGAL_triangulation_precondition( pos != NULL );
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
    return pos;
  }

  bool operator==(const Cell_circulator & ccir) const
  {
    return _tds == ccir._tds &&
	   _c->vertex(_s) == ccir._c->vertex(ccir._s) &&
	   _c->vertex(_t) == ccir._c->vertex(ccir._t) &&
	   pos == ccir.pos;
  }

  bool operator!=(const Cell_circulator & ccir) const
  {
    return ! (*this == ccir);
  }

private:
  Tds* _tds;
  Cell* _c;  // cell containing the considered edge
  int _s;    // index of the source vertex of the edge in _c
  int _t;    // index of the target vertex of the edge in _c
  Cell* pos; // current cell
};

template < class Tds_ >
class Triangulation_ds_facet_circulator_3
  : public Bidirectional_circulator_base<typename Tds_::Facet,
    ptrdiff_t, size_t>,
    public Triangulation_utils_3
{
  // circulates on facets around a given edge
public:

  typedef Tds_                  Tds;
  typedef typename Tds::Cell    Cell;
  typedef typename Tds::Facet   Facet;
  typedef typename Tds::Edge    Edge;
  typedef typename Tds::Vertex  Vertex;

  typedef Triangulation_ds_facet_circulator_3<Tds> Facet_circulator;

  Triangulation_ds_facet_circulator_3()
    : _tds(NULL), _c(NULL), _s(0), _t(0), pos(NULL)
  {}

  Triangulation_ds_facet_circulator_3(Tds * tds, Cell* c, int s, int t)
    : _tds(tds), _c(c), _s(s), _t(t)
  {
    CGAL_triangulation_precondition( c != NULL &&
				     s >= 0 && s < 4 &&
				     t >= 0 && t < 4 );
    //     if ( _tds->dimension() <3 ) useless since precondition in tds
    //     incident_facets
    pos = c;
  }

  Triangulation_ds_facet_circulator_3(Tds * tds, const Edge & e)
    : _tds(tds), _c(e.first), _s(e.second), _t(e.third)
  {
    CGAL_triangulation_precondition( e.first != NULL &&
				     e.second >= 0 && e.second < 4 &&
				     e.third  >= 0 && e.third  < 4);
    //     if ( _tds->dimension() <3 ) useless since precondition in tds
    //     incident_facets
      pos = e.first;
  }

  Triangulation_ds_facet_circulator_3(Tds * tds, Cell* c, int s, int t,
				      Cell* start, int f)
    : _tds(tds), _c(c), _s(s), _t(t)
  {
    CGAL_triangulation_precondition( c != NULL &&
				     s >= 0 && s < 4 &&
				     t >= 0 && t < 4 &&
				     f >= 0 && f < 4 &&
                                     start->has_vertex( c->vertex(s) ) &&
	                             start->has_vertex( c->vertex(t) ) );

    int i = start->index( c->vertex(s) );
    int j = start->index( c->vertex(t) );

    CGAL_triangulation_precondition( f!=i && f!=j );

    if ( f == (int) next_around_edge(i,j) )
	pos = start;
    else
	pos = start->neighbor(f); // other cell with same facet
  }

  Triangulation_ds_facet_circulator_3(Tds * tds, Cell* c, int s, int t,
				      const Facet & start)
    : _tds(tds), _c(c), _s(s), _t(t)
  {
    CGAL_triangulation_precondition( c != NULL &&
				     s >= 0 && s < 4 &&
				     t >= 0 && t < 4 &&
                                     start.first->has_vertex( c->vertex(s) ) &&
	                             start.first->has_vertex( c->vertex(t) ) );

    int i = start.first->index( c->vertex(s) );
    int j = start.first->index( c->vertex(t) );

    CGAL_triangulation_precondition( start.second !=i && start.second !=j );

    if ( start.second == (int) next_around_edge(i,j) )
	pos = start.first;
    else
      pos = start.first->neighbor(start.second); // other cell with same facet
  }

  Triangulation_ds_facet_circulator_3(Tds * tds, const Edge & e,
				       Cell* start, int f)
    : _tds(tds), _c(e.first), _s(e.second), _t(e.third)
  {
    CGAL_triangulation_precondition( e.first != NULL &&
				     e.second >= 0 && e.second < 4 &&
				     e.third  >= 0 && e.third  < 4 &&
				     f >= 0 && f < 4 &&
                     start.first->has_vertex( e.first->vertex(e.second) ) &&
	             start.first->has_vertex( e.first->vertex(e.third) ) );

    int i = start.first->index( e.first->vertex(e.second) );
    int j = start.first->index( e.first->vertex(e.third) );

    CGAL_triangulation_precondition( f!=i && f!=j );

    if ( f == (int) next_around_edge(i,j) )
	pos = start.first;
    else
	pos = start.first->neighbor(f); // other cell with same facet
  }

  Triangulation_ds_facet_circulator_3(Tds * tds, const Edge & e,
				     const Facet & start)
    : _tds(tds), _c(e.first), _s(e.second), _t(e.third)
  {
    CGAL_triangulation_precondition( e.first != NULL &&
				     e.second >= 0 && e.second < 4 &&
				     e.third  >= 0 && e.third  < 4 &&
                       start.first->has_vertex( e.first->vertex(e.second) ) &&
	               start.first->has_vertex( e.first->vertex(e.third) ) );

    int i = start.first->index( e.first->vertex(e.second) );
    int j = start.first->index( e.first->vertex(e.third) );

    if ( start.second == (int) next_around_edge(i,j) )
	pos = start.first;
    else
	pos = start.first->neighbor(start.second);
  }

  Triangulation_ds_facet_circulator_3(const Facet_circulator & ccir)
    : _tds(ccir._tds), _c(ccir._c), _s(ccir._s), _t(ccir._t), pos(ccir.pos)
  {}

  Facet_circulator & operator++()
  {
    CGAL_triangulation_precondition( pos != NULL );
    //then dimension() cannot be < 3

    pos = pos->neighbor( next_around_edge( pos->index(_c->vertex(_s)),
					   pos->index(_c->vertex(_t)) ) );
    return *this;
  }

  Facet_circulator operator++(int)
  {
    CGAL_triangulation_precondition( pos != NULL );
    Facet_circulator tmp(*this);
    ++(*this);
    return tmp;
  }

  Facet_circulator & operator--()
  {
    CGAL_triangulation_precondition( pos != NULL );

    pos = pos->neighbor( next_around_edge( pos->index(_c->vertex(_t)),
					   pos->index(_c->vertex(_s)) ) );
    return *this;
  }

  Facet_circulator operator--(int)
  {
    CGAL_triangulation_precondition( pos != NULL );
    Facet_circulator tmp(*this);
    --(*this);
    return tmp;
  }

  Facet operator*() const
  {
    return std::make_pair(pos,
			  next_around_edge( pos->index(_c->vertex(_s)),
					    pos->index(_c->vertex(_t)) ) );
  }

  bool operator==(const Facet_circulator & ccir) const
  {
    return _tds == ccir._tds &&
	   _c->vertex(_s) == ccir._c->vertex(ccir._s) &&
	   _c->vertex(_t) == ccir._c->vertex(ccir._t) &&
	   pos == ccir.pos;
  }

  bool operator!=(const Facet_circulator & ccir) const
  {
    return ! (*this == ccir);
  }

private:
  Tds* _tds;
  Cell* _c;  // cell containing the considered edge
  int _s;    // index of the source vertex of the edge in _c
  int _t;    // index of the target vertex of the edge in _c
  Cell* pos; // current cell
  // the current facet is the facet of pos numbered
  // next_around_edge( pos->index(_c->vertex(_s)),
  //                   pos->index(_c->vertex(_t)) )
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_DS_CIRCULATORS_3_H
