// Copyright (c) 1999-2004,2006-2009,2014-2016   INRIA Nancy - Grand Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Iordan Iordanov  <Iordan.Iordanov@loria.fr>
//

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_2_H

#include <CGAL/license/Periodic_4_hyperbolic_triangulation_2.h>

#include <CGAL/basic.h>

#include <CGAL/Triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_face_base_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>

#include <CGAL/NT_converter.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/utility.h>
#include <CGAL/use.h>

#include <boost/tuple/tuple.hpp>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <list>

#if defined PROFILING_MODE
#include <CGAL/Timer.h>
extern long calls_predicate_identity;
extern long calls_predicate_non_identity;
extern double time_predicate_identity;
extern double time_predicate_non_identity;
#endif

namespace CGAL {

template <class GT,
          class TDS = Triangulation_data_structure_2<
                        Periodic_4_hyperbolic_triangulation_vertex_base_2<GT>,
                        Periodic_4_hyperbolic_triangulation_face_base_2<GT> > >
class Periodic_4_hyperbolic_triangulation_2
{
  typedef Periodic_4_hyperbolic_triangulation_2<GT, TDS>          Self;

private:
  typedef typename GT::FT                                         FT;

public:
  typedef GT                                                      Geom_traits;
  typedef TDS                                                     Triangulation_data_structure;
  typedef typename GT::Hyperbolic_translation                     Hyperbolic_translation;
  typedef typename TDS::Vertex::Point                             Point;
  typedef typename GT::Hyperbolic_segment_2                       Hyperbolic_segment;
  typedef typename GT::Hyperbolic_triangle_2                      Hyperbolic_triangle;

  typedef std::pair<Point, Hyperbolic_translation>                     Periodic_point;
  typedef std::array< std::pair<Point, Hyperbolic_translation>, 2 >    Periodic_segment;
  typedef std::array< std::pair<Point, Hyperbolic_translation>, 3 >    Periodic_triangle;

  typedef typename TDS::Vertex                                    Vertex;
  typedef typename TDS::Edge                                      Edge;
  typedef typename TDS::Face                                      Face;

  typedef typename TDS::Vertex_handle                             Vertex_handle;
  typedef typename TDS::Face_handle                               Face_handle;

  typedef typename TDS::size_type                                 size_type;
  typedef typename TDS::difference_type                           difference_type;

  typedef typename TDS::Face_iterator                             Face_iterator;
  typedef typename TDS::Edge_iterator                             Edge_iterator;
  typedef typename TDS::Vertex_iterator                           Vertex_iterator;
  typedef typename TDS::Face_circulator                           Face_circulator;
  typedef typename TDS::Edge_circulator                           Edge_circulator;
  typedef typename TDS::Vertex_circulator                         Vertex_circulator;

  typedef Point                                                   value_type;
  typedef const value_type&                                       const_reference;

  enum Locate_type {
    VERTEX = 0,
    EDGE,
    FACE
  };

protected:
  GT _gt;
  TDS _tds;

public:
  Periodic_4_hyperbolic_triangulation_2(const Geom_traits &gt = Geom_traits())
    : _gt(gt), _tds()
  {
    init_tds();
  }

  Periodic_4_hyperbolic_triangulation_2(const Periodic_4_hyperbolic_triangulation_2& tr)
    : _gt(tr.geom_traits())
  {
    CGAL_triangulation_expensive_postcondition(*this == tr);
  }

  Periodic_4_hyperbolic_triangulation_2& operator=(Periodic_4_hyperbolic_triangulation_2 tr)
  {
    swap(tr);
    return *this;
  }

  void swap(Periodic_4_hyperbolic_triangulation_2& tr)
  {
    _tds.swap(tr._tds);
    std::swap(tr._gt, _gt);
  }

  void clear()
  {
    _tds.clear();
    init_tds();
  }

  int cw(int i) const { return _tds.cw(i); }
  int ccw(int i) const { return _tds.ccw(i); }

private:
  void init_tds() { _tds.set_dimension(-2); }

public:
  const Geom_traits& geom_traits() const { return _gt; }

  /// Returns the data structure storing the triangulation.
  const TDS & tds() const { return _tds; }

  /// Returns the data structure storing the triangulation.
  TDS & tds() { return _tds; }

  size_type number_of_vertices() const { return _tds.number_of_vertices(); }
  size_type number_of_edges() const { return _tds.number_of_edges(); }
  size_type number_of_faces() const { return _tds.number_of_faces(); }

  Orientation orientation(const Point& p1, const Point& p2, const Point& p3) const
  {
    return geom_traits().orientation_2_object()(p1, p2, p3);
  }

  Orientation orientation(const Point& p1, const Point& p2, const Point& p3,
                          const Hyperbolic_translation& o1, const Hyperbolic_translation& o2, const Hyperbolic_translation& o3) const
  {
    return geom_traits().orientation_2_object()(p1, p2, p3, o1, o2, o3);
  }

  Oriented_side side_of_oriented_circle(const Point& p0, const Point& p1, const Point& p2, const Point &q) const
  {
    return geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, q);
  }

  Oriented_side side_of_oriented_circle(const Point& p0, const Point& p1, const Point& p2, const Point &q,
                                        const Hyperbolic_translation &o0, const Hyperbolic_translation &o1,
                                        const Hyperbolic_translation &o2, const Hyperbolic_translation &oq) const
  {
    return geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, q, o0, o1, o2, oq);
  }

  Periodic_point periodic_point(const Face_handle c, int i) const
  {
    CGAL_triangulation_precondition(i >= 0 && i <= 2);
    return std::make_pair(c->vertex(i)->point(), c->translation(i));
  }

  Periodic_segment periodic_segment(const Point& p1, const Point& p2,
                                    const Hyperbolic_translation& tr1, const Hyperbolic_translation& tr2) const
  {
    return make_array(std::make_pair(p1, tr1), std::make_pair(p2, tr2));
  }


  Periodic_segment periodic_segment(const Point& p1, const Point& p2) const
  {
    return make_array(std::make_pair(p1, Hyperbolic_translation()),
                      std::make_pair(p2, Hyperbolic_translation()));
  }

  Periodic_segment periodic_segment(const Face_handle c, int i, int j) const
  {
    CGAL_triangulation_precondition(i != j);
    CGAL_triangulation_precondition((i >= 0 && i <= 2) && (j >= 0 && j <= 2));
    return periodic_segment(c->vertex(i)->point(), c->vertex(j)->point(),
                            c->translation(i), c->translation(j));
  }

  Periodic_segment periodic_segment(const Edge & e) const
  {
    return periodic_segment(e.first, e.second, e.third);
  }

  Periodic_segment periodic_segment(const Edge & e, const Hyperbolic_translation& o) const
  {
    Face_handle f = e.first;
    int i = e.second;
    return periodic_segment(f->vertex(cw(i))->point(), f->vertex(ccw(i))->point(),
                            o * f->translation(cw(i)), o * f->translation(ccw(i)));
  }

  Periodic_triangle periodic_triangle(const Point& p1, const Point& p2, const Point& p3,
                                      const Hyperbolic_translation tr1, const Hyperbolic_translation& tr2,
                                      const Hyperbolic_translation& tr3) const
  {
    return make_array(std::make_pair(p1, tr1),
                      std::make_pair(p2, tr2),
                      std::make_pair(p3, tr3));
  }

  Periodic_triangle periodic_triangle(const Point& p1, const Point& p2, const Point& p3) const
  {
    return periodic_triangle(p1, p2, p3,
                             Hyperbolic_translation(), Hyperbolic_translation(), Hyperbolic_translation());
  }

  Periodic_triangle periodic_triangle(const Face & f) const
  {
    return periodic_triangle(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point(),
                             f->translation(0), f->translation(0), f->translation(0));
  }

  Periodic_triangle periodic_triangle(const Face & f, const Hyperbolic_translation& tr) const
  {
    return periodic_triangle(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point(),
                             tr*f->translation(0), tr*f->translation(0), tr*f->translation(0));
  }

  Point construct_point(const Point& p, const Hyperbolic_translation& tr) const
  {
    return geom_traits().construct_hyperbolic_point_2_object()(p, tr);
  }

  Point construct_point(const Periodic_point & pp) const
  {
    return construct_point(pp.first, pp.second);
  }

  Hyperbolic_segment construct_hyperbolic_segment(const Point& src, const Point& tgt) const
  {
    return geom_traits().construct_hyperbolic_segment_2_object()(src, tgt);
  }

  Hyperbolic_segment construct_hyperbolic_segment(const Point& src, const Point& tgt,
                                                  const Hyperbolic_translation& tr1, const Hyperbolic_translation& tr2) const
  {
    return geom_traits().construct_hyperbolic_segment_2_object()(src, tgt, tr1, tr2);
  }

  Hyperbolic_segment construct_hyperbolic_segment(const Face_handle & fh, int idx) const
  {
    CGAL_triangulation_precondition(idx >= 0 && idx <= 2);
    return construct_hyperbolic_segment(fh->vertex(ccw(idx))->point(), fh->vertex(cw(idx))->point(),
                                        fh->translation(ccw(idx)), fh->translation(cw(idx)));
  }

  Hyperbolic_segment construct_hyperbolic_segment(const std::pair<Face_handle, int> & edge) const
  {
    CGAL_triangulation_precondition(edge.second >= 0 && edge.second <= 2);
    return construct_hyperbolic_segment(edge.first, ccw(edge.second));
  }

  Hyperbolic_segment construct_hyperbolic_segment(const Periodic_segment & ps) const
  {
    return construct_hyperbolic_segment(ps[0].first,  ps[1].first, ps[0].second, ps[1].second);
  }

  Hyperbolic_triangle construct_hyperbolic_triangle(const Face_handle & fh) const
  {
    return geom_traits().construct_triangle_2_object()(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point(),
                                                       fh->translation(0), fh->translation(1), fh->translation(2));
  }

  Hyperbolic_triangle construct_hyperbolic_triangle(const Periodic_triangle & pt) const
  {
    return construct_triangle(pt[0].first, pt[1].first, pt[2].first,
                              pt[0].second,pt[1].second,pt[2].second);
  }

  bool is_vertex(Vertex_handle v) const { return _tds.is_vertex(v); }
  bool is_edge(Vertex_handle u, Vertex_handle v, Face_handle& fh, int& i) const
  {
    return _tds.is_edge(u, v, fh, i);
  }

  bool is_face(Vertex_handle u, Vertex_handle v, Vertex_handle w, Face_handle & fh) const
  {
    return _tds.is_face(u, v, w, fh);
  }

  bool has_vertex(const Face_handle f, const Vertex_handle v, int& i) const
  {
    if(f->vertex(0) == v) {
      i = 0;
      return true;
    } else {
      if(f->vertex(1) == v) {
        i = 1;
        return true;
      } else {
        if(f->vertex(2) == v) {
          i = 2;
          return true;
        } else {
          return false;
        }
      }
    }
  }

protected:
  bool has_inexact_negative_orientation(const Point& p, const Point& q, const Point& r,
                                        const Hyperbolic_translation& ofp, const Hyperbolic_translation& ofq, const Hyperbolic_translation& ofr) const
  {
    return has_inexact_negative_orientation(construct_point(p, ofp),
                                            construct_point(q, ofq),
                                            construct_point(r, ofr));
  }

  bool has_inexact_negative_orientation(const Point &p, const Point &q, const Point &r) const
  {
    const double px = to_double(p.x());
    const double py = to_double(p.y());
    const double qx = to_double(q.x());
    const double qy = to_double(q.y());
    const double rx = to_double(r.x());
    const double ry = to_double(r.y());

    const double pqx = qx - px;
    const double pqy = qy - py;
    const double prx = rx - px;
    const double pry = ry - py;

    return (determinant(pqx, pqy, prx, pry) < 0);
  }

  Face_handle inexact_euclidean_locate(const Point& p, Hyperbolic_translation& o, Face_handle fh = Face_handle()) const;

  Face_handle euclidean_locate(const Point& p,
                               Locate_type& lt,
                               int& li,
                               Hyperbolic_translation& lo,
                               Face_handle f = Face_handle()) const;

  Face_handle euclidean_locate(const Point& p,
                               Hyperbolic_translation& lo,
                               Face_handle f = Face_handle()) const
  {
    Locate_type lt;
    int li;
    return euclidean_locate(p, lt, li, lo, f);
  }

  Oriented_side side_of_hyperbolic_triangle(const Point p, const Point q, const Point r,
                                            const Point query, Locate_type &lt, int& li) const
  {
    // Point p is assumed to be at index 0, q at index 1 and r at index 2 in the face.
    li = -1;

    if(query == p)
    {
      lt = VERTEX;
      li = 0;
      return ON_ORIENTED_BOUNDARY;
    }

    if(query == q)
    {
      lt = VERTEX;
      li = 1;
      return ON_ORIENTED_BOUNDARY;
    }

    if(query == r)
    {
      lt = VERTEX;
      li = 2;
      return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side cp1 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(p, q, query);
    if(cp1 == ON_ORIENTED_BOUNDARY)
    {
      lt = EDGE;
      li = 2;
      return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side cp2 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(q, r, query);
    if(cp2 == ON_ORIENTED_BOUNDARY)
    {
      lt = EDGE;
      li = 0;
      return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side cp3 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(r, p, query);
    if(cp3 == ON_ORIENTED_BOUNDARY)
    {
      lt = EDGE;
      li = 1;
      return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side cs1 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(p, q, r);
    Oriented_side cs2 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(q, r, p);
    Oriented_side cs3 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(r, p, q);

    // Cannot be on the boundary here.
    lt = FACE;
    if(cs1 != cp1) {
      li = 2;
      return ON_NEGATIVE_SIDE;
    }
    if(cs2 != cp2){
      li = 0;
      return ON_NEGATIVE_SIDE;
    }
    if(cs3 != cp3) {
      li = 1;
      return ON_NEGATIVE_SIDE;
    }

    return ON_POSITIVE_SIDE;
  }

public:
  Face_handle hyperbolic_periodic_locate(const Point& p, Locate_type& lt, int& li,
                                         Hyperbolic_translation& lo,
                                         const Face_handle fh = Face_handle()) const;
  Face_handle hyperbolic_periodic_locate(const Point& p, Hyperbolic_translation& lo,
                                         const Face_handle fh = Face_handle()) const
  {
    Locate_type lt;
    int li;
    return hyperbolic_periodic_locate(p, lt, li, lo, fh);
  }

  Face_handle hyperbolic_locate(const Point& p,
                                Locate_type& lt,
                                int& li,
                                Face_handle start = Face_handle()) const
  {
    Hyperbolic_translation lo;
    return hyperbolic_periodic_locate(p, lt, li, lo, start);
  }

  Face_handle hyperbolic_locate(const Point& p,
                                Face_handle start = Face_handle()) const
  {
    Hyperbolic_translation lo;
    Locate_type lt;
    int li;
    return hyperbolic_periodic_locate(p, lt, li, lo, start);
  }

public:
  Vertex_iterator vertices_begin() const { return _tds.vertices_begin(); }
  Vertex_iterator vertices_end() const { return _tds.vertices_end(); }
  Edge_iterator edges_begin() const { return _tds.edges_begin(); }
  Edge_iterator edges_end() const { return _tds.edges_end(); }
  Face_iterator faces_begin() const { return _tds.faces_begin(); }
  Face_iterator faces_end() const { return _tds.faces_end(); }

  // Circulators
  Vertex_circulator adjacent_vertices(Vertex_handle v) const {
    return _tds.incident_vertices(v, v->face());
  }
  Vertex_circulator adjacent_vertices(Vertex_handle v, Face_handle f) const {
    return _tds.incident_vertices(v, f);
  }

  Edge_circulator incident_edges(Vertex_handle v) const {
    return _tds.incident_edges(v, v->face());
  }
  Edge_circulator incident_edges(Vertex_handle v, Face_handle f) const {
    return _tds.incident_edges(v, f);
  }

  Face_circulator incident_faces(Vertex_handle v) const {
    return _tds.incident_faces(v, v->face());
  }
  Face_circulator incident_faces(Vertex_handle v, Face_handle f) const {
    return _tds.incident_faces(v, f);
  }

  // around a vertex
  template <class OutputIterator>
  OutputIterator incident_faces(Vertex_handle v, OutputIterator faces) const {
    return _tds.incident_faces(v, faces);
  }

  template <class OutputIterator>
  OutputIterator incident_edges(Vertex_handle v, OutputIterator edges) const {
    return _tds.incident_edges(v, edges);
  }

  template <class OutputIterator>
  OutputIterator adjacent_vertices(Vertex_handle v, OutputIterator vertices) const {
    return _tds.incident_vertices(v, vertices);
  }

  size_type degree(Vertex_handle v) const { return _tds.degree(v); }

  // Functions forwarded from TDS.
  int mirror_index(Face_handle c, int i) const { return _tds.mirror_index(c, i); }
  Vertex_handle mirror_vertex(Face_handle c, int i) const { return _tds.mirror_vertex(c, i); }
  Edge mirror_edge(Edge e) const { return _tds.mirror_edge(e); }

  Hyperbolic_translation neighbor_translation(const Face_handle fh, int i) const
  {
    CGAL_triangulation_precondition(i >= 0 && i <= 2);
    int myi = ccw(i);
    Hyperbolic_translation myof = fh->translation(myi);

    Hyperbolic_translation nbof;
    for(int c=0; c<3; ++c)
    {
      if(fh->neighbor(i)->vertex(c) == fh->vertex(myi))
      {
        nbof = fh->neighbor(i)->translation(c);
        break;
      }
    }

    return (myof - nbof);
  }

private:
  bool has_self_edges() const
  {
    Face_iterator it;
    for(it=faces_begin(); it!=faces_end(); ++it)
      if(has_self_edges(it))
        return true;

    return false;
  }

  bool has_self_edges(Face_handle c) const;

  bool has_cycles_length_2() const
  {
    Vertex_iterator it;
    for(it = vertices_begin(); it!= vertices_end(); ++it)
      if(has_cycles_length_2(it))
        return true;

    return false;
  }

  bool has_cycles_length_2(Vertex_handle v) const;

protected:
  void make_canonical(Face_handle fh) const
  {
    // If all translation are the same, there is an image of the face
    // inside the original octagon. We store that one. This covers the
    // simplest case.
    if(fh->translation(0) == fh->translation(1) && fh->translation(1) == fh->translation(2))
    {
      fh->set_translation(0, Hyperbolic_translation());
      fh->set_translation(1, Hyperbolic_translation());
      fh->set_translation(2, Hyperbolic_translation());
      return;
    }

    // This covers the cases in which two vertices lie in the same domain.
    for(int i=0; i<3; ++i)
    {
      int j = (i + 1) % 3;

      if(fh->translation(i) == fh->translation(j))
      {
        int k = (i + 2) % 3;
        if((fh->translation(i).inverse()*fh->translation(k)) < (fh->translation(k).inverse()*fh->translation(i))) {
          fh->set_translation(k, fh->translation(i).inverse() * fh->translation(k));
          fh->set_translation(i, Hyperbolic_translation());
          fh->set_translation(j, Hyperbolic_translation());
          return;
        }
        else
        {
          fh->set_translation(i, fh->translation(k).inverse() * fh->translation(i));
          fh->set_translation(j, fh->translation(k).inverse() * fh->translation(j));
          fh->set_translation(k, Hyperbolic_translation());
          return;
        }
      }
      else
      {
        continue;
      }
    }

    // Now we know that all vertices lie in different regions.
    Hyperbolic_translation vmin(7, 2, 5);
    Hyperbolic_translation trans;
    for(int i=0; i<3; ++i)
    {
      int j = (i + 1) % 3; // the index of the 'next' vertex
      Hyperbolic_translation tmp = fh->translation(i).inverse() * fh->translation(j);
      if(tmp < vmin)
      {
        vmin = tmp;
        trans = fh->translation(i).inverse();
      }
    }

    if(!trans.is_identity())
    {
      fh->set_translation(0, trans * fh->translation(0));
      fh->set_translation(1, trans * fh->translation(1));
      fh->set_translation(2, trans * fh->translation(2));
    }
  }

public:
  bool is_valid(bool verbose = false) const;
  bool is_valid(Face_handle c, bool verbose = false) const;

  int dimension() const { return _tds.dimension(); }

}; // class Periodic_4_hyperbolic_triangulation_2

/*********** FUNCTION IMPLEMENTATIONS *************/

/// tests if two vertices of one cell are just periodic copies of each other
template < class GT, class TDS >
inline bool
Periodic_4_hyperbolic_triangulation_2<GT, TDS>::
has_self_edges(typename TDS::Face_handle c) const
{
  CGAL_triangulation_assertion((c->vertex(0) != c->vertex(1)) || (c->translation(0) != c->translation(1)));
  CGAL_triangulation_assertion((c->vertex(0) != c->vertex(2)) || (c->translation(0) != c->translation(2)));
  CGAL_triangulation_assertion((c->vertex(1) != c->vertex(2)) || (c->translation(1) != c->translation(2)));

  return ((c->vertex(0) == c->vertex(1)) ||
          (c->vertex(0) == c->vertex(2)) ||
          (c->vertex(1) == c->vertex(2)));
}

template < class GT, class TDS >
inline bool
Periodic_4_hyperbolic_triangulation_2<GT, TDS>::
has_cycles_length_2(typename TDS::Vertex_handle v) const
{
  typename TDS::Vertex_circulator vc = adjacent_vertices(v), done(vc);
  std::set<Vertex_handle> check;
  do
  {
    std::pair<typename std::set<Vertex_handle>::iterator, bool> res = check.insert(vc);
    if(!res.second)
      return true;
  } while(++vc != done);

  return false;
}

/*! \brief Tests if the triangulation is valid.
 *
 * A triangulation is valid if
 * - A cell is not its own neighbor.
 * - A cell has no two equal neighbors
 * - A cell has no two equal vertex-translation pairs
 * - A cell is positively oriented.
 * - The point of a neighbor of cell c that does not belong to c is not inside
 *   the circumcircle of c.
 */
template < class GT, class TDS >
bool
Periodic_4_hyperbolic_triangulation_2<GT,TDS>::
is_valid(bool verbose) const {

  bool error = false;
  for(Face_iterator cit = faces_begin(); cit != faces_end(); ++cit)
  {
    for(int i=0; i<3; ++i)
    {
      CGAL_triangulation_assertion(cit != cit->neighbor(i));
      for(int j=i+1; j<3; ++j)
      {
        CGAL_triangulation_assertion(cit->neighbor(i) != cit->neighbor(j));
        CGAL_triangulation_assertion(cit->vertex(i) != cit->vertex(j));
      }
    }

#if defined PROFILING_MODE
    if(cit->translation(0).is_identity() && cit->translation(1).is_identity() && cit->translation(2).is_identity()) {
      calls_predicate_identity++;
    } else {
      calls_predicate_non_identity++;
    }
#endif

#if defined PROFILING_MODE
    CGAL::Timer tmr;
    tmr.start();
#endif

    Orientation ori = orientation(cit->vertex(0)->point(), cit->vertex(1)->point(), cit->vertex(2)->point(),
                                  cit->translation(0), cit->translation(1), cit->translation(2));

#if defined PROFILING_MODE
    tmr.stop();
    if(cit->translation(0).is_identity() && cit->translation(1).is_identity() && cit->translation(2).is_identity()) {
      time_predicate_identity += tmr.time();
    } else {
      time_predicate_non_identity += tmr.time();
    }
#endif

    if(ori != POSITIVE)
    {
      if(verbose)
        std::cerr << "Orientation test failed!" << std::endl;

      error = true;
    }
  }

  if(!error)
  {
    CGAL_triangulation_assertion(!has_self_edges());
    CGAL_triangulation_assertion(!has_cycles_length_2());
    CGAL_triangulation_assertion(_tds.number_of_vertices() + _tds.number_of_faces() + 2 == _tds.number_of_edges());
    return true;
  }
  else
  {
    return false;
  }
}

template < class GT, class TDS >
bool
Periodic_4_hyperbolic_triangulation_2<GT,TDS>::
is_valid(Face_handle ch, bool verbose) const
{
  if(!_tds.is_valid(ch,verbose))
    return false;

  bool error = false;
  const Point *p[3]; Hyperbolic_translation off[3];

  for(int i=0; i<3; ++i)
  {
    p[i] = &ch->vertex(i)->point();
    off[i] = ch->translation(i);
  }

#if defined PROFILING_MODE
  if(off[0].is_identity() && off[1].is_identity() && off[2].is_identity()) {
    calls_predicate_identity++;
  } else {
    calls_predicate_non_identity++;
  }
#endif

#if defined PROFILING_MODE
  CGAL::Timer tmr;
  tmr.start();
#endif

  if(orientation(*p[0], *p[1], *p[2], off[0], off[1], off[2]) != POSITIVE)
    error = true;

#if defined PROFILING_MODE
  tmr.stop();
  if(off[0].is_identity() && off[1].is_identity() && off[2].is_identity()) {
    time_predicate_identity += tmr.time();
  } else {
    time_predicate_non_identity += tmr.time();
  }
#endif

  return !error;
}

/*********************************************************************************/

template <class GT, class TDS>
typename TDS::Face_handle
Periodic_4_hyperbolic_triangulation_2<GT, TDS>::
inexact_euclidean_locate(const Point& p, Hyperbolic_translation& loff, Face_handle f) const
{
  typename GT::Side_of_original_octagon check;

  if(check(p) != ON_BOUNDED_SIDE)
    return Face_handle();

  if(f == Face_handle())
    f = _tds.faces().begin();

  int curr = 0;
  int succ = ccw(curr);
  int counter = 0;

  for(;;)
  {
#if defined PROFILING_MODE
    if((loff * f->translation(curr)).is_identity() &&
       (loff * f->translation(succ)).is_identity()) {
      calls_predicate_identity++;
    } else {
      calls_predicate_non_identity++;
    }
#endif

#if defined PROFILING_MODE
    CGAL::Timer tmr;
    tmr.start();
#endif

    bool on_negative_side = has_inexact_negative_orientation(
                              f->vertex(curr)->point(), f->vertex(succ)->point(), p,
                              loff * f->translation(curr), loff * f->translation(succ), Hyperbolic_translation());

#if defined PROFILING_MODE
    tmr.stop();
    if((loff * f->translation(curr)).is_identity() &&
       (loff * f->translation(succ)).is_identity()) {
      time_predicate_identity += tmr.time();
    } else {
      time_predicate_non_identity += tmr.time();
    }
#endif

    if(on_negative_side)
    {
      loff = loff * f->neighbor_translation(cw(curr));
      f = f->neighbor(cw(curr));
      curr = ccw(curr);
      succ = ccw(curr);
      counter = 0;
    }
    else
    {
      curr = succ;
      succ = ccw(curr);
      counter++;
    }

    if(counter > 2)
      break;
  }

  return f;
}


template <class GT, class TDS>
typename TDS::Face_handle Periodic_4_hyperbolic_triangulation_2<GT, TDS>::
euclidean_locate(const Point& p,
                 Locate_type& lt,
                 int& li,
                 Hyperbolic_translation& loff, Face_handle f) const
{

  CGAL::Bounded_side side = geom_traits().side_of_original_octagon_object()(p);
  if(side != ON_BOUNDED_SIDE)
    return Face_handle();

  // Handle the case where an initial Face_handle is not given
  if(f == Face_handle())
    f = _tds.faces().begin();

  int curr = 0;
  int succ = ccw(curr);
  int counter = 0;

  for(;;)
  {
#if defined PROFILING_MODE
    if((loff * f->translation(curr)).is_identity() &&
       (loff * f->translation(succ)).is_identity()) {
      calls_predicate_identity++;
    } else {
      calls_predicate_non_identity++;
    }
#endif

#if defined PROFILING_MODE
    CGAL::Timer tmr;
    tmr.start();
#endif

    Orientation o = orientation(f->vertex(curr)->point(), f->vertex(succ)->point(), p,
                                loff * f->translation(curr), loff * f->translation(succ), Hyperbolic_translation());

#if defined PROFILING_MODE
    tmr.stop();
    if((loff * f->translation(curr)).is_identity() &&
       (loff * f->translation(succ)).is_identity()) {
      time_predicate_identity += tmr.time();
    } else {
      time_predicate_non_identity += tmr.time();
    }
#endif

    if(o == NEGATIVE)
    {
      loff = loff * neighbor_translation(f, cw(curr));
      f = f->neighbor(cw(curr));
      curr = ccw(curr);
      succ = ccw(curr);
      counter = 0;
    }
    else
    {
      curr = succ;
      succ = ccw(curr);
      counter++;
    }

    if(counter > 2)
      break;
  }

  const Point p0 = construct_point(f->vertex(0)->point(), loff * f->translation(0));
  const Point p1 = construct_point(f->vertex(1)->point(), loff * f->translation(1));
  const Point p2 = construct_point(f->vertex(2)->point(), loff * f->translation(2));

  lt = FACE;
  li = -1;

  if(p == p0) {
    lt = VERTEX;
    li = 0;
  } else {
    if(p == p1) {
      lt = VERTEX;
      li = 1;
    } else {
      if(p == p2) {
        lt = VERTEX;
        li = 2;
      }
    }
  }

  return f;
}

template <class GT, class TDS>
typename TDS::Face_handle
Periodic_4_hyperbolic_triangulation_2<GT, TDS>::
hyperbolic_periodic_locate(const Point& p,
                           Locate_type& lt,
                           int& li,
                           Hyperbolic_translation& lo,
                           Face_handle start) const
{
  // Get a hint of where the point is located. It's either in lf or in one of its neighbors.
  Face_handle lf = euclidean_locate(p, lt, li, lo, start);
  if(lf == Face_handle())
    return lf;

  // The input point has been located in a vertex, so we can just return here, nothing more to do.
  if(lt == VERTEX)
    return lf;

  //typedef typename GT::Side_of_hyperbolic_triangle_2 Side_of_hyperbolic_triangle_2;
  //Side_of_hyperbolic_triangle_2 sf;

  Point p0 = construct_point(construct_point(lf->vertex(0)->point(), lf->translation(0)) , lo);
  Point p1 = construct_point(construct_point(lf->vertex(1)->point(), lf->translation(1)) , lo);
  Point p2 = construct_point(construct_point(lf->vertex(2)->point(), lf->translation(2)) , lo);

  Oriented_side os = side_of_hyperbolic_triangle(p0, p1, p2, p, lt, li);

  if(os == ON_NEGATIVE_SIDE)
  {
    Hyperbolic_translation tr = lo * neighbor_translation(lf, li);
    Face_handle nf = lf->neighbor(li);
    Point np0 = construct_point(construct_point(nf->vertex(0)->point(), nf->translation(0)) , tr);
    Point np1 = construct_point(construct_point(nf->vertex(1)->point(), nf->translation(1)) , tr);
    Point np2 = construct_point(construct_point(nf->vertex(2)->point(), nf->translation(2)) , tr);
    CGAL_triangulation_assertion_code(Oriented_side os1 = side_of_hyperbolic_triangle(np0, np1, np2, p, lt, li));
    CGAL_triangulation_assertion(os1 == ON_POSITIVE_SIDE);
    lo = tr;
    return nf;
  }

  return lf;
}

} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_2_H
