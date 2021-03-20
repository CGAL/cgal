// Copyright (c) 1997, 2O12, 2019  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Claudia Werner,  Mariette Yvinec

#ifndef CGAL_TRIANGULATION_SPHERE_LINE_FACE_CIRCULATOR_2_H
#define CGAL_TRIANGULATION_SPHERE_LINE_FACE_CIRCULATOR_2_H

#include <CGAL/license/Triangulation_on_sphere_2.h>

#include <CGAL/triangulation_assertions.h>

namespace CGAL {

template <class Triangulation_>
class Triangulation_sphere_line_face_circulator_2
  : public Bidirectional_circulator_base<typename Triangulation_::Triangulation_data_structure::Face,
                                         std::size_t>,
    public Triangulation_cw_ccw_2
{
public:
  typedef Triangulation_sphere_line_face_circulator_2<Triangulation_> Line_face_circulator;
  typedef Triangulation_                                              Triangulation;

  typedef typename Triangulation::Geom_traits                         Gt;
  typedef typename Triangulation_::Triangulation_data_structure       Tds;

  typedef typename Tds::Vertex                                        Vertex;
  typedef typename Tds::Edge                                          Edge;
  typedef typename Tds::Face                                          Face;
  typedef typename Tds::Vertex_handle                                 Vertex_handle;
  typedef typename Tds::Face_handle                                   Face_handle;
  typedef typename Tds::Face_circulator                               Face_circulator;

  typedef typename Gt::Point_2                                        Point;
  typedef typename Triangulation::Locate_type                         Locate_type;

  enum State {undefined = -1,
              vertex_vertex,
              vertex_edge,
              edge_vertex,
              edge_edge};

private:
  Face_handle pos;
  const Triangulation* _tr;
  State s;
  int i;
  Point p, q;

public:
  Triangulation_sphere_line_face_circulator_2()
    : pos(), _tr(NULL), s(undefined), i(-1)
  { }

  Triangulation_sphere_line_face_circulator_2(Vertex_handle v,
                                              const Triangulation* tr,
                                              const Point& dir);

  Line_face_circulator& operator++();
  Line_face_circulator& operator--();
  Line_face_circulator operator++(int);
  Line_face_circulator operator--(int);
  Face* operator->() { return &*pos; }
  Face& operator*() { return *pos; }
  Face_handle handle() { return pos; }

  operator const Face_handle() const { return pos; }
  bool operator==(const Line_face_circulator& lfc) const;
  bool operator!=(const Line_face_circulator& lfc) const;

  bool operator==(const Face_handle& fh) const { return fh == pos; }
  bool operator!=(const Face_handle& fh) const { return fh != pos; }

  bool operator==(Nullptr_t  CGAL_assertion_code(n)) const;
  bool operator!=(Nullptr_t n) const;
  bool  is_empty() const;
  bool locate(const Point& t, Locate_type &lt,  int &li);

  //private:
  Triangulation_sphere_line_face_circulator_2(const Face_handle& face,
                                              int index,
                                              State state,
                                              const Triangulation * t,
                                              const Point& pp,
                                              const Point& qq);
private:
  void increment();
  void decrement();
};

template <typename Triangulation>
inline
bool
operator==(typename Triangulation::Triangulation_data_structure::Face_handle fh,
           Triangulation_sphere_line_face_circulator_2<Triangulation> fc)
{
  return (fc == fh);
}

template <typename Triangulation>
inline
bool
operator!=(typename Triangulation::Triangulation_data_structure::Face_handle fh,
           Triangulation_sphere_line_face_circulator_2<Triangulation> fc)
{
  return (fc != fh);
}

template <typename Triangulation>
Triangulation_sphere_line_face_circulator_2<Triangulation>::
Triangulation_sphere_line_face_circulator_2(Vertex_handle v,
                                            const Triangulation* tr,
                                            const Point& dir)
  : pos(), _tr(tr), s(undefined)
  // begin at the face incident to v, traversed by the ray from v to  dir
  // or null iterator
{
  CGAL_precondition((_tr->dimension() == 2) &&
                                  (! _tr->xy_equal(v->point(), dir)));
  p = v->point();
  q = dir;

  // find a  vertex to the left of pq
  // if there is no, the line_face_circulator is null
  Face_circulator fc = _tr->incident_faces(v);
  Face_circulator done(fc);
  int ic = fc->index(v);
  Vertex_handle vt = fc->vertex(cw(ic));
  //Orientation o = _tr->orientation(p, q, vt->point());
  while( _tr->orientation(p, q, vt->point()) != LEFT_TURN)
  {
    ++fc;
    ic = fc->index(v);
    vt= fc->vertex(cw(ic));

    if(fc == done)
    {
      *this = Line_face_circulator();
      return;
    }
  }

  Vertex_handle vr = fc-> vertex(ccw(ic));
  Orientation pqr = RIGHT_TURN; // warning "pqr might be used uninitialized"
  while((pqr = _tr->orientation(p, q, vr->point())) == LEFT_TURN)
  {
    --fc;
    ic = fc->index(v);
    vr = fc-> vertex(ccw(ic));
  }

  // first intersected face found
  //  [pqr] is COLLINEAR or RIGHT_TURN
  ic = fc->index(v);
  vt= fc->vertex(cw(ic));
  CGAL_assertion (_tr->orientation(p, q, vt->point()) == LEFT_TURN);

  if(pqr == COLLINEAR)
  {
    pos = fc;
    s = vertex_vertex;
    i = ccw(ic);
  }
  else // pqr == RIGHT_TURN
  {
    pos = fc;
    s = vertex_edge;
    i = ic ;
  }
}

template <typename Triangulation>
inline
void
Triangulation_sphere_line_face_circulator_2<Triangulation>::
increment()
{
  CGAL_precondition(pos != Face_handle());
  if(s == vertex_vertex || s == edge_vertex)
  {
    Orientation o;
    do
    {
      Face_handle n = pos->neighbor(cw(i));
      i = n->index(pos);
      pos = n;

      o = _tr->orientation(p, q, pos->vertex(i)->point());
      i = cw(i);
    } while(o == LEFT_TURN);

    if(o == COLLINEAR)
    {
      s = vertex_vertex;
      i = ccw(i);
    }
    else
    {
      s = vertex_edge;
    }
  }
  else
  {
    Face_handle n = pos->neighbor(i);
    int ni = n->index(pos);
    pos = n ;
    Orientation o = _tr->orientation(p, q, pos->vertex(ni)->point());

    switch(o)
    {
      case LEFT_TURN:
        s = edge_edge;
        i = ccw(ni);
        break;
      case RIGHT_TURN:
        s = edge_edge;
        i = cw(ni);
        break;
      default:
        s = edge_vertex;
        i = ni;
    }
  }
}

template <typename Triangulation>
void
Triangulation_sphere_line_face_circulator_2<Triangulation>::
decrement()
{
  CGAL_precondition(pos != Face_handle());
  if(s == vertex_vertex || s == vertex_edge)
  {
    if(s == vertex_vertex)
      i = cw(i);

    Orientation o;
    do
    {
      Face_handle n = pos->neighbor(ccw(i));
      i = n->index(pos);
      pos = n;

      o = _tr->orientation(p, q, pos->vertex(i)->point());
      i = ccw(i);
    }
    while(o == LEFT_TURN);

    s = (o == COLLINEAR) ? vertex_vertex : edge_vertex;

  }
  else // s == edge_edge || s == edge_vertex
  {
    // the following is not nice. A better solution is to say
    // that index i is at the vertex that is alone on one side of l(p, q)
    if(s == edge_edge)
      i = (_tr->orientation(p, q, pos->vertex(i)->point()) == LEFT_TURN) ? cw(i) : ccw(i);

    Face_handle n = pos->neighbor(i);
    i = n->index(pos);
    pos = n;
    Orientation o = _tr->orientation(p, q, pos->vertex(i)->point());

    s = (o == COLLINEAR) ? vertex_edge : edge_edge;
  }
}

template <typename Triangulation>
inline
Triangulation_sphere_line_face_circulator_2<Triangulation>&
Triangulation_sphere_line_face_circulator_2<Triangulation>::
operator++()
{
  CGAL_precondition(pos != Face_handle()) ;
  increment();
  return *this;
}

template <typename Triangulation>
inline
Triangulation_sphere_line_face_circulator_2<Triangulation>&
Triangulation_sphere_line_face_circulator_2<Triangulation>::
operator--()
{
  CGAL_precondition(pos != Face_handle()) ;
  decrement();
  return *this;
}

template <typename Triangulation>
inline
Triangulation_sphere_line_face_circulator_2<Triangulation>
Triangulation_sphere_line_face_circulator_2<Triangulation>::
operator++(int)
{
  Line_face_circulator tmp(*this);
  ++(*this);
  return tmp;
}

template <typename Triangulation>
inline
Triangulation_sphere_line_face_circulator_2<Triangulation>
Triangulation_sphere_line_face_circulator_2<Triangulation>::
operator--(int)
{
  Line_face_circulator tmp(*this);
  --(*this);
  return tmp;
}

template <typename Triangulation>
inline bool
Triangulation_sphere_line_face_circulator_2<Triangulation>::
operator==(const Line_face_circulator& lfc) const
{
  CGAL_precondition(pos != Face_handle() && lfc.pos != Face_handle());
  return (pos == lfc.pos &&  _tr == lfc._tr &&
          s == lfc.s && p==lfc.p && q==lfc.q);
}

template <typename Triangulation>
inline bool
Triangulation_sphere_line_face_circulator_2<Triangulation>::
operator!=(const Line_face_circulator& lfc) const
{
  return !(*this == lfc);
}

template <typename Triangulation>
inline bool
Triangulation_sphere_line_face_circulator_2<Triangulation>::
is_empty() const
{
  return pos == Face_handle();
}

template <typename Triangulation>
inline bool
Triangulation_sphere_line_face_circulator_2<Triangulation>::
operator==(Nullptr_t CGAL_assertion_code(n)) const
{
  CGAL_triangulation_assertion(n == NULL);
  return pos == Face_handle();
}

template <typename Triangulation>
inline bool
Triangulation_sphere_line_face_circulator_2<Triangulation>::
operator!=(Nullptr_t n) const
{
  CGAL_triangulation_assertion(n == NULL);
  return !(*this == n);
}

} // end namespace CGAL

#endif // CGAL_TRIANGULATION_SPHERE_LINE_FACE_CIRCULATOR_2_H
