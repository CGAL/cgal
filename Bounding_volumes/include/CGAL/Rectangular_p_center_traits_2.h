// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
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
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_RECTANGULAR_P_CENTER_TRAITS_2_H
#define CGAL_RECTANGULAR_P_CENTER_TRAITS_2_H 1

#include <CGAL/license/Bounding_volumes.h>


#include <CGAL/Point_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/pierce_rectangles_2.h>
#include <CGAL/Optimisation/assertions.h>

namespace CGAL {

template < class A, class S >
struct Select : public std::binary_function< A, A, A > {
  Select() {}
  Select(const S& s) : s_(s) {}
  A operator()(const A& a, const A& b) const
  { return s_(a, b) ? a : b; }
  A operator()(const A& a, const A& b)
  { return s_(a, b) ? a : b; }
protected:
  S s_;
};

template < class R >
struct I_Signed_x_distance_2
: public std::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const
  { return q1.x() - q2.x(); }
};
template < class R >
struct I_Signed_y_distance_2
: public std::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const
  { return q1.y() - q2.y(); }
};
template < class R >
struct I_Infinity_distance_2
: public std::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const {
    return (std::max)(CGAL_NTS abs(q1.x() - q2.x()),
                      CGAL_NTS abs(q1.y() - q2.y()));
  }
};

template < class R >
struct I_Signed_infinity_distance_2
: public std::binary_function<
  Point_2< R >, Point_2< R >, typename R::FT >
{
  typename R::FT
  operator()(const Point_2< R >& q1, const Point_2< R >& q2) const
  { 
    return (std::max)(q1.x() - q2.x(), q1.y() - q2.y()); 
  }
};

template < class R >
struct I_Construct_point_2_above_right_implicit_point_2 {
  // (p, q, r) |--> (p.x() + r, q.y() + r)
  typedef typename R::FT                     FT;
  typedef typename R::Point_2                P;

  inline P
  operator()(const P& p, const P& q, FT r) const
  { return P(p.x() + r, q.y() + r); }
};

template < class R >
struct I_Construct_point_2_above_left_implicit_point_2 {
  // (p, q, r) |--> (p.x() - r, q.y() + r)
  typedef typename R::FT                     FT;
  typedef typename R::Point_2                P;

  inline P
  operator()(const P& p, const P& q, FT r) const
  { return P(p.x() - r, q.y() + r); }
};

template < class R >
struct I_Construct_point_2_below_left_implicit_point_2 {
  // (p, q, r) |--> (p.x() - r, q.y() - r)
  typedef typename R::FT                     FT;
  typedef typename R::Point_2                P;

  inline P
  operator()(const P& p, const P& q, FT r) const
  { return P(p.x() - r, q.y() - r); }
};

template < class R >
struct I_Construct_point_2_below_right_implicit_point_2 {
  // (p, q, r) |--> (p.x() + r, q.y() - r)
  typedef typename R::FT                     FT;
  typedef typename R::Point_2                P;

  inline P
  operator()(const P& p, const P& q, FT r) const
  { return P(p.x() + r, q.y() - r); }
};


template < class R >
struct Rectangular_p_center_default_traits_2 : public R
{
  // -----------------------------------------------------------------
  // types:
  //
  typedef typename R::FT               FT;
  typedef typename R::Point_2          Point_2;
  typedef typename R::Iso_rectangle_2  Iso_rectangle_2;

  // -----------------------------------------------------------------
  // predicates:
  //

  typedef typename R::Less_x_2         Less_x_2;
  typedef typename R::Less_y_2         Less_y_2;

  // -----------------------------------------------------------------
  // constructions:
  //

  // from the kernel
  typedef typename R::Construct_iso_rectangle_2 Construct_iso_rectangle_2;
  typedef typename R::Construct_vertex_2        Construct_vertex_2;

  // additions
  typedef I_Signed_x_distance_2< R >  Signed_x_distance_2;
  typedef I_Signed_y_distance_2< R >  Signed_y_distance_2;

  typedef I_Infinity_distance_2< R >           Infinity_distance_2;
  typedef I_Signed_infinity_distance_2< R >    Signed_infinity_distance_2;

  typedef I_Construct_point_2_above_right_implicit_point_2< R >
    Construct_point_2_above_right_implicit_point_2;
  typedef I_Construct_point_2_above_left_implicit_point_2< R >
    Construct_point_2_above_left_implicit_point_2;
  typedef I_Construct_point_2_below_right_implicit_point_2< R >
    Construct_point_2_below_right_implicit_point_2;
  typedef I_Construct_point_2_below_left_implicit_point_2< R >
    Construct_point_2_below_left_implicit_point_2;

  // get object methods:
  Signed_x_distance_2
  signed_x_distance_2_object() const
  { return Signed_x_distance_2(); }
  Signed_y_distance_2
  signed_y_distance_2_object() const
  { return Signed_y_distance_2(); }
  Infinity_distance_2
  infinity_distance_2_object() const
  { return Infinity_distance_2(); }
  Signed_infinity_distance_2
  signed_infinity_distance_2_object() const
  { return Signed_infinity_distance_2(); }
  Construct_point_2_above_right_implicit_point_2
  construct_point_2_above_right_implicit_point_2_object() const
  { return Construct_point_2_above_right_implicit_point_2(); }
  Construct_point_2_above_left_implicit_point_2
  construct_point_2_above_left_implicit_point_2_object() const
  { return Construct_point_2_above_left_implicit_point_2(); }
  Construct_point_2_below_left_implicit_point_2
  construct_point_2_below_left_implicit_point_2_object() const
  { return Construct_point_2_below_left_implicit_point_2(); }
  Construct_point_2_below_right_implicit_point_2
  construct_point_2_below_right_implicit_point_2_object() const
  { return Construct_point_2_below_right_implicit_point_2(); }

}; // Rectangular_p_center_default_traits_2

template < class Traits_, class PiercingFunction_ >
struct Rectangular_p_center_matrix_search_traits_2 {
  typedef Traits_                        Traits;
  typedef typename Traits::FT            FT;
  typedef typename Traits::Point_2       Point_2;
  typedef PiercingFunction_              PiercingFunction;
  typedef Staircases< Traits >           LD;
  typedef typename LD::size_type         size_type;

  template < class InputIC >
  Rectangular_p_center_matrix_search_traits_2(
    InputIC f, InputIC l, Traits t, const PiercingFunction& p)
  : ld(f, l, t), pf(p)
#if (!defined(CGAL_OPTIMISATION_NO_ASSERTIONS) && \
  !defined(CGAL_NO_ASSERTIONS) && !defined(NDEBUG))
    , ld_size(ld.size())
#endif // optimisation_assertion_code
  {}

  size_type number_of_points() const { return ld.pts.size(); }

  bool operator()(FT v)
  {
    CGAL_optimisation_assertion(ld.size() == ld_size);
    ld.r = v / FT(2);
    bool ok;
    pf(ld, Wastebasket< Point_2 >(), ok);
    CGAL_optimisation_assertion(ld.size() == ld_size);
    return ok;
  }

  template < class OutputIterator >
  OutputIterator operator()(FT v, OutputIterator o, bool& ok)
  {
    CGAL_optimisation_assertion(ld.size() == ld_size);
    ld.r = v / FT(2);
    OutputIterator n = pf(ld, o, ok);
    CGAL_optimisation_assertion(ld.size() == ld_size);
    return n; //pf(ld, o, ok);
  }

protected:
  // data members:
  LD                 ld;
  PiercingFunction   pf;
  CGAL_optimisation_assertion_code(typename LD::size_type ld_size;)

  // copying this would be too inefficient
  Rectangular_p_center_matrix_search_traits_2(
    const Rectangular_p_center_matrix_search_traits_2&)
  {}

}; // Rectangular_p_center_matrix_search_traits_2< ... >

template < class ForwardIterator, class Traits >
typename Traits::Iso_rectangle_2
bounding_box_2(ForwardIterator f, ForwardIterator l, const Traits& t)
// PRE: f != l.
{
  CGAL_precondition(f != l);
  typedef typename Traits::Less_x_2                  Less_x_2;
  typedef typename Traits::Less_y_2                  Less_y_2;
  typedef typename Traits::Construct_iso_rectangle_2 Rect;
  typedef typename Traits::Construct_vertex_2        CVertex;

  Less_x_2 lessx = t.less_x_2_object();
  Less_y_2 lessy = t.less_y_2_object();
  Rect     rect  = t.construct_iso_rectangle_2_object();
  CVertex  v     = t.construct_vertex_2_object();

  ForwardIterator xmin = f;
  ForwardIterator xmax = f;
  ForwardIterator ymin = f;
  ForwardIterator ymax = f;

  while (++f != l) {
    if (lessx(*f, *xmin)) xmin = f;
    if (lessx(*xmax, *f)) xmax = f;
    if (lessy(*f, *ymin)) ymin = f;
    if (lessy(*ymax, *f)) ymax = f;
  }

  return rect(v(rect(*xmin, *ymin), 0), v(rect(*xmax, *ymax), 2));
} // bounding_box_2(f, l, t)
template < class ForwardIterator >
inline typename
std::iterator_traits< ForwardIterator >::value_type::R::Iso_rectangle_2
bounding_box_2(ForwardIterator f, ForwardIterator l)
// PRE: f != l.
{
  CGAL_precondition(f != l);
  // that is how it is supposed to be ...
  //typedef typename std::iterator_traits< ForwardIterator >::value_type::R
  //  Traits;
  typedef typename std::iterator_traits< ForwardIterator >::value_type::R R;
  typedef Rectangular_p_center_default_traits_2< R > Traits;
  Traits t;
  return bounding_box_2(f, l, t);
} // bounding_box_2(f, l)

template < class Rectangle, class Traits >
inline Rectangle
construct_bounding_box_union_2(const Rectangle& r1,
                               const Rectangle& r2,
                               const Traits& t)
{
  typedef typename Traits::Construct_iso_rectangle_2  Rect;
  typedef typename Traits::Construct_vertex_2         CVertex;
  typedef typename Traits::Less_x_2                   Less_x_2;
  typedef typename Traits::Less_y_2                   Less_y_2;

  Less_x_2 lessx = t.less_x_2_object();
  Less_y_2 lessy = t.less_y_2_object();
  Rect     rect  = t.construct_iso_rectangle_2_object();
  CVertex  v     = t.construct_vertex_2_object();

  return rect(
    v(rect(lessx(v(r1, 0), v(r2, 0)) ? v(r1, 0) : v(r2, 0),
           lessy(v(r1, 0), v(r2, 0)) ? v(r1, 0) : v(r2, 0)),
      0),
    v(rect(lessx(v(r2, 2), v(r1, 2)) ? v(r1, 2) : v(r2, 2),
           lessy(v(r2, 2), v(r1, 2)) ? v(r1, 2) : v(r2, 2)),
      2));
} // construct_bounding_box_union_2(r1, r2, t)

template < class Rectangle >
inline Rectangle
construct_bounding_box_union_2(const Rectangle& r1, const Rectangle& r2)
{
  typename Rectangle::R t;
  return construct_bounding_box_union_2(r1, r2, t);
} // construct_bounding_box_union_2(r1, r2)

} //namespace CGAL

#endif // ! (CGAL_RECTANGULAR_P_CENTER_TRAITS_2_H)
