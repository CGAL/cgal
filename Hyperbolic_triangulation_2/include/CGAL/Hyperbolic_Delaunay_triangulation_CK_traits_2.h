// Copyright (c) 2010-2018   INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mikhail Bogdanov
//                 Monique Teillaud <Monique.Teillaud@inria.fr>

#ifndef CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_CK_TRAITS_2_H
#define CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_CK_TRAITS_2_H

#include <CGAL/license/Hyperbolic_triangulation_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2/Intersection_traits.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/determinant.h>

#include <boost/tuple/tuple.hpp>
#include <boost/variant.hpp>

#include <CGAL/internal/Hyperbolic_Delaunay_triangulation_traits_2_functions.h>

namespace CGAL {

namespace internal {

  template <typename Traits>
  class Construct_hyperbolic_circumcenter_CK_2
  {
    typedef typename Traits::Hyperbolic_Voronoi_point_2     Hyperbolic_Voronoi_point_2;
    typedef typename Traits::Hyperbolic_point_2             Hyperbolic_point_2;
    typedef typename Traits::Euclidean_circle_or_line_2     Euclidean_circle_or_line_2;
    typedef typename Traits::Euclidean_line_2               Euclidean_line_2;
    typedef typename Traits::Circle_2                       Circle_2;
    typedef typename Traits::FT                             FT;
    typedef typename Traits::Circular_arc_point_2           Circular_arc_point_2;

  public:
    Construct_hyperbolic_circumcenter_CK_2(const Traits gt = Traits()): _gt(gt) {}

    Hyperbolic_Voronoi_point_2 operator()(const Hyperbolic_point_2& p,
                                          const Hyperbolic_point_2& q,
                                          const Hyperbolic_point_2& r) const
    {

      Construct_circle_or_line_supporting_bisector<Traits> cclsb(_gt);

      Hyperbolic_point_2 po(CGAL::ORIGIN);
      Circle_2 l_inf(po, FT(1));

      if( _gt.compare_distance_2_object()(po, p, q) == EQUAL &&
          _gt.compare_distance_2_object()(po, p, r) == EQUAL)
        return po;

      // now supporting objects cannot both be Euclidean lines

      Euclidean_circle_or_line_2 bis_pq = cclsb(p, q);
      Euclidean_circle_or_line_2 bis_qr = cclsb(q, r);

      std::pair<Circular_arc_point_2, unsigned> pair;
      Euclidean_line_2* l;
      Circle_2* c;

      if(Circle_2* c_pq = boost::get<Circle_2>(&bis_pq))
      {
        if(Circle_2* c_qr = boost::get<Circle_2>(&bis_qr))
        {
          typedef typename CK2_Intersection_traits<Traits, Circle_2, Circle_2>::type Intersection_result;
          std::vector< Intersection_result > inters;
          intersection(*c_pq, *c_qr, std::back_inserter(inters));

          CGAL_triangulation_assertion(assign(pair, inters[0]));
          if(pair.second == 1)
          {
            if(_gt.has_on_bounded_side_2_object()(l_inf, pair.first))
              return pair.first;

            CGAL_triangulation_assertion(assign(pair, inters[1]));
            return pair.first;
          }
          return pair.first;
        }

        // here bis_qr is a line
        l = boost::get<Euclidean_line_2>(&bis_qr);
        c = c_pq;
      }
      else
      {
        // here bis_pq is a line, and bis_qr is necessarily a circle
        l = boost::get<Euclidean_line_2>(&bis_pq);
        c = boost::get<Circle_2>(&bis_qr);
      }

      typedef typename CK2_Intersection_traits<Traits, Euclidean_line_2, Circle_2>::type Intersection_result;
      std::vector< Intersection_result > inters;
      intersection(*l, *c, std::back_inserter(inters));

      CGAL_triangulation_assertion(assign(pair,inters[0]));
      if(pair.second == 1)
      {
        if(_gt.has_on_bounded_side_2_object()(l_inf, pair.first))
          return pair.first;

        CGAL_triangulation_assertion(assign(pair, inters[1]));
        return pair.first;
      }
      return pair.first;
    }

  private:
    const Traits& _gt;
  }; // end Construct_hyperbolic_circumcenter_2


  template <typename Traits>
  class Construct_hyperbolic_bisector_CK_2
  {
    typedef typename Traits::Hyperbolic_segment_2         Hyperbolic_segment_2;
    typedef typename Traits::Hyperbolic_point_2           Hyperbolic_point_2;
    typedef typename Traits::Circle_2                     Circle_2;
    typedef typename Traits::Euclidean_line_2             Euclidean_line_2;
    typedef typename Traits::Euclidean_circle_or_line_2   Euclidean_circle_or_line_2;
    typedef typename Traits::Circular_arc_2               Circular_arc_2;
    typedef typename Traits::Line_arc_2                   Line_arc_2;
    typedef typename Traits::Circular_arc_point_2         Circular_arc_point_2;
    typedef typename Traits::FT                           FT;

  public:
    Construct_hyperbolic_bisector_CK_2(const Traits& gt = Traits()) : _gt(gt) {}

    // constructs a hyperbolic line
    Hyperbolic_segment_2 operator()(const Hyperbolic_point_2& p,
                                    const Hyperbolic_point_2& q) const
    {

      Construct_circle_or_line_supporting_bisector<Traits>  cclsb(_gt);

      Hyperbolic_point_2 po(CGAL::ORIGIN);
      Circle_2 l_inf = Circle_2(po, FT(1));

      if(_gt.compare_distance_2_object()(po, p, q) == EQUAL)
      {
        Euclidean_line_2 l = _gt.construct_Euclidean_bisector_2_object()(p, q);
        if(_gt.less_y_2_object()(p, q))
          return Line_arc_2(l, l_inf, false, l_inf, true);

        return Line_arc_2(l, l_inf, true, l_inf, false);
      }

      Euclidean_circle_or_line_2 bis_pq = cclsb(p, q);
      Circle_2* c = boost::get<Circle_2>(&bis_pq);

      if(_gt.less_y_2_object()(po, c->center()))
        return Circular_arc_2(*c, l_inf, true, l_inf, false);
      else if(_gt.less_y_2_object()(c->center(), po))
        return Circular_arc_2(*c, l_inf, false, l_inf, true);

      // the center of the circle is on the x-axis
      if(_gt.less_x_2_object()(po, c->center())) { return Circular_arc_2(*c, l_inf, true, l_inf, false); }
      return Circular_arc_2(*c, l_inf, false, l_inf, true);
    }

    // constructs the hyperbolic bisector of segment [p, q] limited by
    // circumcenter(p, q, r) on one side
    // and circumcenter(p, s, q) on the other side
    Hyperbolic_segment_2 operator()(const Hyperbolic_point_2& p,
                                    const Hyperbolic_point_2& q,
                                    const Hyperbolic_point_2& r,
                                    const Hyperbolic_point_2& s) const
    {
      CGAL_triangulation_precondition((_gt.orientation_2_object()(p, q, r) == ON_POSITIVE_SIDE) &&
                                      (_gt.orientation_2_object()(p, s, q) == ON_POSITIVE_SIDE));
      CGAL_triangulation_precondition((_gt.side_of_oriented_circle_2_object()(p, q, r,s) == ON_NEGATIVE_SIDE) &&
                                      (_gt.side_of_oriented_circle_2_object()(p, s, q, r) == ON_NEGATIVE_SIDE));

      Construct_circle_or_line_supporting_bisector<Traits>  cclsb(_gt);
      Construct_hyperbolic_circumcenter_CK_2<Traits>        chc(_gt);
      Hyperbolic_point_2 po(CGAL::ORIGIN);

      // TODO MT this is non-optimal...
      // the bisector is already computed here
      // and it will be recomputed below
      Circular_arc_point_2 a = chc(p, q, r);
      Circular_arc_point_2 b = chc(p, s, q);

      if(_gt.compare_distance_2_object()(po, p, q) == EQUAL)
      {
        Euclidean_line_2 l = _gt.construct_Euclidean_bisector_2_object()(p, q);
        return Line_arc_2(l, a, b);
      }

      Euclidean_circle_or_line_2
          bis_pq = cclsb(p, q);
      Circle_2* c_pq = boost::get<Circle_2>(&bis_pq);

      if(_gt.compare_distance_2_object()(po, p, q) == POSITIVE)
      {
        // then p is inside the supporting circle
        return Circular_arc_2(*c_pq, b, a);
      }

      return Circular_arc_2(*c_pq, a, b);
    }

    // constructs the hyperbolic bisector of segment [p, q]
    // limited by hyperbolic circumcenter(p, q, r) on one side
    // and going to the infinite line on the other side
    Hyperbolic_segment_2 operator()(const Hyperbolic_point_2& p,
                                    const Hyperbolic_point_2& q,
                                    const Hyperbolic_point_2& r) const
    {
      CGAL_triangulation_precondition(_gt.orientation_2_object()(p, q, r) == POSITIVE);

      Construct_circle_or_line_supporting_bisector<Traits>  cclsb(_gt);
      Construct_hyperbolic_circumcenter_CK_2<Traits>        chc(_gt);

      Hyperbolic_point_2 po(CGAL::ORIGIN);
      Circle_2 l_inf(po, FT(1));

      // TODO MT this is non-optimal...
      // the bisector is computed (at least) twice
      Circular_arc_point_2 a = chc(p, q, r);

      if(_gt.compare_distance_2_object()(po, p, q) == EQUAL)
      {
        Euclidean_line_2 bis_pq = _gt.construct_Euclidean_bisector_2_object()(p, q);
        typedef typename
        CK2_Intersection_traits<Traits, Euclidean_line_2, Circle_2>::type
            Intersection_result;
        std::vector< Intersection_result > inters;
        intersection(bis_pq, l_inf, std::back_inserter(inters));
        std::pair<Circular_arc_point_2, unsigned> pair;

        CGAL_triangulation_assertion(assign(pair,inters[0]));
        CGAL_triangulation_assertion(pair.second == 1);
        if(_gt.less_y_2_object()(p, q))
          return Line_arc_2(bis_pq,a,pair.first);

        CGAL_triangulation_assertion(assign(pair,inters[1]));
        CGAL_triangulation_assertion(pair.second == 1);
        return Line_arc_2(bis_pq,a,pair.first);
      }

      Euclidean_circle_or_line_2 bis_pq = cclsb(p, q);
      Circle_2* c_pq = boost::get<Circle_2>(&bis_pq);

      Hyperbolic_point_2 approx_a(to_double(a.x()),to_double(a.y()));

      typedef typename
      CK2_Intersection_traits<Traits, Circle_2, Circle_2>::type            Intersection_result;
      std::vector< Intersection_result > inters;
      intersection(*c_pq, l_inf, std::back_inserter(inters));
      std::pair<Circular_arc_point_2, unsigned> pair;

      CGAL_triangulation_assertion(assign(pair,inters[0]));
      CGAL_triangulation_assertion(pair.second == 1);

      Hyperbolic_point_2 approx_pinf(to_double(pair.first.x()), to_double(pair.first.y()));
      Hyperbolic_point_2 approx_c(to_double(c_pq->center().x()),
                                  to_double(c_pq->center().y()));
      if(_gt.orientation_2_object()(p, q,approx_pinf) == NEGATIVE)
      {
        if(_gt.orientation_2_object()(approx_c,approx_a,approx_pinf) == POSITIVE)
          return Circular_arc_2(*c_pq, a, pair.first);

        return Circular_arc_2(*c_pq, pair.first, a);
      }

      CGAL_triangulation_assertion(assign(pair,inters[1]));
      if(_gt.orientation_2_object()(approx_c,approx_a,approx_pinf) == POSITIVE)
        return Circular_arc_2(*c_pq, pair.first, a);

      return Circular_arc_2(*c_pq, a, pair.first);
    }

  private:
    const Traits& _gt;
  }; // end Construct_hyperbolic_bisector_2


} // end namespace internal





template<class R = CGAL::Circular_kernel_2<
           CGAL::Exact_predicates_inexact_constructions_kernel,
           CGAL::Algebraic_kernel_for_circles_2_2<
             CGAL::Exact_predicates_inexact_constructions_kernel::RT >
           >
         >
class Hyperbolic_Delaunay_triangulation_CK_traits_2
    : public R // R is supposed to be a model of CircularKernel2
{
  typedef Hyperbolic_Delaunay_triangulation_CK_traits_2<R>    Self;
  typedef R                                                   Base;

public:
  typedef typename R::FT                                      FT;

  typedef typename R::Point_2                                 Hyperbolic_point_2;
  typedef typename R::Circle_2                                Circle_2;
  typedef typename R::Line_2                                  Euclidean_line_2;
  typedef boost::variant<Circle_2, Euclidean_line_2>          Euclidean_circle_or_line_2;

  typedef typename R::Circular_arc_2                          Circular_arc_2;
  typedef typename R::Line_arc_2                              Line_arc_2;
  typedef typename R::Circular_arc_point_2                    Circular_arc_point_2;
  typedef Circular_arc_point_2                                Hyperbolic_Voronoi_point_2;
  typedef typename R::Segment_2                               Euclidean_segment_2; // only used internally here
  typedef boost::variant<Circular_arc_2, Line_arc_2>          Hyperbolic_segment_2;

  typedef typename R::Triangle_2                              Hyperbolic_triangle_2;

  typedef typename R::Compare_x_2                             Compare_x_2;
  typedef typename R::Compare_y_2                             Compare_y_2;
  typedef typename R::Compare_distance_2                      Compare_distance_2;
  typedef typename R::Orientation_2                           Orientation_2;
  typedef typename R::Side_of_oriented_circle_2               Side_of_oriented_circle_2;

  // the following types are only used internally in this traits class,
  // so they need not be documented, and they don't need _object()
  typedef typename R::Collinear_2                             Euclidean_collinear_2;
  typedef typename R::Construct_bisector_2                    Construct_Euclidean_bisector_2;
  typedef typename R::Construct_midpoint_2                    Construct_Euclidean_midpoint_2;
  typedef typename R::Compute_squared_distance_2              Compute_squared_Euclidean_distance_2;
  typedef typename R::Has_on_bounded_side_2                   Has_on_bounded_side_2;

  typedef typename R::Less_x_2                                Less_x_2;
  typedef typename R::Less_y_2                                Less_y_2;

  // only kept for demo to please T2graphicsitems
  typedef Euclidean_segment_2                                 Line_segment_2;
  typedef Hyperbolic_segment_2                                Segment_2;

  typedef internal::Construct_hyperbolic_circumcenter_CK_2<Self>   Construct_hyperbolic_circumcenter_2;
  typedef internal::Construct_hyperbolic_bisector_CK_2<Self>       Construct_hyperbolic_bisector_2;
  typedef internal::Is_Delaunay_hyperbolic<Self>                Is_Delaunay_hyperbolic;
  typedef internal::Side_of_oriented_hyperbolic_segment_2<Self> Side_of_oriented_hyperbolic_segment_2;
  typedef typename internal::Construct_circle_or_line_supporting_bisector<Self> Construct_circle_or_line_supporting_bisector;
  typedef internal::Construct_hyperbolic_segment_2<Self>        Construct_hyperbolic_segment_2;
  typedef typename Base::Construct_segment_2                    Construct_segment_2;

public:
  Hyperbolic_Delaunay_triangulation_CK_traits_2(const Base& kernel = Base()) : Base(kernel) {}


  Construct_hyperbolic_segment_2
  construct_hyperbolic_segment_2_object() const
  { return Construct_hyperbolic_segment_2(*this); }

  Construct_segment_2
  construct_segment_2_object() const
  { return this->Base::construct_segment_2_object(); }

  Construct_hyperbolic_circumcenter_2
  construct_hyperbolic_circumcenter_2_object() const
  { return Construct_hyperbolic_circumcenter_2(*this); }

  Construct_hyperbolic_bisector_2
  construct_hyperbolic_bisector_2_object() const
  { return Construct_hyperbolic_bisector_2(*this); }

  Is_Delaunay_hyperbolic
  is_Delaunay_hyperbolic_2_object() const
  { return Is_Delaunay_hyperbolic(*this); }

  Side_of_oriented_hyperbolic_segment_2
  side_of_oriented_hyperbolic_segment_2_object() const
  { return Side_of_oriented_hyperbolic_segment_2(*this); }

  Construct_Euclidean_bisector_2
  construct_Euclidean_bisector_2_object() const
  { return this->Base::construct_bisector_2_object(); }

  Construct_circle_or_line_supporting_bisector
  construct_circle_or_line_supporting_bisector_2_object() const
  { return Construct_circle_or_line_supporting_bisector(*this); }

  Euclidean_collinear_2
  euclidean_collinear_2_object() const
  { return this->Base::collinear_2_object(); }

  Compute_squared_Euclidean_distance_2
  compute_squared_Euclidean_distance_2_object() const
  { return this->Base::compute_squared_distance_2_object(); }

};

// Take out the code below to some separate file

#ifdef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
template  <>
struct Triangulation_structural_filtering_traits< Hyperbolic_Delaunay_triangulation_CK_traits_2<Epeck> > {
  typedef Tag_true Use_structural_filtering_tag;
};
#endif // CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H

#ifdef CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
template <>
struct Triangulation_structural_filtering_traits< Hyperbolic_Delaunay_triangulation_CK_traits_2<Epick> > {
  typedef Tag_true Use_structural_filtering_tag;
};
#endif // CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H

} //namespace CGAL

#endif // CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_CK_TRAITS_2_H
