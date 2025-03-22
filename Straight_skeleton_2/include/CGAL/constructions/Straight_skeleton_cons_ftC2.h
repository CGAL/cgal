// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_CONS_FTC2_H
#define CGAL_STRAIGHT_SKELETON_CONS_FTC2_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>

#include <CGAL/Lazy.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/FPU.h>

#include <boost/optional/optional.hpp>

#include <cstddef>

namespace CGAL {

namespace CGAL_SS_i {

template<class S>
bool are_edges_collinear( S const& e0, S const& e1 )
{
  return   ((e1.source() == e0.source()) || (e1.source() == e0.target()) || collinear(e0.source(),e0.target(),e1.source()))
    && ( (e1.target() == e0.source()) || (e1.target() == e0.target()) || (collinear(e0.source(),e0.target(),e1.target()))) ;
}

template<class S>
inline
bool are_parallel_edges_equally_oriented( S const& e0, S const& e1 )
{
  return angle(e0.source(), e0.target(),
               e1.source(), e1.target()) == ACUTE;
}

template<class S>
bool are_edges_orderly_collinear( S const& e0, S const& e1 )
{
  return are_edges_collinear(e0,e1) && are_parallel_edges_equally_oriented(e0,e1);
}


template<class K, class Caches>
Trisegment_collinearity trisegment_collinearity_no_exact_constructions(const Segment_2_with_ID<K>& e0,
                                                                       const Segment_2_with_ID<K>& e1,
                                                                       const Segment_2_with_ID<K>& e2,
                                                                       Caches& caches)
{
  // 'are_edges_orderly_collinear()' is used to harmonize coefficients, but if the kernel is inexact
  // we could also have that are_edges_orderly_collinear() returns false, but the computed coefficients
  // are identical. In that case, we want to return that there is a collinearity, otherwise the internal
  // computations (even the exact ones) will fail.
  boost::optional<typename K::Line_2> l0 = compute_normalized_line_coeffC2(e0, caches);
  boost::optional<typename K::Line_2> l1 = compute_normalized_line_coeffC2(e1, caches);
  boost::optional<typename K::Line_2> l2 = compute_normalized_line_coeffC2(e2, caches);

  bool is_01 = (l0->a() == l1->a()) && (l0->b() == l1->b()) && (l0->c() == l1->c());
  bool is_02 = (l0->a() == l2->a()) && (l0->b() == l2->b()) && (l0->c() == l2->c());
  bool is_12 = (l1->a() == l2->a()) && (l1->b() == l2->b()) && (l1->c() == l2->c());

  CGAL_STSKEL_TRAITS_TRACE("collinearity: [" << e0.id() << ", " << e1.id() << ", " << e2.id()
                        << "]: " << is_01 << " " << is_02 << " " << is_12);

  if ( is_01 & !is_02 & !is_12 )
    return TRISEGMENT_COLLINEARITY_01;
  else if ( is_02 & !is_01 & !is_12 )
    return TRISEGMENT_COLLINEARITY_02;
  else if ( is_12 & !is_01 & !is_02 )
    return TRISEGMENT_COLLINEARITY_12;
  else if ( !is_01 & !is_02 & !is_12  )
    return TRISEGMENT_COLLINEARITY_NONE;
  else
    return TRISEGMENT_COLLINEARITY_ALL;
}

///

// Attempted to use std::hypot (https://github.com/CGAL/cgal/commit/a1845691d5d8055978662cd95059c6d3f94c17a2)
// but did not notice any gain, and even observed some regressions in the tests.

template <typename NT>
typename Coercion_traits<double, NT>::Type
inexact_sqrt_implementation(const NT& n, CGAL::Null_functor)
{
  typedef CGAL::Interval_nt<false> IFT;
  typename IFT::Protector protector;

  CGAL::NT_converter<NT, IFT> to_ift;
  IFT sqrt_ift = sqrt(to_ift(n));
  CGAL_STSKEL_TRAITS_TRACE("sqrt's interval " << sqrt_ift.inf() << " " << sqrt_ift.sup() ) ;
  CGAL_STSKEL_TRAITS_TRACE("interval delta " << sqrt_ift.sup() - sqrt_ift.inf() ) ;

  return NT(to_double(sqrt_ift));
}

template <typename NT, typename Sqrt>
typename Sqrt::result_type
inexact_sqrt_implementation(const NT& nt, Sqrt sqrt)
{
  CGAL_STSKEL_TRAITS_TRACE("sqrt(" << typeid(NT).name() << ")");
  return sqrt(nt);
}

template <typename NT>
decltype(auto) inexact_sqrt(const NT& nt)
{
  // the initial version of this function was using Algebraic_category
  // for the dispatch but some ring type (like Gmpz) provides a Sqrt
  // functor even if not being Field_with_sqrt.
  typedef CGAL::Algebraic_structure_traits<NT> AST;
  typedef typename AST::Sqrt Sqrt;
  return inexact_sqrt_implementation(nt, Sqrt());
}

template <typename NT>
Quotient<NT>
inexact_sqrt(const Quotient<NT>& q)
{
  return { inexact_sqrt(q.numerator()*q.denominator()), abs(q.denominator()) };
}

template <typename NT>
Lazy_exact_nt<NT>
inexact_sqrt(const Lazy_exact_nt<NT>& lz)
{
  return inexact_sqrt(exact(lz));
}

// Given an oriented 2D straight line segment 'e', computes the normalized coefficients (a,b,c)
// of the supporting line, and weights them with 'aWeight'.
//
// POSTCONDITION: [a,b] is the leftward normal vector.
// POSTCONDITION: In case of overflow, an empty optional<> is returned.
template<class K>
boost::optional< typename K::Line_2> compute_normalized_line_coeffC2( Segment_2<K> const& e )
{
  typedef typename K::FT FT ;

  bool finite = true ;
  FT a(0), b(0), c(0) ;

  if(e.source().y() == e.target().y())
  {
    a = 0 ;
    if(e.target().x() > e.source().x())
    {
      b = 1;
      c = -e.source().y();
    }
    else if(e.target().x() == e.source().x())
    {
      b = 0;
      c = 0;
    }
    else
    {
      b = -1;
      c = e.source().y();
    }

    CGAL_STSKEL_TRAITS_TRACE("HORIZONTAL line; a="<< n2str(a) << ", b=" << n2str(b) << ", c=" << n2str(c) ) ;
  }
  else if(e.target().x() == e.source().x())
  {
    b = 0;
    if(e.target().y() > e.source().y())
    {
      a = -1;
      c = e.source().x();
    }
    else if (e.target().y() == e.source().y())
    {
      a = 0;
      c = 0;
    }
    else
    {
      a = 1;
      c = -e.source().x();
    }

    CGAL_STSKEL_TRAITS_TRACE("VERTICAL line; a="<< n2str(a) << ", b=" << n2str(b) << ", c=" << n2str(c) ) ;
  }
  else
  {
    FT sa = e.source().y() - e.target().y();
    FT sb = e.target().x() - e.source().x();
    FT l2 = square(sa) + square(sb) ;

    if ( CGAL_NTS is_finite(l2) )
    {
      FT l = CGAL_SS_i::inexact_sqrt(l2);
      a = sa / l ;
      b = sb / l ;

      c = -e.source().x()*a - e.source().y()*b;

      CGAL_STSKEL_TRAITS_TRACE("GENERIC line;\nsa="<< n2str(sa) << "\nsb=" << n2str(sb)
                               << "\nnorm²=" << n2str(l2) << "\nnorm=" << n2str(l)
                               << "\na="<< n2str(a) << "\nb=" << n2str(b) << "\nc=" << n2str(c) ) ;
    }
    else
    {
      finite = false ;
    }
  }

  if ( finite )
    if ( !CGAL_NTS is_finite(a) || !CGAL_NTS is_finite(b) || !CGAL_NTS is_finite(c) )
      finite = false ;

  return cgal_make_optional( finite, K().construct_line_2_object()(a,b,c) ) ;
}

template<class K>
boost::optional< typename K::Line_2> compute_normalized_line_coeffC2(const Segment_2_with_ID<K>& e)
{
  typedef typename K::Segment_2 Segment_2 ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~ Unweighted line coefficients for E" << e.id() << " [" << typeid(typename K::FT).name() << "]" );

  return compute_normalized_line_coeffC2(static_cast<const Segment_2&>(e));
}

template<class K, class Caches>
boost::optional< typename K::Line_2 >
compute_normalized_line_coeffC2( Segment_2_with_ID<K> const& e,
                                 Caches& aCaches )
{
  typedef typename K::Line_2 Line_2 ;

  if(aCaches.mCoeff_cache.IsCached(e.mID) )
    return aCaches.mCoeff_cache.Get(e.mID) ;

  boost::optional< Line_2 > rRes = compute_normalized_line_coeffC2 ( e ) ;

  aCaches.mCoeff_cache.Set(e.mID, rRes) ;

  return rRes ;
}

template<class FT>
Rational<FT> squared_distance_from_point_to_lineC2( FT const& px, FT const& py, FT const& sx, FT const& sy, FT const& tx, FT const& ty )
{
  FT ldx = tx - sx ;
  FT ldy = ty - sy ;
  FT rdx = sx - px ;
  FT rdy = sy - py ;

  FT n = CGAL_NTS square(ldx * rdy - rdx * ldy);
  FT d = CGAL_NTS square(ldx) + CGAL_NTS square(ldy);

  return Rational<FT>(n,d) ;
}

//
// Constructs a Trisegment_2 which stores 3 oriented straight line segments e0,e1,e2 along with their collinearity.
//
// NOTE: If the collinearity cannot be determined reliably, a null trisegment is returned.
//
template<class K, class Caches>
Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > >
construct_trisegment ( Segment_2_with_ID<K> const& e0,
                       typename K::FT const& w0,
                       Segment_2_with_ID<K> const& e1,
                       typename K::FT const& w1,
                       Segment_2_with_ID<K> const& e2,
                       typename K::FT const& w2,
                       std::size_t id,
                       Caches& aCaches )
{
  typedef Trisegment_2<K, Segment_2_with_ID<K> >Trisegment_2 ;
  typedef typename Trisegment_2::Self_ptr Trisegment_2_ptr ;

  Trisegment_collinearity lCollinearity = trisegment_collinearity_no_exact_constructions(e0,e1,e2,aCaches);

  return Trisegment_2_ptr( new Trisegment_2(e0, w0, e1, w1, e2, w2, lCollinearity, id) ) ;
}

// Given 3 oriented straight line segments: e0, e1, e2
// returns the OFFSET DISTANCE (n/d) at which the offsetted lines
// intersect at a single point, IFF such intersection exist.
// If the lines intersect to the left, the returned distance is positive.
// If the lines intersect to the right, the returned distance is negative.
// If the lines do not intersect, for example, for collinear edges, or parallel edges but with the same orientation,
// returns 0 (the actual distance is undefined in this case, but 0 is a useful return)
//
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
//
// PRECONDITION: None of e0, e1 and e2 are collinear (but two of them can be parallel)
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
// NOTE: The segments (e0,e1,e2) are stored in the argument as the trisegment st.event()
//
template <class K, class Caches>
boost::optional< Rational< typename K::FT> >
compute_normal_offset_lines_isec_timeC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                          Caches& aCaches )
{
  typedef typename K::FT  FT ;
  typedef typename K::Line_2 Line_2 ;

  typedef boost::optional<Line_2> Optional_line_2 ;

  CGAL_STSKEL_TRAITS_TRACE("\ncompute_normal_offset_lines_isec_timeC2(" << tri->id() << ") [" << typeid(FT).name() << "]") ;

  FT num(0), den(0) ;

  // DETAILS:
  //
  // An offset line is given by:
  //
  //   a*x(t) + b*y(t) + c - w*t = 0
  //
  // with 't > 0' being to the left of the line.
  // If 3 such offset lines intersect at the same offset distance, the intersection 't',
  // or 'time', can be computed solving for 't' in the linear system formed by 3 such equations.
  //
  //   [[ t == ((b2*c1 - b1*c2)*a0 - (a2*c1 - a1*c2)*b0 + (a2*b1 - a1*b2)*c0)/((b2*w1 - b1*w2)*a0 - (a2*w1 - a1*w2)*b0 + (a2*b1 - a1*b2)*w0) ]]

  bool ok = false ;

  Optional_line_2 l0 = compute_normalized_line_coeffC2(tri->e0(), aCaches) ;
  Optional_line_2 l1 = compute_normalized_line_coeffC2(tri->e1(), aCaches) ;
  Optional_line_2 l2 = compute_normalized_line_coeffC2(tri->e2(), aCaches) ;

  if ( l0 && l1 && l2 )
  {
    const FT& a0 = l0->a(); const FT& b0 = l0->b(); const FT& c0 = l0->c(); const FT& w0 = tri->w0();
    const FT& a1 = l1->a(); const FT& b1 = l1->b(); const FT& c1 = l1->c(); const FT& w1 = tri->w1();
    const FT& a2 = l2->a(); const FT& b2 = l2->b(); const FT& c2 = l2->c(); const FT& w2 = tri->w2();

    CGAL_STSKEL_TRAITS_TRACE("coeffs E" << tri->e0().id() << " [" << n2str(a0) << "; " << n2str(b0) << "; " << n2str(c0) << "]"
                        << "\ncoeffs E" << tri->e1().id() << " [" << n2str(a1) << "; " << n2str(b1) << "; " << n2str(c1) << "]"
                        << "\ncoeffs E" << tri->e2().id() << " [" << n2str(a2) << "; " << n2str(b2) << "; " << n2str(c2) << "]") ;
    CGAL_STSKEL_TRAITS_TRACE("weight E" << tri->e0().id() << " [" << n2str(w0) << "] "
                          << "weight E" << tri->e1().id() << " [" << n2str(w1) << "] "
                          << "weight E" << tri->e2().id() << " [" << n2str(w2) << "]") ;

    num = (b2*c1 - b1*c2)*a0 - (a2*c1 - a1*c2)*b0 + (a2*b1 - a1*b2)*c0 ;
    den = (b2*w1 - b1*w2)*a0 - (a2*w1 - a1*w2)*b0 + (a2*b1 - a1*b2)*w0 ;

    // Here we are in a pickle: lines are not collinear but with weights, we could
    // get an infinity of solutions still (particularly easy to create with null weights)
    if (CGAL_NTS certified_is_zero(den)) {
      // If the numerator (which is a 2nd determinant within the augmented matrix) is also 0,
      // then we have a consistent system, so an infinity of solution for a fixed value of 't'
      if (CGAL_NTS certified_is_zero(num)) {
        CGAL_STSKEL_TRAITS_TRACE("Infinite solutions");

        const auto& e0 = tri->e0();
        const auto& e1 = tri->e1();
        const auto& e2 = tri->e2();

        CGAL_assertion(are_edges_parallelC2(e0, e1));
        CGAL_assertion(are_edges_parallelC2(e1, e2));

        // System is consistent, extract t
        bool w0_is_zero = CGAL_NTS certified_is_zero(w0);
        bool w1_is_zero = CGAL_NTS certified_is_zero(w1);
        bool w2_is_zero = CGAL_NTS certified_is_zero(w2);

        // If the intersection is a line, that would mean that all lines are collinear
        // and there would be two lines in the same direction, so edges would be collinear
        // and we would not be here
        CGAL_assertion(!(w0_is_zero && w1_is_zero && w2_is_zero));

        // If at least one line is not moving, evaluate a moving line at a point of the static line.
        if (w0_is_zero && w1_is_zero) {
          CGAL_assertion(a0 * a1 < 0);
          num = a2*e0.source().x() + b2*e0.source().y() + c2;
          den = w0;
        } else if (w0_is_zero && w2_is_zero) {
          CGAL_assertion(a0 * a2 < 0);
          num = a1*e0.source().x() + b1*e0.source().y() + c1;
          den = w1;
        } else if (w1_is_zero && w2_is_zero) {
          CGAL_assertion(a1 * a2 < 0);
          num = a0*e1.source().x() + b0*e1.source().y() + c0;
          den = w0;
        } else if (w0_is_zero) {
          num = a2*e0.source().x() + b2*e0.source().y() + c2;
          den = w2;
        } else if (w1_is_zero) {
          num = a1*e0.source().x() + b1*e0.source().y() + c1;
          den = w1;
        } else if (w2_is_zero) {
          num = a0*e1.source().x() + b0*e1.source().y() + c0;
          den = w0;
        } else {
          // We can't have two lines moving in the same direction with the same speeds,
          // otherwise they are either identical and we would be collinear (contradiction),
          // or there is no solution to the base linear system (also contradiction).
          // So if there is the same speed, they move in opposite directions and the only possibility
          // is at t=0
          if (CGAL_NTS certified_is_zero(w0 - w1)) {
            num = 0;
            den = 1;
          } else {
            num = c1 - c0;
            den = w0 - w1;
          }
        }
      } else {
          CGAL_assertion_msg(false, "Unexpected case: lines are not collinear, but no valid 't' could be computed.");
      }
    } else {
      CGAL_STSKEL_TRAITS_TRACE("No solution exists (inconsistent system).");
    }

    ok = CGAL_NTS is_finite(num) && CGAL_NTS is_finite(den);
  }

  CGAL_STSKEL_TRAITS_TRACE("Event time (normal): n=" << num << " d=" << den << " n/d=" << Rational<FT>(num,den));

  return cgal_make_optional(ok,Rational<FT>(num,den)) ;
}

// Given two oriented straight line segments e0 and e1 such that e-next follows e-prev, returns
// the coordinates of the midpoint of the segment between e-prev and e-next.
// NOTE: the edges can be oriented e0->e1 or e1->e0
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
boost::optional< typename K::Point_2 >
compute_oriented_midpoint ( Segment_2_with_ID<K> const& e0,
                            Segment_2_with_ID<K> const& e1 )
{
  typedef typename K::FT FT ;

  CGAL_STSKEL_TRAITS_TRACE("Computing oriented midpoint between:"
                            << "\n\te0: " << s2str(e0)
                            << "\n\te1: " << s2str(e1)
                          );

  FT delta01 = CGAL::squared_distance(e0.target(),e1.source());
  if( CGAL_NTS is_finite(delta01) && CGAL_NTS is_zero(delta01))
    return cgal_make_optional(e0.target());

  FT delta10 = CGAL::squared_distance(e1.target(),e0.source());
  if( CGAL_NTS is_finite(delta10) && CGAL_NTS is_zero(delta10))
    return cgal_make_optional(e1.target());

  bool ok = false ;
  typename K::Point_2 mp ;

  if ( CGAL_NTS is_finite(delta01) &&  CGAL_NTS is_finite(delta10) )
  {
    if ( delta01 <= delta10 )
      mp = CGAL::midpoint(e0.target(),e1.source());
    else
      mp = CGAL::midpoint(e1.target(),e0.source());

    CGAL_STSKEL_TRAITS_TRACE("\nmp=" << p2str(mp) );

    ok = CGAL_NTS is_finite(mp.x()) && CGAL_NTS is_finite(mp.y());
  }

  return cgal_make_optional(ok,mp);
}


//
// Given 3 oriented straight line segments: e0, e1, e2 and the corresponding offsetted segments: e0*, e1* and e2*,
// returns the point of the left or right seed (offset vertex) (e0*,e1*) or (e1*,e2*)
//
// If the current event (defined by e0,e1,e2) is a propagated event, that is, it follows from a previous event,
// the seeds are skeleten nodes and are given by non-null trisegments.
// If the current event is an initial event the seeds are contour vertices and are given by null trisegmets.
//
// If a seed is a skeleton node, its point has to be computed from the trisegment that defines it.
// That trisegment is exactly the trisegment tree that defined the previous event which produced the skeleton node
// (so the trisegment tree is basically a lazy representation of the seed point).
//
// If a seed is a contour vertex, its point is then simply the target endpoint of e0 or e1 (for the left/right seed).
//
// This method returns the specified seed point (left or right)
//
// NOTE: Split events involve 3 edges but only one seed, the left (that is, only e0*,e1* is connected before the event).
// The trisegment tree for a split event has always a null right child even if the event is not an initial event
// (in which case its left child won't be null).
// If you ask for the right child point for a trisegment tree corresponding to a split event you will just get e1.target()
// which is nonsensical for a non initial split event.
//
// NOTE: There is an abnormal collinearity case which occurs when e0 and e2 are collinear.
// In this case, these lines do not correspond to an offset vertex (because e0* and e2* are never consecutive before the event),
// so the degenerate seed is neither the left or the right seed. In this case, the SEED ID for the degenerate pseudo seed is UNKNOWN.
// If you request the point of such degenerate pseudo seed the oriented midpoint between e0 and e2 is returned.
//
template <class K, class Caches>
boost::optional< typename K::Point_2 >
compute_seed_pointC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                       typename Trisegment_2<K, Segment_2_with_ID<K> >::SEED_ID sid,
                       Caches& aCaches)
{
  boost::optional< typename K::Point_2 > p ;

  typedef Trisegment_2<K, Segment_2_with_ID<K> > Trisegment_2 ;

  switch ( sid )
  {
    case Trisegment_2::LEFT :

       p = tri->child_l() ? construct_offset_lines_isecC2(tri->child_l(), aCaches) // this can recurse
                          : compute_oriented_midpoint(tri->e0(),tri->e1()) ;
       break ;

    case Trisegment_2::RIGHT :

      p = tri->child_r() ? construct_offset_lines_isecC2(tri->child_r(), aCaches) // this can recurse
                         : compute_oriented_midpoint(tri->e1(),tri->e2()) ;
      break ;

    case Trisegment_2::THIRD :

      p = tri->child_t() ? construct_offset_lines_isecC2(tri->child_t(), aCaches) // this can recurse
                         : compute_oriented_midpoint(tri->e0(),tri->e2()) ;

      break ;
  }

  return p ;
}

//
// Given the trisegment tree for an event which is known to have a normal collinearity returns the seed point
// of the degenerate seed.
// A normal collinearity occurs when e0,e1 or e1,e2 are collinear.
template <class K, class Caches>
boost::optional< typename K::Point_2 >
construct_degenerate_seed_pointC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                    Caches& aCaches )
{
  return compute_seed_pointC2( tri, tri->degenerate_seed_id(), aCaches ) ;
}

template <class K, class Caches>
boost::optional< Rational< typename K::FT> >
compute_artifical_isec_timeC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                Caches& aCaches )
{
  typedef typename K::Boolean Boolean ;
  typedef typename K::FT FT ;
  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Segment_2 Segment_2 ;
  typedef typename K::Line_2 Line_2 ;
  typedef typename K::Direction_2 Direction_2 ;
  typedef typename K::Ray_2 Ray_2 ;

  typedef boost::optional<Line_2> Optional_line_2 ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~  Computing artificial isec time [" << typeid(FT).name() << "]");
  CGAL_STSKEL_TRAITS_TRACE("Event:\n" << tri);

  CGAL_precondition(tri->e0() == tri->e1());
  CGAL_precondition(bool(tri->child_l()));

  Optional_line_2 l0 = compute_normalized_line_coeffC2(tri->e0(), aCaches) ;
  if( !l0 )
    return boost::none ;

  boost::optional< typename K::Point_2 > seed = construct_offset_lines_isecC2(tri->child_l(), aCaches ) ;
  if(!seed)
    return boost::none;

  const Segment_2& contour_seg = tri->e0();
  const Direction_2 perp_dir(contour_seg.source().y() - contour_seg.target().y(),
                             contour_seg.target().x() - contour_seg.source().x());
  const Ray_2 ray(*seed, perp_dir);
  const Segment_2& opp_seg = tri->e2();
  Boolean inter_exist = K().do_intersect_2_object()(ray, opp_seg);
  if (!inter_exist) // no intersection
    return cgal_make_optional(Rational<FT>(FT(0),FT(0))) ; // does not exist

  // Compute the intersection point and evalute the time from the line equation of the contour edge
  auto inter_res = K().intersect_2_object()(ray, opp_seg);

  FT t;
  if(const Segment_2* seg = boost::get<Segment_2>(&*inter_res))
  {
    // get the segment extremity closest to the seed
    Boolean res = (K().compare_distance_2_object()(*seed, seg->source(), seg->target()) == CGAL::SMALLER);
    t = res ? l0->a() * seg->source().x() + l0->b() * seg->source().y() + l0->c()
            : l0->a() * seg->target().x() + l0->b() * seg->target().y() + l0->c();
  }
  else
  {
    const Point_2* inter_pt = boost::get<const Point_2>(&*inter_res);
    if(!CGAL_NTS is_finite(inter_pt->x()) || !CGAL_NTS is_finite(inter_pt->y()))
      return boost::none;
    t = l0->a() * inter_pt->x() + l0->b() * inter_pt->y() + l0->c() ;
  }

  bool ok = CGAL_NTS is_finite(t);
  return cgal_make_optional(ok, Rational<FT>(t, tri->w0())) ;
}

// Given 3 oriented straight line segments: e0, e1, e2
// such that two and only two of these edges are collinear, not necessarily consecutive but with the same orientaton;
// returns the OFFSET DISTANCE (n/d) at which a line perpendicular to the collinear edge passing through
// the degenerate seed point intersects the offset line of the non collinear edge
//
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
// DETAILS:
//   The offset lines are defined as:
//     a*x + b*y + c - w*t = 0
//   where (a,b,c) are the normalized line coefficients and w is the weight.
//
template <class K, class Caches>
boost::optional< Rational< typename K::FT> >
compute_degenerate_offset_lines_isec_timeC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                              Caches& aCaches )
{
  typedef typename K::FT FT ;

  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Line_2 Line_2 ;

  typedef boost::optional<Point_2> Optional_point_2 ;
  typedef boost::optional<Line_2>  Optional_line_2 ;

  if(tri->e0() == tri->e1()) // marker for artificial bisectors: they have the same face on both sides
    return compute_artifical_isec_timeC2(tri, aCaches) ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~  Computing degenerate offset lines isec time [" << typeid(FT).name() << "]");
  CGAL_STSKEL_TRAITS_TRACE("Event:\n" << tri);

  // DETAILS:
  //
  // For simplicity, assume e0,e1 are the collinear edges.
  //
  //   (1)
  //   The bisecting line of e0 and e1 is a line perpendicular to e0 (and e1)
  //   which passes through 'q', the degenerate offset vertex (e0*,e1*).
  //   This "degenerate" bisecting line is given by:
  //
  //     B0(t) = p + w0*t*[l0.a,l0.b]
  //
  //   where p is the projection of q along l0 and l0.a,l0.b are the _normalized_ line coefficients for e0 (or e1 which is the same)
  //   Since [a,b] is a _unit_ vector pointing perpendicularly to the left of e0 (and e1);
  //   any point B0(t) is at a scaled distance u from the line supporting e0 and e1 (scaled by w0).
  //
  //   (2)
  //   The bisecting line of e0 and e2 is given by the following SEL
  //
  //    l0.a*x(t) + l0.b*y(t) + l0.c - w0*t = 0
  //    l2.a*x(t) + l2.b*y(t) + l2.c - w2*t = 0
  //
  //   where (l0.a,l0.b,l0.c) and (l2.a,l2.b,l2.c) are the normalized line coefficientes of e0 and e2, resp.
  //
  //     B1(t) = [x(t),y(t)]
  //
  //   (3)
  //   These two bisecting lines B0(t) and B1(t) intersect (if they do) in a single point 'r' whose distance
  //   to the lines supporting the 3 edges is exactly 't' (since those expressions are precisely parametrized in a distance)
  //   Solving the following vectorial equation:
  //
  //     [x(y),y(t)] = p + w0*t*[l0.a,l0.b]
  //
  //   for t gives the result we want.
  //
  //   (4)
  //   To solve this system we can reformulate it as the intersection of three non collinear lines.
  //   The two equations from the non collinear edges:
  //
  //    l0.a*x(t) + l0.b*y(t) + l0.c - w0*t = 0
  //    l2.a*x(t) + l2.b*y(t) + l2.c - w2*t = 0
  //
  //   and if l0 is not vertical (b0 != 0), we add the vertical line:
  //
  //     x(t) = px + w0*t*l0.a
  //
  //   else, because l0 is vertical (==> not horizontal), we add the horizontal line:
  //
  //     y(t) = py + w0*t*l0.b
  //
  // sage: sage: var('a0 b0 c0 a2 b2 c2 x y t w0 w2 px py')
  // sage: eq0 = a0*x + b0*y + c0 - w0*t == 0
  // sage: eq2 = a2*x + b2*y + c2 - w2*t == 0
  // sage: eqv = px + w0 * t * a0 == x
  // sage: solve([eq0,eq2,eqv], x, y, t)
  //    [[ x == -(b0*px*w2 - (a0*b2*c0 - a0*b0*c2 + b2*px)*w0)/((a0*a2*b0 - a0^2*b2 + b2)*w0 - b0*w2),
  //       y == (a0*px*w2 - (a0*a2*c0 - a0^2*c2 + a2*px + c2)*w0 + c0*w2)/((a0*a2*b0 - a0^2*b2 + b2)*w0 - b0*w2),
  //       t == (a0*b2*px - (a2*px + c2)*b0 + b2*c0)/((a0*a2*b0 - a0^2*b2 + b2)*w0 - b0*w2) ]]
  // sage: eqv = py + w0 * t * b0 == y
  // sage: solve([eq0,eq2,eqv], x, y, t)
  //    [[ x == -(b0*py*w2 - (b0*b2*c0 - b0^2*c2 + b2*py + c2)*w0 + c0*w2)/((a2*b0^2 - a0*b0*b2 - a2)*w0 + a0*w2),
  //       y == (a0*py*w2 - (a2*b0*c0 - a0*b0*c2 + a2*py)*w0)/((a2*b0^2 - a0*b0*b2 - a2)*w0 + a0*w2),
  //       t == -(a2*b0*py - (b2*py + c2)*a0 + a2*c0)/((a2*b0^2 - a0*b0*b2 - a2)*w0 + a0*w2) ]]


  Optional_line_2 l0 = compute_normalized_line_coeffC2(tri->collinear_edge(), aCaches) ;
  Optional_line_2 l1 = compute_normalized_line_coeffC2(tri->other_collinear_edge(), aCaches) ;
  Optional_line_2 l2 = compute_normalized_line_coeffC2(tri->non_collinear_edge(), aCaches) ;

  Optional_point_2 q = construct_degenerate_seed_pointC2(tri, aCaches);

  bool ok = false ;

  if ( l0 && l1 && l2 && q )
  {
    CGAL_STSKEL_TRAITS_TRACE("\tCE ID: " << tri->collinear_edge().mID << " w: " << n2str(tri->collinear_edge_weight()) ) ;
    CGAL_STSKEL_TRAITS_TRACE("\tOCE ID: " << tri->other_collinear_edge().mID << " w: " << n2str(tri->other_collinear_edge_weight()) );
    CGAL_STSKEL_TRAITS_TRACE("\tNCE ID: " << tri->non_collinear_edge().mID << " w: " << n2str(tri->non_collinear_edge_weight()) ) ;
    CGAL_STSKEL_TRAITS_TRACE("\tLabc [" << n2str(l0->a()) << "; " << n2str(l0->b()) << "; " << n2str(l0->c()) << "]"
                                 << "[" << n2str(l1->a()) << "; " << n2str(l1->b()) << "; " << n2str(l1->c()) << "]"
                                 << "[" << n2str(l2->a()) << "; " << n2str(l2->b()) << "; " << n2str(l2->c()) << "]") ;

    FT px, py ;
    line_project_pointC2(l0->a(),l0->b(),l0->c(),q->x(),q->y(), px,py);
    CGAL_STSKEL_TRAITS_TRACE("Seed point: " << p2str(*q) << ".\nProjected seed point: (" << n2str(px) << "," << n2str(py) << ")" ) ;

    if ( tri->collinear_edge_weight() == tri->other_collinear_edge_weight() )
    {
      const FT& w0 = tri->collinear_edge_weight();
      const FT& w2 = tri->non_collinear_edge_weight();
      const FT& l0a = l0->a() ;
      const FT& l0b = l0->b() ;
      const FT& l0c = l0->c() ;
      const FT& l2a = l2->a() ;
      const FT& l2b = l2->b() ;
      const FT& l2c = l2->c() ;

      FT num(0), den(0) ;
      if ( ! CGAL_NTS is_zero(l0b) ) // l0 is not vertical, add the vertical component
      {
        num = l0a*l2b*px - (l2a*px + l2c)*l0b + l2b*l0c ;
        den = (l0a*l2a*l0b - l0a*l0a*l2b + l2b)*w0 - l0b*w2 ;

        CGAL_STSKEL_TRAITS_TRACE("Event time (degenerate, non-vertical) n=" << n2str(num) << " d=" << n2str(den) << " n/d=" << Rational<FT>(num,den) )
      }
      else // l0 is vertical, add the horizontal component
      {
        num = -(l2a*l0b*py - (l2b*py + l2c)*l0a + l2a*l0c) ;
        den = (l2a*l0b*l0b - l0a*l0b*l2b - l2a)*w0 + l0a*w2 ;

        CGAL_STSKEL_TRAITS_TRACE("Event time (degenerate, vertical) n=" << n2str(num) << " d=" << n2str(den) << " n/d=" << Rational<FT>(num,den) )
      }

      ok = CGAL_NTS is_finite(num) && CGAL_NTS is_finite(den) ;
     return cgal_make_optional(ok, Rational<FT>(num,den)) ;
    }
    else
    {
      // @todo we can't be there, but understand how (where (if)) is parallel + non-collinear handled

      // l0 and l1 are collinear but with different speeds, so there cannot be an event.
      CGAL_STSKEL_TRAITS_TRACE("Event times (degenerate, inequal norms)")
      CGAL_STSKEL_TRAITS_TRACE("--> Returning 0/0 (no event)");
      // if we return boost::none, exist_offset_lines_isec2() will think it's a numerical error
      return cgal_make_optional(Rational<FT>(FT(0),FT(0))) ;
    }
  }

  return boost::none;
}

//
// Calls the appropriate function depending on the collinearity of the edges.
//
template<class K, class Caches>
boost::optional< Rational< typename K::FT > >
compute_offset_lines_isec_timeC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                   Caches& aCaches )
{
  typedef typename K::FT FT ;
  CGAL_STSKEL_TRAITS_TRACE("compute_offset_lines_isec_timeC2(" << tri->id() << ") [" << typeid(FT).name() << "]" );

  if ( aCaches.mTime_cache.IsCached(tri->id()) )
    return aCaches.mTime_cache.Get(tri->id()) ;

  CGAL_precondition ( tri->collinearity() != TRISEGMENT_COLLINEARITY_ALL ) ;

  boost::optional< Rational<FT> > rRes =
      tri->collinearity() == TRISEGMENT_COLLINEARITY_NONE ? compute_normal_offset_lines_isec_timeC2    (tri, aCaches)
                                                          : compute_degenerate_offset_lines_isec_timeC2(tri, aCaches);

  aCaches.mTime_cache.Set(tri->id(), rRes) ;

  return rRes ;
}

// Given 3 oriented line segments e0, e1 and e2
// such that their offsets at a certain distance intersect in a single point,
// returns the coordinates (x,y) of such a point.
//
// PRECONDITIONS:
// None of e0, e1 and e2 are collinear (but two of them can be parallel)
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K, class Caches>
boost::optional< typename K::Point_2 >
construct_normal_offset_lines_isecC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                       Caches& aCaches)
{
  typedef typename K::FT  FT ;

  typedef typename K::Line_2  Line_2 ;

  typedef boost::optional<Line_2>  Optional_line_2 ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~ Computing normal offset lines isec point [" << typeid(FT).name() << "]");
  CGAL_STSKEL_TRAITS_TRACE("Event:\n" << tri);

  FT x(0), y(0) ;

  Optional_line_2 l0 = compute_normalized_line_coeffC2(tri->e0(), aCaches) ;
  Optional_line_2 l1 = compute_normalized_line_coeffC2(tri->e1(), aCaches) ;
  Optional_line_2 l2 = compute_normalized_line_coeffC2(tri->e2(), aCaches) ;

  bool ok = false ;

  if ( l0 && l1 && l2 )
  {
    // [[ x == ((c2*w1 - c1*w2)*b0 - (b2*w1 - b1*w2)*c0 + (b2*c1 - b1*c2)*w0)/((b2*w1 - b1*w2)*a0 - (a2*w1 - a1*w2)*b0 + (a2*b1 - a1*b2)*w0),
    //    y == -((c2*w1 - c1*w2)*a0 - (a2*w1 - a1*w2)*c0 + (a2*c1 - a1*c2)*w0)/((b2*w1 - b1*w2)*a0 - (a2*w1 - a1*w2)*b0 + (a2*b1 - a1*b2)*w0) ]]

    const FT& a0 = l0->a(); const FT& b0 = l0->b(); const FT& c0 = l0->c(); const FT& w0 = tri->w0();
    const FT& a1 = l1->a(); const FT& b1 = l1->b(); const FT& c1 = l1->c(); const FT& w1 = tri->w1();
    const FT& a2 = l2->a(); const FT& b2 = l2->b(); const FT& c2 = l2->c(); const FT& w2 = tri->w2();

    CGAL_STSKEL_TRAITS_TRACE("coeffs E" << tri->e0().id() << " [" << n2str(a0) << "; " << n2str(b0) << "; " << n2str(c0) << "] weight [" << n2str(w0) << "]"
                        << "\ncoeffs E" << tri->e1().id() << " [" << n2str(a1) << "; " << n2str(b1) << "; " << n2str(c1) << "] weight [" << n2str(w1) << "]"
                        << "\ncoeffs E" << tri->e2().id() << " [" << n2str(a2) << "; " << n2str(b2) << "; " << n2str(c2) << "] weight [" << n2str(w2) << "]") ;

    FT den = (b2*w1 - b1*w2)*a0 - (a2*w1 - a1*w2)*b0 + (a2*b1 - a1*b2)*w0 ;

    CGAL_STSKEL_TRAITS_TRACE("\tden=" << n2str(den) )

    if ( ! CGAL_NTS certified_is_zero(den) )
    {
      FT numX = (c2*w1 - c1*w2)*b0 - (b2*w1 - b1*w2)*c0 + (b2*c1 - b1*c2)*w0 ;
      FT numY = -((c2*w1 - c1*w2)*a0 - (a2*w1 - a1*w2)*c0 + (a2*c1 - a1*c2)*w0) ;

      CGAL_STSKEL_TRAITS_TRACE("\tnumX=" << n2str(numX) << "\n\tnumY=" << n2str(numY) ) ;

      if ( CGAL_NTS is_finite(den) && CGAL_NTS is_finite(numX) && CGAL_NTS is_finite(numY)  )
      {
        ok = true ;

        x = numX / den ;
        y = numY / den ;

       CGAL_STSKEL_TRAITS_TRACE("\n\tx=" << n2str(x) << "\n\ty=" << n2str(y) ) ;
      }
    }
  }

  return cgal_make_optional(ok,K().construct_point_2_object()(x,y)) ;
}

// Given a contour halfedge and a bisector halfedge, constructs the intersection
// between the line orthogonal to the contour halfedge through a given seed and the bisector halfedge.
// This is an artificial vertex added to recover simply-connectedness of a skeleton face
// in weighted skeletons of polygons with holes.
template <class K, class Caches>
boost::optional< typename K::Point_2 >
construct_artifical_isecC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                             Caches& aCaches )
{
  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Segment_2 Segment_2 ;
  typedef typename K::Ray_2 Ray_2 ;
  typedef typename K::Direction_2 Direction_2 ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~  Computing artificial isec point [" << typeid(typename K::FT).name() << "]");
  CGAL_STSKEL_TRAITS_TRACE("Event:\n" << tri);

  CGAL_precondition(tri->e0() == tri->e1());
  CGAL_precondition(!is_zero(tri->w0()));
  CGAL_precondition(bool(tri->child_l()));

  const Segment_2& contour_seg = tri->e0();
  Direction_2 perp_dir ( contour_seg.source().y() - contour_seg.target().y() ,
                         contour_seg.target().x() - contour_seg.source().x() ) ;
  boost::optional< typename K::Point_2 > seed = construct_offset_lines_isecC2(tri->child_l(), aCaches) ;

  if(!seed)
    return boost::none;

  const Ray_2 ray(*seed, perp_dir);
  const Segment_2& opp_seg = tri->e2();
  auto inter_res = K().intersect_2_object()(ray, opp_seg);
  if (!inter_res) // shouldn't be here if there is no intersection
    return boost::none;

  if(const Point_2* inter_pt = boost::get<Point_2>(&*inter_res))
  {
    bool ok = CGAL_NTS is_finite(inter_pt->x()) && CGAL_NTS is_finite(inter_pt->y()) ;
    return cgal_make_optional(ok, *inter_pt) ;
  }
  else if(const Segment_2* seg = boost::get<Segment_2>(&*inter_res))
  {
    // get the segment extremity closest to the seed
    const Point_2& pt = (K().compare_distance_2_object()(*seed,
                                                         seg->source(),
                                                         seg->target()) == CGAL::SMALLER) ? seg->source()
                                                                                          : seg->target() ;
    return cgal_make_optional(pt);
  }

  return boost::none;
}

// Given 3 oriented line segments e0, e1 and e2
// such that their offsets at a certain distance intersect in a single point,
// returns the coordinates (x,y) of such a point.
// two and only two of the edges are collinear, not necessarily consecutive but with the same orientaton
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at a certain distance do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
// See detailed computations in compute_degenerate_offset_lines_isec_timeC2()
template <class K, class Caches>
boost::optional< typename K::Point_2 >
construct_degenerate_offset_lines_isecC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                           Caches& aCaches)
{
  typedef typename K::FT FT ;

  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Line_2  Line_2 ;

  typedef boost::optional<Point_2> Optional_point_2 ;
  typedef boost::optional<Line_2>  Optional_line_2 ;

  if(tri->e0() == tri->e1()) // marker for artificial bisectors: they have the same face on both sides
    return construct_artifical_isecC2(tri, aCaches) ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~ Computing degenerate offset lines isec point [" << typeid(FT).name() << "]");
  CGAL_STSKEL_TRAITS_TRACE("Event:\n" << tri);

  FT x(0),y(0) ;

  Optional_line_2 l0 = compute_normalized_line_coeffC2(tri->collinear_edge(), aCaches) ;
  Optional_line_2 l2 = compute_normalized_line_coeffC2(tri->non_collinear_edge(), aCaches) ;

  Optional_point_2 q = construct_degenerate_seed_pointC2(tri, aCaches);

  bool ok = false ;

  if ( l0 && l2 && q )
  {
    const FT& w0 = tri->collinear_edge_weight();
    const FT& w2 = tri->non_collinear_edge_weight();
    const FT& l0a = l0->a() ;
    const FT& l0b = l0->b() ;
    const FT& l0c = l0->c() ;
    const FT& l2a = l2->a() ;
    const FT& l2b = l2->b() ;
    const FT& l2c = l2->c() ;

    if ( tri->collinear_edge_weight() == tri->other_collinear_edge_weight() )
    {
      FT px, py ;
      line_project_pointC2(l0->a(),l0->b(),l0->c(),q->x(),q->y(), px,py);

      CGAL_STSKEL_TRAITS_TRACE("Degenerate, equal weights " << tri->collinear_edge_weight() ) ;
      CGAL_STSKEL_TRAITS_TRACE("Seed point: " << p2str(*q) << ". Projected seed point: (" << n2str(px) << "," << n2str(py) << ")" ) ;

      FT num_x, num_y, den ;
      if ( ! CGAL_NTS is_zero(l0->b()) ) // Non-vertical
      {
        num_x = -(l0b*px*w2 - (l0a*l2b*l0c - l0a*l0b*l2c + l2b*px)*w0) ;
        num_y = (l0a*px*w2 - (l0a*l2a*l0c - l0a*l0a*l2c + l2a*px + l2c)*w0 + l0c*w2) ;
        den = ((l0a*l2a*l0b - l0a*l0a*l2b + l2b)*w0 - l0b*w2) ;
      }
      else
      {
        num_x = -(l0b*py*w2 - (l0b*l2b*l0c - l0b*l0b*l2c + l2b*py + l2c)*w0 + l0c*w2) ;
        num_y = (l0a*py*w2 - (l2a*l0b*l0c - l0a*l0b*l2c + l2a*py)*w0) ;
        den = ((l2a*l0b*l0b - l0a*l0b*l2b - l2a)*w0 + l0a*w2) ;
      }

      if ( ! CGAL_NTS certified_is_zero(den) && CGAL_NTS is_finite(den) &&
             CGAL_NTS is_finite(num_x) && CGAL_NTS is_finite(num_y) )
      {
        x = num_x / den ;
        y = num_y / den ;
        ok = CGAL_NTS is_finite(x) && CGAL_NTS is_finite(y) ;
      }
    }
    else
    {
      CGAL_assertion(false); // this is COLLINEAR, not parallel so we can never be here

      CGAL_STSKEL_TRAITS_TRACE("Degenerate, different weights " << n2str(tri->collinear_edge_weight())
                                                     << " and " << n2str(tri->other_collinear_edge_weight()));

      // The line parallel to l0 (and l1) passing through q is: l0a*x + l0b*y + lambda = 0, with
      const FT lambda = -l0a*q->x() - l0b*q->y();

      // The bisector between l0 (l1) and l2 is:
      //  l0a*x + l0b*y + l0c - w0*t = 0
      //  l2a*x + l2b*y + l2c - w2*t = 0

      // The intersection point is thus:
      //  l0a*x + l0b*y + l0c - w0*t = 0
      //  l2a*x + l2b*y + l2c - w2*t = 0
      //  l0a*x + l0b*y + lambda = 0

      // const FT t = l0c - lambda ; // (3) - (1)
      const FT den = l2a*l0b - l0a*l2b;

      if ( ! CGAL_NTS certified_is_zero(den) && CGAL_NTS is_finite(den) )

      x =  (l0b*l0c - l0b*(l2c + lambda) + l2b*lambda) / den;
      y = -(l0a*l0c - l0a*(l2c + lambda) + l2a*lambda) / den;
    }
  }

  CGAL_STSKEL_TRAITS_TRACE("Degenerate" << (CGAL_NTS is_zero(l0->b()) ? " (vertical)" : "") << " event point:  x=" << n2str(x) << " y=" << n2str(y) )

  return cgal_make_optional(ok,K().construct_point_2_object()(x,y)) ;
}

// Calls the appropriate function depending on the collinearity of the edges.
template <class K, class Caches>
boost::optional< typename K::Point_2 >
construct_offset_lines_isecC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                Caches& aCaches )
{
  typedef typename K::Point_2 Point_2 ;

  CGAL_STSKEL_TRAITS_TRACE("construct_offset_lines_isecC2(" << tri->id() << ") [" << typeid(typename K::FT).name() << "]" );

  if(aCaches.mPoint_cache.IsCached(tri->id()) )
    return aCaches.mPoint_cache.Get(tri->id()) ;

  CGAL_precondition ( tri->collinearity() != TRISEGMENT_COLLINEARITY_ALL ) ;

  boost::optional< Point_2 > rRes =
      tri->collinearity() == TRISEGMENT_COLLINEARITY_NONE ? construct_normal_offset_lines_isecC2    (tri, aCaches)
                                                          : construct_degenerate_offset_lines_isecC2(tri, aCaches);


  CGAL_STSKEL_TRAITS_TRACE("isec = " << (rRes ? p2str(*rRes) : "none"));

  aCaches.mPoint_cache.Set(tri->id(), rRes) ;

  return rRes ;
}

} // namespace CGAL_SS_i
} // namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_CONS_FTC2_H
