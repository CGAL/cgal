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


template<class K>
Trisegment_collinearity trisegment_collinearity_no_exact_constructions ( Segment_2_with_ID<K> const& e0
                                                                       , Segment_2_with_ID<K> const& e1
                                                                       , Segment_2_with_ID<K> const& e2)
{
  bool is_01 = are_edges_orderly_collinear(e0,e1);
  bool is_02 = are_edges_orderly_collinear(e0,e2);
  bool is_12 = are_edges_orderly_collinear(e1,e2);

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

template <typename NT>
typename Coercion_traits<double, NT>::Type
inexact_sqrt(const NT& n, CGAL::Null_functor)
{
  typedef CGAL::Interval_nt<false> IFT;
  typename IFT::Protector protector;

  CGAL::NT_converter<NT, IFT> to_ift;
  IFT sqrt_ift = sqrt(to_ift(n));
  CGAL_STSKEL_TRAITS_TRACE("interval " << sqrt_ift.inf() << " " << sqrt_ift.sup() ) ;
  CGAL_STSKEL_TRAITS_TRACE("delta " << sqrt_ift.sup() - sqrt_ift.inf() ) ;

  return NT(to_double(sqrt_ift));
}

template <typename NT, typename Sqrt>
typename Sqrt::result_type
inexact_sqrt(const NT& nt, Sqrt sqrt)
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
  return inexact_sqrt(nt, Sqrt());
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

// Currently the norm can be inexact even when we are in the exact pipeline of the traits
// if we use something like K = EPICK because there is no exact sqrt.
//
// @todo Ideally, we could compute how much precision is required for the sqrt
// given all possible operations that are performed in the SLS and the input values
// and compute (in the exact pipeline) an approximate sqrt with such sufficient precision.

template <typename FT>
FT inexact_norm (const FT& x, const FT& y,
                 CGAL::Null_functor /*no_sqrt*/)
{
  CGAL_STSKEL_TRAITS_TRACE("inexact_norm(" << typeid(FT).name() << "," << typeid(FT).name() << ",Null_functor)");
  return std::hypot(CGAL::to_double(x), CGAL::to_double(y));
}

template <typename FT, class Sqrt>
FT inexact_norm (const FT& x, const FT& y,
                 Sqrt sqrt_f)
{
  CGAL_STSKEL_TRAITS_TRACE("inexact_norm(" << typeid(FT).name() << "," << typeid(FT).name() << "," << typeid(Sqrt).name() << ")");
  const FT n = square(x) + square(y);
  return sqrt_f(n);
}

template <typename FT>
FT inexact_norm (const FT& x, const FT& y)
{
  typedef CGAL::Algebraic_structure_traits<FT> AST;
  typedef typename AST::Sqrt Sqrt;
  return inexact_norm(x,y, Sqrt());
}

template <typename NT>
inline Lazy_exact_nt<NT> inexact_norm( Lazy_exact_nt<NT> const& lx,
                                       Lazy_exact_nt<NT> const& ly )
{
  return inexact_norm( exact(lx), exact(ly) ) ;
}

template <typename NT>
inline Quotient<NT> inexact_norm( Quotient<NT> const& x,
                                  Quotient<NT> const& y )
{
  CGAL_STSKEL_TRAITS_TRACE("inexact_norm(Quotient,Quotient)");
  return { inexact_norm(x.numerator()*y.denominator(), y.numerator()*x.denominator()),
           abs(x.denominator()*y.denominator()) } ;
}

template <bool Protected>
inline
Interval_nt<Protected>
inexact_norm (const Interval_nt<Protected>& x, const Interval_nt<Protected>& y)
{
  typename Interval_nt<Protected>::Internal_protector P;

  CGAL_STSKEL_TRAITS_TRACE("inexact_norm(Interval_nt,Interval_nt)");

#if 0 // CGAL_USE_SSE2  _mm_hypot_pd is documented but does not exist...?
  __m128d xx = IA_opacify128(x.simd());
  __m128d yy = IA_opacify128(y.simd());
  __m128d r = _mm_hypot_pd(xx, yy);
  return Interval_nt<Protected>(IA_opacify128(r));
#else

  // @fixme is below correct?

#ifdef CGAL_ALWAYS_ROUND_TO_NEAREST
  double i = std::nextafter(std::hypot(x.inf(), y.inf()), 0.) ;
#else
  FPU_set_cw(CGAL_FE_DOWNWARD);
  double i = CGAL_IA_FORCE_TO_DOUBLE(std::hypot(CGAL_IA_STOP_CPROP(x.inf()),
                                                CGAL_IA_STOP_CPROP(y.inf())));
  FPU_set_cw(CGAL_FE_UPWARD);
#endif

  return Interval_nt<Protected>(i, IA_up(std::hypot(x.sup(), y.sup())));
#endif
}

///

template <typename ET>
CGAL::Lazy_exact_nt<ET> ceil(const CGAL::Lazy_exact_nt<ET>& n)
{
  return { ceil(exact(n)) };
}

///


// Given an oriented 2D straight line segment 'e', computes the normalized coefficients (a,b,c)
// of the supporting line, and weights them with 'aWeight'.
//
// POSTCONDITION: [a,b] is the leftward normal vector.
// POSTCONDITION: In case of overflow, an empty optional<> is returned.
template<class K>
boost::optional< typename K::Line_2> compute_normalized_line_coeffC2( Segment_2<K> const& e )
{
  typedef typename K::FT FT ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~ Unweighted line coefficients for " << s2str(e) );

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
    FT l = inexact_norm(sa, sb);

    if ( CGAL_NTS is_finite(l) )
    {
      a = sa / l ;
      b = sb / l ;

      c = -e.source().x()*a - e.source().y()*b;

      CGAL_STSKEL_TRAITS_TRACE("GENERIC line; sa="<< n2str(sa) << " sb=" << n2str(sb)
                               << "\nl=" << n2str(l)
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

template<class K, class CoeffCache>
boost::optional< typename K::Line_2 >
compute_normalized_line_coeffC2( Segment_2_with_ID<K> const& e,
                                 CoeffCache& aCoeff_cache )
{
  typedef typename K::Segment_2 Segment_2 ;
  typedef typename K::Line_2 Line_2 ;

  if ( aCoeff_cache.IsCached(e.mID) )
    return aCoeff_cache.Get(e.mID) ;

  boost::optional< Line_2 > rRes = compute_normalized_line_coeffC2 ( static_cast<const Segment_2&>(e) ) ;

  aCoeff_cache.Set(e.mID, rRes) ;

  return rRes ;
}

// @todo weightless coefficients are stored because we use them sometimes weighted, and sometimes
// inversely weighted (filtering bound). Should we store them weighted also?
template<class K, class CoeffCache>
boost::optional< typename K::Line_2 > compute_weighted_line_coeffC2( Segment_2_with_ID<K> const& e,
                                                                     typename K::FT const& aWeight,
                                                                     CoeffCache& aCoeff_cache )
{
  typedef typename K::FT FT ;
  typedef typename K::Line_2 Line_2 ;

  CGAL_precondition( CGAL_NTS is_finite(aWeight) && CGAL_NTS is_positive(aWeight) ) ;

  boost::optional< Line_2 > l = compute_normalized_line_coeffC2(e, aCoeff_cache);
  if( ! l )
    return boost::none ;

  FT a = l->a() * aWeight ;
  FT b = l->b() * aWeight ;
  FT c = l->c() * aWeight ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~ Weighted line coefficients for line: " << s2str(e)
                            << "\nweight=" << n2str(aWeight)
                            << "\na="<< n2str(a) << "\nb=" << n2str(b) << "\nc=" << n2str(c)
                            ) ;

  if ( !CGAL_NTS is_finite(a) || !CGAL_NTS is_finite(b) || !CGAL_NTS is_finite(c) )
    return boost::none;
  else
    return cgal_make_optional( true, K().construct_line_2_object()(a,b,c) ) ;
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
template<class K>
Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > >
construct_trisegment ( Segment_2_with_ID<K> const& e0,
                       typename K::FT const& w0,
                       Segment_2_with_ID<K> const& e1,
                       typename K::FT const& w1,
                       Segment_2_with_ID<K> const& e2,
                       typename K::FT const& w2,
                       std::size_t id )
{
  typedef Trisegment_2<K, Segment_2_with_ID<K> >Trisegment_2 ;
  typedef typename Trisegment_2::Self_ptr Trisegment_2_ptr ;

  Trisegment_collinearity lCollinearity = trisegment_collinearity_no_exact_constructions(e0,e1,e2);

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
template <class K, class CoeffCache>
boost::optional< Rational< typename K::FT> >
compute_normal_offset_lines_isec_timeC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                          CoeffCache& aCoeff_cache )
{
  typedef typename K::FT  FT ;
  typedef typename K::Line_2 Line_2 ;

  typedef boost::optional<Line_2> Optional_line_2 ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~ Computing normal offset lines isec time [" << typeid(FT).name() << "]") ;
  CGAL_STSKEL_TRAITS_TRACE("Event" << tri );

  FT num(0), den(0) ;

  // DETAILS:
  //
  // An offset line is given by:
  //
  //   a*x(t) + b*y(t) + c - t = 0
  //
  // were 't > 0' being to the left of the line.
  // If 3 such offset lines intersect at the same offset distance, the intersection 't',
  // or 'time', can be computed solving for 't' in the linear system formed by 3 such equations.
  // The result is :
  //
  // sage: var('a0 b0 c0 a1 b1 c1 a2 b2 c2 x y t w0 w1 w2')
  // (a0, b0, c0, a1, b1, c1, a2, b2, c2, x, y, t, w0, w1, w2)
  // sage:
  // sage: eqw0 = w0*a0*x + w0*b0*y + w0*c0 - t == 0
  // sage: eqw1 = w1*a1*x + w1*b1*y + w1*c1 - t == 0
  // sage: eqw2 = w2*a2*x + w2*b2*y + w2*c2 - t == 0
  // sage:
  // sage: solve([eqw0,eqw1,eqw2], x,y,t)
  //   x ==  (((c1*w1 - c2*w2)*b0 - (b1*w1 - b2*w2)*c0)*w0 - (b2*c1*w2 - b1*c2*w2)*w1) / (((b1*w1 - b2*w2)*a0 - (a1*w1 - a2*w2)*b0)*w0 - (a2*b1*w2 - a1*b2*w2)*w1),
  //   y == -(((c1*w1 - c2*w2)*a0 - (a1*w1 - a2*w2)*c0)*w0 - (a2*c1*w2 - a1*c2*w2)*w1) / (((b1*w1 - b2*w2)*a0 - (a1*w1 - a2*w2)*b0)*w0 - (a2*b1*w2 - a1*b2*w2)*w1),
  //   t == -((b2*c1*w2 - b1*c2*w2)*a0*w1 - (a2*c1*w2 - a1*c2*w2)*b0*w1 + (a2*b1*w2 - a1*b2*w2)*c0*w1)*w0/(((b1*w1 - b2*w2)*a0 - (a1*w1 - a2*w2)*b0)*w0 - (a2*b1*w2 - a1*b2*w2)*w1)

  bool ok = false ;

  Optional_line_2 l0 = compute_weighted_line_coeffC2(tri->e0(), tri->w0(), aCoeff_cache) ;
  Optional_line_2 l1 = compute_weighted_line_coeffC2(tri->e1(), tri->w1(), aCoeff_cache) ;
  Optional_line_2 l2 = compute_weighted_line_coeffC2(tri->e2(), tri->w2(), aCoeff_cache) ;

  if ( l0 && l1 && l2 )
  {
    CGAL_STSKEL_TRAITS_TRACE("coeffs 0 [" << n2str(l0->a()) << "; " << n2str(l0->b()) << "; " << n2str(l0->c()) << "]"
                          << "\ncoeffs 1 [" << n2str(l1->a()) << "; " << n2str(l1->b()) << "; " << n2str(l1->c()) << "]"
                          << "\ncoeffs 2 [" << n2str(l2->a()) << "; " << n2str(l2->b()) << "; " << n2str(l2->c()) << "]") ;

    num = (l2->a()*l0->b()*l1->c())
         -(l2->a()*l1->b()*l0->c())
         -(l2->b()*l0->a()*l1->c())
         +(l2->b()*l1->a()*l0->c())
         +(l1->b()*l0->a()*l2->c())
         -(l0->b()*l1->a()*l2->c());

    den = (-l2->a()*l1->b())
         +( l2->a()*l0->b())
         +( l2->b()*l1->a())
         -( l2->b()*l0->a())
         +( l1->b()*l0->a())
         -( l0->b()*l1->a());

    ok = CGAL_NTS is_finite(num) && CGAL_NTS is_finite(den);
  }

  CGAL_STSKEL_TRAITS_TRACE("Event time (normal): n=" << num << " d=" << den << " n/d=" << Rational<FT>(num,den)  )

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
                            << "\ne0: " << s2str(e0)
                            << "\ne1: " << s2str(e1)
                          );

  FT delta01 = CGAL::squared_distance(e0.target(),e1.source());
  if( CGAL_NTS is_finite(delta01) && CGAL_NTS is_zero(delta01))
    return (true /*ok*/, e0.target());

  FT delta10 = CGAL::squared_distance(e1.target(),e0.source());
  if( CGAL_NTS is_finite(delta10) && CGAL_NTS is_zero(delta10))
    return (true /*ok*/, e1.target());

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
template <class K, class CoeffCache>
boost::optional< typename K::Point_2 >
compute_seed_pointC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                       typename Trisegment_2<K, Segment_2_with_ID<K> >::SEED_ID sid,
                       CoeffCache& aCoeff_cache)
{
  boost::optional< typename K::Point_2 > p ;

  typedef Trisegment_2<K, Segment_2_with_ID<K> > Trisegment_2 ;

  switch ( sid )
  {
    case Trisegment_2::LEFT :

       p = tri->child_l() ? construct_offset_lines_isecC2(tri->child_l(), aCoeff_cache)  // this can recurse
                          : compute_oriented_midpoint(tri->e0(),tri->e1()) ;
       break ;

    case Trisegment_2::RIGHT :

      p = tri->child_r() ? construct_offset_lines_isecC2(tri->child_r(), aCoeff_cache) // this can recurse
                         : compute_oriented_midpoint(tri->e1(),tri->e2()) ;
      break ;

    case Trisegment_2::THIRD :

      p = tri->child_t() ? construct_offset_lines_isecC2(tri->child_t(), aCoeff_cache) // this can recurse
                         : compute_oriented_midpoint(tri->e0(),tri->e2()) ;

      break ;
  }

  return p ;
}

//
// Given the trisegment tree for an event which is known to have a normal collinearity returns the seed point
// of the degenerate seed.
// A normal collinearity occurs when e0,e1 or e1,e2 are collinear.
template <class K, class CoeffCache>
boost::optional< typename K::Point_2 >
compute_degenerate_seed_pointC2 ( boost::intrusive_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                  CoeffCache& aCoeff_cache )
{
  return compute_seed_pointC2( tri, tri->degenerate_seed_id(), aCoeff_cache ) ;
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
template <class K, class CoeffCache>
boost::optional< Rational< typename K::FT> >
compute_degenerate_offset_lines_isec_timeC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                              CoeffCache& aCoeff_cache )
{
  typedef typename K::FT FT ;

  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Line_2 Line_2 ;

  typedef boost::optional<Point_2> Optional_point_2 ;
  typedef boost::optional<Line_2>  Optional_line_2 ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~  Computing degenerate offset lines isec time for: " << tri ) ;

  // DETAILS:
  //
  // For simplicity, assume e0,e1 are the collinear edges.
  //
  //   (1)
  //   The bisecting line of e0 and e1 is a line perpendicular to e0 (and e1)
  //   which passes through 'q', the degenerate offset vertex (e0*,e1*).
  //   This "degenerate" bisecting line is given by:
  //
  //     B0(t) = p + t*[l0.a,l0.b]
  //
  //   where p is the projection of q along l0 and l0.a,l0.b are the _normalized_ line coefficients for e0 (or e1 which is the same)
  //   Since [a,b] is a _unit_ vector pointing perpendicularly to the left of e0 (and e1);
  //   any point B0(k) is at a distance k from the line supporting e0 and e1.
  //
  //   (2)
  //   The bisecting line of e0 and e2 is given by the following SEL
  //
  //    l0.a*x(t) + l0.b*y(t) + l0.c - t = 0
  //    l2.a*x(t) + l2.b*y(t) + l2.c - t = 0
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
  //     [x(y),y(t)] = p + t*[l0.a,l0.b]
  //
  //   for t gives the result we want.
  //
  //
  //   (4)
  //   With weights, the above equations become:
  //
  //   sage: eq0 = w0*a0*x + w0*b0*y + w0*c0 - t == 0
  //   sage: eq2 = w2*a2*x + w2*b2*y + w2*c2 - t == 0
  //
  //   sage: solve([eq0,eq2], x,y)
  //     [[x == -(b2*t*w2 - (b2*c0*w2 - b0*c2*w2 + b0*t)*w0)/((a2*b0*w2 - a0*b2*w2)*w0),
  //       y ==  (a2*t*w2 - (a2*c0*w2 - a0*c2*w2 + a0*t)*w0)/((a2*b0*w2 - a0*b2*w2)*w0) ]]
  //
  //   sage: x0 = -(b2*t*w2 - (b2*c0*w2 - b0*c2*w2 + b0*t)*w0)/((a2*b0*w2 - a0*b2*w2)*w0)
  //   sage: eqb0 = qx + t * a0 / w0 - x0 == 0
  //   sage: solve(eqb0, t)
  //     [t == -(b2*c0 - b0*c2 - (a2*b0 - a0*b2)*qx)*w0*w2/(b0*w0 - (a0*a2*b0 - (a0^2 - 1)*b2)*w2) ]
  //
  //   sage: y0 = (a2*t*w2 - (a2*c0*w2 - a0*c2*w2 + a0*t)*w0)/((a2*b0*w2 - a0*b2*w2)*w0)
  //   sage: eqb1 = qy + t * b0 / w0 - y0 == 0
  //   sage: solve(eqb1, t)
  //     [t == -(a2*c0 - a0*c2 + (a2*b0 - a0*b2)*qy)*w0*w2/(a0*w0 + (a2*b0^2 - a0*b0*b2 - a2)*w2)]

  Optional_line_2 l0 = compute_weighted_line_coeffC2(tri->collinear_edge(), tri->collinear_edge_weight(), aCoeff_cache) ;
  Optional_line_2 l1 = compute_weighted_line_coeffC2(tri->other_collinear_edge(), tri->other_collinear_edge_weight(), aCoeff_cache) ;
  Optional_line_2 l2 = compute_weighted_line_coeffC2(tri->non_collinear_edge(), tri->non_collinear_edge_weight(), aCoeff_cache) ;

  Optional_point_2 q = compute_degenerate_seed_pointC2(tri, aCoeff_cache);

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

    const Comparison_result res = compare(tri->collinear_edge_weight(), tri->other_collinear_edge_weight());
    if ( res == EQUAL )
    {
      const FT& l0a = l0->a() ;
      const FT& l0b = l0->b() ;
      const FT& l0c = l0->c() ;
      const FT& l2a = l2->a() ;
      const FT& l2b = l2->b() ;
      const FT& l2c = l2->c() ;

      // Since l0 and l1 are parallel, we cannot solve the system using:
      //   l0a*x + l0b*y + l0c = 0
      //   l1a*x + l1b*y + l1c = 0
      // Instead, we use the equation of the line orthogonal to l0 (and l1).
      // However, rephrasing
      //   l0a*x + l0b*y + l0c = 0
      // to
      //   [x, y] = projected_seed + t * N
      // requires the norm (l0a² + l0b²) to be exactly '1', which likely isn't the case
      // if we are using inexact square roots. In that case, the norm behaves similarly
      // to a weight, and likewise needs to be inverted. So (now with weights), the equation
      //   w*l0a*x + w*l0b*y + w*l0c = 0
      // with l0a² + l0b² ~= 1 is:
      //   w'*l0a'*x + w'*l0b'*y + w'*l0c' = 0,
      // with l0a'² + l0b'² = 1. The orthogonal displacement is rephrased to:
      //   [x, y] = projected_seed + t / w' * N,
      // with w' = weight * (l0a² + l0b²).
      //
      // @todo further robustification by storing independently the norm (and not the weighted,
      // normalized coeffs?)
      const FT sq_w0 = square(l0a) + square(l0b); // l0a and l0b are *weighted* coefficients

      FT num(0), den(0) ;
      if ( ! CGAL_NTS is_zero(l0->b()) ) // Non-vertical
      {
        num = ((l2a*l0b - l0a*l2b) * px - l2b*l0c + l0b*l2c) * sq_w0 ;
        den = l0a*l0a*l2b - l2b*sq_w0 + l0b*sq_w0 - l0a*l2a*l0b ;

        CGAL_STSKEL_TRAITS_TRACE("Event time (degenerate, non-vertical) n=" << n2str(num) << " d=" << n2str(den) << " n/d=" << Rational<FT>(num,den) )
      }
      else
      {
        num = ((l2a*l0b - l0a*l2b) * py - l0a*l2c + l2a*l0c) * sq_w0 ;
        den = l0a*l0b*l2b - l0b*l0b*l2a + l2a*sq_w0 - l0a*sq_w0;

        CGAL_STSKEL_TRAITS_TRACE("Event time (degenerate, vertical) n=" << n2str(num) << " d=" << n2str(den) << " n/d=" << Rational<FT>(num,den) )
      }

      ok = CGAL_NTS is_finite(num) && CGAL_NTS is_finite(den) && !CGAL_NTS certified_is_zero(den);
      return cgal_make_optional(ok, Rational<FT>(num,den)) ;
    }
    else
    {
      // l0 and l1 are collinear and are meeting up, the time is simply the time at Q,
      // and the event only exists if the two times are identical
      //
      // @todo, which direction the horizontal line uses is likely dictated by the order
      // of the input segments in the triedge, maybe we would like more control over this.
      const FT t0 = l0->a() * q->x() + l0->b() * q->y() + l0->c();
      const FT t1 = l1->a() * q->x() + l1->b() * q->y() + l1->c();

      if( CGAL_NTS is_finite(t0) && CGAL_NTS is_finite(t1) && t0 == t1 )
      {
        CGAL_STSKEL_TRAITS_TRACE("Event time (degenerate, inequal norms) t:" << t0 )
        return cgal_make_optional(true, Rational<FT>(t0,FT(1))) ;
      }
      else
      {
        CGAL_STSKEL_TRAITS_TRACE("Event times (degenerate, inequal norms) t0:" << t0 << " != t1:" << t1 )
        CGAL_STSKEL_TRAITS_TRACE("--> Returning 0/0 (no event)");
        // if we return boost::none, exist_offset_lines_isec2() will think it's a numerical error
        return cgal_make_optional(true, Rational<FT>(FT(0),FT(0))) ;
      }
    }
  }

  return boost::none;
}

//
// Calls the appropriate function depending on the collinearity of the edges.
//
template<class K, class TimeCache, class CoeffCache>
boost::optional< Rational< typename K::FT > >
compute_offset_lines_isec_timeC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                   TimeCache& aTime_cache,
                                   CoeffCache& aCoeff_cache)
{
  if ( aTime_cache.IsCached(tri->id()) )
    return aTime_cache.Get(tri->id()) ;

  CGAL_precondition ( tri->collinearity() != TRISEGMENT_COLLINEARITY_ALL ) ;

  boost::optional< Rational< typename K::FT > > rRes =
      tri->collinearity() == TRISEGMENT_COLLINEARITY_NONE ? compute_normal_offset_lines_isec_timeC2    (tri, aCoeff_cache)
                                                          : compute_degenerate_offset_lines_isec_timeC2(tri, aCoeff_cache);

  aTime_cache.Set(tri->id(), rRes) ;

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
template<class K, class CoeffCache>
boost::optional< typename K::Point_2 >
construct_normal_offset_lines_isecC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                       CoeffCache& aCoeff_cache)
{
  typedef typename K::FT  FT ;

  typedef typename K::Line_2  Line_2 ;

  typedef boost::optional<Line_2>  Optional_line_2 ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~ Computing normal offset lines isec point [" << typeid(FT).name() << "]");
  CGAL_STSKEL_TRAITS_TRACE("Event:" << tri ) ;

  FT x(0), y(0) ;

  Optional_line_2 l0 = compute_weighted_line_coeffC2(tri->e0(), tri->w0(), aCoeff_cache) ;
  Optional_line_2 l1 = compute_weighted_line_coeffC2(tri->e1(), tri->w1(), aCoeff_cache) ;
  Optional_line_2 l2 = compute_weighted_line_coeffC2(tri->e2(), tri->w2(), aCoeff_cache) ;

  bool ok = false ;

  if ( l0 && l1 && l2 )
  {
    FT den = l0->a()*l2->b() - l0->a()*l1->b() - l1->a()*l2->b() + l2->a()*l1->b() + l0->b()*l1->a() - l0->b()*l2->a();

    CGAL_STSKEL_TRAITS_TRACE("\tden=" << n2str(den) )

    if ( ! CGAL_NTS certified_is_zero(den) )
    {
      FT numX = l0->b()*l2->c() - l0->b()*l1->c() - l1->b()*l2->c() + l2->b()*l1->c() + l1->b()*l0->c() - l2->b()*l0->c();
      FT numY = l0->a()*l2->c() - l0->a()*l1->c() - l1->a()*l2->c() + l2->a()*l1->c() + l1->a()*l0->c() - l2->a()*l0->c();

      CGAL_STSKEL_TRAITS_TRACE("\tnumX=" << n2str(numX) << "\n\tnumY=" << n2str(numY) ) ;

      if ( CGAL_NTS is_finite(den) && CGAL_NTS is_finite(numX) && CGAL_NTS is_finite(numY)  )
      {
        ok = true ;

        x =  numX / den ;
        y = -numY / den ;

       CGAL_STSKEL_TRAITS_TRACE("\n\tx=" << n2str(x) << "\n\ty=" << n2str(y) ) ;
      }
    }
  }

  return cgal_make_optional(ok,K().construct_point_2_object()(x,y)) ;
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
template <class K, class CoeffCache>
boost::optional< typename K::Point_2 >
construct_degenerate_offset_lines_isecC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                           CoeffCache& aCoeff_cache)
{
  typedef typename K::FT FT ;

  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Line_2  Line_2 ;

  typedef boost::optional<Point_2> Optional_point_2 ;
  typedef boost::optional<Line_2>  Optional_line_2 ;

  CGAL_STSKEL_TRAITS_TRACE("\n~~ Computing degenerate offset lines isec point [" << typeid(FT).name() << "]");
  CGAL_STSKEL_TRAITS_TRACE("Event: " << tri ) ;

  FT x(0),y(0) ;

  Optional_line_2 l0 = compute_weighted_line_coeffC2(tri->collinear_edge(), tri->collinear_edge_weight(), aCoeff_cache) ;
  Optional_line_2 l2 = compute_weighted_line_coeffC2(tri->non_collinear_edge(), tri->non_collinear_edge_weight(), aCoeff_cache) ;

  Optional_point_2 q = compute_degenerate_seed_pointC2(tri, aCoeff_cache);

  bool ok = false ;

  if ( l0 && l2 && q )
  {
    const Comparison_result res = compare(tri->collinear_edge_weight(), tri->other_collinear_edge_weight());
    if ( res == EQUAL )
    {
      FT px, py ;
      line_project_pointC2(l0->a(),l0->b(),l0->c(),q->x(),q->y(), px,py);

      CGAL_STSKEL_TRAITS_TRACE("Degenerate, equal weights " << tri->collinear_edge_weight() ) ;
      CGAL_STSKEL_TRAITS_TRACE("Seed point: " << p2str(*q) << ". Projected seed point: (" << n2str(px) << "," << n2str(py) << ")" ) ;
      const FT& l0a = l0->a() ;
      const FT& l0b = l0->b() ;
      const FT& l0c = l0->c() ;
      const FT& l2a = l2->a() ;
      const FT& l2b = l2->b() ;
      const FT& l2c = l2->c() ;

      // See details in compute_degenerate_offset_lines_isec_timeC2()
      const FT sq_w0 = square(l0a) + square(l0b);

      // Note that "* sq_w0" is removed from the numerator expression.
      //
      // This is because the speed is inverted while representing the front
      // progression using the orthogonal normalized vector [l0a, l0b]: P = Q + t/w * V with V normalized.
      // However, here l0a & l0b are not normalized but *weighted* coeff, so we need to divide by w0².
      // Hence we can just avoid multiplying by w0² in the numerator in the first place.
      FT num, den ;
      if ( ! CGAL_NTS is_zero(l0->b()) ) // Non-vertical
      {
        num = ((l2a*l0b - l0a*l2b) * px - l2b*l0c + l0b*l2c) /* * sq_w0 */ ;
        den = l0a*l0a*l2b - l2b*sq_w0 + l0b*sq_w0 - l0a*l2a*l0b ;
      }
      else
      {
        num = ((l2a*l0b - l0a*l2b) * py - l0a*l2c + l2a*l0c) /* * sq_w0 */ ;
        den = l0a*l0b*l2b - l0b*l0b*l2a + l2a*sq_w0 - l0a*sq_w0;
      }

      if ( ! CGAL_NTS certified_is_zero(den) && CGAL_NTS is_finite(den) && CGAL_NTS is_finite(num) )
      {
        x = px + l0a * num / den ;
        y = py + l0b * num / den ;

        ok = CGAL_NTS is_finite(x) && CGAL_NTS is_finite(y) ;
      }
    }
    else
    {
      CGAL_STSKEL_TRAITS_TRACE("Degenerate, different weights " << n2str(tri->collinear_edge_weight())
                                                     << " and " << n2str(tri->other_collinear_edge_weight()));

      const FT& l0a = l0->a() ; const FT& l0b = l0->b() ; const FT& l0c = l0->c() ;
      const FT& l2a = l2->a() ; const FT& l2b = l2->b() ; const FT& l2c = l2->c() ;

      // The line parallel to l0 (and l1) passing through q is: l0a*x + l0b*y + lambda = 0, with
      const FT lambda = -l0a*q->x() - l0b*q->y();

      // The bisector between l0 (l1) and l2 is:
      //  l0a*x + l0b*y + l0c - t = 0
      //  l2a*x + l2b*y + l2c - t = 0

      // The intersection point is thus:
      //  l0a*x + l0b*y + l0c - t = 0
      //  l2a*x + l2b*y + l2c - t = 0
      //  l0a*x + l0b*y + lambda = 0

      const FT t = l0c - lambda ; // (3) - (1)
      const FT den = l2a*l0b - l0a*l2b;

      if ( ! CGAL_NTS certified_is_zero(den) && CGAL_NTS is_finite(den) )

      x =  (l0b*l0c - l0b*(l2c + lambda) + l2b*lambda) / den;
      y = -(l0a*l0c - l0a*(l2c + lambda) + l2a*lambda) / den;
    }
  }

  CGAL_STSKEL_TRAITS_TRACE("Degenerate" << (CGAL_NTS is_zero(l0->b()) ? " (vertical)" : "") << " event point:  x=" << n2str(x) << " y=" << n2str(y) )

  return cgal_make_optional(ok,K().construct_point_2_object()(x,y)) ;
}

//
// Calls the appropriate function depending on the collinearity of the edges.
//
template <class K, class CoeffCache>
boost::optional< typename K::Point_2 >
construct_offset_lines_isecC2 ( Trisegment_2_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                CoeffCache& aCoeff_cache)
{
  CGAL_precondition ( tri->collinearity() != TRISEGMENT_COLLINEARITY_ALL ) ;

  return tri->collinearity() == TRISEGMENT_COLLINEARITY_NONE ? construct_normal_offset_lines_isecC2    (tri, aCoeff_cache)
                                                             : construct_degenerate_offset_lines_isecC2(tri, aCoeff_cache) ;
}

} // namespace CGAL_SS_i
} // namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_CONS_FTC2_H
