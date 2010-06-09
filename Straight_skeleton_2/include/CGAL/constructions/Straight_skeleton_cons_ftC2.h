// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_CONS_FTC2_H
#define CGAL_STRAIGHT_SKELETON_CONS_FTC2_H 1

namespace CGAL { 

namespace CGAL_SS_i
{

#ifdef CGAL_USE_CORE  

template<class NT>
inline CORE::BigFloat to_BigFloat( NT const& n )
{
  return CORE::BigFloat( CGAL::to_double(n) ) ;
}

template<>
inline CORE::BigFloat to_BigFloat<MP_Float>( MP_Float const& b )
{
  if (b.is_zero())
    return CORE::BigFloat::getZero();
    
  typedef MP_Float::exponent_type exponent_type;

  const int                    log_limb         = 8 * sizeof(MP_Float::limb);
  const MP_Float::V::size_type limbs_per_double = 2 + 53/log_limb;
    
  exponent_type exp = b.max_exp();
  int steps = (std::min)(limbs_per_double, b.v.size());
  
  CORE::BigFloat d_exp_1 = CORE::BigFloat::exp2(-log_limb);
  
  CORE::BigFloat d_exp   = CORE::BigFloat::getOne() ;
  
  CORE::BigFloat d       = CORE::BigFloat::getZero();

  for ( exponent_type i = exp - 1; i > exp - 1 - steps; i--) 
  {
    d_exp *= d_exp_1;
    d += d_exp * CORE::BigFloat(b.of_exp(i));
  }

  return d * CORE::BigFloat::exp2( static_cast<int>(exp * log_limb) );
}


#endif



template<class NT> 
inline NT inexact_sqrt_implementation( NT const& n, CGAL::Null_functor no_sqrt )
{

#ifdef CGAL_USE_CORE

  CORE::BigFloat nn = to_BigFloat(n) ;
  CORE::BigFloat s  = CORE::sqrt(nn);
  return NT(s.doubleValue());
  
#else

  double nn = CGAL::to_double(n) ;
  
  if ( !CGAL_NTS is_valid(nn) || ! CGAL_NTS is_finite(nn) )
    nn = std::numeric_limits<double>::max BOOST_PREVENT_MACRO_SUBSTITUTION () ;
    
  CGAL_precondition(nn > 0);
  
  double s = CGAL_NTS sqrt(nn);
  
  return NT(s);
  
#endif
}

template<class NT, class Sqrt>
inline NT inexact_sqrt_implementation( NT const& n, Sqrt sqrt_f )
{
  return sqrt_f(n);
}

template<class NT>
inline NT inexact_sqrt( NT const& n )
{
  typedef CGAL::Algebraic_structure_traits<NT> AST;
  typedef typename AST::Sqrt Sqrt; 
  return inexact_sqrt_implementation(n,Sqrt());
}

inline Quotient<MP_Float> inexact_sqrt( Quotient<MP_Float> const& q )
{
  CGAL_precondition(q > 0);
  return Quotient<MP_Float>(CGAL_SS_i::inexact_sqrt(q.numerator()*q.denominator()), q.denominator() );
}


// Given an oriented 2D straight line segment 'e', computes the normalized coefficients (a,b,c) of the
// supporting line.
// POSTCONDITION: [a,b] is the leftward normal _unit_ (a^2+b^2=1) vector.
// POSTCONDITION: In case of overflow, an empty optional<> is returned.
template<class Segment_2>
optional< typename Kernel_traits<Segment_2>::Kernel::Line_2 > compute_normalized_line_ceoffC2( Segment_2 const& e )
{
  bool finite = true ;
  
  typedef typename Kernel_traits<Segment_2>::Kernel K;
  
  typedef typename K::FT FT ;
  
  FT a (0.0),b (0.0) ,c(0.0)  ;

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

    CGAL_STSKEL_TRAITS_TRACE("Line coefficients for HORIZONTAL line:\n" 
                            << s2str(e) 
                            << "\na="<< n2str(a) << ", b=" << n2str(b) << ", c=" << n2str(c)
                           ) ;
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

    CGAL_STSKEL_TRAITS_TRACE("Line coefficients for VERTICAL line:\n"
                            << s2str(e) 
                            << "\na="<< n2str(a) << ", b=" << n2str(b) << ", c=" << n2str(c)
                           ) ;
  }
  else
  {
    FT sa = e.source().y() - e.target().y();
    FT sb = e.target().x() - e.source().x();
    FT l2 = (sa*sa) + (sb*sb) ;

    if ( CGAL_NTS is_finite(l2) )
    {
      FT l = CGAL_SS_i :: inexact_sqrt(l2);
  
      a = sa / l ;
      b = sb / l ;
  
      c = -e.source().x()*a - e.source().y()*b;
      
      CGAL_STSKEL_TRAITS_TRACE("Line coefficients for line:\n"
                               << s2str(e) 
                               << "\nsa="<< n2str(sa) << "\nsb=" << n2str(sb) << "\nl2=" << n2str(l2) << "\nl=" << n2str(l)
                               << "\na="<< n2str(a) << "\nb=" << n2str(b) << "\nc=" << n2str(c)
                               ) ;
    }
    else finite = false ;
    
  }
  
  if ( finite )
    if ( !CGAL_NTS is_finite(a) || !CGAL_NTS is_finite(b) || !CGAL_NTS is_finite(c) ) 
      finite = false ;

  return cgal_make_optional( finite, K().construct_line_2_object()(a,b,c) ) ;
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
template<class Trisegment_2,class Segment_2, class FT>
typename Trisegment_2::Self_ptr construct_trisegment ( Segment_2 const& e0
                                                     , FT        const& w0
                                                     , Segment_2 const& e1
                                                     , FT        const& w1
                                                     , Segment_2 const& e2
                                                     , FT        const& w2
                                                     )
{
  typedef typename Trisegment_2::Self_ptr Trisegment_2_ptr ;
  
  Uncertain<Trisegment_collinearity> lCollinearity = certified_trisegment_collinearity(e0,w0,e1,w1,e2,w2);
  
  if (is_certain(lCollinearity) ) 
       return Trisegment_2_ptr( new Trisegment_2(e0, w0, e1, w1, e2, w2, lCollinearity) ) ;
  else return Trisegment_2_ptr();
}

// Given 3 oriented moving straight line equations: e0, e1, e2
// en : an*x(t) + bn*y(t) + cn - wn*t = 0
// returns the "time" (n/d) at which the weighted offsetted lines
// intersect at a single point, IFF such intersection exist.
// If the lines intersect on the common offsetting side, the returned distance is positive, otherwise it is negative.
// The offsetting side of a line is given by the sign of its assigned weight: positive weight is to the left, negative to the right.
//
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
//
// If the lines do not intersect, for example, for collinear edges, or parallel edges but with the same orientation,
// returns den=0
//
// PRECONDITION: None of e0, e1 and e2 are collinear (but two of them can be parallel)
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//

template<class FT >
optional< Rational< FT > > linsolve_timeC2( FT const & a0, FT const & b0, FT const & c0, FT const & w0
                                          , FT const & a1, FT const & b1, FT const & c1, FT const & w1
                                          , FT const & a2, FT const & b2, FT const & c2, FT const & w2
                                          )
{ 
  // DETAILS:
  //
  // A weighted offset line is given by:
  //
  //   a*x(t) + b*y(t) + c - w*t = 0
  //
  // were 'w*t > 0' being to the left of the line.
  // If 3 such offset lines intersect at the same offset distance, the intersection 't',
  // or 'time', can be computed solving for 't' in the linear system formed by 3 such equations.
  // The result is :
  //
  //      (b2*a1 - a2*b1)*c0 + (a2*b0- b2*a0)*c1 + (b1*a0 - b0*a1)*c2
  //  t = ------------------------------------------------------------
  //      (b2*a1 - a2*b1)*w0 + (a2*b0- b2*a0)*w1 + (b1*a0 - b0*a1)*w2 
 
  FT num = ( b2 * a1 - a2 * b1 ) * c0
          +( a2 * b0 - b2 * a0 ) * c1
          +( b1 * a0 - b0 * a1 ) * c2;
  
  FT den = ( b2 * a1 - a2 * b1 ) * w0
          +( a2 * b0 - b2 * a0 ) * w1
          +( b1 * a0 - b0 * a1 ) * w2;

  bool ok = CGAL_NTS is_finite(num) && CGAL_NTS is_finite(den);     

  return cgal_make_optional(ok,Rational<FT>(num,den)) ;
}


// Given 3 oriented moving straight line equations: e0, e1, e2
// en : an*x(t) + bn*y(t) + cn - wn*t = 0
// returns the "point" pair(x,y) at which the weighted offsetted lines
// intersect, IFF such intersection exist.
// If the lines intersect on the common offsetting side, the returned distance is positive, otherwise it is negative.
// The offsetting side of a line is given by the sign of its assigned weight: positive weight is to the left, negative to the right.
//
// If the lines do not intersect, for example, for collinear edges, or parallel edges but with the same orientation,
// returns den=0
//
// PRECONDITION: None of e0, e1 and e2 are collinear (but two of them can be parallel)
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//

template<class FT >
optional< std::pair< FT, FT > >
  linsolve_pointC2( FT const & a0, FT const & b0, FT const & c0, FT const & w0
                  , FT const & a1, FT const & b1, FT const & c1, FT const & w1
                  , FT const & a2, FT const & b2, FT const & c2, FT const & w2
                  )
{
  bool ok = false;
  FT x(0.0), y(0.0);

  FT den = ( b2 * a1 - a2 * b1 ) * w0
          +( a2 * b0 - b2 * a0 ) * w1
          +( b1 * a0 - b0 * a1 ) * w2 ;  
  
  CGAL_STSKEL_TRAITS_TRACE("den=" << n2str(den) )
  
  if ( ! CGAL_NTS certified_is_zero(den) ) 
  {
    FT numX = ( b2 * w1 - b1 * w2 ) * c0
             +( b0 * w2 - b2 * w0 ) * c1
             +( b1 * w0 - b0 * w1 ) * c2 ;
    
    FT numY = ( w1 * a2 - w2 * a1 ) * c0
             +( w2 * a0 - w0 * a2 ) * c1
             +( w0 * a1 - w1 * a0 ) * c2 ;

    if ( CGAL_NTS is_finite(den) && CGAL_NTS is_finite(numX) && CGAL_NTS is_finite(numY)  )
    {
      ok = true ;
        
      x =  numX / den ;
      y = -numY / den ;
        
      CGAL_STSKEL_TRAITS_TRACE("numX=" << n2str(numX) << "\nnumY=" << n2str(numY)
                               << "\nx=" << n2str(x) << "\ny=" << n2str(y)
                               )
    }
  }   

  return cgal_make_optional(ok,std::make_pair(x,y)) ;
}


// Given 3 oriented straight line segments: e0, e1, e2 
// returns the "time" (n/d) at which the weighted offsetted lines
// intersect at a single point, IFF such intersection exist.
// If the lines intersect on the common offsetting side, the returned distance is positive, otherwise it is negative.
// The offsetting side of a line is given by the sign of its assigned weight: positive weight is to the left, negative to the right.
//
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
//
// If the lines do not intersect, for example, for collinear edges, or parallel edges but with the same orientation,
// returns den=0
//
// PRECONDITION: None of e0, e1 and e2 are collinear (but two of them can be parallel)
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//

template<class Trisegment_2 >
optional< Rational< typename Kernel_traits<Trisegment_2>::Kernel::FT > > compute_normal_offset_lines_isec_timeC2 ( intrusive_ptr< Trisegment_2 > const& tri )
{
  typedef typename Kernel_traits<Trisegment_2>::Kernel K ;
  
  typedef typename K::FT  FT ;
  
  typedef typename K::Line_2 Line_2 ;
  
  typedef optional<Line_2> Optional_line_2 ;
  
  CGAL_STSKEL_TRAITS_TRACE("Computing normal offset lines isec time for: " << tri ) ;
   
  // DETAILS:
  //
  // A weighted offset line is given by:
  //
  //   a*x(t) + b*y(t) + c - w*t = 0
  //
  // were 'w*t > 0' being to the left of the line.
  // If 3 such offset lines intersect at the same offset distance, the intersection 't',
  // or 'time', can be computed solving for 't' in the linear system formed by 3 such equations.
 
  optional< Rational<FT> > result;

  Optional_line_2 l0 = compute_normalized_line_ceoffC2(tri->e0()) ;
  Optional_line_2 l1 = compute_normalized_line_ceoffC2(tri->e1()) ;
  Optional_line_2 l2 = compute_normalized_line_ceoffC2(tri->e2()) ;

  if ( l0 && l1 && l2 ) 
  {
    result = linsolve_timeC2( l0->a(), l0->b(), l0->c(), tri->w0()
                            , l1->a(), l1->b(), l1->c(), tri->w1()
                            , l2->a(), l2->b(), l2->c(), tri->w2()
                            );
  }

  CGAL_STSKEL_TRAITS_TRACE("Event time (normal): " << 
    ( result ? ("n=" << result->n() << " d=" << result->d() << " n/d=" << *result) : "none" ) );

  return result ;
}

// Given two oriented straight line segments e0 and e1, which are known to be along the same line,
// but which can be oriented as e0->e1 or e1->e0, return the the midpoint of the segment between 
// the two.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class Segment_2>
optional< typename Kernel_traits<Segment_2>::Kernel::Point_2 > compute_oriented_midpoint ( Segment_2 const& e0, Segment_2 const& e1 )
{
  bool ok = false ;
  
  typedef typename Kernel_traits<Segment_2>::Kernel K ;
  
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;
  
  FT delta01 = CGAL::squared_distance(e0.target(),e1.source());
  FT delta10 = CGAL::squared_distance(e1.target(),e0.source());
  
  Point_2 mp ;
   
  if ( CGAL_NTS is_finite(delta01) &&  CGAL_NTS is_finite(delta10) )
  {
    if ( delta01 <= delta10 )
         mp = CGAL::midpoint(e0.target(),e1.source());
    else mp = CGAL::midpoint(e1.target(),e0.source());
    
    CGAL_STSKEL_TRAITS_TRACE("Computing oriented midpoint between:"
                             << "\ne0: " << s2str(e0)
                             << "\ne1: " << s2str(e1)
                             << "\nmp=" << p2str(mp)
                            );
                            
    ok = CGAL_NTS is_finite(mp.x()) && CGAL_NTS is_finite(mp.y());
  }
  
  return cgal_make_optional(ok,mp);
}


//
// Given 3 oriented straight line segments: e0, e1, e2 and the corresponding offseted segments: e0*, e1* and e2*,
// returns the point of the left or right seed (offset vertex). That is: (e0*,e1*) or (e1*,e2*) 
// 
// If the current event (defined by e0,e1,e2) is a propagated event, that is, it follows from a previous event,
// the seeds are skeleten nodes and are given by non-null trisegments.
// If the current event is an initial event the seeds are contour vertices and are given by null trisegmets.
//
// If a seed is a skeleton node, its point has to be computed from the trisegment that defines it.
// That trisegment is exactly the trisegment tree that defined the previous event which produced the skeleton node
// (so the trisegment tree is basically a lazy representation of the seed point).
//
// If a seed is a contour vertex, its point is then simply the target endoint of e0 or e1 (for the left/right seed).
//
// This method returns the specified seed point (left or right)
//
// NOTE: Split events involve 3 edges but only one seed, the left (that is, only e0*,e1* is connected before the event).
// The trisegment tree for a split event has always a null right child even if the event is not an initial event
// (in which case its left child won't be null).
// If you ask for the right child point for a trisegment tree corresponding to a split event you will just get e1.target()
// which is nonsensical for a non initial split event.
//
// NOTE: There is an abnormal collinearity case which ocurrs when e0 and e2 are collinear.
// In this case, these lines do not correspond to an offset vertex (because e0* and e2* are never consecutive before the event),
// so the degenerate seed is neither the left or the right seed. In this case, the SEED ID for the degenerate pseudo seed is UNKOWN.
// If you request the point of such degenerate pseudo seed the oriented midpoint bettwen e0 and e2 is returned.
//
// NOTE 2: Terminal vertices are given a redundant trisegment: [eA,eA,eB] for a left terminal, [eA,eB,eB] for a right terminal
// ands require the special handling seend below.
template<class Trisegment_2>
optional< typename Kernel_traits<Trisegment_2>::Kernel::Point_2 > compute_seed_pointC2 ( intrusive_ptr< Trisegment_2 > const& tri, typename Trisegment_2::SEED_ID sid )
{
  typedef typename Kernel_traits<Trisegment_2>::Kernel K ;
  
  typedef typename K::Point_2 Point_2 ;
  
  optional< Point_2 > p ;

  switch ( sid )
  {
    case Trisegment_2::LEFT :
         
       p = tri->child_l() ? construct_offset_lines_isecC2(tri->child_l())  // this can recurse
                          : ( are_edges_coincident(tri->e0(),tri->e1()) ? tri->e0().source() // left-terminal (degree-1) vertex 
                                                                        : compute_oriented_midpoint(tri->e0(),tri->e1()) // degree-2 colinear vertex
                            ) ;
       break ;                     
             
    case Trisegment_2::RIGHT : 
    
      p = tri->child_r() ? construct_offset_lines_isecC2(tri->child_r()) // this can recurse
                         : ( are_edges_coincident(tri->e1(),tri->e2()) ? tri->e1().target() // right-terminal (degree-1) vertex 
                                                                       : compute_oriented_midpoint(tri->e1(),tri->e2())  // degree-2 colinear vertex
                           ) ;
      break ;        
             
    case Trisegment_2::UNKNOWN : 
    
      p = compute_oriented_midpoint(tri->e0(),tri->e2());
      
      break ;
  }
  
  return p ;  
}

//
// Given the trisegment tree for an event which is known to have a normal collinearity, returns the seed point
// of the degenerate seed.
// A normal collinearity occurs when e0,e1 or e1,e2 are collinear
// (but not e0,e2 which ocurrs when edges *become* collinear after an event)
template<class Trisegment_2>
optional< typename Kernel_traits<Trisegment_2>::Kernel::Point_2 > compute_degenerate_seed_pointC2 ( intrusive_ptr< Trisegment_2 > const& tri )
{
  return compute_seed_pointC2( tri, tri->degenerate_seed_id() ) ;
}

// Given 3 oriented straight line segments: e0, e1, e2 
// such that two and only two of these edges are collinear, not neccesarily consecutive but with the same orientaton;
// returns the "time" (n/d) at which a line perpendicular to the collinear edge passing through
// the degenerate seed point intersects the weighted bisector between the lines supporting the collinear and none collinear edges.
//
// NOTE: The result is a explicit rational number returned as a tuple (num,den); the caller must check that den!=0 manually
// (a predicate for instance should return indeterminate in this case)
//
// PRECONDITION: The collinear edges have the same weight.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class Trisegment_2>
optional< Rational< typename Kernel_traits<Trisegment_2>::Kernel::FT > > compute_degenerate_offset_lines_isec_timeC2 ( intrusive_ptr< Trisegment_2 > const& tri )
{
  typedef typename Kernel_traits<Trisegment_2>::Kernel K ;
  
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Line_2  Line_2 ;
  
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  CGAL_STSKEL_TRAITS_TRACE("Computing degenerate offset lines isec time for: " << tri ) ;

  // DETAILS:
  //
  // For simplicity, assume e0,e1 are the collinear edges.
  //
  //   (1)
  //   The bisecting line of e0 and e1 is a line perpendicular to e0 (and e1)
  //   which passes through the degenerate offset vertex q=(e0*,e1*)
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
  //   where (l0.a,l0.b,l0.c) and (l2.a,l2.b,l0.c) are the normalized line coefficientes of e0 and e2 resp.
  //
  //     B1(t)=[x(t),y(t)]
  //
  //   (3)
  //   These two bisecting lines B0(t) and B1(t) intersect (if they do) in a single point 'x' whose distance
  //   to the lines supporting the 3 edges is exactly 't' (since those expressions are precisely parametrized in a distance)
  //   Solving the following vectorial equation:
  //
  //     [x(t),y(t)] = p + w0*t*[l0.a,l0.b]
  //
  //   for t gives the result we want.
  //
  //   (4)
  //   To solve this system we can reformulate it as the intersection of three non colinear lines.
  //   The two equations from the non colinear edges:
  //
  //    l0.a*x(t) + l0.b*y(t) + l0.c - w0*t = 0
  //    l2.a*x(t) + l2.b*y(t) + l2.c - w2*t = 0
  //
  //   and if l0 is not vertical, we add:
  //
  //     x(t) = px + w0*l0.a* t   ==>   (1.0)*x(t) +  (0.0)*y(t) + (-px) - (l0.a*w0)*t = 0
  //
  //   else, because l0 is not horizontal:
  //
  //     y(t) = py + w0*l0.b* t   ==>   (0.0)*x(t) +  (1.0)*y(t) + (-py) - (l0.b*w0)*t = 0
  //
  
  optional< Rational<FT> > result;

  Optional_line_2 l0 = compute_normalized_line_ceoffC2(tri->collinear_edge_A   ()) ;
  Optional_line_2 l2 = compute_normalized_line_ceoffC2(tri->non_collinear_edge()) ;

  Optional_point_2 q = compute_degenerate_seed_pointC2(tri);

  if ( l0 && l2 && q )
  {
    FT px, py ;
    line_project_pointC2(l0->a(),l0->b(),l0->c(),q->x(),q->y(),px,py); 
    
    CGAL_STSKEL_TRAITS_TRACE("Seed point: " << p2str(*q) << ".\nProjected seed point: (" << n2str(px) << "," << n2str(py) << ")" ) ;
    
    FT abs_a  = CGAL_NTS abs (l0->a()) ;
    FT abs_b  = CGAL_NTS abs (l0->b()) ;
    
    result = ( abs_a < abs_b ) ? linsolve_timeC2( l0->a()         , l0->b(), l0->c(),         tri->w0() // Not vertical
                                                , FT(1.0), FT(0.0),     -px, l0->a() * tri->w0()
                                                , l2->a()         , l2->b(), l2->c(),         tri->w2()
                                                )
                               : linsolve_timeC2( l0->a(), l0->b()         , l0->c(),         tri->w0() // Not horizontal
                                                , FT(0.0), FT(1.0),     -py, l0->b() * tri->w0()
                                                , l2->a(), l2->b()         , l2->c(),         tri->w2()
                                                );
  }

  CGAL_STSKEL_TRAITS_TRACE("Event time (normal): " << 
    ( result ? ("n=" << result->n() << " d=" << result->d() << " n/d=" << *result) : "none" ) );
 
  return result ;
}

//
// Calls the appropiate function depending on the collinearity of the edges.
//
template<class Trisegment_2>
optional< Rational< typename Kernel_traits<Trisegment_2>::Kernel::FT > > compute_offset_lines_isec_timeC2 ( intrusive_ptr< Trisegment_2 > const& tri )
{
  CGAL_precondition ( tri->collinearity() != TRISEGMENT_COLLINEARITY_ALL ) ;
 
  return tri->collinearity() == TRISEGMENT_COLLINEARITY_NONE ? compute_normal_offset_lines_isec_timeC2    (tri)
                                                             : compute_degenerate_offset_lines_isec_timeC2(tri);
}


// Given 3 oriented line segments e0, e1 and e2 
// such that their offsets at a certian time intersect in a single point, 
// returns the coordinates (x,y) of such a point.
//
// PRECONDITIONS:
// None of e0, e1 and e2 are collinear (but two of them can be parallel)
// The line coefficients must be normalized: a^2+b^2==1 and (a,b) being the leftward normal vector
// The offsets at a certain time do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class Trisegment_2>
optional< typename Kernel_traits<Trisegment_2>::Kernel::Point_2 > construct_normal_offset_lines_isecC2 ( intrusive_ptr< Trisegment_2 > const& tri )
{
  typedef typename Kernel_traits<Trisegment_2>::Kernel K ;
  
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Line_2  Line_2 ;
  
  typedef optional<Line_2>  Optional_line_2 ;
  
  CGAL_STSKEL_TRAITS_TRACE("Computing normal offset lines isec point for: " << tri ) ;
   
  Optional_line_2 l0 = compute_normalized_line_ceoffC2(tri->e0()) ;
  Optional_line_2 l1 = compute_normalized_line_ceoffC2(tri->e1()) ;
  Optional_line_2 l2 = compute_normalized_line_ceoffC2(tri->e2()) ;
 
  optional< std::pair<FT,FT> > ip;

  if ( l0 && l1 && l2 )
  {
    ip = linsolve_pointC2( l0->a(), l0->b(), l0->c(), tri->w0()
                         , l1->a(), l1->b(), l1->c(), tri->w1()
                         , l2->a(), l2->b(), l2->c(), tri->w2()
                         );
                         
                         
  }
 
  return cgal_make_optional(ip,K().construct_point_2_object()(ip->first,ip->second)) ;
}

// Given 3 oriented line segments e0, e1 and e2
// such that their offsets at a certian time intersect in a single point, 
// and two and only two of the edges are collinear, not neccesarily consecutive but with the same orientaton,
// returns the coordinates (x,y) of such a point.
//
// PRECONDITIONS:
// The line coefficients must be normalized: a^2+b^2==1 and (a,b) being the leftward normal vector
// The offsets at a certain time do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class Trisegment_2>
optional< typename Kernel_traits<Trisegment_2>::Kernel::Point_2 > construct_degenerate_offset_lines_isecC2 ( intrusive_ptr< Trisegment_2 > const& tri )
{
  typedef typename Kernel_traits<Trisegment_2>::Kernel K ;
  
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Line_2  Line_2 ;
  
  typedef optional<Point_2>  Optional_point_2 ;
  typedef optional<Line_2>   Optional_line_2 ;
  
  CGAL_STSKEL_TRAITS_TRACE("Computing degenerate offset lines isec point for: " << tri )  ;
  
  Optional_line_2 l0 = compute_normalized_line_ceoffC2(tri->collinear_edge_A   ()) ;
  Optional_line_2 l2 = compute_normalized_line_ceoffC2(tri->non_collinear_edge()) ;

  Optional_point_2 q = compute_degenerate_seed_pointC2(tri);

  optional< std::pair<FT,FT> > ip;

  if ( l0 && l2 && q )
  {
    FT px, py ;
    line_project_pointC2(l0->a(),l0->b(),l0->c(),q->x(),q->y(),px,py); 
    
    CGAL_STSKEL_TRAITS_TRACE("Seed point: " << p2str(*q) << ".\nProjected seed point: (" << n2str(px) << "," << n2str(py) << ")" ) ;
    
    FT abs_a  = CGAL_NTS abs (l0->a()) ;
    FT abs_b  = CGAL_NTS abs (l0->b()) ;
    
    ip = ( abs_a < abs_b ) ? linsolve_pointC2( l0->a()         , l0->b(), l0->c(),         tri->w0() // Not vertical
                                             , FT(1.0), FT(0.0),     -px, l0->a() * tri->w0()
                                             , l2->a()         , l2->b(), l2->c(),         tri->w2()
                                             )
                           : linsolve_pointC2( l0->a(), l0->b()         , l0->c(),         tri->w0() // Not horizontal
                                             , FT(0.0), FT(1.0),     -py, l0->b() * tri->w0()
                                             , l2->a(), l2->b()         , l2->c(),         tri->w2()
                                             );
                                             
  }
  return cgal_make_optional(!!ip,K().construct_point_2_object()(ip->first,ip->second)) ;
}

//
// Calls the appropiate function depending on the collinearity of the edges.
//
template<class Trisegment_2>
optional< typename Kernel_traits<Trisegment_2>::Kernel::Point_2 > construct_offset_lines_isecC2 ( intrusive_ptr< Trisegment_2 > const& tri )
{
  CGAL_precondition ( tri->collinearity() != TRISEGMENT_COLLINEARITY_ALL ) ;
  
  return tri->collinearity() == TRISEGMENT_COLLINEARITY_NONE ? construct_normal_offset_lines_isecC2    (tri)
                                                             : construct_degenerate_offset_lines_isecC2(tri) ;
}

} // namnepsace CGAIL_SS_i

} //namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_CONS_FTC2_H //
// EOF //
