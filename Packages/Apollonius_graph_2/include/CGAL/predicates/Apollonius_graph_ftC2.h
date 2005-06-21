// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_APOLLONIUS_GRAPH_FTC2_H
#define CGAL_APOLLONIUS_GRAPH_FTC2_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/predicates/check_filter.h>
#include <CGAL/functions_on_signs.h>
#include <CGAL/predicates/Apollonius_graph_predicates_C2.h>

#include <CGAL/Apollonius_graph_kernel_wrapper_2.h>

CGAL_BEGIN_NAMESPACE

template< class RT >
inline
Orientation
ag2_orientation_test_C2(const RT &x1, const RT &y1, const RT &w1,
			const RT &x2, const RT &y2, const RT &w2,
			const RT &x3, const RT &y3, const RT &w3)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                    Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>  Kernel;
  typedef typename Kernel::Point_2                Point_2;
  typedef typename Kernel::Site_2                 Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);

  Ag2_orientation_2<Kernel> f;
  return f(p1, p2, p3);
}

//--------------------------------------------------------------------

template< class RT >
inline
bool
ag_is_hidden_test_ring_C2(const RT &x1, const RT &y1, const RT &w1,
			  const RT &x2, const RT &y2, const RT &w2)
{
#ifdef AG2_PROFILE_PREDICATES
  ag2_predicate_profiler::is_trivial_counter++;
#endif
  must_be_filtered(x1);

  Sign s = CGAL::sign( CGAL::square(x1 - x2) + CGAL::square(y1 - y2)
		       - CGAL::square(w1 - w2)
		       );
  if ( s == POSITIVE ) { return false; }
  return (CGAL::compare(w1, w2) != SMALLER);
}


template< class RT >
inline
bool
ag_is_hidden_test_sqrtf_C2(const RT &x1, const RT &y1, const RT &w1,
			   const RT &x2, const RT &y2, const RT &w2)
{
#ifdef AG2_PROFILE_PREDICATES
  ag2_predicate_profiler::is_trivial_counter++;
#endif
  RT d = CGAL::sqrt(CGAL::square(x1 - x2) + CGAL::square(y1 - y2));
  Sign s = CGAL::sign(d - w1 + w2);

  return ( s != POSITIVE );
}

//--------------------------------------------------------------------



template< class RT >
CGAL_MEDIUM_INLINE
Comparison_result
compare_ag_distances_test_ring_C2(const RT &x1, const RT &y1, const RT &w1,
				  const RT &x2, const RT &y2, const RT &w2,
				  const RT & x, const RT & y)
{
#ifdef AG2_PROFILE_PREDICATES
  ag2_predicate_profiler::side_of_bisector_counter++;
#endif
  must_be_filtered(x1);

  // this function compares the distances of the point(x, y) from the 
  // disks {(x1, y1), w1} and {(x2, y2), w2}
  RT D1 = CGAL::square(x1 - x) + CGAL::square(y1 - y);
  RT D2 = CGAL::square(x2 - x) + CGAL::square(y2 - y);
  RT Dw = w2 - w1;

  Sign sign_of_Dw = CGAL::sign(Dw);
  Comparison_result R = CGAL::compare(D1, D2);

  if ( sign_of_Dw == ZERO ) {
    return R;
  }
  if ( sign_of_Dw == POSITIVE ) {
    if ( R != SMALLER )  return LARGER;

    Sign s = sign_a_plus_b_x_sqrt_c(D1 - D2 + CGAL::square(Dw),
				    RT(2) * Dw, D1);
    return ((s == POSITIVE) ? LARGER : ((s == ZERO) ? EQUAL : SMALLER));
  }

  if ( R != LARGER )  return SMALLER;
  Sign s = sign_a_plus_b_x_sqrt_c(D1 - D2 - CGAL::square(Dw),
				  RT(2) * Dw, D2);

  return ((s == POSITIVE) ? LARGER : ((s == ZERO) ? EQUAL : SMALLER));
}

//--------------------------------------------------------------------


template< class RT >
CGAL_MEDIUM_INLINE
Comparison_result
compare_ag_distances_test_sqrtf_C2(const RT &x1, const RT &y1, const RT &w1,
				   const RT &x2, const RT &y2, const RT &w2,
				   const RT & x, const RT & y)
{
#ifdef AG2_PROFILE_PREDICATES
  ag2_predicate_profiler::side_of_bisector_counter++;
#endif
  // this function compares the distances of the point(x, y) from the 
  // disks {(x1, y1), w1} and {(x2, y2), w2}

  RT D1 = CGAL::square(x1 - x) + CGAL::square(y1 - y);
  RT D2 = CGAL::square(x2 - x) + CGAL::square(y2 - y);

  RT d1 = CGAL::sqrt(D1) - w1;
  RT d2 = CGAL::sqrt(D2) - w2;

  return CGAL::compare(d1, d2);
}

//--------------------------------------------------------------------

template< class RT >
/*CGAL_NO_FILTER*/
Bounded_side
bounded_side_of_segment(const RT& x1, const RT& y1,
			const RT& x2, const RT& y2,
			const RT&  x, const RT&  y)
{
  // this function tests whether the point (x, y) is in the interior
  // of the segment defined by (x1, y1) and (x2, y2).
  // It is assumed that the three points are collinear.
  // The result is ON_BOUNDED_SIDE if the point (x,y) is in the
  // interior of the segment, ON_BOUNDARY if (x,y) coincides with
  // either (x1, y1) or (x1, y1) and ON_UNBOUNDED_SIDE otherwise.

  CGAL_precondition( sign_of_determinant2x2(x1-x,y1-y,x2-x,y2-y) == ZERO );

  Comparison_result rx1 = CGAL::compare(x, x1);
  Comparison_result ry1 = CGAL::compare(y, y1);

  if ( rx1 == EQUAL && ry1 == EQUAL ) { return ON_BOUNDARY; }

  Comparison_result rx2 = CGAL::compare(x, x2);
  Comparison_result ry2 = CGAL::compare(y, y2);

  if ( rx2 == EQUAL && ry2 == EQUAL ) { return ON_BOUNDARY; }

  Comparison_result rx12 = CGAL::compare(x1, x2);

  if ( rx12 == SMALLER ) {
    CGAL_assertion( rx1 != EQUAL && rx2 != EQUAL );
    return ( rx1 == LARGER && rx2 == SMALLER ) ? ON_BOUNDED_SIDE :
      ON_UNBOUNDED_SIDE;
  } else if ( rx12 == LARGER ) {
    CGAL_assertion( rx1 != EQUAL && rx2 != EQUAL );
    return ( rx1 == SMALLER && rx2 == LARGER ) ? ON_BOUNDED_SIDE :
      ON_UNBOUNDED_SIDE;
  }

  Comparison_result ry12 = CGAL::compare(y1, y2);
  CGAL_assertion( ry12 != EQUAL );
  CGAL_assertion( ry1 != EQUAL && ry2 != EQUAL );

  if ( ry12 == SMALLER ) {
    return ( ry1 == LARGER && ry2 == SMALLER ) ? ON_BOUNDED_SIDE :
      ON_UNBOUNDED_SIDE;
  }

  CGAL_assertion( ry12 == LARGER );
  return ( ry1 == SMALLER && ry2 == LARGER ) ? ON_BOUNDED_SIDE :
    ON_UNBOUNDED_SIDE;
}


//--------------------------------------------------------------------

template < class RT >
Sign
ag_incircle_test_sqrtf_C2(const RT &x1, const RT &y1,
			  const RT &w1,
			  const RT &x2, const RT &y2,
			  const RT &w2,
			  const RT &qx, const RT &qy,
			  const RT &qw)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                    Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>  Kernel;
  typedef typename Kernel::Point_2                Point_2;
  typedef typename Kernel::Site_2                 Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2  q(Point_2(qx, qy), qw);

  Incircle_test<Kernel,Sqrt_field_tag> f;
  return f(p1, p2, q);
}


template < class RT >
Sign
ag_incircle_test_ring_C2(const RT &x1, const RT &y1,
			 const RT &w1,
			 const RT &x2, const RT &y2,
			 const RT &w2,
			 const RT &qx, const RT &qy,
			 const RT &qw)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2  q(Point_2(qx, qy), qw);

  Incircle_test<Kernel,Ring_tag> f;
  return f(p1, p2, q);
}




//--------------------------------------------------------------------

template < class RT >
Sign
ag_incircle_test_sqrtf_C2(const RT &x1, const RT &y1,
			  const RT &w1,
			  const RT &x2, const RT &y2,
			  const RT &w2,
			  const RT &x3, const RT &y3,
			  const RT &w3,
			  const RT &qx, const RT &qy,
			  const RT &qw)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                    Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>  Kernel;
  typedef typename Kernel::Point_2                Point_2;
  typedef typename Kernel::Site_2                 Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);
  Site_2  q(Point_2(qx, qy), qw);

  Incircle_test<Kernel,Sqrt_field_tag> f;
  return f(p1, p2, p3, q);
}


template < class RT >
Sign
ag_incircle_test_ring_C2(const RT &x1, const RT &y1,
			 const RT &w1,
			 const RT &x2, const RT &y2,
			 const RT &w2,
			 const RT &x3, const RT &y3,
			 const RT &w3,
			 const RT &qx, const RT &qy,
			 const RT &qw)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);
  Site_2  q(Point_2(qx, qy), qw);

  Incircle_test<Kernel,Ring_tag> f;
  return f(p1, p2, p3, q);
}



//--------------------------------------------------------------------


template < class RT >
bool
ag_finite_edge_test_sqrtf_C2(const RT &x1, const RT &y1,
			     const RT &w1,
			     const RT &x2, const RT &y2,
			     const RT &w2,
			     const RT &x3, const RT &y3,
			     const RT &w3,
			     const RT &x4, const RT &y4,
			     const RT &w4,
			     const RT &qx, const RT &qy,
			     const RT &qw, bool b)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);
  Site_2 p4(Point_2(x4, y4), w4);
  Site_2  q(Point_2(qx, qy), qw);

  Ag2_finite_edge_test_C2<Kernel,Sqrt_field_tag> f;
  return f(p1, p2, p3, p4, q, b);
}

template < class RT >
bool
ag_finite_edge_test_ring_C2(const RT &x1, const RT &y1,
			    const RT &w1,
			    const RT &x2, const RT &y2,
			    const RT &w2,
			    const RT &x3, const RT &y3,
			    const RT &w3,
			    const RT &x4, const RT &y4,
			    const RT &w4,
			    const RT &qx, const RT &qy,
			    const RT &qw, bool b)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);
  Site_2 p4(Point_2(x4, y4), w4);
  Site_2  q(Point_2(qx, qy), qw);

  Ag2_finite_edge_test_C2<Kernel,Ring_tag> f;
  return f(p1, p2, p3, p4, q, b);
}


//--------------------------------------------------------------------

template < class RT >
bool
ag_finite_edge_test_degenerated_sqrtf_C2(const RT &x1, const RT &y1,
					 const RT &w1,
					 const RT &x2, const RT &y2,
					 const RT &w2,
					 const RT &qx, const RT &qy,
					 const RT &qw, bool b)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2  q(Point_2(qx, qy), qw);

  Ag2_finite_edge_test_C2<Kernel,Sqrt_field_tag> f;
  return f(p1, p2, q, b);
}

template < class RT >
bool
ag_finite_edge_test_degenerated_ring_C2(const RT &x1, const RT &y1,
					const RT &w1,
					const RT &x2, const RT &y2,
					const RT &w2,
					const RT &qx, const RT &qy,
					const RT &qw, bool b)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2  q(Point_2(qx, qy), qw);

  Ag2_finite_edge_test_C2<Kernel,Ring_tag> f;
  return f(p1, p2, q, b);
}


template < class RT >
bool
ag_finite_edge_test_degenerated_sqrtf_C2(const RT &x1, const RT &y1,
					 const RT &w1,
					 const RT &x2, const RT &y2,
					 const RT &w2,
					 const RT &x3, const RT &y3,
					 const RT &w3,
					 const RT &qx, const RT &qy,
					 const RT &qw, bool b)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);
  Site_2  q(Point_2(qx, qy), qw);

  Ag2_finite_edge_test_C2<Kernel,Sqrt_field_tag> f;
  return f(p1, p2, p3, q, b);
}


template < class RT >
bool
ag_finite_edge_test_degenerated_ring_C2(const RT &x1, const RT &y1,
					const RT &w1,
					const RT &x2, const RT &y2,
					const RT &w2,
					const RT &x3, const RT &y3,
					const RT &w3,
					const RT &qx, const RT &qy,
					const RT &qw, bool b)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);
  Site_2  q(Point_2(qx, qy), qw);

  Ag2_finite_edge_test_C2<Kernel,Ring_tag> f;
  return f(p1, p2, p3, q, b);
}

//--------------------------------------------------------------------

template < class RT >
bool
ag_infinite_edge_test_sqrtf_C2(const RT &x2, const RT &y2,
			       const RT &w2,
			       const RT &x3, const RT &y3,
			       const RT &w3,
			       const RT &x4, const RT &y4,
			       const RT &w4,
			       const RT &qx, const RT &qy,
			       const RT &qw, bool b)
{ 
  must_be_filtered(x2);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);
  Site_2 p4(Point_2(x4, y4), w4);
  Site_2  q(Point_2(qx, qy), qw);

  Infinite_edge_test<Kernel,Sqrt_field_tag> f;
  return f(p2, p3, p4, q, b);
}

template < class RT >
bool
ag_infinite_edge_test_ring_C2(const RT &x2, const RT &y2,
			      const RT &w2,
			      const RT &x3, const RT &y3,
			      const RT &w3,
			      const RT &x4, const RT &y4,
			      const RT &w4,
			      const RT &qx, const RT &qy,
			      const RT &qw, bool b)
{
  must_be_filtered(x2);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);
  Site_2 p4(Point_2(x4, y4), w4);
  Site_2  q(Point_2(qx, qy), qw);

  Infinite_edge_test<Kernel,Ring_tag> f;
  return f(p2, p3, p4, q, b);
}



//--------------------------------------------------------------------

template < class RT >
bool
ag_is_degenerate_edge_test_sqrtf_C2(const RT &x1, const RT &y1,
				    const RT &w1,
				    const RT &x2, const RT &y2,
				    const RT &w2,
				    const RT &x3, const RT &y3,
				    const RT &w3,
				    const RT &x4, const RT &y4,
				    const RT &w4)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);
  Site_2 p4(Point_2(x4, y4), w4);

  Is_degenerate_edge_test<Kernel,Sqrt_field_tag> f;
  return f(p1, p2, p3, p4);
}

template < class RT >
bool
ag_is_degenerate_edge_test_ring_C2(const RT &x1, const RT &y1,
				   const RT &w1,
				   const RT &x2, const RT &y2,
				   const RT &w2,
				   const RT &x3, const RT &y3,
				   const RT &w3,
				   const RT &x4, const RT &y4,
				   const RT &w4)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>                     Rep;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>   Kernel;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Kernel::Site_2                  Site_2;

  Site_2 p1(Point_2(x1, y1), w1);
  Site_2 p2(Point_2(x2, y2), w2);
  Site_2 p3(Point_2(x3, y3), w3);
  Site_2 p4(Point_2(x4, y4), w4);

  Is_degenerate_edge_test<Kernel,Ring_tag> f;
  return f(p1, p2, p3, p4);
}



//--------------------------------------------------------------------

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#ifndef CGAL_ARITHMETIC_FILTER_APOLLONIUS_GRAPH_FTC2_H
#include <CGAL/Arithmetic_filter/predicates/Apollonius_graph_ftC2.h>
#endif // CGAL_ARITHMETIC_FILTER_APOLLONIUS_GRAPH_FTC2_H
#endif

#endif // CGAL_APOLLONIUS_GRAPH_FTC2_H
