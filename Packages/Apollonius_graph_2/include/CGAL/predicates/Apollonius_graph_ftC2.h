#ifndef CGAL_APOLLONIUS_GRAPH_FTC2_H
#define CGAL_APOLLONIUS_GRAPH_FTC2_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/predicates/check_filter.h>
#include <CGAL/functions_on_signs.h>
#include <CGAL/predicates/Apollonius_graph_predicates_C2.h>

#include <CGAL/Apollonius_graph_method_tags.h>

CGAL_BEGIN_NAMESPACE

//--------------------------------------------------------------------

template< class RT >
inline
bool
ad_is_trivial_test_alg_C2(const RT &x1, const RT &y1, const RT &w1,
			  const RT &x2, const RT &y2, const RT &w2)
{
#ifdef AG2_PROFILE_PREDICATES
  ag2_predicate_profiler::is_trivial_counter++;
#endif
  must_be_filtered(x1);

  Sign s = CGAL_NTS sign(
			 CGAL_NTS square(x1 - x2)
			 + CGAL_NTS square(y1 - y2)
			 - CGAL_NTS square(w1 - w2)
			 );
  if ( s == POSITIVE ) { return false; }
  return (CGAL_NTS compare(w1, w2) != SMALLER);
}


template< class RT >
inline
bool
ad_is_trivial_test_naive_C2(const RT &x1, const RT &y1, const RT &w1,
			    const RT &x2, const RT &y2, const RT &w2)
{
#ifdef AG2_PROFILE_PREDICATES
  ag2_predicate_profiler::is_trivial_counter++;
#endif
  RT d = CGAL_NTS sqrt(CGAL_NTS square(x1 - x2)
		       + CGAL_NTS square(y1 - y2));
  Sign s = CGAL_NTS sign(d - w1 + w2);

  return ( s != POSITIVE );
}

//--------------------------------------------------------------------



template< class RT >
CGAL_MEDIUM_INLINE
Comparison_result
compare_ad_distances_test_alg1_C2(const RT &x1, const RT &y1, const RT &w1,
				  const RT &x2, const RT &y2, const RT &w2,
				  const RT & x, const RT & y)
{
#ifdef AG2_PROFILE_PREDICATES
  ag2_predicate_profiler::side_of_bisector_counter++;
#endif
  must_be_filtered(x1);

  // this function compares the distances of the point(x, y) from the 
  // disks {(x1, y1), w1} and {(x2, y2), w2}
  RT D1 = CGAL_NTS square(x1 - x) + CGAL_NTS square(y1 - y);
  RT D2 = CGAL_NTS square(x2 - x) + CGAL_NTS square(y2 - y);
  RT Dw = w2 - w1;

  Sign sign_of_Dw = CGAL_NTS sign(Dw);
  Comparison_result R = CGAL_NTS compare(D1, D2);

  if ( sign_of_Dw == ZERO ) {
    return R;
  }
  if ( sign_of_Dw == POSITIVE ) {
    if ( R != SMALLER )  return LARGER;

    Sign s = sign_a_plus_b_x_sqrt_c(D1 - D2 + CGAL_NTS square(Dw),
				    RT(2) * Dw, D1);
    return ((s == POSITIVE) ? LARGER : ((s == ZERO) ? EQUAL : SMALLER));
  }

  if ( R != LARGER )  return SMALLER;
  Sign s = sign_a_plus_b_x_sqrt_c(D1 - D2 - CGAL_NTS square(Dw),
				  RT(2) * Dw, D2);

  return ((s == POSITIVE) ? LARGER : ((s == ZERO) ? EQUAL : SMALLER));
}

template< class RT >
/*CGAL_NO_FILTER*/
inline
Comparison_result
compare_ad_distances_test_alg2_C2(const RT &x1, const RT &y1, const RT &w1,
				  const RT &x2, const RT &y2, const RT &w2,
				  const RT & x, const RT & y)
{
  return
    compare_ad_distances_test_alg2_C2(x1, y1, w1, x2, y2, w2, x, y);
}

//--------------------------------------------------------------------


template< class RT >
CGAL_MEDIUM_INLINE
Comparison_result
compare_ad_distances_test_naive_C2(const RT &x1, const RT &y1, const RT &w1,
				   const RT &x2, const RT &y2, const RT &w2,
				   const RT & x, const RT & y)
{
#ifdef AG2_PROFILE_PREDICATES
  ag2_predicate_profiler::side_of_bisector_counter++;
#endif
  // this function compares the distances of the point(x, y) from the 
  // disks {(x1, y1), w1} and {(x2, y2), w2}

  RT D1 = CGAL_NTS square(x1 - x) + CGAL_NTS square(y1 - y);
  RT D2 = CGAL_NTS square(x2 - x) + CGAL_NTS square(y2 - y);

  RT d1 = CGAL_NTS sqrt(D1) - w1;
  RT d2 = CGAL_NTS sqrt(D2) - w2;

  return CGAL_NTS compare(d1, d2);
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

  Comparison_result rx1 = CGAL_NTS compare(x, x1);
  Comparison_result ry1 = CGAL_NTS compare(y, y1);

  if ( rx1 == EQUAL && ry1 == EQUAL ) { return ON_BOUNDARY; }

  Comparison_result rx2 = CGAL_NTS compare(x, x2);
  Comparison_result ry2 = CGAL_NTS compare(y, y2);

  if ( rx2 == EQUAL && ry2 == EQUAL ) { return ON_BOUNDARY; }

  Comparison_result rx12 = CGAL_NTS compare(x1, x2);

  if ( rx12 == SMALLER ) {
    CGAL_assertion( rx1 != EQUAL && rx2 != EQUAL );
    return ( rx1 == LARGER && rx2 == SMALLER ) ? ON_BOUNDED_SIDE :
      ON_UNBOUNDED_SIDE;
  } else if ( rx12 == LARGER ) {
    CGAL_assertion( rx1 != EQUAL && rx2 != EQUAL );
    return ( rx1 == SMALLER && rx2 == LARGER ) ? ON_BOUNDED_SIDE :
      ON_UNBOUNDED_SIDE;
  }

  Comparison_result ry12 = CGAL_NTS compare(y1, y2);
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
ad_incircle_test_naive_C2(const RT &x1, const RT &y1,
			  const RT &w1,
			  const RT &x2, const RT &y2,
			  const RT &w2,
			  const RT &qx, const RT &qy,
			  const RT &qw)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point  q(Point(qx, qy), qw);

  Incircle_test<R> f;
  return f(p1, p2, q, Naive_tag() );
}


template < class RT >
Sign
ad_incircle_test_alg1_C2(const RT &x1, const RT &y1,
			 const RT &w1,
			 const RT &x2, const RT &y2,
			 const RT &w2,
			 const RT &qx, const RT &qy,
			 const RT &qw)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point  q(Point(qx, qy), qw);

  Incircle_test<R> f;
  return f(p1, p2, q, Algebraic1_tag() );
}


template < class RT >
Sign
ad_incircle_test_alg2_C2(const RT &x1, const RT &y1,
			 const RT &w1,
			 const RT &x2, const RT &y2,
			 const RT &w2,
			 const RT &qx, const RT &qy,
			 const RT &qw)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point  q(Point(qx, qy), qw);

  Incircle_test<R> f;
  return f(p1, p2, q, Algebraic2_tag() );
}




//--------------------------------------------------------------------

template < class RT >
Sign
ad_incircle_test_naive_C2(const RT &x1, const RT &y1,
			  const RT &w1,
			  const RT &x2, const RT &y2,
			  const RT &w2,
			  const RT &x3, const RT &y3,
			  const RT &w3,
			  const RT &qx, const RT &qy,
			  const RT &qw)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point  q(Point(qx, qy), qw);

  Incircle_test<R> f;
  return f(p1, p2, p3, q, Naive_tag() );
}


template < class RT >
Sign
ad_incircle_test_alg1_C2(const RT &x1, const RT &y1,
			 const RT &w1,
			 const RT &x2, const RT &y2,
			 const RT &w2,
			 const RT &x3, const RT &y3,
			 const RT &w3,
			 const RT &qx, const RT &qy,
			 const RT &qw)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point  q(Point(qx, qy), qw);

  Incircle_test<R> f;
  return f(p1, p2, p3, q, Algebraic1_tag() );
}


template < class RT >
Sign
ad_incircle_test_alg2_C2(const RT &x1, const RT &y1,
			 const RT &w1,
			 const RT &x2, const RT &y2,
			 const RT &w2,
			 const RT &x3, const RT &y3,
			 const RT &w3,
			 const RT &qx, const RT &qy,
			 const RT &qw)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point  q(Point(qx, qy), qw);

  Incircle_test<R> f;
  return f(p1, p2, p3, q, Algebraic2_tag() );
}




//--------------------------------------------------------------------


template < class RT >
bool
ad_finite_edge_test_naive_C2(const RT &x1, const RT &y1,
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
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point p4(Point(x4, y4), w4);
  Weighted_point  q(Point(qx, qy), qw);

  Finite_edge_test<R> f;
  return f(p1, p2, p3, p4, q, b, Naive_tag() );
}

template < class RT >
bool
ad_finite_edge_test_alg1_C2(const RT &x1, const RT &y1,
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

  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point p4(Point(x4, y4), w4);
  Weighted_point  q(Point(qx, qy), qw);

  Finite_edge_test<R> f;
  return f(p1, p2, p3, p4, q, b, Algebraic1_tag() );
}

template < class RT >
bool
ad_finite_edge_test_alg2_C2(const RT &x1, const RT &y1,
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
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point p4(Point(x4, y4), w4);
  Weighted_point  q(Point(qx, qy), qw);

  Finite_edge_test<R> f;
  return f(p1, p2, p3, p4, q, b, Algebraic2_tag() );
}

//--------------------------------------------------------------------

template < class RT >
bool
ad_finite_edge_test_degenerated_naive_C2(const RT &x1, const RT &y1,
					 const RT &w1,
					 const RT &x2, const RT &y2,
					 const RT &w2,
					 const RT &qx, const RT &qy,
					 const RT &qw, bool b)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point  q(Point(qx, qy), qw);

  Finite_edge_test_degenerated<R> f;
  return f(p1, p2, q, b, Naive_tag() );
}

template < class RT >
bool
ad_finite_edge_test_degenerated_alg1_C2(const RT &x1, const RT &y1,
					const RT &w1,
					const RT &x2, const RT &y2,
					const RT &w2,
					const RT &qx, const RT &qy,
					const RT &qw, bool b)
{
  must_be_filtered(x1);

  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point  q(Point(qx, qy), qw);

  Finite_edge_test_degenerated<R> f;
  return f(p1, p2, q, b, Algebraic1_tag() );
}

template < class RT >
bool
ad_finite_edge_test_degenerated_alg2_C2(const RT &x1, const RT &y1,
					const RT &w1,
					const RT &x2, const RT &y2,
					const RT &w2,
					const RT &qx, const RT &qy,
					const RT &qw, bool b)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point  q(Point(qx, qy), qw);

  Finite_edge_test_degenerated<R> f;
  return f(p1, p2, q, b, Algebraic2_tag() );
}



template < class RT >
bool
ad_finite_edge_test_degenerated_naive_C2(const RT &x1, const RT &y1,
					 const RT &w1,
					 const RT &x2, const RT &y2,
					 const RT &w2,
					 const RT &x3, const RT &y3,
					 const RT &w3,
					 const RT &qx, const RT &qy,
					 const RT &qw, bool b)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point  q(Point(qx, qy), qw);

  Finite_edge_test_degenerated<R> f;
  return f(p1, p2, p3, q, b, Naive_tag() );
}


template < class RT >
bool
ad_finite_edge_test_degenerated_alg1_C2(const RT &x1, const RT &y1,
					const RT &w1,
					const RT &x2, const RT &y2,
					const RT &w2,
					const RT &x3, const RT &y3,
					const RT &w3,
					const RT &qx, const RT &qy,
					const RT &qw, bool b)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point  q(Point(qx, qy), qw);

  Finite_edge_test_degenerated<R> f;
  return f(p1, p2, p3, q, b, Algebraic1_tag() );
}

template < class RT >
bool
ad_finite_edge_test_degenerated_alg2_C2(const RT &x1, const RT &y1,
					const RT &w1,
					const RT &x2, const RT &y2,
					const RT &w2,
					const RT &x3, const RT &y3,
					const RT &w3,
					const RT &qx, const RT &qy,
					const RT &qw, bool b)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point  q(Point(qx, qy), qw);

  Finite_edge_test_degenerated<R> f;
  return f(p1, p2, p3, q, b, Algebraic2_tag() );
}


//--------------------------------------------------------------------

template < class RT >
bool
ad_infinite_edge_test_naive_C2(const RT &x2, const RT &y2,
			       const RT &w2,
			       const RT &x3, const RT &y3,
			       const RT &w3,
			       const RT &x4, const RT &y4,
			       const RT &w4,
			       const RT &qx, const RT &qy,
			       const RT &qw, bool b)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point p4(Point(x4, y4), w4);
  Weighted_point  q(Point(qx, qy), qw);

  Infinite_edge_test<R> f;
  return f(p2, p3, p4, q, b, Naive_tag() );
}

template < class RT >
bool
ad_infinite_edge_test_alg1_C2(const RT &x2, const RT &y2,
			      const RT &w2,
			      const RT &x3, const RT &y3,
			      const RT &w3,
			      const RT &x4, const RT &y4,
			      const RT &w4,
			      const RT &qx, const RT &qy,
			      const RT &qw, bool b)
{
  must_be_filtered(x2);

  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point p4(Point(x4, y4), w4);
  Weighted_point  q(Point(qx, qy), qw);

  Infinite_edge_test<R> f;
  return f(p2, p3, p4, q, b, Algebraic1_tag() );
}

template < class RT >
bool
ad_infinite_edge_test_alg2_C2(const RT &x2, const RT &y2,
			      const RT &w2,
			      const RT &x3, const RT &y3,
			      const RT &w3,
			      const RT &x4, const RT &y4,
			      const RT &w4,
			      const RT &qx, const RT &qy,
			      const RT &qw, bool b)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point p4(Point(x4, y4), w4);
  Weighted_point  q(Point(qx, qy), qw);

  Infinite_edge_test<R> f;
  return f(p2, p3, p4, q, b, Algebraic2_tag() );
}



//--------------------------------------------------------------------

template < class RT >
bool
ad_is_degenerate_edge_test_naive_C2(const RT &x1, const RT &y1,
				    const RT &w1,
				    const RT &x2, const RT &y2,
				    const RT &w2,
				    const RT &x3, const RT &y3,
				    const RT &w3,
				    const RT &x4, const RT &y4,
				    const RT &w4)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point p4(Point(x4, y4), w4);

  Is_degenerate_edge_test<R> f;
  return f(p1, p2, p3, p4, Naive_tag() );
}

template < class RT >
bool
ad_is_degenerate_edge_test_alg1_C2(const RT &x1, const RT &y1,
				   const RT &w1,
				   const RT &x2, const RT &y2,
				   const RT &w2,
				   const RT &x3, const RT &y3,
				   const RT &w3,
				   const RT &x4, const RT &y4,
				   const RT &w4)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point p4(Point(x4, y4), w4);

  Is_degenerate_edge_test<R> f;
  return f(p1, p2, p3, p4, Algebraic1_tag() );
}


template < class RT >
bool
ad_is_degenerate_edge_test_alg2_C2(const RT &x1, const RT &y1,
				   const RT &w1,
				   const RT &x2, const RT &y2,
				   const RT &w2,
				   const RT &x3, const RT &y3,
				   const RT &w3,
				   const RT &x4, const RT &y4,
				   const RT &w4)
{
  typedef Simple_cartesian<RT>             R;
  typedef Point_2<R>                       Point;
  typedef RT                               Weight;
  typedef Weighted_point<Point, Weight>    Weighted_point;

  Weighted_point p1(Point(x1, y1), w1);
  Weighted_point p2(Point(x2, y2), w2);
  Weighted_point p3(Point(x3, y3), w3);
  Weighted_point p4(Point(x4, y4), w4);

  Is_degenerate_edge_test<R> f;
  return f(p1, p2, p3, p4, Algebraic2_tag() );
}



//--------------------------------------------------------------------

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#ifndef CGAL_ARITHMETIC_FILTER_APOLLONIUS_GRAPH_FTC2_H
#include <CGAL/Arithmetic_filter/predicates/Apollonius_graph_ftC2.h>
#endif // CGAL_ARITHMETIC_FILTER_APOLLONIUS_GRAPH_FTC2_H
#endif

#endif // CGAL_APOLLONIUS_GRAPH_FTC2_H
