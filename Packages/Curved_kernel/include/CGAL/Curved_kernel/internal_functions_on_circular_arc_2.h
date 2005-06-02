// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Curved_kernel/predicates_on_circular_arc_2.h

#ifndef CGAL_CURVED_KERNEL_PREDICATES_ON_CIRCULAR_ARC_2_H
#define CGAL_CURVED_KERNEL_PREDICATES_ON_CIRCULAR_ARC_2_H

//#include <CGAL/Curved_kernel/enum.h>

namespace CGAL {
namespace CircularFunctors {
  
  template < class CK >
  inline
  Comparison_result 
  compare_x(const typename CK::Circular_arc_endpoint_2 &p0,
            const typename CK::Circular_arc_endpoint_2 &p1)
  {
    return CGAL::compare(p0.x(), p1.x());
  }

  template < class CK >
  inline
  Comparison_result 
  compare_y(const typename CK::Circular_arc_endpoint_2 &p0,
            const typename CK::Circular_arc_endpoint_2 &p1)
  {
    return CGAL::compare(p0.y(), p1.y());
  }

  template < class CK >
  Comparison_result 
  compare_xy(const typename CK::Circular_arc_endpoint_2 &p0,
             const typename CK::Circular_arc_endpoint_2 &p1)
  {
    Comparison_result compx = compare_x<CK>(p0, p1);
    if (compx != 0)
      return compx;
    return compare_y<CK>(p0, p1);
  }

  template < class CK >
  bool
  point_in_range(const typename CK::Circular_arc_2 &A,
                 const typename CK::Circular_arc_endpoint_2 &p) 
  {
    assert (A.is_x_monotone());
    // range includes endpoints here
    return compare_x<CK>(p, A.begin()) !=
           compare_x<CK>(p, A.end());
  }

  template < class CK >
  Comparison_result
  compare_y_at_x(const typename CK::Circular_arc_endpoint_2 &p,
                 const typename CK::Circular_arc_2 &A1)
  {
    assert (A1.is_x_monotone());
    assert (point_in_range<CK>(A1, p));

    // Compare the ordinate of p with the ordinate of the center.
    Comparison_result sgn =
                  CGAL::compare(p.y(), A1.supporting_circle().center().y());
    // Is the arc on the lower or upper part of the circle ?
    // I.e. it's the comparison of the "ordinate" of the arc with the center.
    Comparison_result cmp = A1.on_upper_part() ? LARGER : SMALLER;
    if (sgn == opposite(cmp))
      return sgn;

    // If not, then we can compute if p is inside the circle or not.
    typedef typename CK::Root_of_2 Root;
    Root dx_sqr = CGAL::square(p.x() - A1.supporting_circle().center().x());
    Root dy_sqr = CGAL::square(p.y() - A1.supporting_circle().center().y());
                  // NB : that one can be factorized with the above...

    // Now we want the comparison of dx_sqr + dy_sqr with squared_radius.
    // It's the same as dx_sqr - squared_radius with -dy_sqr.

    Comparison_result distance_to_center = 
          CGAL::compare(dx_sqr,
                           A1.supporting_circle().squared_radius() - dy_sqr);

    if (cmp > 0)
      return distance_to_center;
    else
      return opposite(distance_to_center);
  }

  template < class CK >
  Comparison_result 
  compare_y_to_right(const typename CK::Circular_arc_2 &A1,
		     const typename CK::Circular_arc_2 &A2, 
		     const typename CK::Circular_arc_endpoint_2 &p)
  {
    // FIXME : add preconditions to check that the 2 arcs are defined at
    // the right of the intersection.
    assert (A1.is_x_monotone());
    assert (A2.is_x_monotone());

    const typename CK::Circle_2 & C1 = A1.supporting_circle();
    const typename CK::Circle_2 & C2 = A2.supporting_circle();

    if (C1 == C2) {
      // The point is either a left vertical tangent point of both,
      // or a normal point (-> EQUAL).
      bool b1 = A1.on_upper_part();
      bool b2 = A2.on_upper_part();
      if (b1 == b2)
        return EQUAL;
      if (b1 == true && b2 == false)
        return LARGER;
      assert (b1 == false && b2 == true);
      return SMALLER;
    }

    typename CK::Root_of_2 b1_y = C1.center().y() - p.y();
    typename CK::Root_of_2 b2_y = C2.center().y() - p.y();

    int s_b1_y = CGAL::sign(b1_y);
    int s_b2_y = CGAL::sign(b2_y);

    if (s_b1_y == 0) {
      // Vertical tangent for A1.
      if (s_b2_y != 0)
        return A1.on_upper_part() ? LARGER : SMALLER;
      // Vertical tangent for A2 also.
      bool b1 = A1.on_upper_part();
      bool b2 = A2.on_upper_part();
      if (b1 == b2)
        return b1 ? compare_x(C1.center(), C2.center())
                  : compare_x(C2.center(), C1.center());
      if (b1 == true && b2 == false)
        return LARGER;
      assert (b1 == false && b2 == true);
      return SMALLER;
    }
    if (s_b2_y == 0) {
      // Vertical tangent for A2.
      return A2.on_upper_part() ? SMALLER : LARGER;
    }

    // No more vertical tangent points.
    assert(s_b1_y != 0);
    assert(s_b2_y != 0);

    typename CK::Root_of_2 b1_x = p.x() - C1.center().x();
    typename CK::Root_of_2 b2_x = p.x() - C2.center().x();

    int s_b1_x = CGAL::sign(b1_x);
    int s_b2_x = CGAL::sign(b2_x);

    // We compute the slope of the 2 tangents, then we compare them.

    Comparison_result cmp = CGAL::compare(s_b1_y * s_b1_x,
                                          s_b2_y * s_b2_x);
    // The slopes have different signs.
    if (cmp != 0)
      return cmp;

    // The slopes have the same signs : we have to square.
    if (CGAL::square(squared_distance(C1.center(), C2.center())
                           - C1.squared_radius() - C2.squared_radius())
                       < 4 * C1.squared_radius() * C2.squared_radius() ) {
      // The two circles are not tangent.
      return static_cast<Comparison_result>(
             CGAL::compare(C1.squared_radius() * CGAL::square(b2_y),
                           C2.squared_radius() * CGAL::square(b1_y))
             * s_b1_y * s_b1_x );
    }

    // tangent circles
    if (s_b1_x * s_b2_x < 0)
      // Circles are on both sides, and the tangent is not horizontal
      return compare_y(C1.center(), C2.center());

    if (s_b1_x * s_b2_x > 0)
      // Circles are on the same side, and the tgt is not horizontal.
      return compare_y(C2.center(), C1.center());

    // The tangent is horizontal.
    assert(s_b1_x == 0 && s_b2_x == 0);
    if (s_b1_y == s_b2_y)
      // The 2 circles are both below or both above the tangent
      return compare_y(C2.center(), C1.center());
    return compare_y(C1.center(), C2.center());
  }

  template < class CK >
  inline
  bool
  equal(const typename CK::Circular_arc_endpoint_2 &p0,
        const typename CK::Circular_arc_endpoint_2 &p1)
  {
#if 0
    CGAL_PROFILER(x1, "point_equal calls");
    // Let's try to see if the 2 points come from the same arc.
    if (p0.circle(0) == p1.circle(0)) {
      CGAL_PROFILER(x2, "point_equal same circle(0)");
      if (p0.circle(1) == p1.circle(1)) {
        CGAL_PROFILER(x2, "point_equal same circle(1) as well");
        if (p0.is_left() == p1.is_left()) {
          CGAL_PROFILER(x2, "point_equal same is_left as well");
        }
      }
    }
    if (p0.circle(0) == p1.circle(1) ||
        p0.circle(1) == p1.circle(0)) {
      CGAL_PROFILER(x2, "point_equal same crossed circles");
    }
    if (p0.circle(0) == p1.circle(0) ||
        p0.circle(1) == p1.circle(1) ||
        p0.circle(1) == p1.circle(0) ||
        p0.circle(0) == p1.circle(1)) {
      CGAL_PROFILER(x2, "point_equal share one common circle");
    }
    if (compare_xy<CK>(p0, p1) == 0) {
      CGAL_PROFILER(x2, "point_equal returns equal");
    }
#endif

    // Filter out some of the stupid cases.
    if (p0.is_left() == p1.is_left() &&
	p0.circle(0) == p1.circle(0) &&
        p0.circle(1) == p1.circle(1) )
      return true;

    // More stupid cases can be caught (or fixed at the source).

    return compare_xy<CK>(p0, p1) == 0;
  }

  template < class CK >
  bool
  equal(const typename CK::Circular_arc_2 &A1,
        const typename CK::Circular_arc_2 &A2)
  {
    assert (A1.is_x_monotone());
    assert (A2.is_x_monotone());

    if ( A1.supporting_circle() != A2.supporting_circle() )
      return false;

    return equal<CK>( A1.begin(), A2.begin() ) &&
           equal<CK>( A1.end(), A2.end() );
  }

  template < class CK >
  bool
  do_overlap(const typename CK::Circular_arc_2 &A1,
	     const typename CK::Circular_arc_2 &A2)
  {
    assert (A1.is_x_monotone());
    assert (A2.is_x_monotone());

    if ( A1.supporting_circle() != A2.supporting_circle() ) return false;
    if ( A1.on_upper_part() != A2.on_upper_part() ) return false;

    return compare_x<CK>(A1.right(), A2.left()) > 0
        && compare_x<CK>(A1.left(), A2.right()) < 0;
  }

  // Small accessory function
  // Tests whether a given point is on an arc, with the precondition that
  // it's (symbolically) on the supporting circle.
  template < class CK >
  bool
  is_on_arc(const typename CK::Circular_arc_2 &a,
	    const typename CK::Circular_arc_endpoint_2 &p)
  {
    assert(a.is_x_monotone());
    assert(a.supporting_circle() == p.circle(0) ||
           a.supporting_circle() == p.circle(1) );

    if (! point_in_range<CK>(a, p) )
      return false;

    int cmp = CGAL::compare(p.y(), a.supporting_circle().center().y());

    return  cmp == 0 || (cmp > 0 &&  a.on_upper_part())
                     || (cmp < 0 && !a.on_upper_part());
  }

  template < class CK >
  void
  split(const typename CK::Circular_arc_2 &A,
	const typename CK::Circular_arc_endpoint_2 &p,
	typename CK::Circular_arc_2 &ca1,
	typename CK::Circular_arc_2 &ca2)
  {
    assert( A.is_x_monotone() );
    assert( point_in_range<CK>( A, p ) );
    assert( A.on_upper_part() == (p.y() >
				  A.supporting_circle().center().y()));
    assert( is_on_arc<CK>(A,p) );

    typedef typename CK::Circular_arc_2  Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2      Circular_arc_endpoint_2;

    if ( p.circle(0) == A.supporting_circle() )
      {
	assert( !(p.circle(1) == A.supporting_circle()) );
	ca1 = Circular_arc_2( A, true,  p.circle(1), p.is_left() );
	ca2 = Circular_arc_2( A, false, p.circle(1), p.is_left() );
	return;
      };

    if ( p.circle(1) == A.supporting_circle() )
      {
	assert( !(p.circle(0) == A.supporting_circle()) );
	ca1 = Circular_arc_2( A, true,  p.circle(0), p.is_left() );
	ca2 = Circular_arc_2( A, false, p.circle(0), p.is_left() );
	return;
      };

    assert( !(p.circle(0) == A.supporting_circle()) );
    assert( !(p.circle(1) == A.supporting_circle()) );
 
    // p is defined by another pair of circles.
    // Then we must determine which intersection it is of the supporting
    // circle of A, with one of p's circles (whichever is ok).
    bool b = (p == Circular_arc_endpoint_2 (A.supporting_circle(), p.circle(0), true));
    assert( b ||
             (p == Circular_arc_endpoint_2 (A.supporting_circle(), p.circle(0), false)));

    ca1 = Circular_arc_2( A, true,  p.circle(0), b );
    ca2 = Circular_arc_2( A, false, p.circle(0), b );
    return;
  }

  template< class CK, class OutputIterator>
  OutputIterator
  construct_intersections_2( const typename CK::Circular_arc_2 &a1,
			     const typename CK::Circular_arc_2 &a2,
			     OutputIterator res )
  {
    typedef typename CK::Circular_arc_endpoint_2  Circular_arc_endpoint_2;
    typedef typename CK::Circular_arc_2           Circular_arc_2;

    assert(a1.is_x_monotone());
    assert(a2.is_x_monotone());
    // todo : implement for general arcs

    // Overlapping curves.
    if (a1.supporting_circle() == a2.supporting_circle()) {
      // The ranges need to overlap in order for the curves to overlap.
      if (compare_x<CK>(a1.left(), a2.right()) > 0 ||
          compare_x<CK>(a2.left(), a1.right()) > 0)
        return res;

      // They both need to be on the same upper/lower part.
      if (a1.on_upper_part() != a2.on_upper_part()) {
        // But they could share the right vertical tangent point.
        if (a1.right() == a2.right()) {
            *res++ = make_object
	      (std::make_pair(a1.right(),static_cast<unsigned>(0)));
        }
        // Or they could share the left vertical tangent point.
        if (a1.left() == a2.left()) {
            *res++ = make_object
	      (std::make_pair(a1.left(),static_cast<unsigned>(0)));
        }
        return res;
      };

      // We know they overlap, determine the extremities of the common subcurve
      // TODO : We should use std::max and std::min, but they require less_x_2.
      const Circular_arc_2 & arctmp = 
	compare_x<CK>(a1.right(), a2.right()) < 0 ? a1 : a2;
      // we know that the right endpoint is correct, let us look for
      // the left now:


      if ( compare_x<CK>(a1.left(), a2.left()) > 0 ) //? a1.left() : a2.left();
	{ //the left endpoint is a1's
	  if (a1.left().circle(0) == arctmp.supporting_circle()) {
	    const Circular_arc_2 & arc =
	      Circular_arc_2(arctmp.supporting_circle(), false,
			     a1.left().circle(1), a1.left().is_left());
	    assert(arc.left()==a1.left());
	    assert(arc.right()==arctmp.right());
	    assert(arc.is_x_monotone());
	    *res++ = make_object(arc);
	  }
	  else {
	    assert(a1.left().circle(1) == arctmp.supporting_circle());
	    const Circular_arc_2 & arc =
	      Circular_arc_2(arctmp.supporting_circle(), false,
			     a1.left().circle(0), a1.left().is_left());
	    assert(arc.left()==a1.left());
	    assert(arc.right()==arctmp.right());
	    assert(arc.is_x_monotone());
	    *res++ = make_object(arc);
	  };
	}
      else { //the left endpoint is a2's
	if (a2.left().circle(0) == arctmp.supporting_circle()) {
	    const Circular_arc_2 & arc =
	      Circular_arc_2(arctmp.supporting_circle(), false,
			     a2.left().circle(1), a2.left().is_left());
	    assert(arc.left()==a2.left());
	    assert(arc.right()==arctmp.right());
	    assert(arc.is_x_monotone());
	    *res++ = make_object(arc);
	  }
	  else {
	    assert(a2.left().circle(1) == arctmp.supporting_circle());
	    const Circular_arc_2 & arc =
	      Circular_arc_2(arctmp.supporting_circle(), false,
			     a2.left().circle(0), a2.left().is_left());
	    assert(arc.left()==a2.left());
	    assert(arc.right()==arctmp.right());
	    assert(arc.is_x_monotone());
	    *res++ = make_object(arc);
	  };
      };

      return res;
    }

    // SHOULD USE INTERSECTIONS ON CIRCLES INSTEAD
    // OR AT LEAST SOLVE...

    // We need to check that the supporting circles
    // do intersect before going further.
    if (! do_intersect(a1.supporting_circle(), a2.supporting_circle())) 
      { return res; }

    // Get the two intersection points of the supporting circles.
    Circular_arc_endpoint_2 
      left (a1.supporting_circle(), a2.supporting_circle(), true);
    Circular_arc_endpoint_2 
      right(a1.supporting_circle(), a2.supporting_circle(), false);

    if ( left != right ) // multiplicity 1
      {
	// We also need to check that these intersection points are on the arc.
	if (is_on_arc<CK>(a1, left) && is_on_arc<CK>(a2, left)) {
	  *res++ = make_object(std::make_pair(left,static_cast<unsigned>(1)));
	}
	if (is_on_arc<CK>(a1, right) && is_on_arc<CK>(a2, right)) {
	  *res++ = make_object(std::make_pair(right,static_cast<unsigned>(1)));
	}
      }
    else // multiplicity 2
      {
	if (is_on_arc<CK>(a1, left) && is_on_arc<CK>(a2, left)) 
	  { *res++ = make_object
	      (std::make_pair(left,static_cast<unsigned>(2))); }
      }
    return res;
  }


  template < class CK >
  bool
  nearest_intersection_to_right(const typename CK::Circular_arc_2 &a1,
                                const typename CK::Circular_arc_2 &a2,
                                const typename CK::Circular_arc_endpoint_2 &pt,
                                      typename CK::Circular_arc_endpoint_2 &p1,
                                      typename CK::Circular_arc_endpoint_2 &p2)
  {
    typedef typename CK::Circular_arc_endpoint_2  Circular_arc_endpoint_2;

    assert(a1.is_x_monotone());
    assert(a2.is_x_monotone());

    // Overlapping curves.
    if (a1.supporting_circle() == a2.supporting_circle()) {
      // The ranges need to overlap in order for the curves to overlap.
      if (compare_x<CK>(a1.left(), a2.right()) > 0 ||
          compare_x<CK>(a2.left(), a1.right()) > 0)
        return false;

      // They both need to be on the same upper/lower part.
      if (a1.on_upper_part() != a2.on_upper_part()) {
        // But they could share the right vertical tangent point.
        if (a1.right() == a2.right() &&
            compare_x<CK>(pt, a1.right()) < 0) {
            p1 = p2 = a1.right();
            return true;
        }
        // Or they could share the left vertical tangent point.
        if (a1.left() == a2.left() &&
            compare_x<CK>(pt, a1.left()) < 0) {
            p1 = p2 = a1.left();
            return true;
        }
        return false;
      }

      // We know they overlap, determine the extremities of the common subcurve
      // TODO : We should use std::max and std::min, but they require less_x_2.
      const Circular_arc_endpoint_2 & tmp2 = compare_x<CK>(a1.right(), a2.right())
                           < 0 ? a1.right() : a2.right();

      // Now we need to compare that with pt.
      if (compare_x<CK>(pt, tmp2) != SMALLER)
        return false;
      p2 = tmp2;
      const Circular_arc_endpoint_2 & tmp1 =
	  compare_x<CK>(a1.left(), a2.left()) > 0 ? a1.left() : a2.left();
      if (compare_x<CK>(pt, tmp1) != SMALLER) {
        p1 = pt;
        return true;
      }
      p1 = tmp1;
      return true;
    }

    // We need to check that the supporting circles
    // do intersect before going further.
    if (! do_intersect(a1.supporting_circle(),
                       a2.supporting_circle())) {
      return false;
    }

    // Get the two intersection points of the supporting circles.
    Circular_arc_endpoint_2 left (a1.supporting_circle(), a2.supporting_circle(), true);
    Circular_arc_endpoint_2 right(a1.supporting_circle(), a2.supporting_circle(), false);

    // We also need to check that these intersection points are on the arc.
    if (is_on_arc<CK>(a1, left) &&
        is_on_arc<CK>(a2, left) &&
        compare_xy<CK>(left, pt) > 0) {
      p1 = p2 = left;
      return true;
    }
    if (is_on_arc<CK>(a1, right) &&
        is_on_arc<CK>(a2, right) &&
        compare_xy<CK>(right, pt) > 0) {
      p1 = p2 = right;
      return true;
    }
    // no intersection.
    return false;
  }

  template < class CK, class OutputIterator >
  OutputIterator
  make_x_monotone( const typename CK::Circular_arc_2 &A,
		   OutputIterator res )
  {
    typedef typename CK::Circular_arc_2           Circular_arc_2;
    typedef typename CK::Circle_2                 Circle_2;
    typedef typename CK::FT                       FT;
    typedef typename CK::Linear_kernel::Point_2   Point_2;

    int cmp_begin = CGAL::compare(A.begin().y(), A.center().y());
    int cmp_end   = CGAL::compare(A.end().y(),   A.center().y());

    int cmp_x = compare_x(A.begin(), A.end());

    // We don't need to split
    if (cmp_begin != opposite(cmp_end) &&
        (((cmp_begin > 0 || cmp_end > 0) && cmp_x > 0) ||
          (cmp_begin < 0 || cmp_end < 0) && cmp_x < 0) ) {
      *res++ = A;
      return res; 
    }

    // Half circles
    if (cmp_begin == 0 && cmp_end == 0 && cmp_x != 0) {
      *res++ = A;
      return res; 
    }

    // We need to split
    assert(!A.is_x_monotone());

    // Define a circle intersecting the supporting circle of A
    // in the 2 vertical tangent points.
    Circle_2 c (Point_2(A.center().x(), A.center().y()-1),
                A.squared_radius()+1);

    if (cmp_begin > 0) {
      *res++ = Circular_arc_2 (A.supporting_circle(),
                               A.begin().circle(1), A.begin().is_left(),
                               c, true);
      if (cmp_end > 0) {
        // We must cut in 3 parts.
        *res++ = Circular_arc_2 (A.supporting_circle(),
                                 c, true, c, false);
        *res++ = Circular_arc_2 (A.supporting_circle(),
                                 c, false,
                                 A.end().circle(1), A.end().is_left());
      } else {
        *res++ = Circular_arc_2 (A.supporting_circle(),
                                 c, true,
                                 A.end().circle(1), A.end().is_left());
      }
    }
    else if (cmp_begin < 0) {
      // Very similar to the previous case.
      *res++ = Circular_arc_2 (A.supporting_circle(),
                               A.begin().circle(1), A.begin().is_left(),
                               c, false);
      if (cmp_end < 0) {
        // We must cut in 3 parts.
        *res++ = Circular_arc_2 (A.supporting_circle(),
                                 c, false, c, true);
        *res++ = Circular_arc_2 (A.supporting_circle(),
                                 c, true,
                                 A.end().circle(1), A.end().is_left());
      } else {
        *res++ = Circular_arc_2 (A.supporting_circle(),
                                 c, false,
                                 A.end().circle(1), A.end().is_left());
      }
    }
    else { // cmp_begin == 0
      if (CGAL::compare(A.begin().x(), A.center().x()) < 0) {
        assert (cmp_end >= 0);
        *res++ = Circular_arc_2 (A.supporting_circle(),
                                 A.begin().circle(1), A.begin().is_left(),
                                 c, false);
        *res++ = Circular_arc_2 (A.supporting_circle(),
                                 c, false,
                                 A.end().circle(1), A.end().is_left());
      }
      else {
        assert (CGAL::compare(A.begin().x(), A.center().x()) > 0);
        assert (cmp_end != LARGER);
        *res++ = Circular_arc_2 (A.supporting_circle(),
                                 A.begin().circle(1), A.begin().is_left(),
                                 c, true);
        *res++ = Circular_arc_2 (A.supporting_circle(),
                                 c, true,
                                 A.end().circle(1), A.end().is_left());
      }
    }

    return res;
  }

} // namespace CircularFunctors 
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_PREDICATES_ON_CIRCULAR_ARC_2_H
