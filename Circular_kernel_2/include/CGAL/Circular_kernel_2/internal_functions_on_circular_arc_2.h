// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_PREDICATES_ON_CIRCULAR_ARC_2_H
#define CGAL_CIRCULAR_KERNEL_PREDICATES_ON_CIRCULAR_ARC_2_H

#include <CGAL/license/Circular_kernel_2.h>


#include <CGAL/Circular_kernel_2/internal_functions_on_circle_2.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Circular_kernel_2/Circular_arc_2.h>
#include <CGAL/int.h>
namespace CGAL {
namespace CircularFunctors {

  template < class CK >
  inline
  Comparison_result
  compare_x(const typename CK::Circular_arc_point_2 &p0,
            const typename CK::Circular_arc_point_2 &p1)
  {
    typedef typename CK::Algebraic_kernel   AK;
    if(p0.equal_ref(p1)) return static_cast<Comparison_result>(0);
    return AK().compare_x_object()(p0.coordinates(), p1.coordinates());
  }

  template < class CK >
  inline
  Comparison_result
  compare_y(const typename CK::Circular_arc_point_2 &p0,
            const typename CK::Circular_arc_point_2 &p1)
  {
    typedef typename CK::Algebraic_kernel   AK;
    if(p0.equal_ref(p1)) return static_cast<Comparison_result>(0);
    return AK().compare_y_object()(p0.coordinates(), p1.coordinates());
  }

  template < class CK >
  Comparison_result
  compare_xy(const typename CK::Circular_arc_point_2 &p0,
             const typename CK::Circular_arc_point_2 &p1)
  {
    if(p0.equal_ref(p1)){
      return EQUAL;
    }
    typedef typename CK::Algebraic_kernel   AK;
    return AK().compare_xy_object()(p0.coordinates(), p1.coordinates());
  }

  // PRE CONDITION:
  // The coordinates of P, Q, R have to have the same
  // delta or (beta == 0 || delta == 0)
  // We cannot code this pre condition because
  // if Root_of_2 is interval_nt "beta", "delta" mean nothing
  template < class CK >
  Orientation
  orientation(const typename CK::Circular_arc_point_2 &p,
              const typename CK::Circular_arc_point_2 &q,
              const typename CK::Circular_arc_point_2 &r)
  {
    typedef typename CK::Root_of_2 Root_of_2;
    const Root_of_2 px = p.x();
    const Root_of_2 py = p.y();
    const Root_of_2 qx = q.x();
    const Root_of_2 qy = q.y();
    const Root_of_2 rx = r.x();
    const Root_of_2 ry = r.y();
    const Root_of_2 a00 = qx-px;
    const Root_of_2 a01 = qy-py;
    const Root_of_2 a10 = rx-px;
    const Root_of_2 a11 = ry-py;
    return CGAL_NTS compare(a00*a11, a10*a01);
  }

//   template < class CK >
//   inline
//   Comparison_result
//   compare_x(const typename CK::Circular_arc_point_2 &p0,
//             const typename CK::Point_2 &p1)
//   {
//     return CGAL::compare(p0.x(), p1.x());
//   }

//   template < class CK >
//   inline
//   Comparison_result
//   compare_x(const typename CK::Point_2 &p0,
//             const typename CK::Circular_arc_point_2 &p1)
//   {
//     return CGAL::compare(p0.x(), p1.x());
//   }


//   template < class CK >
//   inline
//   Comparison_result
//   compare_y(const typename CK::Circular_arc_point_2 &p0,
//             const typename CK::Point_2 &p1)
//   {
//     return CGAL::compare(p0.y(), p1.y());
//   }

//   template < class CK >
//   inline
//   Comparison_result
//   compare_y(const typename CK::Point_2 &p0,
//             const typename CK::Circular_arc_point_2 &p1)
//   {
//     return CGAL::compare(p0.y(), p1.y());
//   }


//   template < class CK >
//   Comparison_result
//   compare_xy(const typename CK::Circular_arc_point_2 &p0,
//              const typename CK::Point_2 &p1)
//   {
//     Comparison_result compx = compare_x<CK>(p0, p1);
//     if (compx != 0)
//       return compx;
//     return compare_y<CK>(p0, p1);
//   }

//   template < class CK >
//   Comparison_result
//   compare_xy(const typename CK::Point_2 &p0,
//              const typename CK::Circular_arc_point_2 &p1)
//   {
//     Comparison_result compx = compare_x<CK>(p0, p1);
//     if (compx != 0)
//       return compx;
//     return compare_y<CK>(p0, p1);
//   }


  template < class CK >
  bool
  point_in_x_range(const typename CK::Circular_arc_point_2 &source,
                   const typename CK::Circular_arc_point_2 &target,
                   const typename CK::Circular_arc_point_2 &p)
  {
    // range includes endpoints here
    return ( (CircularFunctors::compare_x<CK>(p, source) != CircularFunctors::compare_x<CK>(p, target))
             || (CircularFunctors::compare_x<CK>(p, source) == CGAL::EQUAL) );
  }

  template < class CK >
  bool
  point_in_x_range(const typename CK::Circular_arc_2 &A,
                   const typename CK::Circular_arc_point_2 &p)
  {
    //CGAL_kernel_precondition (A.is_x_monotone());
    // range includes endpoints here
    return CircularFunctors::compare_x<CK>( p, A.source()) != CircularFunctors::compare_x<CK>(p, A.target() );
  }

  template < class CK >
  Comparison_result
  compare_y_at_x(const typename CK::Circular_arc_point_2 &p,
                 const typename CK::Circular_arc_2 &A1)
  {
    //CGAL_kernel_precondition (A1.is_x_monotone());
    //CGAL_kernel_precondition (CircularFunctors::point_in_x_range<CK>(A1, p));

    if((p.equal_ref(A1.source())) || (p.equal_ref(A1.target()))){
      return EQUAL;
    }

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
      CGAL::compare(dx_sqr, A1.supporting_circle().squared_radius() - dy_sqr);

    if (cmp > 0)
      return distance_to_center;
    else
      return opposite(distance_to_center);
  }



  template < class CK >
  Comparison_result
  compare_y_to_right(const typename CK::Circular_arc_2 &A1,
                     const typename CK::Circular_arc_2 &A2,
                     const typename CK::Circular_arc_point_2 &p)
  {
    // FIXME : add preconditions to check that the 2 arcs are defined at
    // the right of the intersection.
    //CGAL_kernel_precondition (A1.is_x_monotone());
    //CGAL_kernel_precondition (A2.is_x_monotone());

#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
    // intersection found on the map
    typedef std::vector<CGAL::Object> solutions_container;
    typedef typename CK::Circular_arc_2 Circular_arc_2;

    solutions_container early_sols;
    if(Circular_arc_2::template find_intersection< solutions_container >
      (A1,A2,early_sols)) {
      if(A1.on_upper_part()) return LARGER;
      return SMALLER;
    }
#endif

    const typename CK::Circle_2 & C1 = A1.supporting_circle();
    const typename CK::Circle_2 & C2 = A2.supporting_circle();

    if (CircularFunctors::non_oriented_equal<CK>(C1,C2)) {
      // The point is either a left vertical tangent point of both,
      // or a normal point (-> EQUAL).
      bool b1 = A1.on_upper_part();
      bool b2 = A2.on_upper_part();
      if (b1 == b2)
        return EQUAL;
      if (b1 == true && b2 == false)
        return LARGER;
      CGAL_kernel_assertion (b1 == false && b2 == true);
      return SMALLER;
    }

    const typename CK::Root_of_2 b1_y = C1.center().y() - p.y();
    const typename CK::Root_of_2 b2_y = C2.center().y() - p.y();

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
      CGAL_kernel_assertion (b1 == false && b2 == true);
      return SMALLER;
    }

    if (s_b2_y == 0) {
      // Vertical tangent for A2.
      return A2.on_upper_part() ? SMALLER : LARGER;
    }

    // No more vertical tangent points.
    CGAL_kernel_assertion(s_b1_y != 0);
    CGAL_kernel_assertion(s_b2_y != 0);

    int s_b1_x = (int) CGAL::compare(p.x(), C1.center().x());
    int s_b2_x = (int) CGAL::compare(p.x(), C2.center().x());

    // We compute the slope of the 2 tangents, then we compare them.
    Comparison_result cmp = CGAL::compare(s_b1_y * s_b1_x,
                                          s_b2_y * s_b2_x);
    // The slopes have different signs.
    if (cmp != 0)
      return cmp;

    // The slopes have the same signs : we have to square.
    if (CGAL::square(squared_distance(C1.center(), C2.center())
                     - C1.squared_radius() - C2.squared_radius())
        < 4 * C1.squared_radius() * C2.squared_radius() )
      {
        // The two circles are not tangent.
        return static_cast<Comparison_result>
          (CGAL::compare(C1.squared_radius() * CGAL::square(b2_y),
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
    CGAL_kernel_assertion(s_b1_x == 0 && s_b2_x == 0);
    if (s_b1_y == s_b2_y)
      // The 2 circles are both below or both above the tangent
      return compare_y(C2.center(), C1.center());

    return compare_y(C1.center(), C2.center());
  }

  template < class CK >
  inline
  bool
  equal(const typename CK::Circular_arc_point_2 &p0,
        const typename CK::Circular_arc_point_2 &p1)
  {
    if(p0.equal_ref(p1)) return static_cast<Comparison_result>(1);
    return CircularFunctors::compare_xy<CK>(p0, p1) == 0;
  }

  template < class CK >
  bool
  equal(const typename CK::Circular_arc_2 &A1,
        const typename CK::Circular_arc_2 &A2)
  {
    /*if ((A1.supporting_circle() != A2.supporting_circle()) &&
        (A1.supporting_circle() != A2.supporting_circle().opposite()))
      return false;*/

    if(!CircularFunctors::non_oriented_equal<CK>(
      A1.supporting_circle(), A2.supporting_circle()))
      return false;

    return (CircularFunctors::equal<CK>(A1.source(), A2.source()) &&
            CircularFunctors::equal<CK>(A1.target(), A2.target()));
  }

//   template < class CK >
//   bool
//   equal(const typename CK::Circular_arc_2 &A1,
//         const typename CK::Circular_arc_2 &A2)
//   {
//     CGAL_kernel_precondition (A1.is_x_monotone());
//     CGAL_kernel_precondition (A2.is_x_monotone());

//     if ( A1.supporting_circle() != A2.supporting_circle() )
//       return false;

//     return equal<CK>( A1.source(), A2.source() )
//         && equal<CK>( A1.target(), A2.target() );
//   }

  // Small accessory function
  // Tests whether a given point is on an arc, with the precondition that
  // it's (symbolically) on the supporting circle.
  /*template < class CK >
  bool
  has_on(const typename CK::Circular_arc_2 &a,
         const typename CK::Circular_arc_point_2 &p,
         const bool has_on_supporting_circle = false)
  {
    CGAL_kernel_precondition(a.is_x_monotone());

//     typedef typename CK::Polynomial_for_circles_2_2 Polynomial_for_circles_2_2;
//     Polynomial_for_circles_2_2
//       equation = get_equation<CK>(a.supporting_circle());

//     if(CGAL::sign_at<typename CK::Algebraic_kernel>
//        (equation,p.coordinates())!= ZERO)
//       return false;
    if(!has_on_supporting_circle) {
      if ( ! CircularFunctors::has_on<CK>(a.supporting_circle(),p) )
        return false;
    }

    if (! CircularFunctors::point_in_x_range<CK>(a, p) )
      return false;

    int cmp = CGAL::compare(p.y(), a.supporting_circle().center().y());

    return  cmp == 0 || (cmp > 0 &&  a.on_upper_part())
                     || (cmp < 0 && !a.on_upper_part());
  }*/

  template < class CK >
  bool
  has_on(const typename CK::Circular_arc_2 &a,
         const typename CK::Circular_arc_point_2 &p,
         const bool has_on_supporting_circle = false)
  {

    if( (p.equal_ref(a.source())) || (p.equal_ref(a.source()))) {
      return true;
    }

    if(!has_on_supporting_circle) {
      if ( ! CircularFunctors::has_on<CK>(a.supporting_circle(),p) )
        return false;
    }

    if(a.is_full()) return true;

    if(a.is_x_monotone()) {
      int cmp_ps = CircularFunctors::compare_x<CK>(p,a.source());
      int cmp_pt = CircularFunctors::compare_x<CK>(p,a.target());
      if(cmp_ps == cmp_pt) return false;
      int cmp = CGAL::compare(p.y(), a.supporting_circle().center().y());
      return  cmp == 0 || (cmp > 0 &&  a.on_upper_part())
                       || (cmp < 0 && !a.on_upper_part());
    } else if(a.is_complementary_x_monotone())  {
      int cmp_ps = CircularFunctors::compare_x<CK>(p,a.source());
      int cmp_pt = CircularFunctors::compare_x<CK>(p,a.target());
      if(cmp_ps == cmp_pt) return true;
      if((!cmp_ps) || (!cmp_pt)) return true;
      int cmp = CGAL::compare(p.y(), a.supporting_circle().center().y());
      return  cmp == 0 || (cmp < 0 &&  a.complementary_on_upper_part())
                     || (cmp > 0 && !a.complementary_on_upper_part());
    } else if(a.is_y_monotone()) {
      int cmp_ps = CircularFunctors::compare_y<CK>(p,a.source());
      int cmp_pt = CircularFunctors::compare_y<CK>(p,a.target());
      if(cmp_ps == cmp_pt) return false;
      int cmp = CGAL::compare(p.x(), a.supporting_circle().center().x());
      return  cmp == 0 || (cmp < 0 &&  a.on_left_part())
                       || (cmp > 0 && !a.on_left_part());
    } else if(a.is_complementary_y_monotone()) {
      int cmp_ps = CircularFunctors::compare_y<CK>(p,a.source());
      int cmp_pt = CircularFunctors::compare_y<CK>(p,a.target());
      if(cmp_ps == cmp_pt) return true;
      if((!cmp_ps) || (!cmp_pt)) return true;
      int cmp = CGAL::compare(p.x(), a.supporting_circle().center().x());
      return  cmp == 0 || (cmp > 0 && a.complementary_on_left_part())
                     || (cmp < 0 && !a.complementary_on_left_part());
    } else {
      int cmp_scy = CGAL::compare(a.source().y(), a.supporting_circle().center().y());
      int cmp = CGAL::compare(p.y(), a.supporting_circle().center().y());
      if(cmp_scy < 0) {
        if(cmp > 0) return CGAL::compare(p.x(), a.target().x()) >= 0;
        else return CGAL::compare(p.x(), a.source().x()) >= 0;
      } else {
        if(cmp > 0) return CGAL::compare(p.x(), a.source().x()) <= 0;
        else return CGAL::compare(p.x(), a.target().x()) <= 0;
      }
    }
  }

  template < class CK >
  bool
  do_overlap(const typename CK::Circular_arc_2 &A1,
             const typename CK::Circular_arc_2 &A2)
  {
    //CGAL_kernel_precondition (A1.is_x_monotone());
    //CGAL_kernel_precondition (A2.is_x_monotone());

    /*if ( (A1.supporting_circle() != A2.supporting_circle()) &&
         (A1.supporting_circle() != A2.supporting_circle().opposite()) )
      return false;*/
    if(!CircularFunctors::non_oriented_equal<CK>(
      A1.supporting_circle(), A2.supporting_circle()))
      return false;

    //if ( A1.on_upper_part() != A2.on_upper_part() ) return false;

    //return CircularFunctors::compare_x<CK>(A1.right(), A2.left()) > 0
    //    && CircularFunctors::compare_x<CK>(A1.left(), A2.right()) < 0;
    if(A1.is_full()) return true;
    if(A2.is_full()) return true;
    if((has_on<CK>(A1,A2.target(),true)) ||
       (has_on<CK>(A1,A2.source(),true))) return true;
    return has_on<CK>(A2,A1.source(),true);
  }


  template < class CK >
  void
  split(const typename CK::Circular_arc_2 &A,
        const typename CK::Circular_arc_point_2 &p,
        typename CK::Circular_arc_2 &ca1,
        typename CK::Circular_arc_2 &ca2)
  {
    CGAL_kernel_precondition( CircularFunctors::has_on<CK>(A, p));

    typedef typename CK::Circular_arc_2  Circular_arc_2;

    const Circular_arc_2 &rc1 =
      Circular_arc_2( A.supporting_circle(), A.source(), p);
    const Circular_arc_2 &rc2 =
      Circular_arc_2( A.supporting_circle(), p, A.target());

    if ( CircularFunctors::compare_x<CK>(rc1.source(), rc2.source()) != SMALLER) {
      ca1 = rc2;
      ca2 = rc1;
    } else {
      ca1 = rc1;
      ca2 = rc2;
    }
#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
    std::vector < CGAL::Object > res;

    if(A.is_full()) {
      res.push_back(make_object(std::make_pair(ca1.source(),1u)));
      res.push_back(make_object(std::make_pair(ca2.source(),1u)));
    } else {
      res.push_back(make_object(std::make_pair(p,1u)));
    }

    Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
      (ca1,ca2,res);
#endif

    /*ca1 = Circular_arc_2( A.supporting_circle(), A.source(), p);
    ca2 = Circular_arc_2( A.supporting_circle(), p, A.target());
    //if ( ca1.right()!=ca2.left() )
    if ( CircularFunctors::compare_x<CK>(ca1.left(), ca2.left()) != SMALLER )
      {
        //std::cout << " SWAP " << std::endl;
        std::swap(ca1,ca2);
      }*/
  }

  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Circular_arc_2 &a1,
               const typename CK::Circular_arc_2 &a2,
               OutputIterator res )
  {
    typedef typename CK2_Intersection_traits<CK, typename CK::Circular_arc_2,
                                                 typename CK::Circular_arc_2>::type result_type;

    typedef std::vector<CGAL::Object> solutions_container;
    typedef typename CK::Circular_arc_2 Circular_arc_2;


#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
    // same curve
    if(a1.number() == a2.number()) {
      *res++ = result_type(a1);
       return res;
    }

    // intersection found on the map
    solutions_container early_sols;
    if(Circular_arc_2::template find_intersection< solutions_container >
      (a1,a2,early_sols)) {
      for (typename solutions_container::iterator it = early_sols.begin();
         it != early_sols.end(); ++it) {
        *res++ = *it;
      }
      return res;
    }
#endif

#ifdef  CGAL_CK_EXPLOIT_IDENTITY
    typedef typename CK::Circular_arc_point_2 Circular_arc_point_2;
    bool a1s_a2s = a1.source().equal_ref(a2.source());
    bool a1s_a2t = a1.source().equal_ref(a2.target());
    bool a1t_a2s = a1.target().equal_ref(a2.source());
    bool a1t_a2t = a1.target().equal_ref(a2.target());

    if((a1s_a2s && a1t_a2t) || (a1s_a2t && a1t_a2s)){ // Case 1
      if( (a1.supporting_circle() == a2.supporting_circle()) && ((a1.on_upper_part() && a2.on_upper_part())|| (! a1.on_upper_part() && (! a2.on_upper_part())))){
        *res++ = result_type(a1);
      } else {
        if(compare_x<CK>(a1.source(), a1.target()) == SMALLER){
          *res++ = result_type(std::make_pair(a1.source(),1u));
          *res++ = result_type(std::make_pair(a1.target(),1u));
        } else {
          *res++ = result_type(std::make_pair(a1.target(),1u));
          *res++ = result_type(std::make_pair(a1.source(),1u));
        }
      }
      return res;
    } else if (a1s_a2s || a1t_a2t || a1s_a2t || a1t_a2s) {
      Circular_arc_point_2 p,q,r;

      // Make that q is the middle vertex
      if(a1s_a2s){
        p = a1.target();
        q = a1.source();
        r = a2.target();
      } else if(a1s_a2t){
        p = a1.target();
        q = a1.source();
        r = a2.source();
      } else if(a1t_a2s){
        p = a1.source();
        q = a1.target();
        r = a2.target();
      } else { // a1t_a2t
        p = a1.source();
        q = a1.target();
        r = a2.source();
      }

      bool return_q = false;
      if(CircularFunctors::compare_x<CK>(r,q) == LARGER){
        if (CircularFunctors::point_in_x_range<CK>(p,r,q)){ // Case 2
          return_q = true;
        }  else if (((a1.on_upper_part() && ! a2.on_upper_part()) && CircularFunctors::compare_y_to_right<CK>(a1,a2,q) == SMALLER)
                    || ((! a1.on_upper_part() && a2.on_upper_part()) && CircularFunctors::compare_y_to_right<CK>(a1,a2,q) == LARGER)){
          return_q = true;
        } else if ((a1.on_upper_part() && ! a2.on_upper_part()) && CircularFunctors::compare_y_to_right<CK>(a1,a2,q)==LARGER){
          typename CK::Linear_kernel::Bounded_side p_a2_bs = CircularFunctors::bounded_side<CK>(a2.supporting_circle(),p);
          typename CK::Linear_kernel::Bounded_side r_a1_bs = CircularFunctors::bounded_side<CK>(a1.supporting_circle(),r);
          if(p_a2_bs || r_a1_bs){
            return_q = true;
          } else {
          }
        }
      } else {
        // TODO: treat the cases where the common endpoint is on the right
      }

      if(return_q){

        *res++ = result_type(std::make_pair(q,1u));
        return res;
      }

    }
#endif // CGAL_CK_EXPLOIT_IDENTITY

    const bool sqr1_eq_sqr2 = (a1.squared_radius() == a2.squared_radius());
    const bool c1_eq_c2 = (a1.center() == a2.center());
    typedef typename CK2_Intersection_traits<CK, typename CK::Circular_arc_2,
                                                 typename CK::Circular_arc_2>::type result_type;

    if(sqr1_eq_sqr2 && c1_eq_c2) {
      if(a1.is_full()) {
        *res++ =CGAL::internal::ck2_intersection_return<result_type>(a2);
        //return res;
      }
      else if(a2.is_full()) {
        *res++ =CGAL::internal::ck2_intersection_return<result_type>(a1);
        //return res;
      } else {
        bool t2_in_a1 = has_on<CK>(a1,a2.target(),true);
        bool s2_in_a1 = has_on<CK>(a1,a2.source(),true);
        if(t2_in_a1 && s2_in_a1) {
          bool t1_in_a2 = has_on<CK>(a2,a1.target(),true);
          bool s1_in_a2 = has_on<CK>(a2,a1.source(),true);
          if(t1_in_a2 && s1_in_a2) {
            const Comparison_result comp =
              CircularFunctors::compare_xy<CK>(a1.source(), a2.source());
            if(comp < 0) {
              if(a1.source() == a2.target()) {
                *res++ =CGAL::internal::ck2_intersection_return<result_type>(std::make_pair(a1.source(),1u));
              } else {
                const Circular_arc_2 & arc =
                Circular_arc_2(a1.supporting_circle(),a1.source(),a2.target());
                *res++ =CGAL::internal::ck2_intersection_return<result_type>(arc);
              }
              if(a2.source() == a1.target()) {
                *res++ =CGAL::internal::ck2_intersection_return<result_type>(std::make_pair(a2.source(),1u));
              } else {
                const Circular_arc_2 & arc =
                Circular_arc_2(a1.supporting_circle(),a2.source(),a1.target());
                *res++ =CGAL::internal::ck2_intersection_return<result_type>(arc);
              }
            } else if (comp > 0) {
              if(a2.source() == a1.target()) {
                *res++ =CGAL::internal::ck2_intersection_return<result_type>(std::make_pair(a2.source(),1u));
              } else {
                const Circular_arc_2 & arc =
                Circular_arc_2(a1.supporting_circle(),a2.source(),a1.target());
                *res++ =CGAL::internal::ck2_intersection_return<result_type>(arc);
              }
              if(a1.source() == a2.target()) {
                *res++ =CGAL::internal::ck2_intersection_return<result_type>(std::make_pair(a1.source(),1u));
              } else {
                const Circular_arc_2 & arc =
                Circular_arc_2(a1.supporting_circle(),a1.source(),a2.target());
                *res++ =CGAL::internal::ck2_intersection_return<result_type>(arc);
              }
            } else {
              *res++ =CGAL::internal::ck2_intersection_return<result_type>(a1);
            }
          } else {
            *res++ =CGAL::internal::ck2_intersection_return<result_type>(a2);
          //return res;
          }
        }
        else if(t2_in_a1) {
          if(a1.source() == a2.target())
            *res++ =CGAL::internal::ck2_intersection_return<result_type>(std::make_pair(a1.source(),1u));
          else {
            const Circular_arc_2 & arc =
              Circular_arc_2(a1.supporting_circle(),a1.source(),a2.target());
            *res++ =CGAL::internal::ck2_intersection_return<result_type>(arc);
          } //return res;
        } else if(s2_in_a1) {
          if(a2.source() == a1.target()) {
            *res++ =CGAL::internal::ck2_intersection_return<result_type>(std::make_pair(a2.source(),1u));
          } else {
            const Circular_arc_2 & arc =
              Circular_arc_2(a1.supporting_circle(),a2.source(),a1.target());
            *res++ =CGAL::internal::ck2_intersection_return<result_type>(arc);
          } //return res;
        } else if(has_on<CK>(a2,a1.source(),true)) {
          *res++ =CGAL::internal::ck2_intersection_return<result_type>(a1);
        //return res;
        }
      //return res;
      }
    } else if(!c1_eq_c2) {
      solutions_container solutions;

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
      if(!Circular_arc_2::template
         find_intersection_circle_circle< solutions_container >
         (a1,a2,solutions)) {
#endif

      intersection( a1.supporting_circle(), a2.supporting_circle(),
        std::back_inserter(solutions) );

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
        Circular_arc_2::template
          put_intersection_circle_circle< std::vector < CGAL::Object > >
          (a1,a2,solutions);
      }
#endif

      if(solutions.size() == 0) return res;
      else {
        // The supporting circles are not the same and intersects
        for (typename solutions_container::iterator it = solutions.begin();
         it != solutions.end(); ++it) {
          const std::pair<typename CK::Circular_arc_point_2, unsigned>
            *result = CGAL::object_cast
              <std::pair<typename CK::Circular_arc_point_2, unsigned> > (&(*it));

#ifdef CGAL_CK_TEST_BBOX_BEFORE_HAS_ON
          Bbox_2 rb = result->first.bbox();
          if(do_overlap(a1.bbox(), rb) && do_overlap(a2.bbox(),rb)){
            if (has_on<CK>(a1,result->first,true) &&
                has_on<CK>(a2,result->first,true)) {
              *res++ =CGAL::internal::ck2_intersection_return<result_type>(*result);
            }
          }
#else
          if (has_on<CK>(a1,result->first,true) &&
              has_on<CK>(a2,result->first,true)) {
            *res++ =CGAL::internal::ck2_intersection_return<result_type>(*result);
          }
#endif
        }
      //return res;
      }
    }
    return res;
  }


  // !!!! a lot of useless assertions for debug
  /*template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Circular_arc_2 &a1,
               const typename CK::Circular_arc_2 &a2,
               OutputIterator res )
  {
    typedef typename CK::Circular_arc_point_2  Circular_arc_point_2;
    typedef typename CK::Circular_arc_2        Circular_arc_2;

    if (a1.is_x_monotone() && a2.is_x_monotone()) {
      // Overlapping curves.
      if ( (a1.supporting_circle() == a2.supporting_circle()) ||
           (a1.supporting_circle() == a2.supporting_circle().opposite()) ) {
        // The ranges need to overlap in order for the curves to overlap.
        if ( CircularFunctors::compare_x<CK>(a1.left(), a2.right()) > 0 ||
            CircularFunctors::compare_x<CK>(a2.left(), a1.right()) > 0)
          return res;

        // They both need to be on the same upper/lower part.
        if (a1.on_upper_part() != a2.on_upper_part()) {
          // But they could share the left vertical tangent point.
          if (a1.left() == a2.left())
            *res++ = make_object(std::make_pair(a1.left(),1u));
          // Or they could share the right vertical tangent point.
          if (a1.right() == a2.right())
            *res++ = make_object(std::make_pair(a1.right(),1u));
          return res;
        }

        // We know they overlap, determine the extremities of
        // the common subcurve
        // TODO : We should use std::max and std::min, but they
        // require less_x_2.
        const Circular_arc_2 & arctmp =
          CircularFunctors::compare_x<CK>(a1.right(), a2.right()) < 0 ? a1 : a2;
        // we know that the right endpoint is correct, let us look for
        // the left now:


        //? a1.left() : a2.left();
        if (CircularFunctors::compare_x<CK>(a1.left(), a2.left()) > 0) {
          //the left endpoint is a1's
          if (CircularFunctors::compare_x<CK>(a1.left(), a2.right()) < 0){
            if (a1.on_upper_part()) {
              const Circular_arc_2 & arc =
                Circular_arc_2(a1.supporting_circle(),a2.right(),a1.left());
              CGAL_kernel_assertion(arc.is_x_monotone());
              *res++ = make_object(arc);
            }
            else {
              const Circular_arc_2 & arc =
                Circular_arc_2(a1.supporting_circle(),a1.left(), a2.right());
              CGAL_kernel_assertion(arc.is_x_monotone());
              *res++ = make_object(arc);
            }
          }
          else
            *res++ = make_object(std::make_pair(arctmp.right(),1u));
        }
        else if( CircularFunctors::compare_x<CK>(a1.left(), a2.left()) < 0 ) {
          //the left endpoint is a2's
          if(CircularFunctors::compare_x<CK>(a1.right(), a2.left()) > 0) {
            if(a1.on_upper_part()){
              const Circular_arc_2 & arc =
                Circular_arc_2(a1.supporting_circle(), a1.right(), a2.left());
              CGAL_kernel_assertion(arc.is_x_monotone());
              *res++ = make_object(arc);
            }
            else {
              const Circular_arc_2 & arc =
                Circular_arc_2(a1.supporting_circle(), a2.left(), a1.right());
              CGAL_kernel_assertion(arc.is_x_monotone());
              *res++ = make_object(arc);
            }
          }
          else
            *res++ = make_object(std::make_pair(arctmp.right(),1u));
        }
        else {
          if(CircularFunctors::compare_x<CK>(a1.right(), a2.right()) >= 0)
            *res++ = make_object(a2);
          else if(CircularFunctors::compare_x<CK>(a1.right(), a2.right()) < 0)
            *res++ = make_object(a1);
          else
            *res++ = make_object(std::make_pair(arctmp.right(),1u));
        }
        return res;
      }

//      // We need to check that the supporting circles
//      // do intersect before going further.
//      if (! do_intersect(a1.supporting_circle(), a2.supporting_circle()))
//        { return res; }
//
//      // Get the two intersection points of the supporting circles.
//
//      std::vector<CGAL::Object > intersection_points;
//      CGAL::intersect_2<CK>
//        ( a1.supporting_circle(), a2.supporting_circle(),
//          std::back_inserter(intersection_points) );

      std::vector<CGAL::Object > intersection_points;
      CGAL::intersect_2<CK> ( a1.supporting_circle(), a2.supporting_circle(),
                              std::back_inserter(intersection_points) );
      if(intersection_points.size() == 0) return res;

      const Circular_arc_point_2 &left =
        (CGAL::object_cast< std::pair<Circular_arc_point_2, unsigned> >
         (&(intersection_points[0])))->first;
      if (intersection_points.size() < 2){// multiplicity 2
        if (CircularFunctors::has_on<CK>(a1, left) && CircularFunctors::has_on<CK>(a2, left))
          *res++ = make_object(std::make_pair(left,2u));
      }
      else {// multiplicity 1
        const Circular_arc_point_2 &right =
          (CGAL::object_cast< std::pair<Circular_arc_point_2, unsigned> >
           (&(intersection_points[1])))->first;
        // We also need to check that these intersection points are on the arc.
        if (CircularFunctors::has_on<CK>(a1, left) && CircularFunctors::has_on<CK>(a2, left))
          *res++ = make_object(std::make_pair(left,1u));
        if (CircularFunctors::has_on<CK>(a1, right) && CircularFunctors::has_on<CK>(a2, right))
          *res++ = make_object(std::make_pair(right,1u));
      }
      return res;
    }
    else {//a1 or a2 are not x_monotone
      std::vector< CGAL::Object > arcs_a1_x_monotone;
      make_x_monotone( a1, std::back_inserter(arcs_a1_x_monotone));
      std::vector< CGAL::Object > arcs_a2_x_monotone;
      make_x_monotone( a2, std::back_inserter(arcs_a2_x_monotone));
      std::vector< Circular_arc_2 > circle_arcs;
      std::vector< Circular_arc_point_2 > circle_arc_endpoints;

      for ( std::vector< CGAL::Object >::iterator it1 =
              arcs_a1_x_monotone.begin();
            it1 != arcs_a1_x_monotone.end(); ++it1 ) {
        //CGAL_kernel_assertion(assign( a1_aux, *it1));
        const Circular_arc_2 *a1_aux =
          CGAL::object_cast< Circular_arc_2 >(&*it1);

        for ( std::vector< CGAL::Object >::iterator it2 =
                arcs_a2_x_monotone.begin();
              it2 != arcs_a2_x_monotone.end(); ++it2 ) {
          //CGAL_kernel_assertion(assign( a2_aux, *it2));
          //assign( a2_aux, *it2);
          const Circular_arc_2 *a2_aux =
            CGAL::object_cast<Circular_arc_2>(&*it2);
          std::vector< CGAL::Object > res_aux;
          CircularFunctors::intersect_2<CK>( *a1_aux, *a2_aux, std::back_inserter(res_aux));
          if(res_aux.size() == 2){
            //it can't be a circular_arc_2
            //CGAL_kernel_assertion(assign(the_pair, res_aux[0]));
            const std::pair<Circular_arc_point_2, unsigned int> *the_pair1 =
              CGAL::object_cast<std::pair<Circular_arc_point_2, unsigned int> >
              (&res_aux[0]);
            Circular_arc_point_2 arc_end1 = the_pair1->first;
            //assign(the_pair, res_aux[1]);
            const std::pair<Circular_arc_point_2, unsigned int> *the_pair2 =
              CGAL::object_cast<std::pair<Circular_arc_point_2, unsigned int> >
              (&res_aux[1]);
            Circular_arc_point_2 arc_end2 = the_pair2->first;
            bool exist = false;
            for (typename std::vector< Circular_arc_point_2 >::iterator it
                   = circle_arc_endpoints.begin();
                 it != circle_arc_endpoints.end(); ++it ) {
              if (arc_end1 == *it) {
                exist = true;
                break;
              }
            }
            if (!exist) {
              circle_arc_endpoints.push_back(arc_end1);
            }
            else exist = false;
            for ( typename std::vector< Circular_arc_point_2 >::iterator it
                    = circle_arc_endpoints.begin();
                  it != circle_arc_endpoints.end(); ++it ) {
              if (arc_end2 == *it) {
                exist = true;
                break;
              }
            }
            if (!exist)
              circle_arc_endpoints.push_back(arc_end2);
          }
          else if( res_aux.size() == 1){
            //it can be a Circular_arc_point_2 or a Circular_arc_2
            if(const Circular_arc_2 *arc =
               CGAL::object_cast<Circular_arc_2>(&res_aux[0])) {
              //if(assign(arc,res_aux[0])){
              circle_arcs.push_back(*arc);
            }
            else {
              //CGAL_kernel_assertion(assign(the_pair, res_aux[0]));
              //assign(the_pair, res_aux[0]);
              const std::pair<Circular_arc_point_2, unsigned int> *the_pair =
                CGAL::object_cast<std::pair<Circular_arc_point_2, unsigned int> >
                (&res_aux[0]);
              Circular_arc_point_2 arc_end = the_pair->first;
              if (the_pair->second == 2u) {//there are only one tangent point
                *res++ = res_aux[0];
                return res;
              }
              bool exist = false;
              for (typename std::vector< Circular_arc_point_2 >::iterator it
                     = circle_arc_endpoints.begin();
                   it != circle_arc_endpoints.end(); ++it ) {
                if (arc_end == *it) {
                  exist = true;
                  break;
                }
              }
              if (!exist)
                circle_arc_endpoints.push_back(arc_end);
            }
          }
        }
      }
      //there are not double
      if (circle_arcs.size() > 0){
        std::size_t i = 1;
        while((i < circle_arcs.size()) &&
              (circle_arcs[i-1].target().x() == circle_arcs[i].source().x()) &&
              (circle_arcs[i-1].target().y() == circle_arcs[i].source().y())
              )
          {i++;}

        *res++ = make_object
          (Circular_arc_2(circle_arcs[0].supporting_circle(),
                          circle_arcs[0].source(),
                          circle_arcs[i-1].target()
                          ));
        if (i < circle_arcs.size()) {//there are 2 circle arcs
          std::size_t j = i;
          i++;
          while((i < circle_arcs.size())
                && (circle_arcs[i-1].target() == circle_arcs[i].source()))
            i++;
          *res++ = make_object
            (Circular_arc_2(circle_arcs[j].supporting_circle(),
                            circle_arcs[j].source(),
                            circle_arcs[i-1].target()
                            ));
          return res;
        }
        else {//There are one circle arc and there can be maximum one endpoint
          for (typename std::vector< Circular_arc_point_2 >::iterator it1
                 = circle_arc_endpoints.begin();
               it1 != circle_arc_endpoints.end(); ++it1 ) {
            bool other_point = true;
            for (typename std::vector< Circular_arc_2 >::iterator it2
                   = circle_arcs.begin();
                 it2 != circle_arcs.end(); ++it2 )
              {
                if (CircularFunctors::has_on<CK>(*it2, *it1)) {
                  other_point = false;
                  break;
                }
              }
            if (other_point) {
              *res++ = make_object(std::make_pair(*it1,1u));
              break;
            }
          }
          return res;
        }
      }
      else{//there are one or two endpoint
        if (circle_arc_endpoints.size() > 1){
          *res++ = make_object(std::make_pair(circle_arc_endpoints[0],1u));
          *res++ = make_object(std::make_pair(circle_arc_endpoints[1],1u));
        }
        else if (circle_arc_endpoints.size() == 1)
          *res++ = make_object(std::make_pair(circle_arc_endpoints[0],1u));
        return res;
      }
    }
  }*/

  template < class CK >
  bool
  is_vertical(const typename CK::Circular_arc_2 &)
  {
    return false;
  }

template < class CK, class OutputIterator >
  OutputIterator
  make_x_monotone( const typename CK::Circular_arc_2 &A,
                   OutputIterator res )
  {
    typedef typename CK::Circular_arc_2           Circular_arc_2;
    typedef typename CK::Circular_arc_point_2     Circular_arc_point_2;
    typedef typename CK::Root_for_circles_2_2     Root_for_circles_2_2;

    CGAL_kernel_precondition(A.supporting_circle().squared_radius() != 0);

    if (A.is_x_monotone()) {

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
      // get a number for its supporting circle
      A.circle_number();
#endif

      *res++ = make_object(A);
      return res;
    }

    std::vector< Root_for_circles_2_2 > vector_x_extremal_points;
    CircularFunctors::x_extremal_points<CK>(A.supporting_circle(),
                          std::back_inserter(vector_x_extremal_points));
    Circular_arc_point_2 x_extremal_point1 = vector_x_extremal_points[0];
    Circular_arc_point_2 x_extremal_point2 = vector_x_extremal_points[1];

#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
    std::vector < CGAL::Object > intersecs1;
    std::vector < CGAL::Object > intersecs2;
    std::vector < CGAL::Object > intersecs3;
#endif

    if (A.is_full()) {
      const Circular_arc_2 &ca1 = Circular_arc_2(A.supporting_circle(),
                                           x_extremal_point1,
                                           x_extremal_point2);
      const Circular_arc_2 &ca2 = Circular_arc_2(A.supporting_circle(),
                                           x_extremal_point2,
                                           x_extremal_point1);
      ca1._setx_info(2,1,0); //setting flags outside
      ca2._setx_info(2,2,0);

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
      // get a number for its supporting circle
      unsigned int cn = ca1.circle_number();
      ca2.set_circle_number(cn);
#endif

      *res++ = make_object(ca1);
      *res++ = make_object(ca2);

#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
      intersecs1.push_back(make_object(std::make_pair(x_extremal_point1,1u)));
      intersecs1.push_back(make_object(std::make_pair(x_extremal_point2,1u)));
      Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
        (ca1,ca2,intersecs1);
#endif

      return res;
    }

    int cmp_begin = CGAL::compare(A.source().y(), A.center().y());
    int cmp_end   = CGAL::compare(A.target().y(), A.center().y());

    // Define the 2 Circular_arc_endpoints
    // in the 2 vertical tangent points


    if (cmp_begin > 0) {
      const Circular_arc_2 &ca1 = Circular_arc_2(A.supporting_circle(),
                                                 A.source(),
                                                 x_extremal_point1);
      ca1._setx_info(2,2,0);

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
      unsigned int cn = ca1.circle_number();
#endif

      *res++ = make_object(ca1);
      if (cmp_end > 0) {
        // We must cut in 3 parts.
        const Circular_arc_2 &ca2 = Circular_arc_2(A.supporting_circle(),
                                                   x_extremal_point1,
                                                   x_extremal_point2);
        const Circular_arc_2 &ca3 = Circular_arc_2(A.supporting_circle(),
                                                   x_extremal_point2,
                                                   A.target());
        ca2._setx_info(2,1,0);
        ca3._setx_info(2,2,0);

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
        ca2.set_circle_number(cn);
        ca3.set_circle_number(cn);
#endif

        *res++ = make_object(ca2);
        *res++ = make_object(ca3);
#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
        intersecs1.push_back(make_object(std::make_pair(x_extremal_point1,1u)));
        intersecs2.push_back(make_object(std::make_pair(x_extremal_point2,1u)));
        Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
          (ca1,ca2,intersecs1);
        Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
          (ca2,ca3,intersecs2);
        Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
          (ca1,ca3,intersecs3); //empty - no intersection
#endif
      }
      else {
        const Circular_arc_2 &ca2 = Circular_arc_2(A.supporting_circle(),
                                                   x_extremal_point1,
                                                   A.target());
        ca2._setx_info(2,1,0);

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
        ca2.set_circle_number(cn);
#endif

        *res++ = make_object(ca2);

#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
        intersecs1.push_back(make_object(std::make_pair(x_extremal_point1,1u)));
        Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
          (ca1,ca2,intersecs1);
#endif

      }
    }
    else if (cmp_begin < 0) {
      const Circular_arc_2 &ca1 = Circular_arc_2(A.supporting_circle(),
                                                 A.source(),
                                                 x_extremal_point2);
      ca1._setx_info(2,1,0);

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
      unsigned int cn = ca1.circle_number();
#endif

      *res++ = make_object(ca1);
      if (cmp_end < 0) {
        const Circular_arc_2 &ca2 = Circular_arc_2(A.supporting_circle(),
                                                   x_extremal_point2,
                                                   x_extremal_point1);
        const Circular_arc_2 &ca3 = Circular_arc_2(A.supporting_circle(),
                                                   x_extremal_point1,
                                                   A.target());
        ca2._setx_info(2,2,0);
        ca3._setx_info(2,1,0);

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
        ca2.set_circle_number(cn);
        ca3.set_circle_number(cn);
#endif

        *res++ = make_object(ca2);
        *res++ = make_object(ca3);

#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
        intersecs1.push_back(make_object(std::make_pair(x_extremal_point2,1u)));
        intersecs2.push_back(make_object(std::make_pair(x_extremal_point1,1u)));
        Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
          (ca1,ca2,intersecs1);
        Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
          (ca2,ca3,intersecs2);
        Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
          (ca1,ca3,intersecs3);
#endif

      }
      else {
        const Circular_arc_2 &ca2 = Circular_arc_2(A.supporting_circle(),
                                                   x_extremal_point2,
                                                   A.target());
        ca2._setx_info(2,2,0);

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
        ca2.set_circle_number(cn);
#endif

        *res++ = make_object(ca2);

#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
        intersecs1.push_back(make_object(std::make_pair(x_extremal_point2,1u)));
        Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
          (ca1,ca2,intersecs1);
#endif

      }
    }
    else { // cmp_begin == 0
      if (CGAL::compare(A.source().x(), A.center().x()) < 0) {
        CGAL_kernel_assertion (cmp_end >= 0);
        const Circular_arc_2 &ca1 = Circular_arc_2(A.supporting_circle(),
                                                   x_extremal_point1,
                                                   x_extremal_point2);
        const Circular_arc_2 &ca2 = Circular_arc_2(A.supporting_circle(),
                                                   x_extremal_point2,
                                                   A.target());
        ca1._setx_info(2,1,0);
        ca2._setx_info(2,2,0);

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
        unsigned int cn = ca1.circle_number();
        ca2.set_circle_number(cn);
#endif

        *res++ = make_object(ca1);
        *res++ = make_object(ca2);

#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
        intersecs1.push_back(make_object(std::make_pair(x_extremal_point2,1u)));
        Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
          (ca1,ca2,intersecs1);
#endif

      }
      else {
        CGAL_kernel_assertion
          (CGAL::compare(A.source().x(), A.center().x()) > 0);
        CGAL_kernel_assertion (cmp_end != LARGER);
        const Circular_arc_2 &ca1 = Circular_arc_2(A.supporting_circle(),
                                                   x_extremal_point2,
                                                   x_extremal_point1);
        const Circular_arc_2 &ca2 = Circular_arc_2(A.supporting_circle(),
                                            x_extremal_point1,
                                            A.target());
        ca1._setx_info(2,2,0);
        ca2._setx_info(2,1,0);

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
        unsigned int cn = ca1.circle_number();
        ca2.set_circle_number(cn);
#endif

        *res++ = make_object(ca1);
        *res++ = make_object(ca2);

#ifdef CGAL_INTERSECTION_MAP_FOR_XMONOTONIC_ARC_WITH_SAME_SUPPORTING_CIRCLE
        intersecs1.push_back(make_object(std::make_pair(x_extremal_point1,1u)));
        Circular_arc_2::template put_intersection< std::vector < CGAL::Object > >
          (ca1,ca2,intersecs1);
#endif

      }
    }
    return res;
  }

// This is the make_x_monotone function returning extra information:
// The ouput iterator refers to pairs, the first part of which is an
// object containing the x-monotone arc and the second part is a
// boolean defining whether the arc is on the upper part of the
// circle or not. This extra information returned by make_x_monotone
// and make_xy_monotone helps us to avoid doing twice the same
// comparisons by the functions which call these two in order to define
// the position of the returned arcs on the circle , like in the
// construct_bounding_hexagons function

  template < class CK, class OutputIterator >
  OutputIterator
  advanced_make_x_monotone( const typename CK::Circular_arc_2 &A,
                            OutputIterator res )
  {
    typedef typename CK::Circular_arc_2           Circular_arc_2;
    typedef std::pair<CGAL::Object,bool >         S_pair;


    int cmp_begin_y = CGAL::compare
      (A.source().y(), A.supporting_circle().center().y());
    int cmp_end_y   = CGAL::compare
      (A.target().y(), A.supporting_circle().center().y());

    int cmp_x=compare_x(A.source(),A.target());

    // We don't need to split
    if ((cmp_begin_y != opposite(cmp_end_y)) &&
        ((((cmp_begin_y > 0) || (cmp_end_y > 0)) && (cmp_x > 0)) ||
         (((cmp_begin_y < 0) || (cmp_end_y < 0)) &&
         (cmp_x < 0)))) {
      *res++ = S_pair(make_object(A),(cmp_begin_y>0 || cmp_end_y>0) );
      return res;
    }

    // Half circles
    if (cmp_begin_y == 0 && cmp_end_y == 0 && cmp_x != 0) {
      *res++ = std::make_pair(make_object(A), cmp_x>0 );
      return res;
    }

    // We need to split
    //CGAL_assertion(!A.is_x_monotone());
    if (cmp_begin_y > 0) {

      *res++ = S_pair
        (make_object(Circular_arc_2(A.supporting_circle(), A.source(),
                                    CircularFunctors::x_extremal_point<CK>
                                    (A.supporting_circle(),true))),
         true);

      if (cmp_end_y > 0) {
        // We must cut in 3 parts.
        *res++ = std::make_pair
          (make_object(Circular_arc_2 (A.supporting_circle(),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),true),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),false))),
           false);

        *res++ = std::make_pair
          (make_object(Circular_arc_2 (A.supporting_circle(),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),false),
                                       A.target())),
           true);
      } else {
        *res++ = std::make_pair
          (make_object(Circular_arc_2 (A.supporting_circle(),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),true),
                                       A.target())),
           false);
      }
    }
    else if (cmp_begin_y < 0) {
      // Very similar to the previous case.
      *res++ = std::make_pair
        (make_object(Circular_arc_2 (A.supporting_circle(),
                                     A.source(),
                                     CircularFunctors::x_extremal_point<CK>
                                     (A.supporting_circle(),false))),
         false);

      if (cmp_end_y < CGAL::EQUAL) {
        // We must cut in 3 parts.
        *res++ = std::make_pair
          (make_object(Circular_arc_2 (A.supporting_circle(),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),false),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),true))) ,
           true );

        *res++ = std::make_pair
          (make_object(Circular_arc_2 (A.supporting_circle(),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),true),
                                       A.target())),
           false);
      } else {
        *res++ = std::make_pair
          (make_object(Circular_arc_2 (A.supporting_circle(),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),false),
                                       A.target())),
           true);
      }
    }
    else { // cmp_begin_y == 0
      if ( compare(A.source().x(),A.supporting_circle().center().x())< 0) {
        CGAL_assertion(cmp_end_y >= 0);
        *res++ = std::make_pair
          (make_object(Circular_arc_2 (A.supporting_circle(),
                                       A.source(),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),false))),
           false);

        *res++ = std::make_pair
          (make_object(Circular_arc_2 (A.supporting_circle(),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),false),
                                       A.target())),
           true);
      }
      else {
        CGAL_assertion( compare(A.source().x(),A.supporting_circle().center().x())< 0);
        CGAL_assertion(cmp_end_y != LARGER);
        *res++ = std::make_pair
          (make_object(Circular_arc_2 (A.supporting_circle(),
                                       A.source(),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),true))),
           true);

        *res++ = std::make_pair
          (make_object(Circular_arc_2 (A.supporting_circle(),
                                       CircularFunctors::x_extremal_point<CK>
                                       (A.supporting_circle(),true),
                                       A.target())),
           false);
      }
    }

    return res;
  }

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


// In the same as the advanced_make_x_monotone works, this make_xy_function
// returns extra information, descriptive of the position of the returned
// xy-monotone arcs on the circle: The output iterator refers to pairs, the
// first part of which is the object containing tha arc and the second part
// is another pair containing 2 booleans which equavalently describe whether the
// returned xy-monotone arc is on the upper part and the left side of the circle

template < typename CK , typename Output_iterator>
Output_iterator
advanced_make_xy_monotone( const typename CK::Circular_arc_2 &a,
                           Output_iterator res)
{

  typedef typename CK::Circular_arc_2            Circular_arc_2;
  typedef std::pair<bool, bool>                  relat_pos;
  typedef std::pair< CGAL::Object, bool>         Obj_descr_1;
  typedef std::pair< CGAL::Object, relat_pos>    Obj_descr_2;
  typedef std::vector<Obj_descr_1>               Obj_vector_1;
  typedef std::vector<Obj_descr_2>               Obj_vector_2;

  Obj_vector_1 vec;
  Obj_vector_2 vec2;
  Obj_descr_2  dscr2;

  advanced_make_x_monotone<CK>(a,std::back_inserter(vec));

  for(unsigned int i=0;i<vec.size();i++) {
    const Circular_arc_2 *tmp_arc =
      CGAL::object_cast<Circular_arc_2>(&vec.at(i).first);

    int cmp_begin_x = CGAL::compare
      (tmp_arc->source().x(), tmp_arc->supporting_circle().center().x());
    int cmp_end_x   = CGAL::compare
      (tmp_arc->target().x(), tmp_arc->supporting_circle().center().x());

    if(cmp_begin_x!=opposite(cmp_end_x) || cmp_begin_x==CGAL::EQUAL) {
      dscr2.first=vec.at(i).first;
      dscr2.second.first=vec.at(i).second;
      dscr2.second.second= (cmp_begin_x==CGAL::SMALLER ||
                            cmp_end_x==CGAL::SMALLER   ) ?
        true : false;
      *res++=dscr2; // The arc is xy_monotone
    }
    else{ //We have to split the x_monotone_arc into 2 y_monotone arcs

      Obj_descr_1 tmp=vec.at(i);
      Obj_descr_2 tmp1,tmp2;
      const Circular_arc_2 *tmp_arc =
        CGAL::object_cast<Circular_arc_2>(&tmp.first);

      tmp1.first = make_object
        (Circular_arc_2(a.supporting_circle(),tmp_arc->source(),
                        CircularFunctors::y_extremal_point<CK>
                        (a.supporting_circle(),!tmp.second)));

      tmp1.second.first=tmp.second;
      tmp1.second.second= (tmp.second)? false : true ;

      tmp2.first = make_object
        (Circular_arc_2(a.supporting_circle(),
                        CircularFunctors::y_extremal_point<CK>
                        (a.supporting_circle(),!tmp.second),
                        tmp_arc->target()));

      tmp2.second.first=tmp.second;
      tmp2.second.second= (tmp.second)? true  : false ;

      *res++=tmp1;
      *res++=tmp2;
    }

  }

  return res;

}

  template <class CK>
  CGAL::Bbox_2 circular_arc_bbox
  ( const typename CK::Kernel_base::Circular_arc_2 & a)
  {
    typedef CGAL::Interval_nt<false>::Protector IntervalProtector;
    typedef CGAL::Interval_nt<false> Interval;

    if(a.is_x_monotone()) {
        // The arc is xy-monotone so we just add the bboxes of the endpoints
      if(a.is_y_monotone())
        return a.left().bbox() + a.right().bbox();

      // Just x-monotone, so we have to find the y-critical point

      bool is_on_upper = a.on_upper_part();

      Bbox_2
        left_bb  = a.left().bbox(),
        right_bb = a.right().bbox();

      IntervalProtector ip;
      Interval cy = to_interval(a.center().y());
      Interval r2 = to_interval(a.squared_radius());
      Interval r = CGAL::sqrt(r2);

      double ymin, ymax;

      if(is_on_upper) {
        ymin = (CGAL::min)(left_bb.ymin(),right_bb.ymin());
        ymax = cy.sup() + r.sup();
      } else {
        ymin = cy.inf() - r.sup();
        ymax = (CGAL::max)(left_bb.ymax(),right_bb.ymax());
      }
      /*
      double ymin = (is_on_upper) ?
        (CGAL::min)(left_bb.ymin(),right_bb.ymin()) :
        to_interval
        ( CircularFunctors::y_extremal_point<CK>(a.supporting_circle(),true).y()).first;
      double ymax = (is_on_upper) ?
        to_interval
        ( CircularFunctors::y_extremal_point<CK>(a.supporting_circle(),false).y() ).second :
        CGAL::max(left_bb.ymax(),right_bb.ymax());
      */
      return Bbox_2(left_bb.xmin(),ymin,right_bb.xmax(),ymax);
    }

    if(a.is_y_monotone()) {
      bool is_on_left = a.on_left_part();
      IntervalProtector ip;
      Bbox_2
            source_bb  = a.source().bbox(),
            target_bb = a.target().bbox();
      Interval cx = to_interval(a.center().x());
      Interval r2 = to_interval(a.squared_radius());
      Interval r = CGAL::sqrt(r2);
      double xmin, xmax;
      if(is_on_left) {
        xmax = (CGAL::max)(source_bb.xmax(), target_bb.xmax());
        xmin = cx.inf() - r.sup();
      } else {
        xmax = cx.sup() + r.sup();
        xmin = (CGAL::min)(source_bb.xmin(), target_bb.xmin());
      }
      return Bbox_2(xmin,
                    (CGAL::min)(source_bb.ymin(),target_bb.ymin()),
                    xmax,
                    (CGAL::max)(source_bb.ymax(),target_bb.ymax()));
    }

    // Else return the bounding box of the circle.
    return a.supporting_circle().bbox();
    /*  More precise version for non-x-monotone arcs.
        double xmin,xmax,ymin,ymax;

        // In this case, we can't avoid doing these heavy comparisons

        Comparison_result cmp_source_x=compare(a.source().x(),a.supporting_circle().center().x()),
        cmp_target_x=compare(a.target().x(),a.supporting_circle().center().x()),
        cmp_source_y=compare(a.source().y(),a.supporting_circle().center().y()),
        cmp_target_y=compare(a.target().y(),a.supporting_circle().center().y());

        //Since it's not x-monotone, it must include at least one x-critical point

        if(cmp_source_y==cmp_target_y || cmp_source_y==0 || cmp_target_y==0)
        {
        if(cmp_source_x==cmp_target_x || cmp_source_x==0 || cmp_target_x==0)
        return a.supporting_circle().bbox();

        xmin=to_interval( x_extremal_points<CK>(a.supporting_circle(),true).x() ).first;
        xmax=to_interval( x_extremal_points<CK>(a.supporting_circle(),false).x() ).second;

        if( cmp_source_y==LARGER || cmp_target_y==LARGER)
        {
        ymin=to_interval( y_extremal_point<CK>(a.supporting_circle(),true).y() ).first;
        ymax=(CGAL::max)(to_interval(a.source().y()).second,to_interval(a.target().y()).second);
        }
        else{
        ymax=to_interval( y_extremal_point<CK>(a.supporting_circle(),false).y() ).second;
        ymin=(CGAL::min)(to_interval(a.source().y()).first,to_interval(a.target().y()).first);
        }

        return Bbox_2(xmin,ymin,xmax,ymax);
        }

        if(cmp_source_y > EQUAL)
        {
        xmin=to_interval(x_extremal_points<CK>(a.supporting_circle(),true).x()).first;
        xmax=(CGAL::max)(to_interval(a.source().x()).second,to_interval(a.target().x()).second);
        }
        else
        {
        xmin=(CGAL::min)(to_interval(a.source().x()).first,to_interval(a.target().x()).first);
        xmax=to_interval(x_extremal_points<CK>(a.supporting_circle(),false).x()).second;
        }


        if( ( cmp_source_y== LARGER && cmp_source_x>= EQUAL) ||
        ( cmp_target_y== LARGER && cmp_target_x<= EQUAL) )
        ymax=to_interval(y_extremal_point<CK>(a.supporting_circle(),false).y()).second;
        else
        ymax=(CGAL::max)(to_interval(a.source().y()).second,to_interval(a.target().y()).second);


        if( ( cmp_source_y== SMALLER && cmp_source_x<= EQUAL) ||
        ( cmp_target_y== SMALLER && cmp_target_x>= EQUAL) )
        ymin=to_interval(y_extremal_point<CK>(a.supporting_circle(),true).y()).first;
        else
        ymin=(CGAL::min)(to_interval(a.source().y()).first,to_interval(a.target().y()).first);

        return Bbox_2(xmin,ymin,xmax,ymax);
    */
  }

} // namespace CircularFunctors
} // namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_PREDICATES_ON_CIRCULAR_ARC_2_H
