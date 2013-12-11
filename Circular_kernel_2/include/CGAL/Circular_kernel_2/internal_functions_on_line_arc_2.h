// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_PREDICATES_ON_LINE_ARC_2_H
#define CGAL_CIRCULAR_KERNEL_PREDICATES_ON_LINE_ARC_2_H

#include <CGAL/Circular_kernel_2/internal_functions_on_line_2.h>
#include <CGAL/Circular_kernel_2/internal_functions_on_circular_arc_2.h>
#include <CGAL/Circular_kernel_2/Intersection_traits.h>

namespace CGAL {
namespace CircularFunctors {

  template < class CK >
  bool
  point_in_x_range(const typename CK::Line_arc_2 &A,
		   const typename CK::Circular_arc_point_2 &p)
  {
    // range includes endpoints here
    return ( (CircularFunctors::compare_x<CK>(p, A.source()) != CircularFunctors::compare_x<CK>(p, A.target()))
	     || (CircularFunctors::compare_x<CK>(p, A.source()) == CGAL::EQUAL) );
  }

  template < class CK >
  bool
  equal(const typename CK::Line_arc_2 &A1,
        const typename CK::Line_arc_2 &A2)
  {
    if (!LinearFunctors::non_oriented_equal<CK>(
      A1.supporting_line(),A2.supporting_line()))
      return false;

    return ( (equal<CK>(A1.source(), A2.source()) &&
	      equal<CK>(A1.target(), A2.target())) ||
	     (equal<CK>(A1.target(), A2.source()) &&
	      equal<CK>(A1.source(), A2.target())) );
  }

  template < class CK >
  bool
  do_overlap(const typename CK::Line_arc_2 &A1,
	     const typename CK::Line_arc_2 &A2)
  {
    if (!LinearFunctors::non_oriented_equal<CK>(
      A1.supporting_line(),A2.supporting_line()))
      return false;

    return CircularFunctors::compare_xy<CK>(A1.right(), A2.left()) >= 0
        && CircularFunctors::compare_xy<CK>(A1.left(), A2.right()) <= 0;
  }


  template < class CK >
  bool
  has_on(const typename CK::Line_arc_2 &a,
	 const typename CK::Circular_arc_point_2 &p,
         const bool has_on_supporting_line = false)
  {
    if(p.equal_ref(a.source()) || p.equal_ref(a.target())){
      return true;
    }

    if(!has_on_supporting_line) {
      if(!CGAL::LinearFunctors::has_on<CK>(a.supporting_line(),p))
        return false;
    }

    return (CircularFunctors::compare_xy<CK>(p, a.source()) != CircularFunctors::compare_xy<CK>(p, a.target()));
  }


  template< class CK>
  bool
  is_vertical(const typename CK::Line_arc_2 &l)
  {
    return l.supporting_line().is_vertical();
  }

  template< class CK>
  bool
  is_x_monotone(const typename CK::Line_arc_2& )
  {
    return true;
  }

  template< class CK>
  bool
  is_y_monotone(const typename CK::Line_arc_2& )
  {
    return true;
  }

  template < class CK >
  Comparison_result
  compare_y_at_x(const typename CK::Circular_arc_point_2 &p,
                 const typename CK::Line_arc_2 &A1)
  {
    if(p.equal_ref(A1.source()) || p.equal_ref(A1.target())){
      return EQUAL;
    }
    //CGAL_kernel_precondition (CircularFunctors::point_in_x_range<CK>(A1, p));
    //vertical case
    if (CircularFunctors::is_vertical<CK>(A1)) {
      if (p.y() <= A1.right().y()) {
	if(A1.left().y() <= p.y()) {
	  return CGAL::EQUAL;
	}
	return CGAL::SMALLER;
      }
      return CGAL::LARGER;
    }
    //general case
    typedef typename CK::Polynomial_1_2    Polynomial_1_2;
    typedef typename CK::Root_of_2         Root_of_2;
    Polynomial_1_2 equation =
      CGAL::LinearFunctors::get_equation<CK>(A1.supporting_line());
    Root_of_2 y((-p.x()*equation.a() - equation.c())/equation.b());
    if (y == p.y())
      return CGAL::EQUAL;
    else if (y < p.y())
      return CGAL::LARGER;
    else return CGAL::SMALLER;
  }

  template < class CK >
  Comparison_result
  compare_y_to_right(const typename CK::Line_arc_2 &A1,
		     const typename CK::Line_arc_2 &A2,
		     const typename CK::Circular_arc_point_2 &/*p*/)
  {
    if(A1.supporting_line().is_vertical()){
      if(A2.supporting_line().is_vertical())
	return CGAL::EQUAL;
      return CGAL::LARGER;
    }
    if(A2.supporting_line().is_vertical())
      return CGAL::SMALLER;

    typedef typename CK::Polynomial_1_2       Polynomial_1_2;
    typedef typename CK::Root_of_2            Root_of_2;

    Polynomial_1_2 equation;
    if(A1.right().x() < A2.right().x()){
      equation = CGAL::LinearFunctors::get_equation<CK>(A2.supporting_line());
      Root_of_2 y((-A1.right().x()*equation.a() - equation.c())/equation.b());
      Root_of_2 A1_right_y = A1.right().y();
      if (y == A1_right_y)
	return CGAL::EQUAL;
      if (y < A1_right_y)
	return CGAL::LARGER;
      return CGAL::SMALLER;
    }
    else{
      equation = CGAL::LinearFunctors::get_equation<CK>(A1.supporting_line());
      Root_of_2 y((-A2.right().x()*equation.a() - equation.c())/equation.b());
      Root_of_2 A2_right_y = A2.right().y();
      if (y == A2_right_y)
	return CGAL::EQUAL;
      if (y < A2_right_y)
	return CGAL::SMALLER;
      return CGAL::LARGER;
    }

  }

  template < class CK >
  Comparison_result
  compare_y_to_right(const typename CK::Line_arc_2 &A1,
		     const typename CK::Circular_arc_2 &A2,
		     const typename CK::Circular_arc_point_2 &p)
  {
    //CGAL_kernel_precondition (A2.is_x_monotone());
    if(A1.supporting_line().is_vertical())
      return CGAL::LARGER;

    typedef typename CK::Polynomial_1_2          Polynomial_1_2;
    typedef typename CK::Root_of_2               Root_of_2;

    const typename CK::Circle_2 & C2 = A2.supporting_circle();

    Root_of_2 b2_y = C2.center().y() - p.y();

    int s_b2_y = CGAL::sign(b2_y);

    if (s_b2_y == 0) {
      // Vertical tangent for A1.
      return A2.on_upper_part() ? CGAL::SMALLER : CGAL::LARGER;
    }

    typename CK::Root_of_2 b2_x = C2.center().x() - p.x();

    Root_of_2 tangent_2_x;
    Root_of_2 tangent_2_y;
    if (b2_y < 0){
       tangent_2_x = -b2_y;
       tangent_2_y = b2_x;
    }
    else{
      tangent_2_x = b2_y;
      tangent_2_y = -b2_x;
    }

    Polynomial_1_2 equation =
      CGAL::LinearFunctors::get_equation<CK>(A1.supporting_line());
    typedef typename CK::FT FT;
    FT tangent_1_x;
    FT tangent_1_y;
    if (equation.b() < 0){
      tangent_1_x = -equation.b();
      tangent_1_y = equation.a();
    }
    else{
      tangent_1_x = equation.b();
      tangent_1_y = -equation.a();
    }



    if (((tangent_1_x < 0) && (tangent_2_x > 0)) ||
	((tangent_1_x > 0) && (tangent_2_x < 0))){
      Root_of_2 prod_left = tangent_1_y * tangent_2_x;
      Root_of_2 prod_right = tangent_2_y * tangent_1_x;
      if (prod_left < prod_right){
	return CGAL::LARGER;
      }
      if (prod_left == prod_right){
	return  A2.on_upper_part() ? CGAL::LARGER : CGAL::SMALLER;
      }
      return CGAL::SMALLER;
    }
    else{
      Root_of_2 prod_left = tangent_1_y * tangent_2_x;
      Root_of_2 prod_right = tangent_2_y * tangent_1_x;
      if (prod_left < prod_right){
	return CGAL::SMALLER;
      }
      if (prod_left == prod_right){
	return A2.on_upper_part() ? CGAL::LARGER : CGAL::SMALLER;
      }
      return CGAL::LARGER;
    }
  }


  template < class CK >
  Comparison_result
  compare_y_to_left(const typename CK::Line_arc_2 &A1,
		    const typename CK::Circular_arc_2 &A2,
		    const typename CK::Circular_arc_point_2 &p)
  {
    //CGAL_kernel_precondition (A2.is_x_monotone());
    if(A1.supporting_line().is_vertical())
      return CGAL::LARGER;

    typedef typename CK::Polynomial_1_2          Polynomial_1_2;
    typedef typename CK::Root_of_2               Root_of_2;

    const typename CK::Circle_2 & C2 = A2.supporting_circle();

    Root_of_2 b2_y = C2.center().y() - p.y();

    int s_b2_y = CGAL::sign(b2_y);

    if (s_b2_y == 0) {
      return A2.on_upper_part() ? CGAL::SMALLER : CGAL::LARGER;
    }

    typename CK::Root_of_2 b2_x = C2.center().x() - p.x();

    Root_of_2 tangent_2_x;
    Root_of_2 tangent_2_y;
    if (b2_y < 0){
       tangent_2_x = -b2_y;
       tangent_2_y = b2_x;
    }
    else{
      tangent_2_x = b2_y;
      tangent_2_y = -b2_x;
    }

    Polynomial_1_2 equation =
      CGAL::LinearFunctors::get_equation<CK>(A1.supporting_line());
    typedef typename CK::FT FT;
    FT tangent_1_x;
    FT tangent_1_y;
    if (equation.b() < 0){
      tangent_1_x = -equation.b();
      tangent_1_y = equation.a();
    }
    else{
      tangent_1_x = equation.b();
      tangent_1_y = -equation.a();
    }



    if (((tangent_1_x < 0) && (tangent_2_x > 0)) ||
	((tangent_1_x > 0) && (tangent_2_x < 0))){
      Root_of_2 prod_left = tangent_1_y * tangent_2_x;
      Root_of_2 prod_right = tangent_2_y * tangent_1_x;
      if (prod_left < prod_right){
	return CGAL::SMALLER;
      }
      if (prod_left == prod_right){
	return  A2.on_upper_part() ? CGAL::LARGER : CGAL::SMALLER;
      }
      return CGAL::LARGER;
    }
    else{
      Root_of_2 prod_left = tangent_1_y * tangent_2_x;
      Root_of_2 prod_right = tangent_2_y * tangent_1_x;
      if (prod_left < prod_right){
	return CGAL::LARGER;
      }
      if (prod_left == prod_right){
	      return A2.on_upper_part() ? CGAL::LARGER : CGAL::SMALLER;
      }
      return CGAL::SMALLER;
    }
  }


  template < class CK >
  Comparison_result
  compare_y_to_right(const typename CK::Circular_arc_2 &A1,
		     const typename CK::Line_arc_2 &A2,
		     const typename CK::Circular_arc_point_2 &p)
  {
    if (compare_y_to_right<CK>(A2, A1, p) == CGAL::LARGER)
      return CGAL::SMALLER;
    return CGAL::LARGER;
  }


  template < class CK >
  void
  split(const typename CK::Line_arc_2 &A,
	const typename CK::Circular_arc_point_2 &p,
	typename CK::Line_arc_2 &ca1,
	typename CK::Line_arc_2 &ca2)
  {
    CGAL_kernel_precondition( has_on<CK>(A, p));

    typedef typename CK::Line_arc_2  Line_arc_2;

    ca1 = Line_arc_2( A.supporting_line(), A.source(), p);
    ca2 = Line_arc_2( A.supporting_line(), p, A.target());

    if ( CircularFunctors::compare_xy<CK>(ca1.left(), ca2.left()) != SMALLER )
        {
          std::swap(ca1,ca2);
        }

    return;
  }

  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Line_2 & l,
	       const typename CK::Circle_2 & c,
	       OutputIterator res )
  {
    typedef typename CK2_Intersection_traits<CK, typename CK::Line_2, typename CK::Circle_2>
      ::type result_type;
    typedef typename CK::Algebraic_kernel            AK;
    typedef typename CK::Polynomial_1_2              Equation_line;
    typedef typename CK::Polynomial_for_circles_2_2  Equation_circle;
    typedef typename CK::Root_for_circles_2_2        Root_for_circles_2_2;

    Equation_line e1 = CK().get_equation_object()(l);
    Equation_circle e2 = CK().get_equation_object()(c);

    typedef std::vector< std::pair < Root_for_circles_2_2, unsigned > >
      solutions_container;
    solutions_container solutions;

    AK().solve_object()(e1, e2, std::back_inserter(solutions));
    // to be optimized

    typedef typename CK::Circular_arc_point_2 Circular_arc_point_2;

    for ( typename solutions_container::iterator it = solutions.begin();
	  it != solutions.end(); ++it )
      {
	*res++ = CGAL::internal::ck2_intersection_return<result_type>
	  (std::make_pair(Circular_arc_point_2(it->first), it->second ));
      }

    return res;
  }

  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Line_arc_2 &a1,
	       const typename CK::Line_arc_2 &a2,
	       OutputIterator res )
  {
    typedef typename CK2_Intersection_traits<CK, typename CK::Line_arc_2, typename CK::Line_arc_2>
      ::type result_type;
    typedef typename CK::Circular_arc_point_2  Circular_arc_point_2;
    typedef typename CK::Line_arc_2               Line_arc_2;
    typedef typename CK::Point_2                  Point_2;
    // typedef typename CK::Root_for_circles_2_2     Root_for_circles_2_2;

#ifdef CGAL_CK_EXPLOIT_IDENTITY
    bool a1s_a2s = a1.source().equal_ref(a2.source());
    bool a1s_a2t = a1.source().equal_ref(a2.target());
    bool a1t_a2s = a1.target().equal_ref(a2.source());
    bool a1t_a2t = a1.target().equal_ref(a2.target());

    if((a1s_a2s && a1t_a2t) || (a1s_a2t && a1t_a2s)){

      *res++ = result_type(a1);
      return res;
    }
    if(a1s_a2s || a1s_a2t || a1t_a2s || a1t_a2t) {
      if(! LinearFunctors::non_oriented_equal<CK>(a1.supporting_line(),a2.supporting_line())){
	if(a1s_a2s || a1s_a2t) *res++ = CGAL::internal::ck2_intersection_return<result_type>(std::make_pair(a1.source(), 1u));
	if(a1t_a2s || a1t_a2t) *res++ = CGAL::internal::ck2_intersection_return<result_type>(std::make_pair(a1.target(), 1u));
      return res;
      }
    }
#endif

    if(LinearFunctors::non_oriented_equal<CK>(
      a1.supporting_line(),a2.supporting_line())) {
      if(compare_xy(a1.left(),a2.left()) < 0) {
	int comparison = compare_xy(a2.left(),a1.right());
	if(comparison < 0){
	  if(compare_xy(a1.right(),a2.right()) <= 0){
	    *res++ = CGAL::internal::ck2_intersection_return<result_type>
	      (Line_arc_2(a1.supporting_line(), a2.left(), a1.right() ));
	  } else{
	    *res++ = CGAL::internal::ck2_intersection_return<result_type>
	      (Line_arc_2(a1.supporting_line(), a2.left(), a2.right() ));
	  }
	}
	else if (comparison == 0){
	  *res++ =CGAL::internal::ck2_intersection_return<result_type>
	    ( std::make_pair(a2.left(),1u));
	}
	return res;
      }
      else{
	int comparison = compare_xy(a1.left(),a2.right());
	if(comparison < 0){
	  if(compare_xy(a1.right(),a2.right()) <= 0){
	    *res++ = CGAL::internal::ck2_intersection_return<result_type>
	      (Line_arc_2(a1.supporting_line(), a1.left(), a1.right() ));
	  }
	  else{
	    *res++ = CGAL::internal::ck2_intersection_return<result_type>
	      (Line_arc_2(a1.supporting_line(), a1.left(), a2.right() ));
	  }
	}
	else if (comparison == 0){
	  *res++ = CGAL::internal::ck2_intersection_return<result_type>
	    ( std::make_pair(a1.left(),1u));
	}
	return res;
      }
    }

    typename Intersection_traits<CK, typename CK::Line_2, typename CK::Line_2>::result_type
      v = CGAL::internal::intersection(a1.supporting_line(), a2.supporting_line(), CK());
    if(!v) return res;

    const Point_2 *pt = CGAL::internal::intersect_get<Point_2>(v);
    if(pt == NULL) return res;
    Circular_arc_point_2 intersect_point = Circular_arc_point_2(*pt);
    //      (Root_for_circles_2_2(Root_of_2(pt->x()),Root_of_2(pt->y())));

    if ((CircularFunctors::compare_xy<CK>(intersect_point, a1.source()) !=
	 CircularFunctors::compare_xy<CK>(intersect_point, a1.target())) &&
	(CircularFunctors::compare_xy<CK>(intersect_point, a2.source()) !=
	 CircularFunctors::compare_xy<CK>(intersect_point, a2.target())))
      *res++ = CGAL::internal::ck2_intersection_return<result_type>(std::make_pair(intersect_point, 1u));

    return res;
  }

  template<typename CK, typename T>
  struct Has_on_visitor : public boost::static_visitor<bool> {
    Has_on_visitor(const T* l) : l(l){}
    const T* l;
    bool operator()(const std::pair<typename CK::Circular_arc_point_2, unsigned>& pair) const {
      return has_on<CK>(*l,pair.first,true);
    }
  };

  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Line_arc_2 &l,
	       const typename CK::Circle_2 &c,
	       OutputIterator res )
  {
    typedef std::vector<typename CK2_Intersection_traits<CK, typename CK::Line_2, typename CK::Circle_2>::type>
      solutions_container;
    solutions_container solutions;

    CircularFunctors::intersect_2<CK>
      ( l.supporting_line(), c, std::back_inserter(solutions) );

    for (typename solutions_container::iterator it = solutions.begin();
	 it != solutions.end(); ++it) {
      #if CGAL_INTERSECTION_VERSION  < 2
      if(const std::pair<typename CK::Circular_arc_point_2, unsigned>* p =
         object_cast< std::pair< typename CK::Circular_arc_point_2, unsigned> >(& (*it))) {
        Has_on_visitor<CK, typename CK::Line_arc_2> vis(&l);
        if(vis(*p)) {
	*res++ = *it;
    }
  }
      #else
      if(boost::apply_visitor(Has_on_visitor<CK, typename CK::Line_arc_2>(&l), *it))
	*res++ = *it;
      #endif
    }

    return res;
  }

  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Circle_2 &c,
	       const typename CK::Line_arc_2 &l,
	       OutputIterator res )
  {
    return intersect_2<CK>(l,c,res);
  }

  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Line_arc_2 &l,
	       const typename CK::Circular_arc_2 &c,
	       OutputIterator res )
  {
    typedef std::vector<CGAL::Object > solutions_container;
    typedef typename CK::Circular_arc_point_2 Circular_arc_point_2;

#if defined(CGAL_CK_EXPLOIT_IDENTITY) || \
  defined(CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES)
    typedef typename CK::Line_arc_2 Line_arc_2;
#endif

#ifdef CGAL_CK_EXPLOIT_IDENTITY
    typedef typename CK::Circular_arc_2 Circular_arc_2;
    typedef typename CK2_Intersection_traits<CK, Line_arc_2, Circular_arc_2 >
      ::type result_type;
    bool ls_cs = l.source().equal_ref(c.source());
    bool ls_ct = l.source().equal_ref(c.target());
    bool lt_cs = l.target().equal_ref(c.source());
    bool lt_ct = l.target().equal_ref(c.target());

    if((ls_cs && lt_ct) || (ls_ct && lt_cs)){ // Case 0
      if (compare_xy<CK>(l.source(), l.target()) == SMALLER){
	*res++ = result_type(std::make_pair(l.source(), 1u));
	*res++ = result_type(std::make_pair(l.target(), 1u));
      } else {
	*res++ = result_type(std::make_pair(l.target(), 1u));
	*res++ = result_type(std::make_pair(l.source(), 1u));
      }
      return res;
    } else if (ls_cs || lt_ct || ls_ct || lt_cs) {
      Circular_arc_point_2 p,q,r;

      if(ls_cs){
	p = l.target();
	q = l.source();
	r = c.target();
      } else if(ls_ct){
	p = l.target();
	q = l.source();
	r = c.source();
      } else if(lt_cs){
	p = l.source();
	q = l.target();
	r = c.target();
      } else { // lt_ct
	p = l.source();
	q = l.target();
	r = c.source();
      }

      if( (CircularFunctors::compare_x<CK>(p,q) == EQUAL) || CircularFunctors::point_in_x_range<CK>(p,r,q)){ // Case 1
	*res++ = result_type(std::make_pair(q,1u));
	return res;
      }
       else if(c.on_upper_part()){ // Case 2
	if(CircularFunctors::compare_x<CK>(r,q) == LARGER){
	  if(CircularFunctors::orientation<CK>(q,r,p) == RIGHT_TURN   // Case 2
	     || CircularFunctors::compare_y_to_right<CK>(l,c,q) == LARGER){ // Case 3a
	    *res++ = result_type(std::make_pair(q,1u));
	    return res;
	  }
	} else {
	  if(CircularFunctors::orientation<CK>(r,q,p) == RIGHT_TURN  // Case 2
	  || CircularFunctors::compare_y_to_left<CK>(l,c,q) == LARGER){ // Case 3c
	    *res++ = result_type(std::make_pair(q,1u));
	    return res;
	  }
	}
       } else { // c on lower part
	if(CircularFunctors::compare_x<CK>(r,q) == LARGER){
	  if (CircularFunctors::orientation<CK>(q,r,p) == LEFT_TURN // Case 2
	    || CircularFunctors::compare_y_to_right<CK>(l,c,q) == SMALLER){  // Case 3b
	    *res++ = result_type(std::make_pair(q,1u));
	    return res;
	  }
	} else {
	  if(CircularFunctors::orientation<CK>(r,q,p) == LEFT_TURN
	     || CircularFunctors::compare_y_to_left<CK>(l,c,q) == SMALLER){ // Case 3d
	    *res++ = result_type(std::make_pair(q,1u));
	    return res;
	  }
	}
       }

      typename CK::Linear_kernel::Bounded_side bs = CircularFunctors::bounded_side<CK>(c.supporting_circle(),p);
      if(bs == ON_BOUNDED_SIDE){
	*res++ = result_type(std::make_pair(q,1u));
	return res;
      } else { //Case 4b
        typedef std::vector<typename CK2_Intersection_traits<CK, typename CK::Line_2, typename CK::Circle_2>::type>
          container;
        container solutions;
	CGAL::CircularFunctors::intersect_2<CK>( l.supporting_line(), c.supporting_circle(),
					       std::back_inserter(solutions) );

	if(CircularFunctors::compare_x<CK>(r,q) == LARGER){
	  *res++ = result_type(std::make_pair(q,1u));
	  *res++ = solutions.back();
	  return res;
	} else {
	  *res++ = solutions.front();
	  *res++ = result_type(std::make_pair(q,1u));
	  return res;
	}

      } // Case 4b


      std::cout << "oops we missed a case" << std::endl;

    }
#endif // CGAL_CK_EXPLOIT_IDENTITY

    typedef std::vector<CGAL::Object > solutions_container;
    solutions_container solutions;

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
    if(!Line_arc_2::template
       find_intersection_circle_line< solutions_container >
       (c,l,solutions)) {
#endif

      CGAL::CircularFunctors::intersect_2<CK>
      ( l.supporting_line(), c.supporting_circle(),
	std::back_inserter(solutions) );

#ifdef CGAL_INTERSECTION_MAP_FOR_SUPPORTING_CIRCLES
        Line_arc_2::template
          put_intersection_circle_line< std::vector < CGAL::Object > >
          (c,l,solutions);
      }
#endif

    for (typename solutions_container::iterator it = solutions.begin();
	 it != solutions.end(); ++it) {
      const std::pair<Circular_arc_point_2, unsigned>
        *result = CGAL::object_cast
	  <std::pair<Circular_arc_point_2, unsigned> > (&(*it));
#ifdef CGAL_CK_TEST_BBOX_BEFORE_HAS_ON
	  Bbox_2 rb = result->first.bbox();
	  if(do_overlap(l.bbox(), rb) && do_overlap(c.bbox(),rb)){
	    if (has_on<CK>(l,result->first,true) &&
		has_on<CK>(c,result->first,true)) {
	      *res++ = *it;
	    }
	  }
#else
      if (has_on<CK>(l,result->first,true) &&
          has_on<CK>(c,result->first,true)) {
	*res++ = *it;
      }
#endif
    }
    return res;
  }




  /*template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Line_arc_2 &l,
	       const typename CK::Circular_arc_2 &c,
	       OutputIterator res )
  {
    typedef typename CK::Circular_arc_2 Circular_arc_2;
    typedef std::vector<CGAL::Object > solutions_container;

    solutions_container solutions;
    CGAL::LinearFunctors::intersect_2<CK>
      ( l.supporting_line(), c.supporting_circle(),
	std::back_inserter(solutions) );

    solutions_container objects_monotone;
    std::vector<const Circular_arc_2*> arcs_x_monotone;
    make_x_monotone( c, std::back_inserter(objects_monotone));
    for(typename solutions_container::iterator it2 = objects_monotone.begin();
	it2 != objects_monotone.end(); ++it2){
      arcs_x_monotone.push_back(CGAL::object_cast<Circular_arc_2>(&(*it2)));
    }

    for (typename solutions_container::iterator it = solutions.begin();
	 it != solutions.end(); ++it){
      const std::pair<typename CK::Circular_arc_point_2, unsigned> *result;
      result = CGAL::object_cast
	<std::pair<typename CK::Circular_arc_point_2, unsigned> > (&(*it));
      if ( has_on<CK>(l,result->first)) {
	bool is_on_arc = false;
	for(typename std::vector<const Circular_arc_2*>::iterator
	      it2 = arcs_x_monotone.begin();
	    it2 != arcs_x_monotone.end(); ++it2){
	  if(has_on<CK>(**it2, result->first)){
	    is_on_arc = true;
	    break;
	  }
	}
	if(is_on_arc)
	  *res++ = *it;
      }
    }
    return res;
  }*/

  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Line_2 &l,
	       const typename CK::Line_arc_2 &la,
	       OutputIterator res )
  {
    typedef typename CK::Circular_arc_point_2  Circular_arc_point_2;
    typedef typename CK::Point_2                  Point_2;
    typedef typename CK::Line_2                   Line_2;
    typedef typename CK::Line_arc_2               Line_arc_2;
    typedef typename CK2_Intersection_traits<CK, Line_2, Line_arc_2>::type result_type;

    if(LinearFunctors::non_oriented_equal<CK>(l, la.supporting_line())) {
      *res++ = result_type(la);
    }

    typename Intersection_traits<CK, Line_2, Line_2>::result_type
      v = intersection(l, la.supporting_line());

    if(!v) return res;
    const Point_2 *pt = boost::get<Point_2>(&*v);
    if(pt == NULL) return res;

    Circular_arc_point_2 intersect_point = Circular_arc_point_2(*pt);

    if (CircularFunctors::compare_xy<CK>(intersect_point, la.source()) !=
	 CircularFunctors::compare_xy<CK>(intersect_point, la.target()))
      *res++ = result_type(std::make_pair(intersect_point, 1u));

    return res;
  }

  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Line_2 &l,
	       const typename CK::Circular_arc_2 &c,
	       OutputIterator res )
  {
    typedef typename CK::Circular_arc_2 Circular_arc_2;

    typedef typename CK::Line_2 Line_2;
    typedef typename CK::Circle_2 Circle_2;
    typedef std::vector< typename CK2_Intersection_traits<CK, Line_2, Circle_2 >::type>
      solutions_container;

    solutions_container solutions;

    CGAL::CircularFunctors::intersect_2<CK>
      ( l, c.supporting_circle(),
	std::back_inserter(solutions) );

    for (typename solutions_container::const_iterator it = solutions.begin();
	 it != solutions.end(); ++it) {
#if CGAL_INTERSECTION_VERSION < 2
      typedef typename CK::Circular_arc_point_2 Circular_arc_point_2;
      const std::pair<Circular_arc_point_2, unsigned>* p =
        object_cast<std::pair<Circular_arc_point_2, unsigned> >(& (*it));
      Has_on_visitor<CK, Circular_arc_2> vis(&c);
      if(vis(*p)) *res++ = *it;
#else
      if(boost::apply_visitor(Has_on_visitor<CK, Circular_arc_2>(&c), *it))
        *res++ = *it;
#endif
    }
    return res;
  }

  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Circular_arc_2 &c,
	       const typename CK::Line_arc_2 &l,
	       OutputIterator res )
  {
    return intersect_2<CK>(l,c,res);
  }

  template < class CK, class OutputIterator >
  OutputIterator
  make_x_monotone( const typename CK::Line_arc_2 &A,
		   OutputIterator res )
  {
    *res++ = make_object(A);
    return res;
  }

} // namespace CircularFunctors
} // namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_PREDICATES_ON_LINE_ARC_2_H
