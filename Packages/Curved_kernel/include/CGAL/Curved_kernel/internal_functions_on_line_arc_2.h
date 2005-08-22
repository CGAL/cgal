// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Julien Hazebrouck
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

#ifndef CGAL_CURVED_KERNEL_PREDICATES_ON_LINE_ARC_2_H
#define CGAL_CURVED_KERNEL_PREDICATES_ON_LINE_ARC_2_H

#include <CGAL/Simple_cartesian.h>

namespace CGAL {
namespace CircularFunctors {

template < class CK >
  bool
  point_in_range(const typename CK::Line_arc_2 &A,
                 const typename CK::Circular_arc_endpoint_2 &p) 
  {
    // range includes endpoints here
    return ((compare_x<CK>(p, A.source()) !=
	     compare_x<CK>(p, A.target())) 
	    || (compare_x(p, A.source()) == 
		CGAL::EQUAL));
  }


 template < class CK >
  bool
  equal(const typename CK::Line_arc_2 &A1,
        const typename CK::Line_arc_2 &A2)
  {
    if ((A1.supporting_line() != A2.supporting_line()) && 
	(A1.supporting_line() != A2.supporting_line().opposite()))
      return false;

      return ((equal<CK>(A1.source(), A2.source()) &&
	       equal<CK>(A1.target(), A2.target())) ||
	      (equal<CK>(A1.target(), A2.source()) &&
	       equal<CK>(A1.source(), A2.target())));
  }

  template < class CK >
  bool
  do_overlap(const typename CK::Line_arc_2 &A1,
	     const typename CK::Line_arc_2 &A2)
  {
    if ( (A1.supporting_line() != A2.supporting_line()) &&
	 (A1.supporting_line() != A2.supporting_line().opposite())) return false;

    return compare_xy<CK>(A1.right(), A2.left()) > 0
      && compare_xy<CK>(A1.left(), A2.right()) < 0;
  }


 template < class CK >
  bool
  has_on(const typename CK::Line_arc_2 &a,
	    const typename CK::Circular_arc_endpoint_2 &p)
  {
    typedef typename CK::Polynomial_1_2 Polynomial_1_2;
    Polynomial_1_2 equation = CGAL::LinearFunctors::get_equation<CK>(a.supporting_line());
    CGAL_kernel_precondition(p.x()*equation.a() ==
			     -equation.c() - p.y()*equation.b());
    
    return (compare_xy<CK>(p, a.source()) !=
	    compare_xy<CK>(p, a.target()));
  }

 template < class CK >
  Comparison_result
  compare_y_at_x(const typename CK::Circular_arc_endpoint_2 &p,
                 const typename CK::Line_arc_2 &A1)
  {
    CGAL_kernel_precondition (point_in_range<CK>(A1, p));
    //vertical case
    if(A1.source().x() == A1.target().x()){
      if (p.y() <= A1.right().y()){
	if(A1.left().y() <= p.y()){
	  return CGAL::EQUAL;
	}
	return CGAL::SMALLER;
      }
      return CGAL::LARGER;
    }
    //general case
    typedef typename CK::Polynomial_1_2 Polynomial_1_2;
    typedef typename CK::Root_of_2                Root_of_2;
    Polynomial_1_2 equation = CGAL::LinearFunctors::get_equation<CK>(A1.supporting_line());
    Root_of_2 y((-p.x()*equation.a() - equation.c())/equation.b());
    if(y == p.y())
      return CGAL::EQUAL;
    else if(y < p.y())
      return CGAL::LARGER;
    else return CGAL::SMALLER;
  }

  template < class CK >
  Comparison_result 
  compare_y_to_right(const typename CK::Line_arc_2 &A1,
		     const typename CK::Line_arc_2 &A2, 
		     const typename CK::Circular_arc_endpoint_2 &p)
  {
    //CGAL_kernel_precondition (A1.source().x() != A1.target().x());
    //CGAL_kernel_precondition (A2.source().x() != A2.target().x());
    if(A1.supporting_line().is_vertical()){
      if(A2.supporting_line().is_vertical())
	return CGAL::EQUAL;
      return CGAL::LARGER;
    }
    if(A2.supporting_line().is_vertical())
      return CGAL::SMALLER;
    
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
    typedef typename CK::Polynomial_1_2          Polynomial_1_2;
    typedef typename CK::Root_of_2               Root_of_2;
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
		     const typename CK::Circular_arc_endpoint_2 &p)
  {
    CGAL_kernel_precondition (A2.is_x_monotone());
    if(A1.supporting_line().is_vertical())
      return CGAL::LARGER;
    typedef typename CK::Polynomial_1_2          Polynomial_1_2;
    typedef typename CK::Root_of_2   Root_of_2;
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

    Polynomial_1_2 equation = CGAL::LinearFunctors::get_equation<CK>(A1.supporting_line());
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
      


    if (((tangent_1_x < 0) && (tangent_2_x > 0)) || ((tangent_1_x > 0) && (tangent_2_x < 0))){
      Root_of_2 prod_left = tangent_1_y * tangent_2_x;
      Root_of_2 prod_right = tangent_2_y * tangent_1_x;
      if (prod_left < prod_right)
	return CGAL::LARGER;
      if (prod_left == prod_right)
	return  A2.on_upper_part() ? CGAL::LARGER : CGAL::SMALLER;
      return CGAL::SMALLER;
    }
    else{
      Root_of_2 prod_left = tangent_1_y * tangent_2_x;
      Root_of_2 prod_right = tangent_2_y * tangent_1_x;
      if (prod_left < prod_right)
	return CGAL::SMALLER;
      if (prod_left == prod_right)
	return A2.on_upper_part() ? CGAL::LARGER : CGAL::SMALLER;
      return CGAL::LARGER;
    }
  }
   
     template < class CK >
  Comparison_result 
  compare_y_to_right(const typename CK::Circular_arc_2 &A1,
		     const typename CK::Line_arc_2 &A2, 
		     const typename CK::Circular_arc_endpoint_2 &p)
  {
    if (compare_y_to_right<CK>(A2, A1, p) == CGAL::LARGER)
      return CGAL::SMALLER;
    return CGAL::LARGER;
  }

   
 template < class CK >
   void
  split(const typename CK::Line_arc_2 &A,
	const typename CK::Circular_arc_endpoint_2 &p,
	typename CK::Line_arc_2 &ca1,
	typename CK::Line_arc_2 &ca2)
  {   
    CGAL_kernel_precondition( has_on<CK>(A, p));
   
    typedef typename CK::Line_arc_2  Line_arc_2;

    ca1 = Line_arc_2( A.supporting_line(), A.source(), p);
    ca2 = Line_arc_2( A.supporting_line(), p, A.target());

   if ( ca1.right()!=ca2.left() )
	    {
	      std::swap(ca1,ca2);
	    }
    
    return;
  }

 template< class CK, class OutputIterator>
   OutputIterator
   construct_intersections_2( const typename CK::Line_arc_2 &a1,
			      const typename CK::Line_arc_2 &a2,
			      OutputIterator res )
  {
    typedef typename CK::Circular_arc_endpoint_2  Circular_arc_endpoint_2;
    typedef typename CK::Line_arc_2               Line_arc_2;
    typedef typename CK::Linear_kernel::Point_2   Point_2;
    typedef typename CK::Root_of_2                Root_of_2;
    typedef typename CGAL::Simple_cartesian<Root_of_2>::Point_2
      Numeric_point_2;
    if ((a1.supporting_line() == a2.supporting_line())
	|| (a1.supporting_line() == a2.supporting_line().opposite())){
      if(compare_xy(a1.left(),a2.left()) < 0){
	int comparison = compare_xy(a2.left(),a1.right());
	if(comparison < 0){
	  if(compare_xy(a1.right(),a2.right()) <= 0){
	    *res++ = make_object
	      (Line_arc_2(a1.supporting_line(),
			  a2.left(),
			  a1.right()
			  ));
	  }
	  else{
	    *res++ = make_object
	      (Line_arc_2(a1.supporting_line(),
			  a2.left(),
			  a2.right()
			  ));
	  }
	}
	else if (comparison == 0){
	  *res++ =make_object
	    ( std::make_pair(a2.left(),1u));
	}
	return res;
      }
      else{
	int comparison = compare_xy(a1.left(),a2.right());
	if(comparison < 0){
	  if(compare_xy(a1.right(),a2.right()) <= 0){
	    *res++ = make_object
	      (Line_arc_2(a1.supporting_line(),
			  a1.left(),
			  a1.right()
			  ));
	  }
	  else{
	    *res++ = make_object
	      (Line_arc_2(a1.supporting_line(),
			  a1.left(),
			  a2.right()
			  ));
	  }
	}
	else if (comparison == 0){
	  *res++ = make_object
	    ( std::make_pair(a1.left(),1u));
	}
	return res;
      }
    }
    
    if(!do_intersect(a1.supporting_line(), a2.supporting_line()))
      return res;
    Object obj = intersection(a1.supporting_line(), a2.supporting_line());
    const Point_2 *pt = CGAL::object_cast<Point_2>(&obj);
    Circular_arc_endpoint_2 intersect_point = Circular_arc_endpoint_2(
					       Numeric_point_2(
							       Root_of_2(pt->x()),
							       Root_of_2(pt->y())
							       )
					       );
    if ((has_on<CK>(a1, intersect_point)) && (has_on<CK>(a2, intersect_point)))
      *res++ = make_object(std::make_pair(intersect_point, 1u));
    
    return res;
  }



 template< class CK, class OutputIterator>
   OutputIterator
   construct_intersections_2( const typename CK::Line_arc_2 &l,
			      const typename CK::Circle_2 &c,
			      OutputIterator res )
  { 
    typedef std::vector<CGAL::Object >
      solutions_container;
    
    solutions_container solutions;
    CGAL::LinearFunctors::construct_intersections_2<CK>
      ( l.supporting_line(), c, std::back_inserter(solutions) );
    
    for (typename solutions_container::iterator it = solutions.begin(); it != solutions.end(); ++it){
      const std::pair<typename CK::Circular_arc_endpoint_2, uint> *result;
      result = CGAL::object_cast< 
	std::pair<typename CK::Circular_arc_endpoint_2, uint> >(&(*it));
      if ( has_on<CK>(l,result->first ))
	*res++ = *it;
    }
    return res;
  }

  template< class CK, class OutputIterator>
   OutputIterator
   construct_intersections_2( const typename CK::Circle_2 &c,
			      const typename CK::Line_arc_2 &l,
			      OutputIterator res )
  { 
    return construct_intersections_2<CK>(l,c,res);
  }


 template< class CK, class OutputIterator>
   OutputIterator
   construct_intersections_2( const typename CK::Line_arc_2 &l,
			      const typename CK::Circular_arc_2 &c,
			      OutputIterator res )
  {
    typedef typename CK::Circular_arc_2 Circular_arc_2;
    typedef std::vector<CGAL::Object >
      solutions_container;
    
    solutions_container solutions;
    CGAL::LinearFunctors::construct_intersections_2<CK>
      ( l.supporting_line(), c.supporting_circle(), std::back_inserter(solutions) );
    
    solutions_container objects_monotone;
    std::vector<const Circular_arc_2*> arcs_x_monotone; 
    make_x_monotone( c, std::back_inserter(objects_monotone));
    for(typename solutions_container::iterator it2 = objects_monotone.begin();
	it2 != objects_monotone.end(); ++it2){
      arcs_x_monotone.push_back(CGAL::object_cast<Circular_arc_2>(&(*it2)));
    }
      
    for (typename solutions_container::iterator it = solutions.begin(); it != solutions.end(); ++it){
      const std::pair<typename CK::Circular_arc_endpoint_2, uint> *result;
      result = CGAL::object_cast< 
	std::pair<typename CK::Circular_arc_endpoint_2, uint> >(&(*it));
      if ( has_on<CK>(l,result->first)){
	bool is_on_arc = false;
	for(typename std::vector<const Circular_arc_2*>::iterator it2 = arcs_x_monotone.begin();
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
  }

  template< class CK, class OutputIterator>
   OutputIterator
   construct_intersections_2( const typename CK::Circular_arc_2 &c,
			      const typename CK::Line_arc_2 &l,
			      OutputIterator res )
  {
    return construct_intersections_2<CK>(l,c,res);
  }
   
  template< class CK>
    bool
    is_vertical(const typename CK::Line_arc_2 &l)
  {
    return l.supporting_line().is_vertical();
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

#endif // CGAL_CURVED_KERNEL_PREDICATES_ON_LINE_ARC_2_H
