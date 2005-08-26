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

// file : include/CGAL/Curved_kernel/predicates_on_circular_arc_2.h

#ifndef CGAL_CURVED_KERNEL_PREDICATES_ON_CIRCULAR_ARC_2_H
#define CGAL_CURVED_KERNEL_PREDICATES_ON_CIRCULAR_ARC_2_H

#include <CGAL/Curved_kernel/internal_functions_on_circle_2.h>
#include <CGAL/global_functions_on_roots_and_polynomials_2_2.h>

namespace CGAL {
namespace CircularFunctors {
  

  template < class CK >
  inline
  Comparison_result 
  compare_x(const typename CK::Circular_arc_point_2 &p0,
            const typename CK::Circular_arc_point_2 &p1)
  {
    return CGAL::compare_x(p0.coordinates(), p1.coordinates());
  }

  template < class CK >
  inline
  Comparison_result 
  compare_y(const typename CK::Circular_arc_point_2 &p0,
            const typename CK::Circular_arc_point_2 &p1)
  {
    return CGAL::compare_y(p0.coordinates(), p1.coordinates());
  }

  template < class CK >
  Comparison_result 
  compare_xy(const typename CK::Circular_arc_point_2 &p0,
             const typename CK::Circular_arc_point_2 &p1)
  {
    return compare_xy(p0.coordinates(), p1.coordinates());
  }


  template < class CK >
  inline
  Comparison_result 
  compare_x(const typename CK::Circular_arc_point_2 &p0,
            const typename CK::Point_2 &p1)
  {
    return CGAL::compare(p0.x(), p1.x());
  }

  template < class CK >
  inline
  Comparison_result 
  compare_x(const typename CK::Point_2 &p0,
            const typename CK::Circular_arc_point_2 &p1)
  {
    return CGAL::compare(p0.x(), p1.x());
  }


  template < class CK >
  inline
  Comparison_result 
  compare_y(const typename CK::Circular_arc_point_2 &p0,
            const typename CK::Point_2 &p1)
  {
    return CGAL::compare(p0.y(), p1.y());
  }

  template < class CK >
  inline
  Comparison_result 
  compare_y(const typename CK::Point_2 &p0,
            const typename CK::Circular_arc_point_2 &p1)
  {
    return CGAL::compare(p0.y(), p1.y());
  }


  template < class CK >
  Comparison_result 
  compare_xy(const typename CK::Circular_arc_point_2 &p0,
             const typename CK::Point_2 &p1)
  {
    Comparison_result compx = compare_x<CK>(p0, p1);
    if (compx != 0)
      return compx;
    return compare_y<CK>(p0, p1);
  }

  template < class CK >
  Comparison_result 
  compare_xy(const typename CK::Point_2 &p0,
             const typename CK::Circular_arc_point_2 &p1)
  {
    Comparison_result compx = compare_x<CK>(p0, p1);
    if (compx != 0)
      return compx;
    return compare_y<CK>(p0, p1);
  }



  template < class CK >
  bool
  point_in_range(const typename CK::Circular_arc_2 &A,
                 const typename CK::Circular_arc_point_2 &p) 
  {
    CGAL_kernel_precondition (A.is_x_monotone());
    // range includes endpoints here
    return compare_x<CK>(p, A.source()) !=
           compare_x<CK>(p, A.target());
  }

  template < class CK >
  Comparison_result
  compare_y_at_x(const typename CK::Circular_arc_point_2 &p,
                 const typename CK::Circular_arc_2 &A1)
  {
    CGAL_kernel_precondition (A1.is_x_monotone());
    CGAL_kernel_precondition (point_in_range<CK>(A1, p)); 

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
		     const typename CK::Circular_arc_point_2 &p)
  {
    // FIXME : add preconditions to check that the 2 arcs are defined at
    // the right of the intersection.
    CGAL_kernel_precondition (A1.is_x_monotone());
    CGAL_kernel_precondition (A2.is_x_monotone());

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
      CGAL_kernel_assertion (b1 == false && b2 == true);
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
    return compare_xy<CK>(p0, p1) == 0;
  }

  template < class CK >
  bool
  equal(const typename CK::Circular_arc_2 &A1,
        const typename CK::Circular_arc_2 &A2)
  {
    CGAL_kernel_precondition (A1.is_x_monotone());
    CGAL_kernel_precondition (A2.is_x_monotone());

    if ( A1.supporting_circle() != A2.supporting_circle() )
      return false;

    return equal<CK>( A1.source(), A2.source() ) &&
           equal<CK>( A1.target(), A2.target() );
  }

  template < class CK >
  bool
  do_overlap(const typename CK::Circular_arc_2 &A1,
	     const typename CK::Circular_arc_2 &A2)
  {
    CGAL_kernel_precondition (A1.is_x_monotone());
    CGAL_kernel_precondition (A2.is_x_monotone());

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
  has_on(const typename CK::Circular_arc_2 &a,
	    const typename CK::Circular_arc_point_2 &p)
  {
    CGAL_kernel_precondition(a.is_x_monotone());

    typedef typename CK::Polynomial_for_circles_2_2 Polynomial_for_circles_2_2;
    Polynomial_for_circles_2_2 equation = get_equation<CK>(a.supporting_circle());
    if(CGAL::sign_at<typename CK::Algebraic_kernel>
       (equation,p.coordinates())!= ZERO)
      return false;
    
    if (! point_in_range<CK>(a, p) )
      return false;
 
    int cmp = CGAL::compare(p.y(), a.supporting_circle().center().y());

    return  cmp == 0 || (cmp > 0 &&  a.on_upper_part())
                     || (cmp < 0 && !a.on_upper_part());
  }

  template < class CK >
  void
  split(const typename CK::Circular_arc_2 &A,
	const typename CK::Circular_arc_point_2 &p,
	typename CK::Circular_arc_2 &ca1,
	typename CK::Circular_arc_2 &ca2)
  {
    CGAL_kernel_precondition( A.is_x_monotone() );
    CGAL_kernel_precondition( point_in_range<CK>( A, p ) );
    CGAL_kernel_precondition( A.on_upper_part() == (p.y() >
				  A.supporting_circle().center().y()));
    CGAL_kernel_precondition( has_on<CK>(A, p));
   
    typedef typename CK::Circular_arc_2  Circular_arc_2;

    ca1 = Circular_arc_2( A.supporting_circle(), A.source(), p);
    ca2 = Circular_arc_2( A.supporting_circle(), p, A.target());
    if ( ca1.right()!=ca2.left() )
	    {
	      //std::cout << " SWAP " << std::endl;
	      std::swap(ca1,ca2);
	    }
    
    return;
  }

  // !!!! a lot of useless assertions for debug
  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Circular_arc_2 &a1,
			     const typename CK::Circular_arc_2 &a2,
			     OutputIterator res )
  {
    typedef typename CK::Circular_arc_point_2  Circular_arc_point_2;
    typedef typename CK::Circular_arc_2           Circular_arc_2;

    if(a1.is_x_monotone() && a2.is_x_monotone()){
      // Overlapping curves.
      if (a1.supporting_circle() == a2.supporting_circle()) {
	
	// The ranges need to overlap in order for the curves to overlap.
	if (compare_x<CK>(a1.left(), a2.right()) > 0 ||
	    compare_x<CK>(a2.left(), a1.right()) > 0)
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
	
	// We know they overlap, determine the extremities of the common subcurve
	// TODO : We should use std::max and std::min, but they require less_x_2.
	const Circular_arc_2 & arctmp = 
	  compare_x<CK>(a1.right(), a2.right()) < 0 ? a1 : a2;
	// we know that the right endpoint is correct, let us look for
	// the left now:
	
	
	if ( compare_x<CK>(a1.left(), a2.left()) > 0 ) //? a1.left() : a2.left();
	  { //the left endpoint is a1's
	    if(compare_x<CK>(a1.left(), a2.right()) < 0){
	      if(a1.on_upper_part()){
		const Circular_arc_2 & arc =
		  Circular_arc_2(a1.supporting_circle(), a2.right(), a1.left());
		CGAL_kernel_assertion(arc.is_x_monotone());
		*res++ = make_object(arc);
	      }
	      else{
		const Circular_arc_2 & arc =
		  Circular_arc_2(a1.supporting_circle(), a1.left(), a2.right());
		CGAL_kernel_assertion(arc.is_x_monotone());
		*res++ = make_object(arc);
	      }
	    }
	    else
	      *res++ = make_object(std::make_pair(arctmp.right(),1u));
	  }
	else if( compare_x<CK>(a1.left(), a2.left()) < 0 ) //the left endpoint is a2's
	  {
	    if(compare_x<CK>(a1.right(), a2.left()) > 0){
	       if(a1.on_upper_part()){
		const Circular_arc_2 & arc =
		  Circular_arc_2(a1.supporting_circle(), a1.right(), a2.left());
		CGAL_kernel_assertion(arc.is_x_monotone());
		*res++ = make_object(arc);
	      }
	      else{
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
	    if(compare_x<CK>(a1.right(), a2.right()) >= 0)
	      *res++ = make_object(a2);
	    else if(compare_x<CK>(a1.right(), a2.right()) < 0)
	      *res++ = make_object(a1);
	    else
	      *res++ = make_object(std::make_pair(arctmp.right(),1u));
	}
	return res;
      }

      // We need to check that the supporting circles
      // do intersect before going further.
      if (! do_intersect(a1.supporting_circle(), a2.supporting_circle())) 
	{ return res; }
      
      // Get the two intersection points of the supporting circles.
      
      std::vector<CGAL::Object > intersection_points;
      CGAL::intersect_2<CK>
	( a1.supporting_circle(), a2.supporting_circle(),
	  std::back_inserter(intersection_points) );
      
      Circular_arc_point_2 left =
	(CGAL::object_cast< std::pair<Circular_arc_point_2, uint> >
	 (&(intersection_points[0])))->first;
      if (intersection_points.size() < 2){// multiplicity 2
	if (has_on<CK>(a1, left) && has_on<CK>(a2, left)) 
	  *res++ = make_object(std::make_pair(left,2u));
      }
      else {// multiplicity 1
	Circular_arc_point_2 right = 
	  (CGAL::object_cast< std::pair<Circular_arc_point_2, uint> >
	   (&(intersection_points[1])))->first;
	// We also need to check that these intersection points are on the arc.
	  if (has_on<CK>(a1, left) && has_on<CK>(a2, left))
	    *res++ = make_object(std::make_pair(left,1u));
	  if (has_on<CK>(a1, right) && has_on<CK>(a2, right))
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
	    it1 != arcs_a1_x_monotone.end(); ++it1 )
	{
	    //CGAL_kernel_assertion(assign( a1_aux, *it1));
	  const Circular_arc_2 *a1_aux = 
	    CGAL::object_cast< Circular_arc_2 >(&*it1);
	  for ( std::vector< CGAL::Object >::iterator it2 = 
		  arcs_a2_x_monotone.begin(); 
	    it2 != arcs_a2_x_monotone.end(); ++it2 )
	    {
	      //CGAL_kernel_assertion(assign( a2_aux, *it2));
	      //assign( a2_aux, *it2);
	      const Circular_arc_2 *a2_aux = 
		CGAL::object_cast<Circular_arc_2>(&*it2);
	      std::vector< CGAL::Object > res_aux;
	      intersect_2<CK>( *a1_aux, *a2_aux,
					     std::back_inserter(res_aux));
	      if(res_aux.size() == 2){
		//it can't be a circular_arc_2
		//CGAL_kernel_assertion(assign(the_pair, res_aux[0]));
		const std::pair<Circular_arc_point_2, std::size_t> *the_pair1 = 
		  CGAL::object_cast<std::pair<Circular_arc_point_2, std::size_t> >(&res_aux[0]);
		Circular_arc_point_2 arc_end1 = the_pair1->first;
		//assign(the_pair, res_aux[1]);
		const std::pair<Circular_arc_point_2, std::size_t> *the_pair2 = 
		  CGAL::object_cast<std::pair<Circular_arc_point_2, std::size_t> >(&res_aux[1]);
		Circular_arc_point_2 arc_end2 = the_pair2->first;
		bool exist = false;
		for (typename std::vector< Circular_arc_point_2 >::iterator it 
		       = circle_arc_endpoints.begin(); 
		      it != circle_arc_endpoints.end(); ++it )
		  {
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
		      it != circle_arc_endpoints.end(); ++it )
		  {
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
		   CGAL::object_cast<Circular_arc_2>(&res_aux[0])){
		//if(assign(arc,res_aux[0])){
		  circle_arcs.push_back(*arc);
		}
		else{
		  //CGAL_kernel_assertion(assign(the_pair, res_aux[0]));
		  //assign(the_pair, res_aux[0]);
		   const std::pair<Circular_arc_point_2, std::size_t> *the_pair = 
		     CGAL::object_cast< std::pair<Circular_arc_point_2, std::size_t> >(&res_aux[0]);
		  Circular_arc_point_2 arc_end = the_pair->first;
		  if (the_pair->second == 2u) {//there are only one tangent point
		    *res++ = res_aux[0];
		    return res;
		  }
		  bool exist = false;
		  for (typename std::vector< Circular_arc_point_2 >::iterator it 
			 = circle_arc_endpoints.begin(); 
		       it != circle_arc_endpoints.end(); ++it )
		    {
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
	      ) {i++;}

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
	else{//There are one circle arc and there can be maximum one endpoint
	   for (typename std::vector< Circular_arc_point_2 >::iterator it1 
		 = circle_arc_endpoints.begin(); 
	       it1 != circle_arc_endpoints.end(); ++it1 )
	    {
	      bool other_point = true;
	      for (typename std::vector< Circular_arc_2 >::iterator it2 
		     = circle_arcs.begin(); 
		   it2 != circle_arcs.end(); ++it2 )
		{
		  if (has_on<CK>(*it2, *it1)) {
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
  }
  
  

  template < class CK >
    bool
    is_vertical(const typename CK::Circular_arc_2 a)
  {
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
    CGAL_kernel_precondition(A.supporting_circle().squared_radius() != 0);
    int cmp_begin = CGAL::compare(A.source().y(), A.center().y());
    int cmp_end   = CGAL::compare(A.target().y(),   A.center().y());

    int cmp_x = compare_x(A.source(), A.target());

    // We don't need to split
    if (cmp_begin != opposite(cmp_end) &&
        (((cmp_begin > 0 || cmp_end > 0) && cmp_x > 0) ||
          (cmp_begin < 0 || cmp_end < 0) && cmp_x < 0) ) {
      *res++ = make_object(A); 
      return res; 
    }

    // Half circles
    if (cmp_begin == 0 && cmp_end == 0 && cmp_x != 0) {
      *res++ = make_object(A);
      return res; 
    }

    // We need to split
    CGAL_kernel_assertion(!A.is_x_monotone());

    // Define the 2 Circular_arc_endpoints 
    // in the 2 vertical tangent points
    Circular_arc_2 half_circle( A.supporting_circle(),
				x_critical_points<CK>(A.supporting_circle(),true),
				x_critical_points<CK>(A.supporting_circle(),false));

    
    if (cmp_begin > 0) {
      *res++ = make_object(Circular_arc_2 (A.supporting_circle(),
					   A.source(),
					   half_circle.source()));
      if (cmp_end > 0) {
        // We must cut in 3 parts.
        *res++ = make_object(Circular_arc_2(A.supporting_circle(),
					    half_circle.source(),
					    half_circle.target()));
        *res++ = make_object(Circular_arc_2 (A.supporting_circle(),
					     half_circle.target(),
					     A.target()));
      } else {
        *res++ = make_object(Circular_arc_2 (A.supporting_circle(),
					     half_circle.source(),
					     A.target()));
      }
    }
    else if (cmp_begin < 0) {
      // Very similar to the previous case.
      *res++ = make_object(Circular_arc_2 (A.supporting_circle(),
					   A.source(),
					   half_circle.target()));
      if (cmp_end < 0) {
        // We must cut in 3 parts.
        *res++ = make_object(Circular_arc_2 (A.supporting_circle(),
					     half_circle.target(),
					     half_circle.source()));
        *res++ = make_object(Circular_arc_2 (A.supporting_circle(),
					     half_circle.source(),
					     A.target()));
      } else {
        *res++ = make_object(Circular_arc_2 (A.supporting_circle(),
					     half_circle.target(),
					     A.target()));
      }
    }
    else { // cmp_begin == 0
      if (CGAL::compare(A.source().x(), A.center().x()) < 0) {
        CGAL_kernel_assertion (cmp_end >= 0);
        *res++ = make_object(half_circle);
        *res++ = make_object(Circular_arc_2 (A.supporting_circle(),
					     half_circle.target(),
					     A.target()));
      }
      else {
        CGAL_kernel_assertion (CGAL::compare(A.source().x(), A.center().x()) > 0);
        CGAL_kernel_assertion (cmp_end != LARGER);
        *res++ = make_object(Circular_arc_2 (A.supporting_circle(),
					     A.source(),
					     half_circle.source()));
        *res++ = make_object(Circular_arc_2 (A.supporting_circle(),
					     half_circle.source(),
					     A.target()));
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
    typedef typename CK::Circle_2                 Circle_2;
    typedef typename CK::FT                       FT;
    typedef typename CK::Linear_kernel::Point_2   Point_2;
    typedef std::pair<CGAL::Object,bool >         S_pair;


    int cmp_begin_y = CGAL::compare(A.source().y(), A.center().y());
    int cmp_end_y   = CGAL::compare(A.target().y(), A.center().y());
    
    int cmp_x=compare_x(A.source(),A.target());
    

    // We don't need to split
    if (cmp_begin_y != opposite(cmp_end_y) &&
        (((cmp_begin_y > 0 || cmp_end_y > 0) && cmp_x > 0) ||
          (cmp_begin_y < 0 || cmp_end_y < 0) && cmp_x < 0) ) {
            
      *res++ = S_pair(make_object(A),(cmp_begin_y>0 || cmp_end_y>0) );
      return res; 
    }

    // Half circles
    if (cmp_begin_y == 0 && cmp_end_y == 0 && cmp_x != 0) {
      *res++ = std::make_pair(make_object(A), cmp_x>0 );
      return res; 
    }

    // We need to split
    //assert(!A.is_x_monotone());



    if (cmp_begin_y > 0) {
    
      *res++ = S_pair(make_object(Circular_arc_2 (A.supporting_circle(),A.source(),
			          x_critical_points<CK>(A.supporting_circle(),true))),
			          true);
			       
      if (cmp_end_y > 0) {
        // We must cut in 3 parts.
        *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
			        x_critical_points<CK>(A.supporting_circle(),true),
			        x_critical_points<CK>(A.supporting_circle(),false))),
			        false);
			   
        *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
                                x_critical_points<CK>(A.supporting_circle(),false),
			        A.target())),true);
      } else {    
        *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
                                x_critical_points<CK>(A.supporting_circle(),true),
			        A.target())) , false);
      }
    }
    else if (cmp_begin_y < 0) {
      // Very similar to the previous case.
            
      *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
			      A.source(),
			      x_critical_points<CK>(A.supporting_circle(),false))),
			      false);
			 
      if (cmp_end_y < CGAL::EQUAL) {
      
        // We must cut in 3 parts.
	
        *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
			        x_critical_points<CK>(A.supporting_circle(),false), 
			        x_critical_points<CK>(A.supporting_circle(),true))) ,
			        true );
				     					     
					     
        *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
                                x_critical_points<CK>(A.supporting_circle(),true),
			        A.target())),false);
      } else {
        *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
                                x_critical_points<CK>(A.supporting_circle(),false),
			        A.target())) ,true);
      }
    }
    else { // cmp_begin_y == 0
      if ( compare(A.source().x(),A.center().x())< 0) {
        assert (cmp_end_y >= 0);
        *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
                                A.source(),
			        x_critical_points<CK>(A.supporting_circle(),false))),
			        false);
	
        *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
                                x_critical_points<CK>(A.supporting_circle(),false),
			        A.target())) ,true);
      }
      else {
        assert ( compare(A.source().x(),A.center().x())< 0);
        assert (cmp_end_y != LARGER);
        *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
			        A.source(),
			        x_critical_points<CK>(A.supporting_circle(),true))),
			        true);		   
			   
        *res++ = std::make_pair(make_object(Circular_arc_2 (A.supporting_circle(),
                                x_critical_points<CK>(A.supporting_circle(),true),
			        A.target())),false);
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
Output_iterator advanced_make_xy_monotone( const typename CK::Circular_arc_2 &a, 
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
	
       	for(unsigned int i=0;i<vec.size();i++)
	{
	
		const Circular_arc_2 *tmp_arc =CGAL::object_cast<Circular_arc_2>(&vec.at(i).first);
		
    		int cmp_begin_x = CGAL::compare(tmp_arc->source().x(), tmp_arc->center().x());
    		int cmp_end_x   = CGAL::compare(tmp_arc->target().x(), tmp_arc->center().x());
		
		if(cmp_begin_x!=opposite(cmp_end_x) || cmp_begin_x==CGAL::EQUAL)
		{
			   dscr2.first=vec.at(i).first;
			   dscr2.second.first=vec.at(i).second;
			   dscr2.second.second= (cmp_begin_x==CGAL::SMALLER ||
			   			 cmp_end_x==CGAL::SMALLER   )? 
			                         true : false;
			   *res++=dscr2; // The arc is xy_monotone
		}	   
		else{ //We have to split the x_monotone_arc into 2 y_monotone arcs
		
		
		
		Obj_descr_1 tmp=vec.at(i);
		Obj_descr_2 tmp1,tmp2;
		const Circular_arc_2 *tmp_arc =CGAL::object_cast<Circular_arc_2>(&tmp.first);
				 
					 			 
		tmp1.first=make_object(Circular_arc_2(a.supporting_circle(),tmp_arc->source(),
			               y_critical_points<CK>(a.supporting_circle(),!tmp.second)));
			    	    
			    
		tmp1.second.first=tmp.second;
		tmp1.second.second= (tmp.second)? false : true ;
		
		
		tmp2.first=make_object(Circular_arc_2(a.supporting_circle(),
                                       y_critical_points<CK>(a.supporting_circle(),!tmp.second),
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
 const CGAL::Bbox_2 circular_arc_bbox( const typename CK::Kernel_base::Circular_arc_2 & a)
{	
	typedef typename CK::Root_of_2 	   Root_of_2;
	typedef typename CK::FT 		   FT;

	if(a.is_x_monotone())
	{
		// The arc is xy-monotone so we just add the bboxes of the endpoints
		if(a.is_y_monotone())
			return a.left().bbox() + a.right().bbox();
					
		// Just x-monotone, so we have to find the y-critical point

		bool is_on_upper=a.on_upper_part();

		Bbox_2 left_bb=a.left().bbox(), 
		       right_bb=a.right().bbox();
		
		double ymin= (is_on_upper) ? CGAL::min(left_bb.ymin(),right_bb.ymin()) :
			   	to_interval( y_critical_points<CK>(a.supporting_circle(),true).y() ).first;
		double ymax= (is_on_upper) ? 
		                to_interval( y_critical_points<CK>(a.supporting_circle(),false).y() ).second:
			   	CGAL::max(left_bb.ymax(),right_bb.ymax()); 
		
				
		return Bbox_2(left_bb.xmin(),ymin,right_bb.xmax(),ymax);
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

		        xmin=to_interval( x_critical_points<CK>(a.supporting_circle(),true).x() ).first;
			xmax=to_interval( x_critical_points<CK>(a.supporting_circle(),false).x() ).second;

			if( cmp_source_y==LARGER || cmp_target_y==LARGER)
			{
				ymin=to_interval( y_critical_points<CK>(a.supporting_circle(),true).y() ).first;
				ymax=CGAL::max(to_interval(a.source().y()).second,to_interval(a.target().y()).second);
			}
			else{
				ymax=to_interval( y_critical_points<CK>(a.supporting_circle(),false).y() ).second;
				ymin=CGAL::min(to_interval(a.source().y()).first,to_interval(a.target().y()).first);
			}

			return Bbox_2(xmin,ymin,xmax,ymax);
		}
		
		if(cmp_source_y > EQUAL)
		{
			xmin=to_interval(x_critical_points<CK>(a.supporting_circle(),true).x()).first;
			xmax=CGAL::max(to_interval(a.source().x()).second,to_interval(a.target().x()).second);
		}
		else
		{
			xmin=CGAL::min(to_interval(a.source().x()).first,to_interval(a.target().x()).first);
			xmax=to_interval(x_critical_points<CK>(a.supporting_circle(),false).x()).second;
		}


		if( ( cmp_source_y== LARGER && cmp_source_x>= EQUAL) ||		
		    ( cmp_target_y== LARGER && cmp_target_x<= EQUAL) )
		    ymax=to_interval(y_critical_points<CK>(a.supporting_circle(),false).y()).second;
		else
		    ymax=CGAL::max(to_interval(a.source().y()).second,to_interval(a.target().y()).second);


		if( ( cmp_source_y== SMALLER && cmp_source_x<= EQUAL) ||		
		    ( cmp_target_y== SMALLER && cmp_target_x>= EQUAL) )
		    ymin=to_interval(y_critical_points<CK>(a.supporting_circle(),true).y()).first;
		else
		    ymin=CGAL::min(to_interval(a.source().y()).first,to_interval(a.target().y()).first);

		return Bbox_2(xmin,ymin,xmax,ymax);
		*/
}
  


} // namespace CircularFunctors 
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_PREDICATES_ON_CIRCULAR_ARC_2_H
