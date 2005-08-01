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
    CGAL_kernel_precondition (A.is_x_monotone());
    // range includes endpoints here
    return compare_x<CK>(p, A.source()) !=
           compare_x<CK>(p, A.target());
  }

  template < class CK >
  Comparison_result
  compare_y_at_x(const typename CK::Circular_arc_endpoint_2 &p,
                 const typename CK::Circular_arc_2 &A1)
  {
    // debug !!! 
    std::cout << "[compare_y_at_x]"
	      << std::endl << "arc : " << A1 << std::endl ; 
    std::cout << CGALi::print(std::cout, A1) << std::endl;
    std::cout << " point : " << p << std::endl;

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
		     const typename CK::Circular_arc_endpoint_2 &p)
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
  equal(const typename CK::Circular_arc_endpoint_2 &p0,
        const typename CK::Circular_arc_endpoint_2 &p1)
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
	    const typename CK::Circular_arc_endpoint_2 &p)
  {
    CGAL_kernel_precondition(a.is_x_monotone());
   // CGAL_kernel_precondition(a.supporting_circle() == p.circle(0) ||
   //        a.supporting_circle() == p.circle(1) );

    typedef typename CK::Polynomial_for_circles_2_2 Polynomial_for_circles_2_2;
    Polynomial_for_circles_2_2 equation = get_equation<CK>(a.supporting_circle());
    CGAL_kernel_precondition(square(p.x() - equation.a()) ==
			     equation.r_sq() - square(p.y() - equation.b()));
    
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
	      std::cout << " SWAP " << std::endl;
	      std::swap(ca1,ca2);
	    }
    
    return;
  }

  // !!!! a lot of useless assertions for debug
  template< class CK, class OutputIterator>
  OutputIterator
  construct_intersections_2( const typename CK::Circular_arc_2 &a1,
			     const typename CK::Circular_arc_2 &a2,
			     OutputIterator res )
  {
    typedef typename CK::Circular_arc_endpoint_2  Circular_arc_endpoint_2;
    typedef typename CK::Circular_arc_2           Circular_arc_2;

    if(a1.is_x_monotone() && a2.is_x_monotone()){
      std::cout << "<construct_intersection_monotone>" << std::endl; 
      // Overlapping curves.
      if (a1.supporting_circle() == a2.supporting_circle()) {
	
	// The ranges need to overlap in order for the curves to overlap.
	if (compare_x<CK>(a1.left(), a2.right()) > 0 ||
	    compare_x<CK>(a2.left(), a1.right()) > 0){
	  std::cout << "</construct_intersection_monotone>" << std::endl;
	  return res;}
	
	// They both need to be on the same upper/lower part.
	if (a1.on_upper_part() != a2.on_upper_part()) {
	  // But they could share the left vertical tangent point.
	  if (a1.left() == a2.left()) {
	    *res++ =make_object
	      ( std::make_pair(a1.left(),1u));
	  }
	  // Or they could share the right vertical tangent point.
	  if (a1.right() == a2.right()) {
	    *res++ = make_object
	      (std::make_pair(a1.right(),1u));
	  }
	  
	  std::cout << "</construct_intersection_monotone>" << std::endl;
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
	    
	    else{
	      *res++ = make_object
		(std::make_pair(arctmp.right(),1u));
	    }
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
	    else{
	      *res++ = make_object
		(std::make_pair(arctmp.right(),1u));
	    }
	  }
	else {
	    if(compare_x<CK>(a1.right(), a2.right()) >= 0){
	      
	      *res++ = make_object(a2);
	    }
	    else if(compare_x<CK>(a1.right(), a2.right()) < 0){
	      
	      *res++ = make_object(a1);
	    }
	    else{
	      *res++ = make_object
		(std::make_pair(arctmp.right(),1u));
	    }
	}
	std::cout << "</construct_intersection_monotone>" << std::endl;
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
	  if (has_on<CK>(a1, left) && has_on<CK>(a2, left)) {
	    *res++ = make_object(std::make_pair(left,1u));
	  }
	  if (has_on<CK>(a1, right) && has_on<CK>(a2, right)) {
	    *res++ = make_object(std::make_pair(right,1u));
	  }
	}
      else // multiplicity 2
	{
	  if (has_on<CK>(a1, left) && has_on<CK>(a2, left)) 
	    { *res++ = make_object
		(std::make_pair(left,2u)); }
	}
      std::cout << "</construct_intersection_monotone>" << std::endl;
      return res;
    }
    else {//a1 or a2 are not x_monotone
      std::cout << "<construct_intersection_no_monotone>" << std::endl; 
      std::vector< CGAL::Object > arcs_a1_x_monotone;
      make_x_monotone( a1, std::back_inserter(arcs_a1_x_monotone));
      std::vector< CGAL::Object > arcs_a2_x_monotone;
      make_x_monotone( a2, std::back_inserter(arcs_a2_x_monotone));
      std::vector< Circular_arc_2 > circle_arcs;
      std::vector< Circular_arc_endpoint_2 > circle_arc_endpoints;

      for ( std::vector< CGAL::Object >::iterator it1 = arcs_a1_x_monotone.begin(); 
	    it1 != arcs_a1_x_monotone.end(); ++it1 )
	{
	    //CGAL_kernel_assertion(assign( a1_aux, *it1));
	  const Circular_arc_2 *a1_aux = CGAL::object_cast< Circular_arc_2 >(&*it1);
	  for ( std::vector< CGAL::Object >::iterator it2 = arcs_a2_x_monotone.begin(); 
	    it2 != arcs_a2_x_monotone.end(); ++it2 )
	    {
	      //CGAL_kernel_assertion(assign( a2_aux, *it2));
	      //assign( a2_aux, *it2);
	      const Circular_arc_2 *a2_aux = CGAL::object_cast<Circular_arc_2>(&*it2);
	      std::vector< CGAL::Object > res_aux;
	      construct_intersections_2<CK>( *a1_aux, *a2_aux,
					     std::back_inserter(res_aux));
	      if(res_aux.size() == 2){
		//it can't be a circular_arc_2
		//CGAL_kernel_assertion(assign(the_pair, res_aux[0]));
		const std::pair<Circular_arc_endpoint_2, std::size_t> *the_pair1 = CGAL::object_cast<std::pair<Circular_arc_endpoint_2, std::size_t> >(&res_aux[0]);
		Circular_arc_endpoint_2 arc_end1 = the_pair1->first;
		//assign(the_pair, res_aux[1]);
		const std::pair<Circular_arc_endpoint_2, std::size_t> *the_pair2 = CGAL::object_cast<std::pair<Circular_arc_endpoint_2, std::size_t> >(&res_aux[1]);
		Circular_arc_endpoint_2 arc_end2 = the_pair2->first;
		bool exist = false;
		for (typename std::vector< Circular_arc_endpoint_2 >::iterator it 
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
		for ( typename std::vector< Circular_arc_endpoint_2 >::iterator it
			= circle_arc_endpoints.begin(); 
		      it != circle_arc_endpoints.end(); ++it )
		  {
		    if (arc_end2 == *it) {
		      exist = true;
		      break;
		    }
		  }
		if (!exist) {
		  circle_arc_endpoints.push_back(arc_end2);
		}
	      }
	      else if( res_aux.size() == 1){
		//it can be a Circular_arc_endpoint_2 or a Circular_arc_2
		if(const Circular_arc_2 *arc = CGAL::object_cast<Circular_arc_2>(&res_aux[0])){
		//if(assign(arc,res_aux[0])){
		  circle_arcs.push_back(*arc);
		  
		}
		else{
		  //CGAL_kernel_assertion(assign(the_pair, res_aux[0]));
		  //assign(the_pair, res_aux[0]);
		   const std::pair<Circular_arc_endpoint_2, std::size_t> *the_pair = CGAL::object_cast< std::pair<Circular_arc_endpoint_2, std::size_t> >(&res_aux[0]);
		  Circular_arc_endpoint_2 arc_end = the_pair->first;
		  if (the_pair->second == 2u) {//there are only one tangent point
		    *res++ = res_aux[0];
		    return res;
		  }
		  bool exist = false;
		  for (typename std::vector< Circular_arc_endpoint_2 >::iterator it 
			 = circle_arc_endpoints.begin(); 
		       it != circle_arc_endpoints.end(); ++it )
		    {
		      if (arc_end == *it) {
			exist = true;
			break;
		      }
		    }
		  if (!exist) {
		    circle_arc_endpoints.push_back(arc_end);
		  }
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
	  std::cout << "</construct_intersection_no_monotone>" << std::endl; 
	  return res;
	}
	else{//There are one circle arc and there can be maximum one endpoint
	   for (typename std::vector< Circular_arc_endpoint_2 >::iterator it1 
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
	   std::cout << "</construct_intersection_no_monotone>" << std::endl; 
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
	std::cout << "</construct_intersection_no_monotone>" << std::endl;
	return res;
      }
    }
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

    CGAL_kernel_precondition(a1.is_x_monotone());
    CGAL_kernel_precondition(a2.is_x_monotone());

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
    if (has_on<CK>(a1, left) &&
        has_on<CK>(a2, left) &&
        compare_xy<CK>(left, pt) > 0) {
      p1 = p2 = left;
      return true;
    }
    if (has_on<CK>(a1, right) &&
        has_on<CK>(a2, right) &&
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
    CGAL_kernel_precondition(A.supporting_circle().squared_radius() != 0);
    std::cout << "<make_x_monotone>" << std::endl; 
    int cmp_begin = CGAL::compare(A.source().y(), A.center().y());
    int cmp_end   = CGAL::compare(A.target().y(),   A.center().y());

    int cmp_x = compare_x(A.source(), A.target());

    // We don't need to split
    if (cmp_begin != opposite(cmp_end) &&
        (((cmp_begin > 0 || cmp_end > 0) && cmp_x > 0) ||
          (cmp_begin < 0 || cmp_end < 0) && cmp_x < 0) ) {
      *res++ = make_object(A);
      
      std::cout << "</make_x_monotone>" << std::endl; 
      return res; 
    }

    // Half circles
    if (cmp_begin == 0 && cmp_end == 0 && cmp_x != 0) {
      *res++ = make_object(A);
      std::cout << "</make_x_monotone>" << std::endl; 
      return res; 
    }

    // We need to split
    CGAL_kernel_assertion(!A.is_x_monotone());

    // Define a circle intersecting the supporting circle of A
    // in the 2 vertical tangent points.
    Circle_2 c (Point_2(A.center().x(), A.center().y()-1),
                A.squared_radius()+1);

    // Define the 2 Circular_arc_endpoints 
    // in the 2 vertical tangent points
    Circular_arc_2 half_circle( A.supporting_circle(),
				c, true,
				c, false);
    
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
    std::cout << "</make_x_monotone>" << std::endl; 
    return res;
  }

} // namespace CircularFunctors 
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_PREDICATES_ON_CIRCULAR_ARC_2_H
