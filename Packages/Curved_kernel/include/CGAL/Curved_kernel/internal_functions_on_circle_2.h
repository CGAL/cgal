#ifndef CGAL_CURVED_KERNEL_FUNCTIONS_ON_CIRCLE_2_H
#define CGAL_CURVED_KERNEL_FUNCTIONS_ON_CIRCLE_2_H

namespace CGAL {
namespace CircularFunctors {

  template < class CK >
  typename CK::Polynomial_for_circles_2_2
  get_equation( const typename CK::Circle_2 & c )
  {
    typedef typename CK::RT RT;
 
    typedef typename CK::Algebraic_kernel   AK;

    return AK().construct_polynomial_circle_2_2_object()
      ( c.center().x(), c.center().y(), c.squared_radius() );
  }
  
  template < class CK >
  typename CK::Circle_2  
  construct_circle_2( const typename CK::Polynomial_for_circles_2_2 &eq )
  {
    return typename CK::Circle_2( typename CK::Point_2(eq.a(), eq.b()), eq.r_sq() ); 
  }
  
  template< class CK, class OutputIterator>
  OutputIterator
  construct_intersections_2( const typename CK::Circle_2 & c1,
			     const typename CK::Circle_2 & c2,
			     OutputIterator res )
  {
    typedef typename CK::Polynomial_for_circles_2_2 Equation; 

    Equation e1 = CGAL::get_equation<CK>(c1);
    Equation e2 = CGAL::get_equation<CK>(c2);
    
    typedef std::vector
      < std::pair 
          < std::pair<typename CK::Root_of_2,typename CK::Root_of_2>, 
          unsigned > > 
      solutions_container;
    solutions_container solutions;
    
    CGAL::solve<typename CK::Algebraic_kernel>
      ( e1,e2, std::back_inserter(solutions) ); // to be optimized
    
    typedef typename CK::Circular_arc_endpoint_2::Numeric_point_2
                                                    Numeric_point_2;

    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
    bool first = true;

    for ( typename solutions_container::iterator it = solutions.begin(); 
	    it != solutions.end(); ++it )
      {
	*res++ = make_object
	  (std::make_pair(Circular_arc_endpoint_2(c1, c2, first, Numeric_point_2( it->first.first, it->first.second ))
			  , it->second ));
	first = false;
      }

    return res;
  }

} // namespace CircularFunctors
} // namespace CGAL
#endif // CGAL_CURVED_KERNEL_FUNCTIONS_ON_CIRCLE_2_H
