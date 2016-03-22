// Copyright (c) 2016  INRIA Nancy - Grand Est (France).
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
// $URL: $
// $Id:  $
//
//
// Author(s)     : Iordan Iordanov


#ifndef CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>

#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Triangulation_hyperbolic_traits_2.h>
#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"

namespace CGAL {

template< class R >
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 {

public:

  typedef R                                               Kernel;

  typedef Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<R>
                                                          Self;

	typedef Triangulation_hyperbolic_traits_2<R> 				    Base;	                                                          

	// Basic types
	typedef typename Base::RT 									            RT;
	typedef typename Base::FT 									            FT;
	typedef typename Base::Point_2     					            Point_2;
	typedef Point_2 											                  Point; 		// Compatibility issues
  typedef typename Base::Vector_2    					            Vector_2;
  typedef typename Base::Triangle_2  					            Triangle_2;
  typedef typename Base::Line_2      					            Line_2;
  typedef typename Base::Ray_2       					            Ray_2;  
  typedef typename Base::Vector_3    					            Vector_3;
  typedef typename Base::Point_3     					            Point_3;
  typedef typename Base::Angle_2                          Angle_2;
  typedef typename Base::Iso_rectangle_2 			            Iso_rectangle_2;
  typedef typename Base::Circle_2 						            Circle_2;
  typedef typename Base::Arc_2 								            Arc_2;
  typedef typename Base::Line_segment_2 			            Line_segment_2;
  typedef typename Base::Segment_2 						            Segment_2;
  typedef typename Base::Euclidean_line_2 		            Euclidean_line_2;

	// Predicates and comparisons
 	typedef typename Base::Less_x_2                   			Less_x_2;
 	typedef typename Base::Less_y_2                   			Less_y_2;
 	typedef typename Base::Compare_x_2                			Compare_x_2;
 	typedef typename Base::Compare_y_2                			Compare_y_2;
 	typedef typename Base::Orientation_2              			Orientation_2;
 	typedef typename Base::Side_of_oriented_circle_2  			Side_of_oriented_circle_2;
 	typedef typename Base::Compare_distance_2         			Compare_distance_2;
 	typedef typename Base::Compute_squared_distance_2    		Compute_squared_distance_2;

 	// Constructions
 	typedef typename Base::Construct_bisector_2       			Construct_bisector_2;
 	typedef typename Base::Construct_hyperbolic_bisector_2 	Construct_hyperbolic_bisector_2;
 	typedef typename Base::Construct_triangle_2       			Construct_triangle_2;
 	typedef typename Base::Construct_direction_2      			Construct_direction_2;
 	typedef typename Base::Construct_midpoint_2          		Construct_midpoint_2;
 	typedef typename Base::Construct_circumcenter_2 			  Construct_circumcenter_2;
 	typedef typename Base::Construct_segment_2 					    Construct_segment_2;
 	typedef typename Base::Construct_ray_2 						      Construct_ray_2;
 	typedef typename Base::Is_hyperbolic 						        Is_hyperbolic;

private:
	Circle_2 _domain;

public:
	Periodic_4_hyperbolic_Delaunay_triangulation_traits_2() : 
  		_domain(Point_2(0, 0), 1*1)
  	{}
  
  	Periodic_4_hyperbolic_Delaunay_triangulation_traits_2(FT r) : 
  		_domain(Point_2(0, 0), r*r)
  	{}    
  
  	Periodic_4_hyperbolic_Delaunay_triangulation_traits_2(const Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 & other) : 
  		_domain(other._domain)
  	{}
  
  	Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 &operator=(const Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 &)
  	{
    	return *this;
  	}


  	Construct_circumcenter_2
  	construct_circumcenter_2_object() const {
  		return Construct_circumcenter_2(_domain);
  	}

  	Construct_segment_2
  	construct_segment_2_object() const {
  		return Construct_segment_2(_domain);
  	}

  	Construct_midpoint_2
  	construct_midpoint_2_object() const {
  		return Construct_midpoint_2(&_domain);
  	}

  	Less_x_2
  	less_x_2_object() const { 
  		return Less_x_2();
  	}
  
  	Less_y_2
  	less_y_2_object() const { 
  		return Less_y_2();
  	}
  
  	Compare_x_2
  	compare_x_2_object() const { 
  		return Compare_x_2();
  	}
  
  	Compare_y_2
  	compare_y_2_object() const { 
  		return Compare_y_2();
  	}
  
  	Orientation_2
  	orientation_2_object() const { 
  		return Orientation_2();
  	}
  
  	Side_of_oriented_circle_2
  	side_of_oriented_circle_2_object() const {
  		return Side_of_oriented_circle_2();
  	}

  	Construct_hyperbolic_bisector_2
  	construct_hyperbolic_bisector_2_object() const { 
  		return Construct_hyperbolic_bisector_2(_domain);
  	}
  
  	Construct_bisector_2
  	construct_bisector_2_object() const {
  		return Construct_bisector_2();
  	}
  
  	Compare_distance_2
  	compare_distance_2_object() const {
  		return Compare_distance_2();
  	}
  
  	Construct_triangle_2  
  	construct_triangle_2_object() const {
  		return Construct_triangle_2();
  	}
  
  	Construct_direction_2  
  	construct_direction_2_object() const {
  		return Construct_direction_2();
  	}

  	Construct_ray_2  
  	construct_ray_2_object() const {
  		return Construct_ray_2(_domain);
  	}

  	Is_hyperbolic 
  	Is_hyperbolic_object() const { 
  		return Is_hyperbolic(); 
  	}


}; // class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2

} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H








