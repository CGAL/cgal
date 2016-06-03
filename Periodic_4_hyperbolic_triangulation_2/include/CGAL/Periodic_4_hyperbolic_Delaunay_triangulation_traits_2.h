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
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include "../../../Hyperbolic_triangulation_2/include/CGAL/Triangulation_hyperbolic_traits_2.h"
#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"

namespace CGAL {


template < class K, class Predicate_ >
class Hyperbolic_traits_with_offsets_2_adaptor
{
  typedef K Kernel;
  typedef Predicate_ Predicate;

  typedef typename Kernel::Point_2                  Point;
  typedef Hyperbolic_word_4<unsigned short int>     Offset;

  // Use the construct_point_2 predicate from the kernel to convert the periodic points to Euclidean points
  typedef typename Kernel::Construct_point_2        Construct_point_2;
public:
  typedef typename Kernel::Circle_2                 Circle_2;

public:
  typedef typename Predicate::result_type           result_type;

  Hyperbolic_traits_with_offsets_2_adaptor(const Circle_2 * dom) : _domain(dom) { }

  result_type operator()(const Point& p0, const Point& p1,
                         const Offset& o0, const Offset& o1) const
  {
    return Predicate()(pp(p0, o0), pp(p1, o1));
  }
  result_type operator()(const Point& p0, const Point& p1, const Point& p2,
                         const Offset& o0, const Offset& o1, const Offset& o2) const
  {
    return Predicate()(pp(p0, o0), pp(p1, o1), pp(p2, o2));
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3,
                         const Offset& o0, const Offset& o1,
                         const Offset& o2, const Offset& o3) const
  {
    return Predicate()(pp(p0, o0), pp(p1, o1), pp(p2, o2), pp(p3, o3));
  }

  result_type operator()(const Point& p0, const Point& p1) const
  {
    return Predicate()(p0, p1);
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2) const
  {
    return Predicate()(p0, p1, p2);
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3) const
  {
    return Predicate()(p0, p1, p2, p3);
  }

private:
  Point pp(const Point &p, const Offset &o) const
  {
    return o.apply(p);
  }
public:
  const Circle_2* _domain;
};



template < typename K, typename Construct_point_base>
class Periodic_4_hyperbolic_construct_point_2 : public Construct_point_base
{
  typedef K Kernel;

public:
  typedef typename Kernel::Point_2         Point;
  typedef typename Kernel::Offset          Offset;
  typedef typename Kernel::Circle_2        Circle_2;

  typedef Point                            result_type;

  Periodic_4_hyperbolic_construct_point_2(const Circle_2 & dom) : _dom(dom) { }

  Point operator() ( const Point& p, const Offset& o ) const
  {
    return o.apply(p);
  }

private:
  Circle_2 _dom;
};




template< class R >
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 {

public:

  typedef R                                               Kernel;
  typedef R                                               K;
  typedef typename R::Kernel_base                         Kernel_base;
  typedef typename R::Direction_2                         Direction_2;

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
          

 	// Constructions

  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Less_x_2>                   Less_x_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Less_y_2>                   Less_y_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Compare_x_2>                Compare_x_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Compare_y_2>                Compare_y_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Compare_xy_2>               Compare_xy_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Orientation_2>              Orientation_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Side_of_oriented_circle_2>  Side_of_oriented_circle_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Compute_squared_radius_2>   Compute_squared_radius_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Construct_center_2>         Construct_center_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Construct_circumcenter_2>   Construct_circumcenter_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Construct_bisector_2>       Construct_bisector_2;
  typedef Periodic_4_hyperbolic_construct_point_2<Self,  typename K::Construct_point_2>          Construct_point_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Construct_midpoint_2>       Construct_midpoint_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Compare_distance_2>         Compare_distance_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Construct_segment_2>        Construct_segment_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Construct_triangle_2>       Construct_triangle_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Construct_direction_2>      Construct_direction_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename K::Construct_ray_2>            Construct_ray_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, 
                                              typename Base::Construct_hyperbolic_bisector_2>    Construct_hyperbolic_bisector_2;
  typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename Base::Is_hyperbolic>           Is_hyperbolic;

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
  		return Orientation_2(&_domain);
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
  	is_hyperbolic_object() const { 
  		return Is_hyperbolic(); 
  	}


    class Side_of_fundamental_octagon {
    public:
      Side_of_fundamental_octagon() {}

      CGAL::Bounded_side operator()(Point_2 p) {

        CGAL::Aff_transformation_2<Kernel> rotate(CGAL::ROTATION, std::sqrt(0.5), std::sqrt(0.5));
        Point_2  CenterA ( FT( sqrt((sqrt(2.) + 1.) / 2.) ), FT(0.) );
        FT       Radius2 ( (sqrt(2.) - 1.) / 2. );

        Circle_2 Poincare( Point(0, 0),       1*1 );
        Circle_2 DomainA ( CenterA,           Radius2 );
        Circle_2 DomainBb( rotate(CenterA),   Radius2 );

        FT x(FT(p.x()) > FT(0) ? p.x() : -p.x());
        FT y(FT(p.y()) > FT(0) ? p.y() : -p.y());

        if (y > x) {
          FT tmp = x;
          x = y;
          y = tmp;
        }

        Point t(x, y);

        CGAL::Bounded_side bs1 = Poincare.bounded_side(t);
        CGAL::Bounded_side bs2 = DomainA.bounded_side(t);
        CGAL::Bounded_side bs3 = DomainBb.bounded_side(t);

        if (  (bs1 != CGAL::ON_UNBOUNDED_SIDE) &&
              (bs2 != CGAL::ON_BOUNDED_SIDE)   &&
              (bs3 != CGAL::ON_BOUNDED_SIDE)      ) {
          return CGAL::ON_BOUNDED_SIDE;
        } 
        else if ( (bs1 == CGAL::ON_UNBOUNDED_SIDE) ||
                  (bs2 == CGAL::ON_BOUNDED_SIDE)   ||
                  (bs3 == CGAL::ON_BOUNDED_SIDE)      ) {
          return CGAL::ON_UNBOUNDED_SIDE;
        }
       
        return CGAL::ON_BOUNDARY;

      }

    };


    Side_of_fundamental_octagon
    side_of_fundamental_octagon_object() const {
      return Side_of_fundamental_octagon();
    }











}; // class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2

} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H








