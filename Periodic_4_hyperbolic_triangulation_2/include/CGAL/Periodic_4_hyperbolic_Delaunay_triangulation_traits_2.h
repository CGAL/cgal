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
// Author(s)     : 	Iordan Iordanov
// 					Monique Teillaud
//


#ifndef CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/Hyperbolic_octagon_word_4.h>
#include <CGAL/exact_complex.h>
#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"


namespace CGAL {


template < class K, class Predicate_ >
class Hyperbolic_traits_with_offsets_2_adaptor
{
	typedef K 				Kernel;
	typedef Predicate_ 		Predicate;

	typedef typename Kernel::Point_2    	Point;
	typedef typename Kernel::FT 			FT;
	typedef typename Kernel::Offset  	 	Offset;

	// Use the construct_point_2 predicate from the kernel to convert the periodic points to Euclidean points
	typedef typename Kernel::Construct_point_2        			Construct_point_2;

public:
	typedef typename Predicate::result_type           			result_type;


	Hyperbolic_traits_with_offsets_2_adaptor() { }

	result_type operator()(	const Point& p0, 	const Point& p1,
							const Offset& o0, 	const Offset& o1) const
	{
		return Predicate()(pp(p0, o0), pp(p1, o1));
	}
	result_type operator()(	const Point& p0, 	const Point& p1, 	const Point& p2,
							const Offset& o0, 	const Offset& o1, 	const Offset& o2) const
	{
		return Predicate()(pp(p0, o0), pp(p1, o1), pp(p2, o2));
	}
	result_type operator()(	const Point& p0, 	const Point& p1,
							const Point& p2, 	const Point& p3,
							const Offset& o0, 	const Offset& o1,
							const Offset& o2, 	const Offset& o3) const
	{
		return Predicate()(pp(p0, o0), pp(p1, o1), pp(p2, o2), pp(p3, o3));
	}

	result_type operator()(	const Point& p0, 	const Point& p1) const
	{
		return Predicate()(p0, p1);
	}
	result_type operator()(	const Point& p0, 	const Point& p1,
							const Point& p2) 	const
	{
		return Predicate()(p0, p1, p2);
	}
	result_type operator()(	const Point& p0, 	const Point& p1,
							const Point& p2, 	const Point& p3) const
	{
		return Predicate()(p0, p1, p2, p3);
	}

private:
	Point pp(const Point &p, const Offset &o) const
	{
		return o.apply(p);
	}

};



template < typename K, typename Construct_point_base>
class Periodic_4_hyperbolic_construct_point_2 : public Construct_point_base
{
	typedef K Kernel;

public:
	typedef typename Kernel::Point_2         Point;
	typedef typename Kernel::Offset          Offset;

	typedef Point                            result_type;

	Periodic_4_hyperbolic_construct_point_2() { }

	Point operator() ( const Point& p, const Offset& o ) const
	{
		return o.apply(p);
	}

};



template <class R>
class Simple_circular_arc_2 {
	typedef typename R::FT 			FT;
	typedef typename R::Point_2 	Point;
	typedef typename R::Circle_2 	Circle;

private:
	Circle _c;
	Point _s, _t;

public:
	Simple_circular_arc_2() :
		_c(Point(FT(0),FT(0)), FT(0)), _s(FT(0),FT(0)), _t(FT(0),FT(0)) {}
	
	Simple_circular_arc_2(Circle c, Point source, Point target) :
		_c(c), _s(source), _t(target) {}

	Circle circle() const {
		return _c;
	}

	Point source() const {
		return _s;
	}

	Point target() const {
		return _t;
	}

	FT squared_radius() const {
		return _c.squared_radius();
	}

	Point center() const {
		return _c.center();
	}

	Bbox_2 bbox(void) const {
    	return typename R::Construct_bbox_2()(*this);
  	}

};



template< class R >
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 : public R {

typedef Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<R>  Self;  

public:

	typedef typename R::FT          								FT;
	typedef Hyperbolic_octagon_word_4<unsigned short int, FT>		Offset;
	typedef typename R::Point_2     								Point_2;
	typedef Point_2                 								Point;
	typedef typename R::Circle_2    								Circle_2;
	typedef typename R::Line_2      								Euclidean_line_2;
	typedef boost::variant<Circle_2,Euclidean_line_2>    			Euclidean_circle_or_line_2; 
	typedef Simple_circular_arc_2<R>         						Circular_arc_2;
	typedef typename R::Segment_2                       			Euclidean_segment_2; //only used internally here
	typedef boost::variant<Circular_arc_2, Euclidean_segment_2>  	Hyperbolic_segment_2;
	

	// Wrappers for the offset adapter
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename R::Compare_x_2>                 Compare_x_2;
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename R::Compare_y_2>                 Compare_y_2;
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename R::Orientation_2>               Orientation_2;
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename R::Compare_distance_2>          Compare_distance_2;
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, typename R::Side_of_oriented_circle_2>   Side_of_oriented_circle_2;

	
	

	// only kept for demo to please T2graphicsitems
	typedef Euclidean_segment_2  							Line_segment_2;
	typedef Hyperbolic_segment_2 							Segment_2;

	// the following types are only used internally in this traits class, 
	// so they need not be documented, and they don't need _object()
	typedef typename R::Collinear_2                			Euclidean_collinear_2;
	typedef typename R::Construct_bisector_2       			Construct_Euclidean_bisector_2;
	typedef typename R::Construct_midpoint_2       			Construct_Euclidean_midpoint_2;
	typedef typename R::Compute_squared_distance_2 			Compute_squared_Euclidean_distance_2;
	typedef typename R::Has_on_bounded_side_2 				Has_on_bounded_side_2;

	typedef typename R::Less_x_2                   			Less_x_2;
	typedef typename R::Less_y_2                   			Less_y_2;
			
public:

	class Construct_hyperbolic_segment_2 {
		typedef exact_complex<FT> 	cplx;
	public:
		typedef Segment_2 result_type;

		Construct_hyperbolic_segment_2() {}

		Segment_2 operator()(const Point_2& p1, const Point_2& p2) const {
			cplx p(p1), q(p2);
			cplx O(0,0);
			cplx inv;
			if (p == O) {
				inv = q.invert_in_unit_circle();
			} else {
				inv = p.invert_in_unit_circle();
			}

			Point ip(inv.real(), inv.imag());

			if (Orientation_2()(p1, p2, ip) == COLLINEAR) {
				Euclidean_segment_2 seg(p1, p2);
				return seg;
			} else {
				Circle_2 c(p1, p2, ip);
				if(Orientation_2()(p1, p2, c.center()) == LEFT_TURN) {
					return Circular_arc_2(c, p1, p2);
				}
				return Circular_arc_2(c, p2, p1);
			}

		}
		
	}; // end Construct_hyperbolic_segment_2
	
	Construct_hyperbolic_segment_2
		construct_hyperbolic_segment_2_object() const
	{ return Construct_hyperbolic_segment_2(); }


	// wrong names kept for demo
	typedef Hyperbolic_traits_with_offsets_2_adaptor<Self, Construct_hyperbolic_segment_2> Construct_segment_2;
	Construct_segment_2
	construct_segment_2_object() const
	{ return Construct_segment_2(); }
	
	
	Compare_x_2 
		compare_x_2_object() const 
	{ return Compare_x_2();} 
	
	Compare_y_2 
		compare_y_2_object() const 
	{ return Compare_y_2();} 
	
	Orientation_2
		orientation_2_object() const
	{ return Orientation_2();}
	
	Side_of_oriented_circle_2
		side_of_oriented_circle_2_object() const
	{ return Side_of_oriented_circle_2(); }
	


	class Construct_hyperbolic_circle_2 {

		typedef exact_complex<FT> 	cplx;

	public: 
		Construct_hyperbolic_circle_2() {}

		Circle_2 operator()(Point_2 hcenter, Point_2 p) {
			
			Origin o;

			if (hcenter == o) {
				return Circle_2(o, squared_distance(o, p));
			} else if (Orientation_2()(hcenter, p, o) != COLLINEAR) {
				Euclidean_line_2 ell(hcenter, o);
				
				cplx p1(hcenter), p2(p);
				cplx inv;
				if (p1 == cplx(0,0)) {
					inv = p2.invert_in_unit_circle();
				} else {
					inv = p1.invert_in_unit_circle();
				}
				Point ip(inv.real(), inv.imag());

				Circle_2 schl(hcenter, p, ip);
				Euclidean_line_2 line_through_p(schl.center(), p);
				Euclidean_line_2 tangent_at_p = line_through_p.perpendicular(p);

				// assume that ell := ax + by + c = 0 and tangent_at_p := dx + ey + f = 0
				FT a = ell.a(), b = ell.b(), c = ell.c();
				FT d = tangent_at_p.a(), e = tangent_at_p.b(), f = tangent_at_p.c();
				FT py = (c*d - a*f)/(a*e - d*b);
				FT px = (-c -b*py)/a;
				Point intersection(px, py);
				return Circle_2(intersection, squared_distance(intersection, p));
			} else {  // if the given points and the origin are collinear, we need to treat them differently
				cplx hcinv = cplx(hcenter).invert_in_unit_circle();
				Point ip(hcinv.real(), hcinv.imag());
				Point mp = midpoint(hcenter, ip);
				Circle_2 tmpc(mp, hcenter);
				cplx res = cplx(p).invert_in_circle(tmpc);
				Point pres(res.real(), res.imag());
				return Circle_2(p, pres);
			}
		}
	};


	class Construct_hyperbolic_bisector_2 {

		typedef exact_complex<FT> 	cplx;

	public:      
		Construct_hyperbolic_bisector_2() 
			{}
		
		// constructs a hyperbolic line 
		Hyperbolic_segment_2 operator()(Point_2 p, Point_2 q) const
		{

			Origin o; 
			Point_2 po = Point_2(o);
			if ( Compare_distance_2()(po, p, q) == EQUAL ){      
				return Construct_Euclidean_bisector_2()(p, q);
			}

			Circle_2 c1 = Construct_hyperbolic_circle_2()(p, q);
			Circle_2 c2 = Construct_hyperbolic_circle_2()(q, p);

			typedef typename 
			CK2_Intersection_traits<R, Circle_2, Circle_2>::type Intersection_result; 
			std::vector< Intersection_result > inters;
			intersection(c1, c2, std::back_inserter(inters));
			std::pair<Point_2, unsigned> p1, p2;

			CGAL_triangulation_assertion(assign(p1,inters[0]));
			CGAL_triangulation_assertion(assign(p2,inters[2]));

			cplx cp1(p1), cp2(p2);
			cplx inv;
			if (cp1 == cplx(0,0)) {
				inv = cp2.invert_in_unit_circle();
			} else {
				inv = cp1.invert_in_unit_circle();
			}
			
			Point cpi(inv.real(), inv.imag());

			Circle_2 c(p1, p2, cpi);

			inters.clear();
			intersection(c, Circle_2(Point(0,0), 1), std::back_inserter(inters));			

			CGAL_triangulation_assertion(assign(p1,inters[0]));
			CGAL_triangulation_assertion(assign(p2,inters[2]));

			if(Orientation_2()(p1, p2, c.center()) == LEFT_TURN) {
				return Circular_arc_2(c, p1, p2);
			}
			return Circular_arc_2(c, p2, p1);

		}
	}; // end Construct_hyperbolic_bisector_2
	
	Construct_hyperbolic_bisector_2
	construct_hyperbolic_bisector_2_object() const
	{ return Construct_hyperbolic_bisector_2(); }
	
	Construct_Euclidean_bisector_2
	construct_Euclidean_bisector_2_object() const
	{ return Construct_Euclidean_bisector_2(); }	


	/****************************************************/
	class Side_of_hyperbolic_face_2 {
		
	public:
		typedef Bounded_side result_type;

		Side_of_hyperbolic_face_2() 
			{}



		template<class Face_handle, class Offset>
		Bounded_side operator()(const Point_2 p, Bounded_side sides[3], const Face_handle fh, const Offset o) const {

			Point_2 p1 = o.append(fh->offset(0)).apply(fh->vertex(0)->point());
			Point_2 p2 = o.append(fh->offset(1)).apply(fh->vertex(1)->point());
			Point_2 p3 = o.append(fh->offset(2)).apply(fh->vertex(2)->point());

			Bounded_side cp1 = side_of_segment_2(p,  p2, p3);
			sides[0] = cp1;
			if (cp1 == ON_BOUNDARY) {
				return ON_BOUNDARY;
			}

			Bounded_side cp2 = side_of_segment_2(p,  p3, p1);
			sides[1] = cp2;
			if (cp2 == ON_BOUNDARY) {
				return ON_BOUNDARY;
			}

			Bounded_side cp3 = side_of_segment_2(p,  p1, p2);
			sides[2] = cp3;
			if (cp3 == ON_BOUNDARY) {
				return ON_BOUNDARY;
			}			

			Bounded_side cs1 = side_of_segment_2(p1, p2, p3);
			Bounded_side cs2 = side_of_segment_2(p2, p3, p1);
			Bounded_side cs3 = side_of_segment_2(p3, p1, p2);


			// Cannot be on the boundary here.
			if (cs1 != cp1 || cs2 != cp2 || cs3 != cp3) {
				return ON_UNBOUNDED_SIDE;
			} else {
				return ON_BOUNDED_SIDE;	
			}


		}



		template<class Face_handle>
		Bounded_side operator()(const Point_2 p, Bounded_side sides[3], const Face_handle fh) const {
			return operator()(p, sides, fh, Offset());
		}


	private:
		Bounded_side side_of_segment_2(const Point_2 query, const Point_2 p, const Point_2 q) const {
			
			Point_2 o(0, 0);

			// Invert p or q through the unit circle.
			// The inversion depends on the distance from the origin, so to increase 
			// numerical stability we choose to invert the point further from (0,0).
			Point_2 inv;
			FT dp = squared_distance(o, p), dq = squared_distance(o, q);
			if (dq < dp) {
				inv = Point_2( p.x()/dp, p.y()/dp );
			} else {
				inv = Point_2( q.x()/dq, q.y()/dq );
			}

			// If iq is on the line defined by p and q, we need to work with a line instead of a circle
			if (orientation(p, q, inv) == COLLINEAR) {
				Orientation oquery = orientation(p, q, query);
				if (oquery == COLLINEAR) {
					return ON_BOUNDARY;
				} else if (oquery == LEFT_TURN) {
					return ON_BOUNDED_SIDE; 	// this is just a convention
				} else {
					return ON_UNBOUNDED_SIDE;
				}
			} else { // this means that we work in the circle
				Circle_2 c(p, q, inv);
				return c.bounded_side(query);
			}

		}

	};


	Side_of_hyperbolic_face_2
	side_of_hyperbolic_face_2_object() const {
		return Side_of_hyperbolic_face_2();
	}

	/****************************************************/




public:
		Periodic_4_hyperbolic_Delaunay_triangulation_traits_2() {}
	
		Periodic_4_hyperbolic_Delaunay_triangulation_traits_2(const Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 & other) {}
	
		Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 &operator=(const Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 &)
		{
			return *this;
		}

		class Side_of_fundamental_octagon {
		public:
			Side_of_fundamental_octagon() {}

			CGAL::Bounded_side operator()(Point_2 p) {

				// Rotation by pi/4
				CGAL::Aff_transformation_2<R> rotate(CGAL::ROTATION, std::sqrt(0.5), std::sqrt(0.5));
				
				// The center of the Euclidean circle corresponding to the side s_1 (east)
				Point_2  CenterA ( FT( sqrt((sqrt(2.) + 1.) / 2.) ), FT(0.) );

				// The squared radius of the Eucliden circle corresponding to the side s_1
				FT       Radius2 ( (sqrt(2.) - 1.) / 2. );

				// Poincare disk (i.e., unit Euclidean disk)
				Circle_2 Poincare    ( Point(0, 0),       1*1 );

				// Euclidean circle corresponding to s_1
				Circle_2 EuclidCircA ( CenterA,           Radius2 );

				// Euclidean circle corresponding to s_2 (just rotate the center, radius is the same)
				Circle_2 EuclidCircBb( rotate(CenterA),   Radius2 );

				// This transformation brings the point in the first quadrant (positive x, positive y)
				FT x(FT(p.x()) > FT(0) ? p.x() : -p.x());
				FT y(FT(p.y()) > FT(0) ? p.y() : -p.y());

				// This brings the point in the first octant (positive x and y < x)
				if (y > x) {
					FT tmp = x;
					x = y;
					y = tmp;
				}

				// This tells us whether the point is in the side of the open boundary
				bool on_open_side = ( ( p.y() + tan(CGAL_PI / 8.) * p.x() ) < 0.0 );

				Point t(x, y);

				CGAL::Bounded_side PoincareSide = Poincare.bounded_side(t);
				CGAL::Bounded_side CircASide    = EuclidCircA.bounded_side(t);
				CGAL::Bounded_side CircBbSide   = EuclidCircBb.bounded_side(t);

				// First off, the point needs to be inside the Poincare disk. if not, there's no hope.
				if ( PoincareSide == CGAL::ON_BOUNDED_SIDE ) {
					
					// Inside the Poincare disk, but still outside the fundamental domain
					if ( CircASide  == CGAL::ON_BOUNDED_SIDE || 
							 CircBbSide == CGAL::ON_BOUNDED_SIDE   ) {
						return CGAL::ON_UNBOUNDED_SIDE;
					}

					// Inside the Poincare disk and inside the fundamental domain
					if ( CircASide  == CGAL::ON_UNBOUNDED_SIDE && 
							 CircBbSide == CGAL::ON_UNBOUNDED_SIDE ) {
						return CGAL::ON_BOUNDED_SIDE;
					} 

					// This is boundary, but we only consider the upper half. The lower half means outside.
					if (on_open_side) {
						return CGAL::ON_UNBOUNDED_SIDE;
					} else {
						return CGAL::ON_BOUNDED_SIDE;
					}

				} else {
					return CGAL::ON_UNBOUNDED_SIDE;
				}

			}

		};


		Side_of_fundamental_octagon
		side_of_fundamental_octagon_object() const {
			return Side_of_fundamental_octagon();
		}



}; // class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2

} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H








