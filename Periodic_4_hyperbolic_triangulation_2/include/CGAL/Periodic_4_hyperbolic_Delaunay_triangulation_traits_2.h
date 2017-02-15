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
#include <CGAL/Ray_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Hyperbolic_octagon_word_4.h>
#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"


namespace CGAL {


template < class K, class Predicate_ >
class Hyperbolic_traits_with_offsets_2_adaptor
{
	typedef K 				Kernel;
	typedef Predicate_ 		Predicate;

	typedef typename Kernel::Point_2                  			Point;
	typedef Hyperbolic_octagon_word_4<unsigned short int, K>  	Offset;

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




template< class R >
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 : public R {

typedef Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<R>  Self;  

public:

	typedef typename R::FT          						FT;

	typedef typename R::Point_2     						Point_2;
	typedef Point_2                 						Point;
	typedef typename R::Circle_2    						Circle_2;
	typedef typename R::Line_2      						Euclidean_line_2;
	typedef boost::variant<Circle_2,Euclidean_line_2>    	Euclidean_circle_or_line_2; 

	typedef typename R::Circular_arc_2         				Circular_arc_2;
	typedef typename R::Line_arc_2             				Line_arc_2; 
	typedef typename R::Circular_arc_point_2   				Circular_arc_point_2;
	typedef typename R::Segment_2                       	Euclidean_segment_2; //only used internally here
	typedef boost::variant<Circular_arc_2, Line_arc_2>  	Hyperbolic_segment_2;
	

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

	class Construct_hyperbolic_segment_2
	{
		typedef typename CGAL::Regular_triangulation_euclidean_traits_2<R> 				Regular_geometric_traits_2;
		typedef typename Regular_geometric_traits_2::Construct_weighted_circumcenter_2 	Construct_weighted_circumcenter_2;
		typedef typename Regular_geometric_traits_2::Weighted_point_2 					Weighted_point_2;
		typedef typename Regular_geometric_traits_2::Bare_point 						Bare_point;

	public:

		typedef Segment_2 result_type;

		Construct_hyperbolic_segment_2() 
			{}
		
		Hyperbolic_segment_2 operator()(const Point_2& p, const Point_2& q) const
		{
			Origin o;
			if(Euclidean_collinear_2()(p, q, Point_2(o))){
				return Euclidean_segment_2(p, q);
			}
			
			Weighted_point_2 wp(p);
			Weighted_point_2 wq(q);
			Weighted_point_2 wo(Point_2(o), FT(1)); // Poincaré circle 
			
			Bare_point center = Construct_weighted_circumcenter_2()(wp, wo, wq);
			FT sq_radius = Compute_squared_Euclidean_distance_2()(p, center);
			
			Circle_2 circle(center, sq_radius);
			// uncomment!!!
			//assert(circle.has_on_boundary(p) && circle.has_on_boundary(q));
			
			if(Orientation_2()(p, q, center) == LEFT_TURN) {
				return Circular_arc_2(circle, p, q);
			}
			return Circular_arc_2(circle, q, p);
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
	


	class Construct_hyperbolic_circumcenter_2
	{
	public:
		
		Circular_arc_point_2 operator()(Point_2 p, Point_2 q, Point_2 r)
		{ 
			Origin o; 
			Point_2 po = Point_2(o);
			Circle_2 l_inf(po, FT(1));
		 
			Euclidean_circle_or_line_2 bis_pq = Construct_circle_or_line_supporting_bisector()(p,q);
			Euclidean_circle_or_line_2 bis_qr = Construct_circle_or_line_supporting_bisector()(q,r);

			if ( Compare_distance_2()(po,p,q) == EQUAL &&
		 Compare_distance_2()(po,p,r) == EQUAL ) 
	return po; 
			// now supporting objects cannot both be Euclidean lines

			std::pair<Circular_arc_point_2, unsigned> pair;

			Euclidean_line_2* l;
			Circle_2* c;

			if ( Circle_2* c_pq = boost::get<Circle_2>(&bis_pq) )
	{
		if ( Circle_2* c_qr = boost::get<Circle_2>(&bis_qr) )
			{
				typedef typename CK2_Intersection_traits<R, Circle_2, Circle_2>::type Intersection_result; 
				std::vector< Intersection_result > inters;
				intersection(*c_pq, *c_qr, std::back_inserter(inters));
				
				CGAL_triangulation_assertion(assign(pair,inters[0]));
				if ( pair.second == 1 )
		{ 
			if ( Has_on_bounded_side_2()( l_inf, pair.first ) )
				return pair.first;
			
			CGAL_triangulation_assertion(assign(pair,inters[1]));
			return pair.first;
		}
				return pair.first;
			}

		// here bis_qr is a line
		l = boost::get<Euclidean_line_2>(&bis_qr);
		c = c_pq;
	}

			// here bis_pq is a line
			l = boost::get<Euclidean_line_2>(&bis_pq);
			c = boost::get<Circle_2>(&bis_qr);

			typedef typename CK2_Intersection_traits<R, Euclidean_line_2, Circle_2>::type Intersection_result; 
			std::vector< Intersection_result > inters;
			intersection(*l, *c, std::back_inserter(inters));

			CGAL_triangulation_assertion(assign(pair,inters[0]));
			if ( pair.second == 1 )
	{ 
		if ( Has_on_bounded_side_2()( l_inf, pair.first ) )
			return pair.first;
		
		CGAL_triangulation_assertion(assign(pair,inters[1]));
		return pair.first;
	}
			return pair.first;
		}

	}; // end Construct_hyperbolic_circumcenter_2
	
	
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
	
	Construct_hyperbolic_circumcenter_2
		construct_hyperbolic_circumcenter_2_object() const
	{ return Construct_hyperbolic_circumcenter_2(); }
	
	class Construct_hyperbolic_bisector_2
	{    
	public:      
		Construct_hyperbolic_bisector_2() 
			{}
		
		// constructs a hyperbolic line 
		Hyperbolic_segment_2 operator()(Point_2 p, Point_2 q) const
		{
			Origin o; 
			Point_2 po = Point_2(o);
			Circle_2 l_inf = Circle_2(po,FT(1));
			
			if ( Compare_distance_2()(po,p,q) == EQUAL ){      
				Euclidean_line_2 l = Construct_Euclidean_bisector_2()(p,q);

	if ( Less_y_2()(p,q) ) 
		{ return Line_arc_2( l, l_inf, false, l_inf, true ); }
	return Line_arc_2( l, l_inf, true, l_inf, false );
			}
			
			Euclidean_circle_or_line_2 
	bis_pq = Construct_circle_or_line_supporting_bisector()(p, q);
			Circle_2* c = boost::get<Circle_2>(&bis_pq);
			
			if ( Less_y_2()(po,c->center()) )
	{ return Circular_arc_2( *c, l_inf, true, l_inf, false ); }
			else if ( Less_y_2()(c->center(),po) )
	{ return Circular_arc_2( *c, l_inf, false, l_inf, true ); }
			// the center of the circle is on the x-axis
			if ( Less_x_2()(po,c->center()) )
	{ return Circular_arc_2( *c, l_inf, true, l_inf, false ); }
			return Circular_arc_2( *c, l_inf, false, l_inf, true);
		}

		// constructs the hyperbolic bisector of segment [p,q] limited by 
		// circumcenter(p,q,r) on one side
		// and circumcenter(p,s,q) on the other side
		Hyperbolic_segment_2 
			operator()(Point_2 p, Point_2 q, Point_2 r, Point_2 s)
		{
			CGAL_triangulation_precondition
	( (Orientation_2()(p,q,r) == ON_POSITIVE_SIDE) 
		&& (Orientation_2()(p,s,q) == ON_POSITIVE_SIDE) );
			CGAL_triangulation_precondition
	( (Side_of_oriented_circle_2()(p,q,r,s) == ON_NEGATIVE_SIDE) 
		&& (Side_of_oriented_circle_2()(p,s,q,r) == ON_NEGATIVE_SIDE) );

			Origin o; 
			Point_2 po = Point_2(o);

			// TODO MT this is non-optimal... 
			// the bisector is already computed here
			// and it will be recomputed below
			Circular_arc_point_2 a =  Construct_hyperbolic_circumcenter_2()(p,q,r);
			Circular_arc_point_2 b =  Construct_hyperbolic_circumcenter_2()(p,s,q);

			if ( Compare_distance_2()(po, p, q) == EQUAL ){      
				Euclidean_line_2 l = Construct_Euclidean_bisector_2()(p,q);
	return Line_arc_2(l,a,b);
			}

			Euclidean_circle_or_line_2 
	bis_pq = Construct_circle_or_line_supporting_bisector()(p, q);
			Circle_2* c_pq = boost::get<Circle_2>(&bis_pq);

			if ( Compare_distance_2()(po,p,q) == POSITIVE )
	// then p is inside the supporting circle
	{ return Circular_arc_2(*c_pq,b,a);}
			return Circular_arc_2(*c_pq,a,b);
		} 

		// constructs the hyperbolic bisector of segment [p,q] 
		// limited by hyperbolic circumcenter(p,q,r) on one side
		// and going to the infinite line on the other side
		Hyperbolic_segment_2 operator()(Point_2 p, Point_2 q, Point_2 r)
		{
			CGAL_triangulation_precondition
	( Orientation_2()(p,q,r) == POSITIVE );

			Origin o; 
			Point_2 po = Point_2(o);
			Circle_2 l_inf(po, FT(1)); 
		
			// TODO MT this is non-optimal... 
			// the bisector is computed (at least) twice
			Circular_arc_point_2 a =  Construct_hyperbolic_circumcenter_2()(p,q,r);

			if ( Compare_distance_2()(po, p, q) == EQUAL ){      
				Euclidean_line_2 bis_pq = Construct_Euclidean_bisector_2()(p,q);
	typedef typename 
		CK2_Intersection_traits<R, Euclidean_line_2, Circle_2>::type 
		Intersection_result; 
	std::vector< Intersection_result > inters;
	intersection(bis_pq, l_inf, std::back_inserter(inters));
	std::pair<Circular_arc_point_2, unsigned> pair;

	CGAL_triangulation_assertion(assign(pair,inters[0]));
	CGAL_triangulation_assertion(pair.second == 1);
	if ( Less_y_2()(p,q) )
		return Line_arc_2(bis_pq,a,pair.first);
	
	CGAL_triangulation_assertion(assign(pair,inters[1]));
	CGAL_triangulation_assertion(pair.second == 1);
	return Line_arc_2(bis_pq,a,pair.first);
			}

			Euclidean_circle_or_line_2 
	bis_pq = Construct_circle_or_line_supporting_bisector()(p, q);
			Circle_2* c_pq = boost::get<Circle_2>(&bis_pq);
			
			Point_2 approx_a(to_double(a.x()),to_double(a.y()));

			typedef typename 
	CK2_Intersection_traits<R, Circle_2, Circle_2>::type 
	Intersection_result; 
			std::vector< Intersection_result > inters;
			intersection(*c_pq, l_inf, std::back_inserter(inters));
			std::pair<Circular_arc_point_2, unsigned> pair;

			CGAL_triangulation_assertion(assign(pair,inters[0]));
			CGAL_triangulation_assertion(pair.second == 1);
			
			Point_2 approx_pinf(to_double(pair.first.x()), to_double(pair.first.y()));
			Point_2 approx_c(to_double(c_pq->center().x()),
					 to_double(c_pq->center().y()));
			 if ( Orientation_2()(p,q,approx_pinf) == NEGATIVE ) {
	if ( Orientation_2()(approx_c,approx_a,approx_pinf) == POSITIVE )
		return Circular_arc_2( *c_pq, a, pair.first );
	return Circular_arc_2( *c_pq, pair.first, a);
			}

			CGAL_triangulation_assertion(assign(pair,inters[1]));
			if ( Orientation_2()(approx_c,approx_a,approx_pinf) == POSITIVE )
	return Circular_arc_2( *c_pq, pair.first, a);
			return Circular_arc_2( *c_pq, a, pair.first);
		}
	}; // end Construct_hyperbolic_bisector_2
	
	Construct_hyperbolic_bisector_2
	construct_hyperbolic_bisector_2_object() const
	{ return Construct_hyperbolic_bisector_2(); }
	
	Construct_Euclidean_bisector_2
	construct_Euclidean_bisector_2_object() const
	{ return Construct_Euclidean_bisector_2(); }	
	

	// do not document
	// constructs the Euclidean circle or line supporting the hyperbolic
	// bisector of two points  
	class Construct_circle_or_line_supporting_bisector
	{
	public:
		Construct_circle_or_line_supporting_bisector()
			{}

		Euclidean_circle_or_line_2 
			operator()(Point_2 p, Point_2 q) const
		{
			Origin o; 
			Point_2 po = Point_2(o);
			typedef typename R::Point_3     Point_3;
		
			if ( Compare_distance_2()(po,p,q) == EQUAL )
	{ return Construct_Euclidean_bisector_2()(p,q); }

			FT dop2 = p.x()*p.x() + p.y()*p.y();
			FT doq2 = q.x()*q.x() + q.y()*q.y();
			Point_3 p3( p.x(), p.y(), dop2 );
			Point_3 q3( q.x(), q.y(), doq2 );
			
			// TODO MT improve 
			
			// The cirle belongs to the pencil with limit points p and q
			// p, q are zero-circles
			// (x, y, xˆ2 + yˆ2 - rˆ2) = alpha*(xp, yp, xpˆ2 + ypˆ2) + (1-alpha)*(xq, yq, xqˆ2 + yqˆ2)
			// xˆ2 + yˆ2 - rˆ2 = 1 (= radius of the Poincare disc)
			FT op = p.x()*p.x() + p.y()*p.y();
			FT oq = q.x()*q.x() + q.y()*q.y();
			FT alpha = (FT(1) - oq) / (op - oq); 
			
			FT x = alpha*p.x() + (1-alpha)*q.x();
			FT y = alpha*p.y() + (1-alpha)*q.y();
			FT sq_radius = x*x + y*y - FT(1);
			
			// TODO MT improve 
			// ?? orientation should depend on
			// Compare_distance(O,p,q)
			// so that p always on positive side
			// ??? 
			// CK does not care about orientation, circular arcs are
			// considered in CCW order in any case

			Euclidean_line_2 l = Construct_Euclidean_bisector_2()(p, q);
			Point_2 middle = Construct_Euclidean_midpoint_2()(p, q);
			Point_2 temp = middle + l.to_vector();

			if (Orientation_2()(middle, temp, Point_2(x, y)) == ON_POSITIVE_SIDE)
	{ return Circle_2(Point_2(x, y), sq_radius, CLOCKWISE); }
			return Circle_2(Point_2(x, y), sq_radius, COUNTERCLOCKWISE);
		}
	}; // end Construct_supporting_circle_of_bisector




	/****************************************************/
	class Side_of_hyperbolic_face_2 {

		typedef Hyperbolic_octagon_word_4<unsigned short int, R>  	Offset;
		
	public:
		typedef Bounded_side result_type;

		Side_of_hyperbolic_face_2() 
			{}



		template<class Face_handle, class Offset>
		Bounded_side operator()(const Point_2 p, const Face_handle fh, const Offset o) const {

			//cout << "Checking face " << fh->get_number() << " with offset " << o << endl;

			Point_2 p1 = o.append(fh->offset(0)).apply(fh->vertex(0)->point());
			Point_2 p2 = o.append(fh->offset(1)).apply(fh->vertex(1)->point());
			Point_2 p3 = o.append(fh->offset(2)).apply(fh->vertex(2)->point());
		
			// Construct_segment_2 bld;
			// Segment_2 s1, s2, s3;
			
			// s1 = bld(p2, p3);
			// s2 = bld(p3, p1);
			// s3 = bld(p1, p2);

			Bounded_side cs1 = side_of_segment_2(p1, p2, p3);
			Bounded_side cp1 = side_of_segment_2(p,  p2, p3);

			Bounded_side cs2 = side_of_segment_2(p2, p3, p1);
			Bounded_side cp2 = side_of_segment_2(p,  p3, p1);

			Bounded_side cs3 = side_of_segment_2(p3, p1, p2);
			Bounded_side cp3 = side_of_segment_2(p,  p1, p2);

			if (cs1 != cp1 || cs2 != cp2 || cs3 != cp3) {
				return ON_UNBOUNDED_SIDE;
			} else {
				if (cp1 == ON_BOUNDARY || cp2 == ON_BOUNDARY || cp3 == ON_BOUNDARY) {
					return ON_BOUNDARY;
				} else {
					return ON_BOUNDED_SIDE;
				}
			}

		}



		template<class Face_handle>
		Bounded_side operator()(const Point_2 p, const Face_handle fh) const {
			Point_2 p1 = fh->offset(0).apply(fh->vertex(0)->point());
			Point_2 p2 = fh->offset(2).apply(fh->vertex(1)->point());
			Point_2 p3 = fh->offset(1).apply(fh->vertex(2)->point());
			
			// Construct_segment_2 bld;
			// Segment_2 s1, s2, s3;

			// s1 = bld(p2, p3);
			// s2 = bld(p3, p1);
			// s3 = bld(p1, p2);

			Bounded_side cs1 = side_of_segment_2(p1, p2, p3);
			Bounded_side cp1 = side_of_segment_2(p,  p2, p3);

			Bounded_side cs2 = side_of_segment_2(p2, p3, p1);
			Bounded_side cp2 = side_of_segment_2(p,  p3, p1);

			Bounded_side cs3 = side_of_segment_2(p3, p1, p2);
			Bounded_side cp3 = side_of_segment_2(p,  p1, p2);

			if (cs1 != cp1 || cs2 != cp2 || cs3 != cp3) {
				return ON_UNBOUNDED_SIDE;
			} else {
				if (cp1 == ON_BOUNDARY || cp2 == ON_BOUNDARY || cp3 == ON_BOUNDARY) {
					return ON_BOUNDARY;
				} else {
					return ON_BOUNDED_SIDE;
				}
			}

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








