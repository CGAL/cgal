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
#include <CGAL/Exact_complex.h>
#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"
#include <CGAL/Hyperbolic_octagon_translation.h>


#if defined PROFILING_MODE
	#include <CGAL/Timer.h>
	extern long calls_predicate_identity;
	extern long calls_predicate_non_identity;
	extern double time_predicate_identity;
	extern double time_predicate_non_identity;
#endif	


using std::pair;
using std::make_pair;

namespace CGAL {


template < class K, class Predicate_ >
class Hyperbolic_traits_with_translations_2_adaptor
{
	typedef K 									Kernel;
	typedef Predicate_ 							Predicate;
	typedef typename Kernel::Construct_point_2 	CP2;
	typedef typename Kernel::Point_2    		Point;
	typedef typename Kernel::FT 				FT;
	typedef typename Kernel::Hyperbolic_translation  	 	Hyperbolic_translation;

	// Use the construct_point_2 predicate from the kernel to convert the periodic points to Euclidean points
	typedef typename Kernel::Construct_point_2        			Construct_point_2;

public:
	typedef typename Predicate::result_type           			result_type;


	Hyperbolic_traits_with_translations_2_adaptor() { }

	result_type operator()(	const Point& p0, 	const Point& p1,
							const Hyperbolic_translation& o0, 	const Hyperbolic_translation& o1) const
	{
		return Predicate()(pp(p0, o0), pp(p1, o1));
	}
	result_type operator()(	const Point& p0, 	const Point& p1, 	const Point& p2,
							const Hyperbolic_translation& o0, 	const Hyperbolic_translation& o1, 	const Hyperbolic_translation& o2) const
	{
		return Predicate()(pp(p0, o0), pp(p1, o1), pp(p2, o2));
	}
	result_type operator()(	const Point& p0, 	const Point& p1,
							const Point& p2, 	const Point& p3,
							const Hyperbolic_translation& o0, 	const Hyperbolic_translation& o1,
							const Hyperbolic_translation& o2, 	const Hyperbolic_translation& o3) const
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
	Point pp(const Point &p, const Hyperbolic_translation &o) const
	{
		return CP2()(p, o);
	}

};



template < typename K, typename Construct_point_base>
class Periodic_4_hyperbolic_construct_point_2 : public Construct_point_base
{

private:
	typedef K 												Kernel;
	typedef typename Kernel::FT 							NT;
	typedef typename Kernel::Point_2         				Point;
	typedef typename Kernel::Hyperbolic_translation         Hyperbolic_translation;
	
public:
	
	typedef Point                            				result_type;

	Periodic_4_hyperbolic_construct_point_2() { }

	Point operator() ( const Point& p, const Hyperbolic_translation& tr ) const
	{
		std::pair<NT, NT> a, b;
		a = tr.alpha();
		b = tr.beta();

		// Prepare variables
		NT ax(a.first);
		NT bx(b.first);
		NT zx(p.x());
		NT ay(a.second);
		NT by(b.second);
		NT zy(p.y());

		// Compute parts of fraction
		NT rn(ax*zx - ay*zy + bx);  // real part of numerator
		NT in(ay*zx + ax*zy + by);  // imaginary part of numerator
		NT rd(bx*zx + by*zy + ax);	// real part of denominator
		NT id(bx*zy - by*zx - ay);	// imaginary part of denominator

		// The denominator cannot be zero
		NT den(rd*rd + id*id);
		CGAL_assertion(den != NT(0));

		// Compute real and imaginary part of result
		NT rx((rn*rd + in*id)); 
		NT ry((in*rd - rn*id));

		return Point(rx/den, ry/den);
	}

};



template <class R>
class Simple_circular_arc_2 {
	typedef typename R::FT 			  FT;
	typedef Exact_complex<FT> 		  Cplx;
	typedef typename R::Point_2 	  Point;
	typedef typename R::Circle_2 	  Circle;
	typedef typename R::Orientation_2 Orientation_2;

private:
	Circle _c;
	Point _s, _t;

public:
	Simple_circular_arc_2() :
		_c(Point(FT(0),FT(0)), FT(0)), _s(FT(0),FT(0)), _t(FT(0),FT(0)) {}
	
	Simple_circular_arc_2(Circle c, Point source, Point target) :
		_c(c), _s(source), _t(target) {}

	Simple_circular_arc_2(Point p1, Point p2) {
		Cplx p(p1), q(p2);
		Cplx O(0,0);
		Cplx inv;
		if (p == O) {
			inv = q.invert_in_unit_circle();
		} else {
			inv = p.invert_in_unit_circle();
		}

		Point ip(inv.real(), inv.imag());

		_c = Circle(p1, p2, ip);
		if (Orientation_2()(p1, p2, _c.center()) == LEFT_TURN) {
			_s = p1;
			_t = p2;
		} else {
			_s = p2;
			_t = p1;
		}

	}

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



template< class Kernel, template<typename> class Translation_type = Hyperbolic_octagon_translation>
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 : public Kernel {

typedef Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel, Translation_type>  Self;  

public:

	typedef typename Kernel::FT          								 	FT;
	typedef Translation_type<FT> 											Hyperbolic_translation;

	typedef typename Kernel::Point_2     									Point_2;
	typedef Point_2                 										Point;
	typedef typename Kernel::Circle_2    									Circle_2;
	typedef typename Kernel::Line_2      									Euclidean_line_2;
	typedef boost::variant<Circle_2,Euclidean_line_2>    					Euclidean_circle_or_line_2; 
	typedef Simple_circular_arc_2<Kernel>				 					Circular_arc_2;
	typedef typename Kernel::Segment_2                       				Euclidean_segment_2; //only used internally here
	typedef boost::variant<Circular_arc_2, Euclidean_segment_2>  			Hyperbolic_segment_2;

	// Wrappers for the translation adapter
	typedef Hyperbolic_traits_with_translations_2_adaptor<Self, typename Kernel::Orientation_2>               Orientation_2;
	typedef Hyperbolic_traits_with_translations_2_adaptor<Self, typename Kernel::Side_of_oriented_circle_2>   Side_of_oriented_circle_2;

	// only kept for demo to please T2graphicsitems
	typedef Euclidean_segment_2  									Line_segment_2;
	typedef Hyperbolic_segment_2 									Segment_2;

	// the following types are only used internally in this traits class, 
	// so they need not be documented, and they don't need _object()
	typedef typename Kernel::Construct_bisector_2       													Construct_Euclidean_bisector_2;
	typedef Hyperbolic_traits_with_translations_2_adaptor<Self, typename Kernel::Compare_distance_2>        Compare_distance_2;
	typedef Periodic_4_hyperbolic_construct_point_2<Self, typename Kernel::Construct_point_2>         		Construct_point_2;
			
private:

	class Compute_circle_in_pencil {
	public:
		typedef Circle_2 result_type;

		Compute_circle_in_pencil() {}

		// Code by Olivier Devillers (CGAL_ipelets)
		result_type operator()(Circle_2 c, Circle_2 c1, Circle_2 c2) {
			Point_2 origin = ORIGIN;
			FT lambda = squared_distance(c.center(), origin);
			lambda -= c.squared_radius();
			FT l1 = squared_distance(c1.center(),origin) - c1.squared_radius() ;
			FT l2 = squared_distance(c2.center(),origin) - c2.squared_radius() ;
			l1 += -2*((c1.center()-origin)*(c.center()-origin));
			l2 += -2*((c2.center()-origin)*(c.center()-origin));
			
			if (l1==l2){ // degenerate case, radical axis
				return Circle_2();
			}

			lambda= -(lambda+l2)/(l1-l2);
			Point_2 center = origin+lambda*(c1.center()-origin)+(1-lambda)*(c2.center()-origin);
			FT sqradius = - lambda*(CGAL::squared_distance(c1.center(),origin)-c1.squared_radius())
						  - (1-lambda)*(CGAL::squared_distance(c2.center(),origin)-c2.squared_radius())
						  + CGAL::squared_distance(center,origin);
			Circle_2 circ(center,sqradius);
			return circ;
		}
	};


	
	class Compute_circle_orthogonal {
	public:
		typedef Circle_2 result_type;

		Compute_circle_orthogonal() {}

		// Code by Olivier Devillers (CGAL_ipelets)
		result_type operator()(Circle_2 c, Circle_2 c1, Circle_2 c2) {
			Point_2 origin = ORIGIN;
			FT z  = squared_distance(c.center() ,origin) -  c.squared_radius();
			FT z1 = squared_distance(c1.center(),origin) - c1.squared_radius();
			FT z2 = squared_distance(c2.center(),origin) - c2.squared_radius();
			FT det=	-(c1.center().x() * c2.center().y() - c1.center().y() * c2.center().x())
					+(c.center().x() * c2.center().y() - c.center().y() * c2.center().x())
					-(c.center().x() * c1.center().y() - c.center().y() * c1.center().x());
			
			if (det == FT(0)){ // degenerate casse, radical axis
				return Circle_2();
			}
			
			FT x =	( -(z1 * c2.center().y() - c1.center().y() * z2)
					+(z * c2.center().y() - c.center().y() * z2)
					-(z * c1.center().y() - c.center().y() * z1))/FT(2)/det;
			
			FT y =  ( -(c1.center().x() * z2 - z1 * c2.center().x())
					+(c.center().x() * z2 - z * c2.center().x())
					-(c.center().x() * z1 - z * c1.center().x()))/FT(2)/det;
			
			FT rr= -(  (c1.center().x() * c2.center().y() - c1.center().y() * c2.center().x())*z
					-(c.center().x() * c2.center().y() - c.center().y() * c2.center().x())*z1
					+(c.center().x() * c1.center().y() - c.center().y() * c1.center().x())*z2)/det+x*x+y*y;
			
			Point_2 center(x,y);
			Circle_2 circ(center,rr);
			return circ;
		}
	};



public:


	class Compute_hyperbolic_diameter {
	public:

		typedef double result_type;

		Compute_hyperbolic_diameter() {}

		result_type operator()(Circle_2 c) {
		
			typedef Euclidean_line_2       				Line;
			typedef Circle_2     						Circle;
			typedef Construct_inexact_intersection_2 	Intersection;

			Point  p0(0, 0);
			Circle c0(p0, 1);
			Line  ell(p0, c.center());

			if (ell.is_degenerate()) {
				return 5.;
			} 

			pair<Point, Point> res1 = Intersection()(c0, ell);
			pair<Point, Point> res2 = Intersection()(c , ell);

			Point a = res1.first;
			Point b = res1.second;

			Point p = res2.first;
			Point q = res2.second;

			double aq = sqrt(to_double(squared_distance(a, q)));
			double pb = sqrt(to_double(squared_distance(p, b)));
			double ap = sqrt(to_double(squared_distance(a, p)));
			double qb = sqrt(to_double(squared_distance(q, b)));

			double hyperdist = fabs(log(to_double((aq*pb)/(ap*qb))));

			return hyperdist;
		}
	};


	Construct_point_2 construct_point_2_object() const {
		return Construct_point_2();
	}

	class Construct_hyperbolic_segment_2 {
	public:
		typedef Segment_2 result_type;

		Construct_hyperbolic_segment_2() {}

		result_type operator()(const Point_2& p1, const Point_2& p2) const {
			
			// Check first if the points are collinear with the origin
			Circle_2 poincare(Point_2(FT(0),FT(0)), FT(1));
			Orientation ori = orientation(poincare.center(), p1, p2);
			if (ori == COLLINEAR) {
				Euclidean_segment_2 seg(p1, p2);
				return seg;
			}

			Compute_circle_orthogonal comp;
			Circle_2 supp = comp(Circle_2(p1, FT(0)), Circle_2(p2, FT(0)), poincare);
 
			if (ori == LEFT_TURN) {
				Circular_arc_2 carc(supp, p2, p1);
				return carc;
			} else {
				Circular_arc_2 carc(supp, p1, p2);
				return carc;
			}
		}
		
	}; // end Construct_hyperbolic_segment_2
	
	Construct_hyperbolic_segment_2
		construct_hyperbolic_segment_2_object() const
	{ return Construct_hyperbolic_segment_2(); }


	// wrong names kept for demo
	typedef Hyperbolic_traits_with_translations_2_adaptor<Self, Construct_hyperbolic_segment_2> Construct_segment_2;
	Construct_segment_2
	construct_segment_2_object() const
	{ return Construct_segment_2(); }
	


	class Construct_hyperbolic_line_2 {
	public:
		typedef Segment_2 result_type;

		Construct_hyperbolic_line_2() {}


		result_type operator()(const Point_2& p1, const Point_2& p2) const {
			
			// Check first if the points are collinear with the origin
			Circle_2 poincare(Point_2(FT(0),FT(0)), FT(1));
			Orientation ori = orientation(poincare.center(), p1, p2);
			if (ori == COLLINEAR) {
				Euclidean_line_2 ell(p1, p2);
				pair<Point, Point> res = Construct_inexact_intersection_2()(ell, poincare);
				return Euclidean_segment_2(res.first, res.second);
			}

			Compute_circle_orthogonal comp;
			Circle_2 supp = comp(Circle_2(p1, FT(0)), Circle_2(p2, FT(0)), poincare);
 			pair<Point, Point> res = Construct_inexact_intersection_2()(supp, poincare);
 			Point pp1 = res.first;
 			Point pp2 = res.second;

			if (ori == LEFT_TURN) {
				Circular_arc_2 carc(supp, pp2, pp1);
				return carc;
			} else {
				Circular_arc_2 carc(supp, pp1, pp2);
				return carc;
			}
		}
	}; // end Construct_hyperbolic_line_2


	Construct_hyperbolic_line_2
		construct_hyperbolic_line_2_object() const { 
			return Construct_hyperbolic_line_2(); 
		}
	
	Orientation_2
		orientation_2_object() const
	{ return Orientation_2();}
	
	Side_of_oriented_circle_2
		side_of_oriented_circle_2_object() const
	{ return Side_of_oriented_circle_2(); }
	


	class Construct_hyperbolic_circle_2 {
	public:
		typedef Circle_2 result_type;

		Construct_hyperbolic_circle_2() {}

		result_type operator()(Point_2 center, Point_2 p) {
			Circle_2 poincare(Point_2(FT(0),FT(0)), FT(1));
			Circle_2 circ = Compute_circle_in_pencil()(Circle_2(p,FT(0)), poincare, Circle_2(center,FT(0)));
			return circ;
		}
	};


	Construct_hyperbolic_circle_2
	construct_hyperbolic_circle_2_object() const 
	{
		return Construct_hyperbolic_circle_2();
	}


	class Construct_inexact_hyperbolic_bisector_2 {
	public:      
		typedef Segment_2 result_type;

		Construct_inexact_hyperbolic_bisector_2() {}
		
		result_type operator()(Point_2 p1, Point_2 p2) {
			Circle_2 poincare(Point_2(FT(0),FT(0)), FT(1));
			
			double np1 = sqrt(to_double(p1.x()*p1.x() + p1.y()*p1.y()));
			double np2 = sqrt(to_double(p2.x()*p2.x() + p2.y()*p2.y()));
			if ( fabs(np1 - np2) < 1e-16 ) {      
				Euclidean_line_2 ell = Construct_Euclidean_bisector_2()(p1, p2);
				pair<Point_2, Point_2> pts = Construct_inexact_intersection_2()(ell, poincare);
				return Euclidean_segment_2(pts.first, pts.second);
			}

			Circle_2 supp = Compute_circle_in_pencil()(poincare, Circle_2(p1,FT(0)), Circle_2(p2,FT(0)));
			pair<Point, Point> res = Construct_inexact_intersection_2()(supp, poincare);
 			Point pp1 = res.first;
 			Point pp2 = res.second;

 			Orientation ori = orientation(poincare.center(), pp1, pp2);

			if (ori == LEFT_TURN) {
				Circular_arc_2 carc(supp, pp2, pp1);
				return carc;
			} else {
				Circular_arc_2 carc(supp, pp1, pp2);
				return carc;
			}

		}

	}; // end Construct_hyperbolic_bisector_2
	
	Construct_inexact_hyperbolic_bisector_2
	construct_inexact_hyperbolic_bisector_2_object() const
	{ return Construct_inexact_hyperbolic_bisector_2(); }



	class Construct_hyperbolic_bisector_2 {
	public:      
		typedef Segment_2 result_type;

		Construct_hyperbolic_bisector_2() {}
		

		result_type operator()(Point_2 p1, Point_2 p2) {
			Circle_2 poincare(Point_2(FT(0),FT(0)), FT(1));
				
			if ( Compare_distance_2()(Point_2(FT(0),FT(0)), p1, p2) == EQUAL ){      
				Euclidean_line_2 ell = Construct_Euclidean_bisector_2()(p1, p2);
				pair<Point_2, Point_2> pts = Construct_intersection_2()(ell, poincare);
				return Euclidean_segment_2(pts.first, pts.second);
			}	

			Circle_2 supp = Compute_circle_in_pencil()(poincare, Circle_2(p1,FT(0)), Circle_2(p2,FT(0)));
			pair<Point, Point> res = Construct_inexact_intersection_2()(supp, poincare);
 			Point pp1 = res.first;
 			Point pp2 = res.second;

 			Orientation ori = orientation(poincare.center(), pp1, pp2);

			if (ori == LEFT_TURN) {
				Circular_arc_2 carc(supp, pp2, pp1);
				return carc;		
			} else {
				Circular_arc_2 carc(supp, pp1, pp2);
				return carc;
			}
		}

	}; // end Construct_hyperbolic_bisector_2
	
	Construct_hyperbolic_bisector_2
	construct_hyperbolic_bisector_2_object() const
	{ return Construct_hyperbolic_bisector_2(); }
	
	Construct_Euclidean_bisector_2
	construct_Euclidean_bisector_2_object() const
	{ return Construct_Euclidean_bisector_2(); }	



	class Construct_intersection_2 {
	public:
		Construct_intersection_2() {}

		Point_2 operator()(Euclidean_line_2 ell1, Euclidean_line_2 ell2) {
			if (ell1.b() == FT(0)) {
				std::swap(ell1, ell2);
			}
			
			CGAL_assertion(ell1.b() != FT(0));
			if (ell2.b() != FT(0)) {
				CGAL_assertion( ell1.a()/ell1.b() != ell2.a()/ell2.b() );
			}

			FT lambda1 = -ell1.a()/ell1.b();
			FT mu1     = -ell1.c()/ell1.b();
			FT x = ( -ell2.c() - mu1*ell2.b() )/( ell2.a() + lambda1*ell2.b() );
			FT y = lambda1*x + mu1;
			return Point_2(x, y);
		}

		std::pair<Point_2, Point_2> operator()(Euclidean_line_2 ell, Circle_2 c) {
			if (ell.b() == FT(0)) {
				FT p 	= c.center().x();
				FT q 	= c.center().y(); 
				FT y1  	= q + CGAL::sqrt(c.squared_radius() - p*p);
				FT y2  	= q - CGAL::sqrt(c.squared_radius() - p*p);
				Point_2 p1(FT(0), y1);
				Point_2 p2(FT(0), y2);
				return make_pair(p1, p2);
			}

			FT lambda = -ell.a()/ell.b();
			FT mu 	  = -ell.c()/ell.b();
			FT p 	  = c.center().x();
			FT q 	  = c.center().y(); 
			FT A = FT(1) + lambda*lambda;
			FT B = FT(2)*( lambda * mu - lambda*q - p);
			FT C = p*p + mu*mu + q*q - c.squared_radius() - FT(2)*q*mu;
			FT Delta = B*B - FT(4)*A*C;
			FT x1 = (-B + CGAL::sqrt(Delta))/(FT(2)*A);
			FT x2 = (-B - CGAL::sqrt(Delta))/(FT(2)*A);
			FT y1 = lambda*x1 + mu;
			FT y2 = lambda*x2 + mu;

			Point_2 sol1(x1, y1);
			Point_2 sol2(x2, y2);
			return make_pair(sol1, sol2);
		}

		std::pair<Point_2, Point_2> operator()(Circle_2 c, Euclidean_line_2 ell) {
			return operator()(ell, c);
		}

		std::pair<Point_2, Point_2> operator()(Circle_2 c1, Circle_2 c2) {
			
			FT xa = c1.center().x(), ya = c1.center().y();
			FT xb = c2.center().x(), yb = c2.center().y();
			FT d2 = squared_distance(c1.center(), c2.center());
			FT ra = CGAL::sqrt(c1.squared_radius());
			FT rb = CGAL::sqrt(c2.squared_radius());
			FT K  = CGAL::sqrt(((ra+rb)*(ra+rb)-d2)*(d2-(ra-rb)*(ra-rb)))/FT(4); 

			FT xbase = (xb + xa)/FT(2) + (xb - xa)*(ra*ra - rb*rb)/d2/FT(2);
			FT xdiff = FT(2)*(yb - ya)*K/d2;
			FT x1 = xbase + xdiff;
			FT x2 = xbase - xdiff;

			FT ybase = (yb + ya)/FT(2) + (yb - ya)*(ra*ra - rb*rb)/d2/FT(2);
			FT ydiff = FT(-2)*(xb - xa)*K/d2;
			FT y1 = ybase + ydiff;
			FT y2 = ybase - ydiff;

			Point_2 res1(x1, y1);
			Point_2 res2(x2, y2);
			return make_pair(res1, res2);
		}


		Point_2 operator()(Hyperbolic_segment_2 s1, Hyperbolic_segment_2 s2) {
			if (Circular_arc_2* c1 = boost::get<Circular_arc_2>(&s1)) {
				if (Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2)) {
					pair<Point_2, Point_2> res = operator()(c1->circle(), c2->circle());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					pair<Point_2, Point_2> res = operator()(c1->circle(), ell2->supporting_line());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;
				}
			} else {
				Euclidean_segment_2* ell1 = boost::get<Euclidean_segment_2>(&s1);
				if (Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2)) {
					pair<Point_2, Point_2> res = operator()(ell1->supporting_line(), c2->circle());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;	
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					Point_2 p1 = operator()(ell1->supporting_line(), ell2->supporting_line());
					CGAL_assertion(CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1));
					return p1;
				}
			}
		}

	};


	Construct_intersection_2
	construct_intersection_2_object() const {
		return Construct_intersection_2();
	}




	class Construct_inexact_intersection_2 {
	public:
		Construct_inexact_intersection_2() {}

		Point_2 operator()(Euclidean_line_2 ell1, Euclidean_line_2 ell2) {

			if (fabs(to_double(ell1.b())) < 1e-16) {
				std::swap(ell1, ell2);
			}
			
			double a1 = to_double(ell1.a()), b1 = to_double(ell1.b()), c1 = to_double(ell1.c());
			double a2 = to_double(ell2.a()), b2 = to_double(ell2.b()), c2 = to_double(ell2.c());

			CGAL_assertion(fabs(b1) > 1e-16);
			if (fabs(b2) > 1e-16) {
				CGAL_assertion( fabs(a1/b1 - a2/b2) > 1e-16 );
			}

			double lambda1 = -a1/b1;
			double mu1     = -c1/b1;
			double x = ( -c2 - mu1*b2 )/( a2 + lambda1*b2 );
			double y = lambda1*x + mu1;
			return Point_2(x, y);
		}

		std::pair<Point_2, Point_2> operator()(Euclidean_line_2 ell, Circle_2 cc) {
			double a = to_double(ell.a()), b = to_double(ell.b()), c = to_double(ell.c());
			double p = to_double(cc.center().x()), q = to_double(cc.center().y()), r2 = to_double(cc.squared_radius());
			
			double A, B, C, D;
			double x1, y1, x2, y2;
			if (fabs(a) < 1e-16) {
				y1 = -c/b;  y2 = -c/b;
				A = b*p;
				D = -b*b*q*q + b*b*r2 - 2.*b*c*q - c*c;
				x1 = (A + sqrt(D))/b;
				x2 = (A - sqrt(D))/b;
			} else if (fabs(b) < 1e-16) {
				x1 = -c/a;  x2 = -c/a;
				A = q*a;
				D = -a*a*p*p + r2*a*a - 2.*a*c*p - c*c;
				y1 = (A + sqrt(D))/a;
				y2 = (A - sqrt(D))/a;
			} else {
				A = a*a*q - a*b*p-b*c;
				C = (-b*q - c)*a*a + b*b*p*a;
				D = -a*a*( b*b*q*q + 2.*q*(p*a + c)*b - b*b*r2 + (p*p - r2)*a*a + 2.*a*c*p + c*c );
				B = a*a + b*b;

				y1 = (A + sqrt(D))/B;
				y2 = (A - sqrt(D))/B;
				x1 = (C - b*sqrt(D))/(a*(a*a + b*b));
				x2 = (C + b*sqrt(D))/(a*(a*a + b*b));
			}

			Point_2 p1(x1, y1);
			Point_2 p2(x2, y2);

			return make_pair(p1, p2);
		}

		std::pair<Point_2, Point_2> operator()(Circle_2 c, Euclidean_line_2 ell) {
			return operator()(ell, c);
		}

		std::pair<Point_2, Point_2> operator()(Circle_2 c1, Circle_2 c2) {
			double xa = to_double(c1.center().x()), ya = to_double(c1.center().y());
			double xb = to_double(c2.center().x()), yb = to_double(c2.center().y());
			double d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);
			double ra = sqrt(to_double(c1.squared_radius()));
			double rb = sqrt(to_double(c2.squared_radius()));
			double K  = sqrt(((ra+rb)*(ra+rb)-d2)*(d2-(ra-rb)*(ra-rb)))/4.; 

			double xbase = (xb + xa)/2. + (xb - xa)*(ra*ra - rb*rb)/d2/2.;
			double xdiff = 2.*(yb - ya)*K/d2;
			double x1 = xbase + xdiff;
			double x2 = xbase - xdiff;

			double ybase = (yb + ya)/2. + (yb - ya)*(ra*ra - rb*rb)/d2/2.;
			double ydiff = -2.*(xb - xa)*K/d2;
			double y1 = ybase + ydiff;
			double y2 = ybase - ydiff;

			Point_2 res1(x1, y1);
			Point_2 res2(x2, y2);
			return make_pair(res1, res2);
		}


		Point_2 operator()(Hyperbolic_segment_2 s1, Hyperbolic_segment_2 s2) {
			if (Circular_arc_2* c1 = boost::get<Circular_arc_2>(&s1)) {
				if (Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2)) {
					pair<Point_2, Point_2> res = operator()(c1->circle(), c2->circle());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					pair<Point_2, Point_2> res = operator()(c1->circle(), ell2->supporting_line());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;
				}
			} else {
				Euclidean_segment_2* ell1 = boost::get<Euclidean_segment_2>(&s1);
				if (Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2)) {
					pair<Point_2, Point_2> res = operator()(ell1->supporting_line(), c2->circle());
					Point_2 p1 = res.first;
					if (CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(CGAL::sqrt(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1));
					return p2;	
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					Point_2 p1 = operator()(ell1->supporting_line(), ell2->supporting_line());
					CGAL_assertion(CGAL::sqrt(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1));
					return p1;
				}
			}
		}

	};



	Construct_inexact_intersection_2
	construct_inexact_intersection_2_object() const {
		return Construct_inexact_intersection_2();
	}




	class Construct_hyperbolic_circumcenter_2_base {
	public:

		typedef Point_2 result_type;

		Construct_hyperbolic_circumcenter_2_base() {}

		Point_2 operator()(Point_2 p, Point_2 q, Point_2 r) {

			Hyperbolic_segment_2 s1 = Construct_hyperbolic_bisector_2()(p, q);
			Hyperbolic_segment_2 s2 = Construct_hyperbolic_bisector_2()(p, r);

			Point_2 rp = Construct_intersection_2()(s1, s2);
			return rp;
		}

	};


	typedef Hyperbolic_traits_with_translations_2_adaptor<Self, Construct_hyperbolic_circumcenter_2_base> Construct_hyperbolic_circumcenter_2;


	Construct_hyperbolic_circumcenter_2
	construct_hyperbolic_circumcenter_2_object() const {
		return Construct_hyperbolic_circumcenter_2();
	}





	class Construct_inexact_hyperbolic_circumcenter_2_base {
	public:

		typedef Point_2 result_type;

		Construct_inexact_hyperbolic_circumcenter_2_base() {}

		Point_2 operator()(Point_2 p, Point_2 q, Point_2 r) {

			Hyperbolic_segment_2 s1 = Construct_inexact_hyperbolic_bisector_2()(p, q);
			Hyperbolic_segment_2 s2 = Construct_inexact_hyperbolic_bisector_2()(p, r);

			Point_2 rp = Construct_inexact_intersection_2()(s1, s2);
			return rp;
		}

	};


	typedef Hyperbolic_traits_with_translations_2_adaptor<Self, Construct_inexact_hyperbolic_circumcenter_2_base> Construct_inexact_hyperbolic_circumcenter_2;


	Construct_inexact_hyperbolic_circumcenter_2
	construct_inexact_hyperbolic_circumcenter_2_object() const {
		return Construct_inexact_hyperbolic_circumcenter_2();
	}


	/****************************************************/
	class Side_of_hyperbolic_face_2 {
		
	public:
		typedef Bounded_side result_type;

		Side_of_hyperbolic_face_2() {}


		template<class Face_handle, class Hyperbolic_translation>
		Bounded_side operator()(const Point_2 p, Bounded_side sides[3], const Face_handle fh, const Hyperbolic_translation o) const {

			Point_2 p1 = Construct_point_2()(fh->vertex(0)->point(), o * fh->translation(0));
			Point_2 p2 = Construct_point_2()(fh->vertex(0)->point(), o * fh->translation(1));
			Point_2 p3 = Construct_point_2()(fh->vertex(0)->point(), o * fh->translation(2));

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
			return operator()(p, sides, fh, Hyperbolic_translation());
		}


	private:
		Bounded_side side_of_segment_2(const Point_2 query, const Point_2 p, const Point_2 q) const {
			
			// Check first if the points are collinear with the origin
			Circle_2 poincare(Point_2(FT(0),FT(0)), FT(1));
			Orientation ori = orientation(poincare.center(), p, q);
			if (ori == COLLINEAR) {
				Euclidean_line_2 seg(p, q);
				Orientation qori = orientation(query, p, q);
				if (qori == COLLINEAR) {
					return ON_BOUNDARY;
				} else {
					// It is sufficient that these are consistent.
					if (qori == LEFT_TURN) {
						return ON_BOUNDED_SIDE;
					} else {
						return ON_UNBOUNDED_SIDE;
					}
				}
			}

			Compute_circle_orthogonal comp;
			Circle_2 supp = comp(Circle_2(p, FT(0)), Circle_2(q, FT(0)), poincare);
 			return supp.bounded_side(query);
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



		class Side_of_original_octagon {
		public:
			Side_of_original_octagon() {}

			template <class Point_2_template>
			CGAL::Bounded_side operator()(Point_2_template p) {

				FT F2(2);
				FT qty = CGAL::sqrt(F2 + F2*CGAL::sqrt(F2));
				// The center of the Euclidean circle corresponding to the side s_1 (east)
				Point_2 CenterA ( qty/F2, FT(0) );
				Point_2 CenterB ( qty*CGAL::sqrt(F2)/FT(4), qty*CGAL::sqrt(F2)/FT(4) );

				// The squared radius of the Euclidean circle corresponding to the side s_1
				FT      Radius2 ( (CGAL::sqrt(F2) - FT(1)) / F2 );

				// Poincare disk (i.e., unit Euclidean disk)
				Circle_2 Poincare    ( Point(0, 0), 1*1 );

				// Euclidean circle corresponding to s_1
				Circle_2 EuclidCircA ( CenterA,     Radius2 );

				// Euclidean circle corresponding to s_2 (just rotate the center, radius is the same)
				Circle_2 EuclidCircBb( CenterB,     Radius2 );

				// This transformation brings the point in the first quadrant (positive x, positive y)
				FT x(FT(p.x()) > FT(0) ? p.x() : -p.x());
				FT y(FT(p.y()) > FT(0) ? p.y() : -p.y());

				// This brings the point in the first octant (positive x and y < x)
				if (y > x) {
					FT tmp = x;
					x = y;
					y = tmp;
				}

				// This tells us whether the point is on the side of the open boundary
				bool on_open_side = ( ( p.y() + tan(CGAL_PI / 8.) * p.x() ) < 0.0 );

				Point t(x, y);

				CGAL::Bounded_side PoincareSide = Poincare.bounded_side(t);
				CGAL::Bounded_side CircASide    = EuclidCircA.bounded_side(t);
				CGAL::Bounded_side CircBbSide   = EuclidCircBb.bounded_side(t);

				// First off, the point needs to be inside the Poincare disk. if not, there's no hope.
				if ( PoincareSide == CGAL::ON_BOUNDED_SIDE ) {
					
					// Inside the Poincare disk, but still outside the original domain
					if ( CircASide  == CGAL::ON_BOUNDED_SIDE || 
							 CircBbSide == CGAL::ON_BOUNDED_SIDE   ) {
						return CGAL::ON_UNBOUNDED_SIDE;
					}

					// Inside the Poincare disk and inside the original domain
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


		Side_of_original_octagon
		side_of_original_octagon_object() const {
			return Side_of_original_octagon();
		}



}; // class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2

} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H








