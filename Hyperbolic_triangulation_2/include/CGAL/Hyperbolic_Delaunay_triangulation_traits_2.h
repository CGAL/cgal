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


#ifndef CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
#define CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/utility.h>
#include <CGAL/Origin.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/internal/Exact_complex.h>
#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"


using std::pair;
using std::make_pair;

namespace CGAL {



template< class Kernel >
class Hyperbolic_Delaunay_triangulation_traits_2 {

typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>  Self;  

private:


	class Circular_arc_2 {
		typedef typename Kernel::FT 	  		FT;
		typedef Exact_complex<FT> 		  		Cplx;
		typedef typename Kernel::Point_2 	  	Point;
		typedef typename Kernel::Circle_2 	  	Circle;
		typedef typename Kernel::Orientation_2 	Orientation_2;

	private:
		Circle _c;
		Point _s, _t;

	public:
		Circular_arc_2() :
			_c(Point(FT(0),FT(0)), FT(0)), _s(FT(0),FT(0)), _t(FT(0),FT(0)) {}
		
		Circular_arc_2(Circle c, Point source, Point target) :
			_c(c), _s(source), _t(target) {}

		Circular_arc_2(Point p1, Point p2) {
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

		Circle supporting_circle() const {
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
	    	return typename Kernel::Construct_bbox_2()(*this);
	  	}

	};

public:

	typedef typename Kernel::FT          								FT;
	typedef typename Kernel::RT 										RT;
	typedef typename Kernel::Kernel_base 								Kernel_base;
	typedef typename Kernel::Point_2     								Point_2;
	typedef Point_2                 									Point;
	typedef Point_2 													Voronoi_point;
	typedef typename Kernel::Circle_2    								Circle_2;
	typedef typename Kernel::Line_2      								Euclidean_line_2;
	typedef boost::variant<Circle_2,Euclidean_line_2>    				Euclidean_circle_or_line_2; 
	typedef Self::Circular_arc_2										Circular_arc_2;
	typedef typename Kernel::Segment_2                       			Euclidean_segment_2; //only used internally here
	typedef boost::variant<Circular_arc_2, Euclidean_segment_2>  		Hyperbolic_segment_2;

	typedef typename Kernel::Triangle_2 								Triangle_2;	
	//typedef typename Kernel::Bbox_2 									Bbox_2;
	typedef typename Kernel::Orientation_2               				Orientation_2;
	typedef typename Kernel::Side_of_oriented_circle_2   				Side_of_oriented_circle_2;

	// only kept for demo to please T2graphicsitems
	typedef Euclidean_segment_2  										Line_segment_2;
	typedef Hyperbolic_segment_2 										Segment_2;

	typedef typename Kernel::Compare_x_2                				Compare_x_2;
	typedef typename Kernel::Compare_y_2                				Compare_y_2;

	typedef typename Kernel::Less_x_2                   				Less_x_2;
  	typedef typename Kernel::Less_y_2                   				Less_y_2;

  	// The objects Ray_2, Iso_rectangle_2 and Line_2 are needed by the CGAL::Qt::PainterOstream
  	typedef typename Kernel::Direction_2 								Direction_2;
  	typedef typename Kernel::Vector_2 									Vector_2;
  	typedef typename Kernel::Ray_2 										Ray_2;
  	typedef typename Kernel::Iso_rectangle_2 							Iso_rectangle_2;
  	typedef Euclidean_line_2 											Line_2;

	// the following types are only used internally in this traits class, 
	// so they need not be documented, and they don't need _object()
	typedef typename Kernel::Construct_bisector_2       				Construct_Euclidean_bisector_2;
	typedef typename Kernel::Construct_triangle_2       				Construct_triangle_2;
	typedef typename Kernel::Compare_distance_2        					Compare_distance_2;
	typedef typename Kernel::Construct_point_2         					Construct_point_2;
		
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
			FT sqradius = - lambda*(squared_distance(c1.center(),origin)-c1.squared_radius())
						  - (1-lambda)*(squared_distance(c2.center(),origin)-c2.squared_radius())
						  + squared_distance(center,origin);
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

	Compare_distance_2 
	compare_distance_2_object() const 
	{ return Compare_distance_2();} 

	Compare_x_2 
	compare_x_2_object() const 
	{ return Compare_x_2();} 

	Compare_y_2 
	compare_y_2_object() const 
	{ return Compare_y_2();} 

	Less_y_2 
	less_y_2_object() const 
	{ return Less_y_2();} 

	Less_x_2 
	less_x_2_object() const 
	{ return Less_x_2();} 

	Construct_triangle_2
	construct_triangle_2_object() const 
	{ return Construct_triangle_2();} 


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
	typedef Construct_hyperbolic_segment_2 Construct_segment_2;
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
		typedef Exact_complex<FT> 		Exact_complex;
	public:      
		typedef Segment_2 result_type;

		Construct_inexact_hyperbolic_bisector_2() {}
		
		result_type operator()(Point_2 p1, Point_2 p2) {
			Circle_2 poincare(Point_2(FT(0),FT(0)), FT(1));

			Circle_2 supp = Compute_circle_in_pencil()(poincare, Circle_2(p1,FT(0)), Circle_2(p2,FT(0)));
			pair<Point, Point> res = Construct_inexact_intersection_2()(supp, poincare);
 			Point pp1 = res.first;
 			Point pp2 = res.second;

 			return Construct_hyperbolic_segment_2()(pp1, pp2);
		}


	}; // end Construct_hyperbolic_bisector_2
	
	Construct_inexact_hyperbolic_bisector_2
	construct_inexact_hyperbolic_bisector_2_object() const
	{ return Construct_inexact_hyperbolic_bisector_2(); }



	class Construct_hyperbolic_inversion_2 {
	public:
		Construct_hyperbolic_inversion_2() {};

		Point_2 operator()(Point_2 p, Hyperbolic_segment_2 seg) {
			typedef Exact_complex<FT> 		Exact_complex;
			Exact_complex cp(p.x(), p.y());
			Circular_arc_2 * supp = boost::get<Circular_arc_2>(&seg);
			if (supp) {
				Exact_complex rp = cp.invert_in_circle(supp->supporting_circle());
				return Point_2(rp.real(), rp.imag());
			} else {
				std::cout << "Inversion FAILED!!" << std::endl;
				return Point(0,0);
			}
		}
	};


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
			pair<Point, Point> res = Construct_intersection_2()(supp, poincare);
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


		// constructs the hyperbolic bisector of segment [p,q] 
	    // limited by hyperbolic circumcenter(p,q,r) on one side
	    // and going to the infinite line on the other side
		result_type operator()(Point_2 p, Point_2 q, Point_2 r) {
			result_type res = this->operator()(p,q);
			Point_2 c = Construct_hyperbolic_circumcenter_2()(p,q,r);
			Orientation ori = orientation(c, p, q);
			Circle_2 poincare(Point_2(0,0), 1);
			if (Euclidean_segment_2* seg = boost::get<Euclidean_segment_2>(&res)) {	
				std::pair<Point_2, Point_2> ip = Construct_intersection_2()(poincare, seg->supporting_line());
				if (ori == LEFT_TURN)
					return Euclidean_segment_2(c, ip.first);
				else
					return Euclidean_segment_2(c, ip.second);
			} else {
				Circular_arc_2* supp = boost::get<Circular_arc_2>(&res);
				std::pair<Point_2, Point_2> ip = Construct_intersection_2()(poincare, supp->supporting_circle());

				if (orientation(ip.first, p, q) == LEFT_TURN) {
					if (orientation(supp->supporting_circle().center(), c, ip.second) == LEFT_TURN) {
						return Circular_arc_2(supp->supporting_circle(), c, ip.second);
					} else {
						return Circular_arc_2(supp->supporting_circle(), ip.second, c);
					}
				} else {
					if (orientation(supp->center(), c, ip.first) == LEFT_TURN) {
						return Circular_arc_2(supp->supporting_circle(), c, ip.first);
					} else {
						return Circular_arc_2(supp->supporting_circle(), ip.first, c);
					}
				}
			}
		}


		// constructs the hyperbolic bisector of segment [p,q] limited by 
	    // circumcenter(p,q,r) on one side
	    // and circumcenter(p,s,q) on the other side
		result_type operator()(Point_2 p, Point_2 q, Point_2 r, Point_2 s) {
			result_type res = this->operator()(p,q);
			Point_2 c1 = Construct_hyperbolic_circumcenter_2()(p,q,r);
			Point_2 c2 = Construct_hyperbolic_circumcenter_2()(p,s,q);

			if (Euclidean_segment_2* seg = boost::get<Euclidean_segment_2>(&res)) {	
				return Euclidean_segment_2(c1, c2);
			} else {
				Circular_arc_2* supp = boost::get<Circular_arc_2>(&res);
				if (orientation(Point(0,0), c1, c2) == LEFT_TURN) {
					return Circular_arc_2(supp->supporting_circle(), c2, c1);
				} else { 
					return Circular_arc_2(supp->supporting_circle(), c1, c2);
				}
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
					pair<Point_2, Point_2> res = operator()(c1->supporting_circle(), c2->supporting_circle());
					Point_2 p1 = res.first;
					if (p1.x()*p1.x() + p1.y()*p1.y() < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
                    CGAL_assertion(p2.x()*p2.x() + p2.y()*p2.y() < FT(1));
					return p2;
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					pair<Point_2, Point_2> res = operator()(c1->supporting_circle(), ell2->supporting_line());
					Point_2 p1 = res.first;
					if (p1.x()*p1.x() + p1.y()*p1.y() < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(p2.x()*p2.x() + p2.y()*p2.y() < FT(1));
					return p2;
				}
			} else {
				Euclidean_segment_2* ell1 = boost::get<Euclidean_segment_2>(&s1);
				if (Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2)) {
					pair<Point_2, Point_2> res = operator()(ell1->supporting_line(), c2->supporting_circle());
					Point_2 p1 = res.first;
					if (p1.x()*p1.x() + p1.y()*p1.y() < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(p2.x()*p2.x() + p2.y()*p2.y() < FT(1));
					return p2;	
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					Point_2 p1 = operator()(ell1->supporting_line(), ell2->supporting_line());
					CGAL_assertion(p1.x()*p1.x() + p1.y()*p1.y() < FT(1));
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
					pair<Point_2, Point_2> res = operator()(c1->supporting_circle(), c2->supporting_circle());
					Point_2 p1 = res.first;
					if (p1.x()*p1.x() + p1.y()*p1.y() < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(p2.x()*p2.x() + p2.y()*p2.y() < FT(1));
					return p2;
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					pair<Point_2, Point_2> res = operator()(c1->supporting_circle(), ell2->supporting_line());
					Point_2 p1 = res.first;
					if (p1.x()*p1.x() + p1.y()*p1.y() < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(p2.x()*p2.x() + p2.y()*p2.y() < FT(1));
					return p2;
				}
			} else {
				Euclidean_segment_2* ell1 = boost::get<Euclidean_segment_2>(&s1);
				if (Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2)) {
					pair<Point_2, Point_2> res = operator()(ell1->supporting_line(), c2->supporting_circle());
					Point_2 p1 = res.first;
					if (p1.x()*p1.x() + p1.y()*p1.y() < FT(1)) {
						return p1;
					}
					Point_2 p2 = res.second;
					CGAL_assertion(p2.x()*p2.x() + p2.y()*p2.y() < FT(1));
					return p2;	
				} else {
					Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
					Point_2 p1 = operator()(ell1->supporting_line(), ell2->supporting_line());
					CGAL_assertion(p1.x()*p1.x() + p1.y()*p1.y() < FT(1));
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

		typedef Voronoi_point result_type;

		Construct_hyperbolic_circumcenter_2_base() {}

		Point_2 operator()(Point_2 p, Point_2 q, Point_2 r) {

			Hyperbolic_segment_2 s1 = Construct_hyperbolic_bisector_2()(p, q);
			Hyperbolic_segment_2 s2 = Construct_hyperbolic_bisector_2()(q, r);
			Hyperbolic_segment_2 s3 = Construct_hyperbolic_bisector_2()(p, r);

			Circular_arc_2* arc1 = boost::get<Circular_arc_2>(&s1);
			Circular_arc_2* arc2 = boost::get<Circular_arc_2>(&s2);
			Circular_arc_2* arc3 = boost::get<Circular_arc_2>(&s3);

			FT r1(0);
			FT r2(0);
			FT r3(0);

			if (arc1) {
				r1 = arc1->squared_radius();
			}
			if (arc2) {
				r2 = arc2->squared_radius();
			}
			if (arc3) {
				r3 = arc3->squared_radius();
			}

			Point_2 rp;
			if (r1 < r2) {
				if (r1 < r3) {
					if (r2 < r3) { 	// r1 < r2 < r3
						rp = Construct_intersection_2()(s1, s2);
					} else { 		// r1 < r3 < r2
						rp = Construct_intersection_2()(s1, s3);
					}
				} else { 			// r3 < r1 < r2
					rp = Construct_intersection_2()(s3, s1);
				}
			} else { 				// r2 < r1
				if (r1 < r3) { 		// r2 < r1 < r3
					rp = Construct_intersection_2()(s2, s1);
				} else { 			// r2 < r1, r3 < r1
					if (r2 < r3) {	// r2 < r3 < r1
						rp = Construct_intersection_2()(s2, s3);
					} else {		// r3 < r2 < r1
						rp = Construct_intersection_2()(s3, s2);
					}
				}
			}

			return rp;
		}

	};


	typedef Construct_hyperbolic_circumcenter_2_base Construct_hyperbolic_circumcenter_2;


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
			Hyperbolic_segment_2 s3 = Construct_inexact_hyperbolic_bisector_2()(q, r);

			Circular_arc_2* arc1 = boost::get<Circular_arc_2>(&s1);
			Circular_arc_2* arc2 = boost::get<Circular_arc_2>(&s2);
			Circular_arc_2* arc3 = boost::get<Circular_arc_2>(&s3);

			double r1(0);
			double r2(0);
			double r3(0);

			if (arc1) {
				r1 = CGAL::to_double(arc1->squared_radius());
			}
			if (arc2) {
				r2 = CGAL::to_double(arc2->squared_radius());
			}
			if (arc3) {
				r3 = CGAL::to_double(arc3->squared_radius());
			}

			Point_2 rp;
			if (r1 < r2) {
				if (r1 < r3) {
					if (r2 < r3) { 	// r1 < r2 < r3
						rp = Construct_inexact_intersection_2()(s1, s2);
					} else { 		// r1 < r3 < r2
						rp = Construct_inexact_intersection_2()(s1, s3);
					}
				} else { 			// r3 < r1 < r2
					rp = Construct_inexact_intersection_2()(s3, s1);
				}
			} else { 				// r2 < r1
				if (r1 < r3) { 		// r2 < r1 < r3
					rp = Construct_inexact_intersection_2()(s2, s1);
				} else { 			// r2 < r1, r3 < r1
					if (r2 < r3) {	// r2 < r3 < r1
						rp = Construct_inexact_intersection_2()(s2, s3);
					} else {		// r3 < r2 < r1
						rp = Construct_inexact_intersection_2()(s3, s2);
					}
				}
			}

			return rp;

		}

	};


	typedef Construct_inexact_hyperbolic_circumcenter_2_base Construct_inexact_hyperbolic_circumcenter_2;


	Construct_inexact_hyperbolic_circumcenter_2
	construct_inexact_hyperbolic_circumcenter_2_object() const {
		return Construct_inexact_hyperbolic_circumcenter_2();
	}


	// For details see the JoCG paper (5:56-85, 2014)
  class Is_hyperbolic
  {
  public:
    typedef typename Kernel::Vector_3    Vector_3;
    typedef typename Kernel::Point_3     Point_3;

    bool operator() (const Point_2& p0, const Point_2& p1, const Point_2& p2) const
    {
      Vector_3 v0 = Vector_3(p0.x()*p0.x() + p0.y()*p0.y(), 
                             p1.x()*p1.x() + p1.y()*p1.y(), 
                             p2.x()*p2.x() + p2.y()*p2.y());
      
      Vector_3 v1 = Vector_3(p0.x(), p1.x(), p2.x());
      Vector_3 v2 = Vector_3(p0.y(), p1.y(), p2.y());
      Vector_3 v3 = Vector_3(FT(1), FT(1), FT(1));
      
      FT dt0 = determinant(v0, v1, v3);
      FT dt1 = determinant(v0, v2, v3);
      FT dt2 = determinant(v0 - v3, v1, v2);
      
      return dt0*dt0 + dt1*dt1 - dt2*dt2 < 0;
    }
    
    bool operator() (const Point_2& p0, const Point_2& p1, const Point_2& p2, int& ind) const
    {
      if (this->operator()(p0, p1, p2) == false) {
        ind = find_non_hyperbolic_edge(p0, p1, p2);
        return false;
      }
      return true;
    }
    
  private:
    
    // assume the face (p0, p1, p2) is non-hyperbolic
    int find_non_hyperbolic_edge(const Point_2& p0, const Point_2& p1, const Point_2& p2) const
    {
      typedef typename Kernel::Direction_2 Direction_2;
      
      Vector_3 v0 = Vector_3(p0.x()*p0.x() + p0.y()*p0.y(), 
                             p1.x()*p1.x() + p1.y()*p1.y(), 
                             p2.x()*p2.x() + p2.y()*p2.y());
      
      Vector_3 v1 = Vector_3(p0.x(), p1.x(), p2.x());
      Vector_3 v2 = Vector_3(p0.y(), p1.y(), p2.y());
      Vector_3 v3 = Vector_3(FT(1), FT(1), FT(1));
      
      FT dt0 = determinant(v0, 2*v2, -v3);
      FT dt1 = determinant(2*v1, v0, -v3);
      FT dt2 = determinant(2*v1, 2*v2, -v3);
      
      Direction_2 d0(p0.x()*dt2 - dt0, p0.y()*dt2 - dt1);
      Direction_2 d1(p1.x()*dt2 - dt0, p1.y()*dt2 - dt1);
      Direction_2 d2(p2.x()*dt2 - dt0, p2.y()*dt2 - dt1);
      
      Direction_2 d(dt0, dt1);
      
      if(d.counterclockwise_in_between(d0, d1)) {
        return 2;
      }
      
      if(d.counterclockwise_in_between(d1, d2)) {
        return 0;
      }
      
      return 1;
    }
  }; // end Is_hyperbolic

  Is_hyperbolic 
    Is_hyperbolic_object() const
  { return Is_hyperbolic(); }


	/****************************************************/
	class Side_of_hyperbolic_face_2 {
		
	public:
		typedef Bounded_side result_type;

		Side_of_hyperbolic_face_2() {}


		template<class Face_handle>
		Bounded_side operator()(const Point_2 p, Bounded_side sides[3], const Face_handle fh) const {

			Point_2 p1 = fh->vertex(0)->point();
			Point_2 p2 = fh->vertex(1)->point();
			Point_2 p3 = fh->vertex(2)->point();

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
		Hyperbolic_Delaunay_triangulation_traits_2() {}
	
		Hyperbolic_Delaunay_triangulation_traits_2(const Hyperbolic_Delaunay_triangulation_traits_2 & other) {}
	
		Hyperbolic_Delaunay_triangulation_traits_2 &operator=(const Hyperbolic_Delaunay_triangulation_traits_2 &)
		{
			return *this;
		}

}; // class Hyperbolic_Delaunay_triangulation_traits_2


} // namespace CGAL

#endif // CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H








