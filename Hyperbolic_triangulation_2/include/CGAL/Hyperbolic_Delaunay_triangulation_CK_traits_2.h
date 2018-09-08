// Copyright (c) 2010-2016   INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: 
// $Id: 
// 
//
// Author(s)     : Mikhail Bogdanov
//                 Monique Teillaud <Monique.Teillaud@inria.fr>

#ifndef CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_CK_TRAITS_2_H
#define CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_CK_TRAITS_2_H

#include <CGAL/Circular_kernel_2/Intersection_traits.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/triangulation_assertions.h>
#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"
#include <CGAL/determinant.h>

namespace CGAL {

template < class R = CGAL::Exact_circular_kernel_2 >
class Hyperbolic_Delaunay_triangulation_CK_traits_2
  : public R
  // R is supposed to be a model of CircularKernel2
  {
public:
  typedef typename R::FT          FT;

  typedef typename R::Point_2     Hyperbolic_point_2;
  typedef typename R::Circle_2    Circle_2;
  typedef typename R::Line_2      Euclidean_line_2;
  typedef boost::variant<Circle_2,Euclidean_line_2>    Euclidean_circle_or_line_2; 

  typedef typename R::Triangle_2             Hyperbolic_triangle_2;

  typedef typename R::Circular_arc_2         Circular_arc_2;
  typedef typename R::Line_arc_2             Line_arc_2; 
  typedef typename R::Circular_arc_point_2   Circular_arc_point_2;
  typedef Circular_arc_point_2               Hyperbolic_Voronoi_point_2;
  typedef typename R::Segment_2                       Euclidean_segment_2; //only used internally here
  typedef boost::variant<Circular_arc_2, Line_arc_2>  Hyperbolic_segment_2;

  typedef typename R::Compare_x_2                Compare_x_2;
  typedef typename R::Compare_y_2                Compare_y_2;
  typedef typename R::Compare_distance_2         Compare_distance_2;
  typedef typename R::Orientation_2              Orientation_2;
  typedef typename R::Side_of_oriented_circle_2  Side_of_oriented_circle_2;

  // only kept for demo to please T2graphicsitems
  typedef Euclidean_segment_2  Line_segment_2;
  typedef Hyperbolic_segment_2 Segment_2;

  // the following types are only used internally in this traits class, 
  // so they need not be documented, and they don't need _object()
  typedef typename R::Collinear_2                Euclidean_collinear_2;
  typedef typename R::Construct_bisector_2       Construct_Euclidean_bisector_2;
  typedef typename R::Construct_midpoint_2       Construct_Euclidean_midpoint_2;
  typedef typename Kernel::Construct_point_2     Construct_hyperbolic_point_2;
  typedef typename R::Compute_squared_distance_2 Compute_squared_Euclidean_distance_2;
  typedef typename R::Has_on_bounded_side_2      Has_on_bounded_side_2;

  typedef typename R::Less_x_2                   Less_x_2;
  typedef typename R::Less_y_2                   Less_y_2;
      
public:

  Construct_hyperbolic_point_2
  construct_hyperbolic_point_2_object() const {
    return Construct_hyperbolic_point_2();
  }

  class Construct_hyperbolic_segment_2
  {
    
    typedef typename R::Construct_weighted_circumcenter_2 Construct_weighted_circumcenter_2;
    typedef typename R::Weighted_point_2 Weighted_point_2;
    typedef typename R::Point_2 Bare_point;

  public:
    Construct_hyperbolic_segment_2() 
      {}
    
    Hyperbolic_segment_2 operator()(const Hyperbolic_point_2& p, const Hyperbolic_point_2& q) const
    {
      Origin o;
      if(Euclidean_collinear_2()(p, q, Hyperbolic_point_2(o))){
        return Euclidean_segment_2(p, q);
      }
      
      Weighted_point_2 wp(p);
      Weighted_point_2 wq(q);
      Weighted_point_2 wo(Hyperbolic_point_2(o), FT(1)); // Poincaré circle 
      
      Bare_point center = Construct_weighted_circumcenter_2()(wp, wo, wq);
      FT sq_radius = Compute_squared_Euclidean_distance_2()(p, center);
      
      Circle_2 circle(center, sq_radius);
      // uncomment!!!
      assert(circle.has_on_boundary(p) && circle.has_on_boundary(q));
      
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
  typedef Construct_hyperbolic_segment_2 Construct_segment_2;
  Construct_segment_2
    construct_segment_2_object() const
  { return Construct_hyperbolic_segment_2(); }
  
  class Construct_hyperbolic_circumcenter_2
  {
  public:
    
    Hyperbolic_Voronoi_point_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 r) { 
      Origin o; 
      Hyperbolic_point_2 po = Hyperbolic_point_2(o);
      Circle_2 l_inf(po, FT(1));
     
      if ( Compare_distance_2()(po,p,q) == EQUAL && Compare_distance_2()(po,p,r) == EQUAL ) 
        return po; 
      
      // now supporting objects cannot both be Euclidean lines

      Euclidean_circle_or_line_2 bis_pq = Construct_circle_or_line_supporting_bisector()(p,q);
      Euclidean_circle_or_line_2 bis_qr = Construct_circle_or_line_supporting_bisector()(q,r);
      
      std::pair<Circular_arc_point_2, unsigned> pair;
      Euclidean_line_2* l;
      Circle_2* c;

      if ( Circle_2* c_pq = boost::get<Circle_2>(&bis_pq) )	{
        if ( Circle_2* c_qr = boost::get<Circle_2>(&bis_qr) ) {
	       typedef typename CK2_Intersection_traits<R, Circle_2, Circle_2>::type Intersection_result; 
	       std::vector< Intersection_result > inters;
	       intersection(*c_pq, *c_qr, std::back_inserter(inters));
	      
	       CGAL_triangulation_assertion(assign(pair,inters[0]));
	       if ( pair.second == 1 ) { 
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
      } else {
        // here bis_pq is a line, and bis_qr is necessarily a circle
        l = boost::get<Euclidean_line_2>(&bis_pq);
        c = boost::get<Circle_2>(&bis_qr);        
      }

      typedef typename CK2_Intersection_traits<R, Euclidean_line_2, Circle_2>::type Intersection_result; 
      std::vector< Intersection_result > inters;
      intersection(*l, *c, std::back_inserter(inters));

      CGAL_triangulation_assertion(assign(pair,inters[0]));
      if ( pair.second == 1 )	{ 
        if ( Has_on_bounded_side_2()( l_inf, pair.first ) )
          return pair.first;
	  
        CGAL_triangulation_assertion(assign(pair,inters[1]));
        return pair.first;
      }
      return pair.first;
    }

  }; // end Construct_hyperbolic_circumcenter_2
  
  Hyperbolic_Delaunay_triangulation_CK_traits_2() 
  {}
  
  Hyperbolic_Delaunay_triangulation_CK_traits_2(const Hyperbolic_Delaunay_triangulation_CK_traits_2 & other)
  {}
  
  Hyperbolic_Delaunay_triangulation_CK_traits_2 &operator=
  (const Hyperbolic_Delaunay_triangulation_CK_traits_2 &)
  {
    return *this;
  }
  
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
    Hyperbolic_segment_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q) const
    {
      Origin o; 
      Hyperbolic_point_2 po = Hyperbolic_point_2(o);
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
      operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 r, Hyperbolic_point_2 s)
    {
      CGAL_triangulation_precondition
	( (Orientation_2()(p,q,r) == ON_POSITIVE_SIDE) 
	  && (Orientation_2()(p,s,q) == ON_POSITIVE_SIDE) );
      CGAL_triangulation_precondition
	( (Side_of_oriented_circle_2()(p,q,r,s) == ON_NEGATIVE_SIDE) 
	  && (Side_of_oriented_circle_2()(p,s,q,r) == ON_NEGATIVE_SIDE) );

      Origin o; 
      Hyperbolic_point_2 po = Hyperbolic_point_2(o);

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
    Hyperbolic_segment_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 r)
    {
      CGAL_triangulation_precondition
	( Orientation_2()(p,q,r) == POSITIVE );

      Origin o; 
      Hyperbolic_point_2 po = Hyperbolic_point_2(o);
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
      
      Hyperbolic_point_2 approx_a(to_double(a.x()),to_double(a.y()));

      typedef typename 
	CK2_Intersection_traits<R, Circle_2, Circle_2>::type 
	Intersection_result; 
      std::vector< Intersection_result > inters;
      intersection(*c_pq, l_inf, std::back_inserter(inters));
      std::pair<Circular_arc_point_2, unsigned> pair;

      CGAL_triangulation_assertion(assign(pair,inters[0]));
      CGAL_triangulation_assertion(pair.second == 1);
      
      Hyperbolic_point_2 approx_pinf(to_double(pair.first.x()), to_double(pair.first.y()));
      Hyperbolic_point_2 approx_c(to_double(c_pq->center().x()),
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
    


  class Side_of_oriented_hyperbolic_segment_2 {
    typedef typename R::Construct_weighted_circumcenter_2  Construct_weighted_circumcenter_2;
    typedef typename R::Weighted_point_2                   Weighted_point_2;
    typedef typename R::Point_2                            Bare_point;
  
  public: 
    Side_of_oriented_hyperbolic_segment_2() {}

    typedef Oriented_side result_type;

    result_type operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 query) const {
      
      // Check first if the points are collinear with the origin
      Circle_2 poincare(Hyperbolic_point_2(FT(0),FT(0)), FT(1));
      Hyperbolic_point_2 O(FT(0), FT(0));
      Orientation ori = orientation(p, q, O);
      if (ori == COLLINEAR) {
        Euclidean_line_2 seg(p, q);
        Orientation qori = orientation(p, q, query);
        if (qori == COLLINEAR) {
          return ON_ORIENTED_BOUNDARY;
        } else {
          // It is sufficient that these are consistent.
          if (qori == LEFT_TURN) {
            return ON_POSITIVE_SIDE;
          } else {
            return ON_NEGATIVE_SIDE;
          }
        }
      }

      Weighted_point_2 wp(p);
      Weighted_point_2 wq(q);
      Weighted_point_2 wo(O, FT(1)); // Poincaré circle 

      Bare_point center = Construct_weighted_circumcenter_2()(wp, wo, wq);
      FT sq_radius = Compute_squared_Euclidean_distance_2()(p, center);

      Circle_2 circle(center, sq_radius);
      Bounded_side bs = circle.bounded_side(query);
      if (bs == ON_BOUNDARY) {
        return ON_ORIENTED_BOUNDARY;
      } else {
        if (bs == ON_BOUNDED_SIDE) {
          return ON_POSITIVE_SIDE;
        } else {
          return ON_NEGATIVE_SIDE;
        }
      }
    }
  };


  Side_of_oriented_hyperbolic_segment_2
  side_of_oriented_hyperbolic_segment_2_object() {
    return Side_of_oriented_hyperbolic_segment_2();
  }


  // For details see the JoCG paper (5:56-85, 2014)
  class Is_Delaunay_hyperbolic
  {
  public:
    typedef typename R::Vector_3    Vector_3;
    typedef typename R::Point_3     Point_3;

    bool operator() (const Hyperbolic_point_2& p0, const Hyperbolic_point_2& p1, const Hyperbolic_point_2& p2) const
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
    
    bool operator() (const Hyperbolic_point_2& p0, const Hyperbolic_point_2& p1, const Hyperbolic_point_2& p2, int& ind) const
    {
      if (this->operator()(p0, p1, p2) == false) {
        ind = find_non_hyperbolic_edge(p0, p1, p2);
        return false;
      }
      return true;
    }
    
  private:
    
    // assume the face (p0, p1, p2) is non-hyperbolic
    int find_non_hyperbolic_edge(const Hyperbolic_point_2& p0, const Hyperbolic_point_2& p1, const Hyperbolic_point_2& p2) const
    {
      typedef typename R::Direction_2 Direction_2;
      
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
  }; // end Is_Delaunay_hyperbolic

  Is_Delaunay_hyperbolic 
    is_Delaunay_hyperbolic_object() const
  { return Is_Delaunay_hyperbolic(); }

  // do not document
  // constructs the Euclidean circle or line supporting the hyperbolic
  // bisector of two points  
  class Construct_circle_or_line_supporting_bisector
  {
  public:
    Construct_circle_or_line_supporting_bisector()
      {}

    Euclidean_circle_or_line_2 
      operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q) const
    {
      Origin o; 
      Hyperbolic_point_2 po = Hyperbolic_point_2(o);
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
      Hyperbolic_point_2 middle = Construct_Euclidean_midpoint_2()(p, q);
      Hyperbolic_point_2 temp = middle + l.to_vector();

      if (Orientation_2()(middle, temp, Hyperbolic_point_2(x, y)) == ON_POSITIVE_SIDE)
	{ return Circle_2(Hyperbolic_point_2(x, y), sq_radius, CLOCKWISE); }
      return Circle_2(Hyperbolic_point_2(x, y), sq_radius, COUNTERCLOCKWISE);
    }
  }; // end Construct_supporting_circle_of_bisector

};


// Take out the code below to some separate file

#ifdef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
template  <>
struct Triangulation_structural_filtering_traits< Hyperbolic_Delaunay_triangulation_CK_traits_2<Epeck> > {
  typedef Tag_true Use_structural_filtering_tag;
};
#endif // CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H

#ifdef CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
template <>
struct Triangulation_structural_filtering_traits< Hyperbolic_Delaunay_triangulation_CK_traits_2<Epick> > {
  typedef Tag_true Use_structural_filtering_tag;
};
#endif // CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H

} //namespace CGAL 

#endif // CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_CK_TRAITS_2_H
