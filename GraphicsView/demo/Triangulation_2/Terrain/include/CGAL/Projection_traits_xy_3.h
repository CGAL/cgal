// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri, Mariette Yvinec



#ifndef CGAL_PROJECTION_TRAITS_XY_3_H
#define CGAL_PROJECTION_TRAITS_XY_3_H



#include <CGAL/intersection_2.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


CGAL_BEGIN_NAMESPACE 


namespace Projection_traits {
template <class R>
class Orientation_xy_3 
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT x(const Point &p) const { return p.x(); }
  typename R::FT y(const Point &p) const { return p.y(); }
  typedef typename R::Point_2 PP;
  CGAL::Orientation operator()(const Point& p,
			       const Point& q,
			       const Point& r)
  {
    typename R::Orientation_2 orient = R().orientation_2_object();
#ifdef CGAL_3_2  
    return orient(p,q,r);
#else
    typename R::Construct_point_2 construct_point_2 = R().construct_point_2_object();
    typename R::Point_2 pp = construct_point_2(x(p), y(p));
    typename R::Point_2 qq = construct_point_2(x(q), y(q));
    typename R::Point_2 rr = construct_point_2(x(r), y(r));
    
    return orient(pp, qq, rr);
#endif   
  }
};




template <class R>
class  Intersect_xy_3
{
public:
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Segment_3 Segment_3;
  typedef typename R::Point_2   Point_2; 
  typedef typename R::Vector_2  Vector_2; 
  typedef typename R::Segment_2 Segment_2;
  
  typename R::FT x(const Point_3 &p) const { return p.x(); }
  typename R::FT y(const Point_3 &p) const { return p.y(); }

  Point_2 project(const Point_3& p)
  {
    return Point_2(x(p),y(p));
  } 

#if 0
  Object operator()(const Segment_3& s1, const Segment_3& s2)
  {
    typedef typename CGAL::Exact_predicates_exact_constructions_kernel::Point_2 ExactPoint_2;
    typedef typename CGAL::Exact_predicates_exact_constructions_kernel::Segment_2 ExactSegment_2;
    
    
    

    ExactPoint_2 s1p(s1.source().x(), s1.source().y());
    ExactPoint_2 t1p(s1.target().x(), s1.target().y());
    ExactSegment_2 s1_2(s1p, t1p);
    ExactSegment_2 s2_2(ExactPoint_2(s2.source().x(), s2.source().y()),
			ExactPoint_2(s2.target().x(), s2.target().y()));
    Object o = intersection(s1_2,s2_2);
    ExactPoint_2 pi;
    if(assign(pi,o)){
      double l1 = std::sqrt(to_double(squared_distance(s1p,t1p)));
      double l2 = std::sqrt(to_double(squared_distance(s1p,pi)));
      double ratio = l2/l1;
      double z = s1.source().z() + ratio * (s1.target().z() - s1.source().z());
      Point_3 res(to_double(pi.x()), to_double(pi.y()), z);
      return make_object(res);
    } else {
      std::cerr << "NOT YET IMPLEMENTED: Intersection is not a point" << std::endl;
      Point_3 res;
      return make_object(res);
    }
  }

#else 

  Object operator()(const Segment_3& s1, const Segment_3& s2)
  {
    Point_2 s1p = project(s1.source());
    Point_2 t1p = project(s1.target());
    Segment_2 s1_2(s1p, t1p);
    Segment_2 s2_2(project(s2.source()), project(s2.target()));
    Object o = intersection(s1_2,s2_2);
    Point_2 pi;
    if(assign(pi,o)){
      double l1 = std::sqrt(to_double(squared_distance(s1p,t1p)));
      double l2 = std::sqrt(to_double(squared_distance(s1p,pi)));
      double ratio = l2/l1;
      Point_3 p = s1.source() + ratio * (s1.target() - s1.source());
      Point_3 res(pi.x(), pi.y(), p.z());
      return make_object(res);
    } else {
      std::cerr << "NOT YET IMPLEMENTED: Intersection is not a point" << std::endl;
      return Object();
    }
  }
#endif




};

template <class R>
class Compare_distance_xy_3
{
public:
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Point_2   Point_2;   
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return p.x(); }
  typename R::FT y(const Point_3 &p) const { return p.y(); }

  Point_2 project(const Point_3& p)
  {
    return Point_2(x(p),y(p));
  }

  Comparison_result operator()(const Point_3& p,const Point_3& q,const Point_3& r)
  {
    Point_2 p2 = project(p);
    Point_2 q2 = project(q);
    Point_2 r2 = project(r);
    return compare_distance_to_point(p2,q2,r2);
  }
};


template <class R>
class Squared_distance_xy_3
{
public:
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Point_2   Point_2; 
  typedef typename R::Line_3    Line_3; 
  typedef typename R::Line_2    Line_2;
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return p.x(); }
  typename R::FT y(const Point_3 &p) const { return p.y(); }

  Point_2 project(const Point_3& p)
  {
    return Point_2(x(p),y(p));
  }

  RT operator()(const Line_3& l, const Point_3& p)
  {
    Point_2 p2 = project(p);
    Line_2 l2 = Line_2(project(l.point(0)), project(l.point(1)));
    return squared_distance(p2, l2);
  }
};


template <class R>
class Side_of_oriented_circle_xy_3 
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT x(const Point &p) const { return p.x(); }
  typename R::FT y(const Point &p) const { return p.y(); }

  CGAL::Oriented_side operator() (const Point &p, 
				  const Point &q,
				  const Point &r, 
				  const Point &s) 
    {
      typename R::Side_of_oriented_circle_2 side_of_oriented_circle_2 = R().side_of_oriented_circle_2_object();
#ifdef CGAL_3_2

      return side_of_oriented_circle_2(p,
				       q,
				       r,
				       s);
#else
      typename R::Point_2 pp(x(p), y(p)), qq(x(q), y(q)), rr(x(r), y(r)), ss(x(s), y(s));
      return side_of_oriented_circle_2(pp,
				       qq,
				       rr,
				       ss);
#endif
    }
};

} // namespace Projection_traits

template < class R >
class Projection_traits_xy_3 {
public:
  typedef Projection_traits_xy_3<R> Traits;
  typedef R Rp;
  typedef typename Rp::FT                   FT;
  typedef typename Rp::Point_3     Point_2;
  typedef typename Rp::Segment_3   Segment_2;
  typedef typename Rp::Triangle_3  Triangle_2;
  typedef typename Rp::Line_3      Line_2;

  typedef typename Rp::Less_x_3             Less_x_2;
  typedef typename Rp::Less_y_3             Less_y_2;
  typedef typename Rp::Compare_x_3          Compare_x_2;
  typedef typename Rp::Compare_y_3          Compare_y_2;


  typedef Projection_traits::Orientation_xy_3<Rp>              Orientation_2;
  //typedef typename Rp::Orientation_2              Orientation_2;

  typedef Projection_traits::Side_of_oriented_circle_xy_3<Rp>  Side_of_oriented_circle_2;
  //typedef typename Rp::Side_of_oriented_circle_2 Side_of_oriented_circle_2;

  typedef Projection_traits::Squared_distance_xy_3<Rp>         Compute_squared_distance_2;
  typedef Projection_traits::Compare_distance_xy_3<Rp>         Compare_distance_2;
  typedef typename Rp::Construct_segment_3  Construct_segment_2;
  typedef typename Rp::Construct_line_3     Construct_line_2;
  typedef typename Rp::Construct_triangle_3 Construct_triangle_2;
  typedef Projection_traits::Intersect_xy_3<Rp>                Intersect_2;
  
  Projection_traits_xy_3()
  {}

  Projection_traits_xy_3(const Projection_traits_xy_3&)
  {}

  Projection_traits_xy_3 &
  operator=(const Projection_traits_xy_3&)
  {
    return *this;
  }

  typename Rp::FT 
  x(const Point_2 &p) const 
  { 
    return p.x(); 
  }

  typename Rp::FT 
  y(const Point_2 &p) const 
  {
    return p.y(); 
  }
    
 
  Less_x_2
  less_x_2_object() const
    { return Less_x_2();}

  Less_y_2
  less_y_2_object() const
    { return Less_y_2();}

  Compare_x_2
  compare_x_2_object() const
    { return Compare_x_2();}

  Compare_y_2
  compare_y_2_object() const
    { return Compare_y_2();}
  
  Compare_distance_2
  compare_distance_2_object() const
  {
    return Compare_distance_2();
  }
  
  Orientation_2
  orientation_2_object() const
    { return Orientation_2();}

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
    {return Side_of_oriented_circle_2();}
  
  Intersect_2
  intersect_2_object () const
  {
    return Intersect_2();
  }

  Compute_squared_distance_2
  compute_squared_distance_2_object () const
  {
    return Compute_squared_distance_2();
  }
  
  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}
  
  Construct_line_2  construct_line_2_object() const
    {return Construct_line_2();}
  
  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}

};
  

CGAL_END_NAMESPACE 

#endif // CGAL_PROJECTION_TRAITS_XY_3_H

