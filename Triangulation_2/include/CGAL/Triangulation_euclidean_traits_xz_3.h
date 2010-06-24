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
// Author(s)     : Mariette Yvinec

#ifndef CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_XZ_3_H
#define CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_XZ_3_H

#include <CGAL/triangulation_assertions.h>

#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Orientation_xz_3 
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT x(const Point &p) const { return p.x(); }
  typename R::FT y(const Point &p) const { return p.z(); }

  CGAL::Orientation operator()(const Point& p,
			       const Point& q,
			       const Point& r)
    {
      return orientationC2(x(p), y(p), x(q), y(q), x(r), y(r));
    }
};

template <class R>
class Side_of_oriented_circle_xz_3 
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT x(const Point &p) const { return p.x(); }
  typename R::FT y(const Point &p) const { return p.z(); }

  CGAL::Oriented_side operator() (const Point &p, 
				  const Point &q,
				  const Point &r, 
				  const Point &s) 
    {
      return side_of_oriented_circleC2(x(p), y(p),
				       x(q), y(q),
				       x(r), y(r),
				       x(s), y(s));
    }
};

template <class R>
class Compare_distance_xz_3
{
public:
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Point_2   Point_2;   
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return p.x(); }
  typename R::FT y(const Point_3 &p) const { return p.z(); }

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


template < class R >
class Triangulation_euclidean_traits_xz_3 {
public:
  typedef Triangulation_euclidean_traits_xz_3<R> Traits;
  typedef R Rp;
  typedef typename Rp::Point_3    Point_2;
  typedef typename Rp::Segment_3  Segment_2;
  typedef typename Rp::Triangle_3 Triangle_2;

  typedef typename Rp::Less_x_3             Less_x_2;
  typedef typename Rp::Less_z_3             Less_y_2;
  typedef typename Rp::Compare_x_3          Compare_x_2;
  typedef typename Rp::Compare_z_3          Compare_y_2;
  typedef Orientation_xz_3<Rp>              Orientation_2;
  typedef Side_of_oriented_circle_xz_3<Rp>  Side_of_oriented_circle_2;
  typedef Compare_distance_xz_3<Rp>         Compare_distance_2;
  typedef typename Rp::Construct_segment_3   Construct_segment_2;
  typedef typename Rp::Construct_triangle_3  Construct_triangle_2;

  // for compatibility with previous versions
  typedef Point_2      Point;
  typedef Segment_2    Segment;
  typedef Triangle_2   Triangle;
  
  Triangulation_euclidean_traits_xz_3(){}
  Triangulation_euclidean_traits_xz_3(
	      const Triangulation_euclidean_traits_xz_3& ){}
  Triangulation_euclidean_traits_xz_3 &operator=(
	 const Triangulation_euclidean_traits_xz_3& ){return *this;}

  typename Rp::FT x(const Point_2 &p) const { return p.x(); }
  typename Rp::FT y(const Point_2 &p) const { return p.z(); }

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

  Orientation_2
  orientation_2_object() const
    { return Orientation_2();}

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
    {return Side_of_oriented_circle_2();}

  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}

  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}

  Compare_distance_2
  compare_distance_2_object() const
  {
    return Compare_distance_2();
  }    
};

CGAL_END_NAMESPACE


#endif // CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_XZ_3_H
