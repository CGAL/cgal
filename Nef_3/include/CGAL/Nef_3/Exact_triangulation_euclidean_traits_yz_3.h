// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
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
// $URL$ 
// $Id$ 
// 
//
// Author(s)     : Ralf Osbild <osbild@mpi-sb.mpg.de>

#ifndef CGAL_NEF3_EXACT_TRIANGULATION_EUCLIDEAN_TRAITS_YZ_3_H
#define CGAL_NEF3_EXACT_TRIANGULATION_EUCLIDEAN_TRAITS_YZ_3_H

#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/intersections.h>

namespace CGAL {

class Homogeneous_tag;
class Cartesian_tag;
template<class R, class Tag> struct Exact_intersect_yz_2;

template <class R>
struct Exact_intersect_yz_2 <R,Cartesian_tag>
{
   typedef typename R::Point_2     Point_2; 
   typedef typename R::Segment_2   Segment_2; 

   typedef typename R::Point_3     Point_3; 
   typedef typename R::Segment_3   Segment_3; 

   CGAL::Object operator() (Segment_3 s3, Segment_3 t3)
   {  Point_2 p2, q2;
      Point_3 p3, q3;

      // convert Segment_3 to Segment_2
      p3 = s3.source();
      q3 = s3.target();
      Segment_2 s2 (Point_2(p3.y(),p3.z()), Point_2(q3.y(),q3.z()));

      p3 = t3.source();
      q3 = t3.target();
      Segment_2 t2 (Point_2(p3.y(),p3.z()), Point_2(q3.y(),q3.z()));

      // convert intersection from Object_2 to Object_3
      // Note: there is not necessarily a spartial intersection,
      //       so all first components are faked!
      CGAL::Object obj = intersection (s2,t2);
      if ( CGAL::assign(p2, obj) )
      {  obj = make_object (Point_3 (0,p2.x(),p2.y()));
      }
      else if ( CGAL::assign(s2, obj) )
      {  p2 = s2.source();
	 q2 = s2.target();
	 obj = make_object( Segment_3(
	       Point_3(0,p2.x(),p2.y()), Point_3(0,q2.x(),q2.y()) ) );
      }
      return obj;
   }
};

template <class R>
struct Exact_intersect_yz_2 <R,Homogeneous_tag>
{  // Homogeneous
   typedef typename R::Point_2     Point_2; 
   typedef typename R::Segment_2   Segment_2; 

   typedef typename R::Point_3     Point_3; 
   typedef typename R::Segment_3   Segment_3; 

   CGAL::Object operator() (Segment_3 s3, Segment_3 t3)
   {  Point_2 p2, q2;
      Point_3 p3, q3;

      // convert Segment_3 to Segment_2
      p3 = s3.source();
      q3 = s3.target();
      Segment_2 s2 (Point_2(p3.hy(),p3.hz(),p3.hw()),
                    Point_2(q3.hy(),q3.hz(),q3.hw()));

      p3 = t3.source();
      q3 = t3.target();
      Segment_2 t2 (Point_2(p3.hy(),p3.hz(),p3.hw()),
                    Point_2(q3.hy(),q3.hz(),q3.hw()));

      // convert intersection from Object_2 to Object_3
      // Note: there is not necessarily a spartial intersection,
      //       so all first components are faked!
      CGAL::Object obj = intersection (s2,t2);
      if ( CGAL::assign(p2, obj) )
      {  obj = make_object (Point_3 (0,p2.hx(),p2.hy(),p2.hw()));
      }
      else if ( CGAL::assign(s2, obj) )
      {  p2 = s2.source();
	 q2 = s2.target();
	 obj = make_object( Segment_3(
	    Point_3 (0,p2.hx(),p2.hy(),p2.hw()),
	    Point_3 (0,q2.hx(),q2.hy(),q2.hw()) ) );
      }
      return obj;
   }
};

template <class R>
struct Exact_triangulation_euclidean_traits_yz_3
 : public CGAL::Projection_traits_yz_3<R>
{
   typedef CGAL::Exact_triangulation_euclidean_traits_yz_3<R> Traits;
   typedef CGAL::Projection_traits_yz_3<R>       Tet;
   typedef R Rp;

   typedef typename Tet::Point_2         Point_2;
   typedef typename Tet::Segment_2       Segment_2;
   typedef typename Tet::Triangle_2      Triangle_2;
   typedef typename Tet::Compare_x_2     Compare_x_2;
   typedef typename Tet::Compare_y_2     Compare_y_2;
   typedef typename Tet::Orientation_2   Orientation_2;
   typedef typename Tet::Side_of_oriented_circle_2 Side_of_oriented_circle_2;
   typedef typename Tet::Construct_segment_2   Construct_segment_2;
   typedef typename Tet::Construct_triangle_2  Construct_triangle_2;

   // for compatibility with previous versions
   typedef typename Tet::Point      Point;
   typedef typename Tet::Segment    Segment;
   typedef typename Tet::Triangle   Triangle;

   Exact_triangulation_euclidean_traits_yz_3(){}
   Exact_triangulation_euclidean_traits_yz_3(
      const Exact_triangulation_euclidean_traits_yz_3&) : Tet() {}
   Exact_triangulation_euclidean_traits_yz_3 &operator=(
      const Exact_triangulation_euclidean_traits_yz_3&) {return *this;}

   typedef Exact_intersect_yz_2<Rp,typename Rp::Kernel_tag> Intersect_2;

   Intersect_2 intersect_2_object () const
      { return Intersect_2(); }
};

} //namespace CGAL
#endif // CGAL_NEF3_EXACT_TRIANGULATION_EUCLIDEAN_TRAITS_YZ_3_H
