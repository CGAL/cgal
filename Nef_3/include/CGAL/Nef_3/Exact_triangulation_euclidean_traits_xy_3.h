// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ralf Osbild <osbild@mpi-sb.mpg.de>

#ifndef CGAL_NEF3_EXACT_TRIANGULATION_EUCLIDEAN_TRAITS_XY_3_H
#define CGAL_NEF3_EXACT_TRIANGULATION_EUCLIDEAN_TRAITS_XY_3_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/intersections.h>

namespace CGAL {

class Homogeneous_tag;
class Cartesian_tag;
template<class R, class Tag> struct Exact_intersect_xy_2;

template <class R>
struct Exact_intersect_xy_2 <R,Cartesian_tag>
{
  typedef typename R::Point_2     Point_2;
  typedef typename R::Segment_2   Segment_2;

  typedef typename R::Point_3     Point_3;
  typedef typename R::Segment_3   Segment_3;

  typedef  boost::variant<Point_3, Segment_3> variant_type;

  boost::optional<variant_type>
  operator() (const Segment_3& s3, const Segment_3& t3)
   {  Point_2 p2, q2;
      Point_3 p3, q3;

      // convert Segment_3 to Segment_2
      p3 = s3.source();
      q3 = s3.target();
      Segment_2 s2 (Point_2(p3.x(),p3.y()), Point_2(q3.x(),q3.y()));

      p3 = t3.source();
      q3 = t3.target();
      Segment_2 t2 (Point_2(p3.x(),p3.y()), Point_2(q3.x(),q3.y()));

      // convert intersection from Object_2 to Object_3
      // Note: there is not necessarily a spartial intersection,
      //       so all third components are faked!
      auto obj = intersection (s2,t2);
      if(! obj){
        return boost::none;
      }
      if (const Point_2* pi =  boost::get<Point_2>(&*obj))
      {
        return boost::make_optional(variant_type(Point_3(p2.x(),p2.y(),0)));
      }

      const Segment_2* si = boost::get<Segment_2>(&*obj);
      p2 = si->source();
      q2 = si->target();

      return boost::make_optional(variant_type(Segment_3(Point_3(p2.x(),p2.y(),0),
                                                         Point_3(q2.x(),q2.y(),0) ) ));
   }
};

template <class R>
struct Exact_intersect_xy_2 <R,Homogeneous_tag>
{  // Homogeneous
   typedef typename R::Point_2     Point_2;
   typedef typename R::Segment_2   Segment_2;

   typedef typename R::Point_3     Point_3;
   typedef typename R::Segment_3   Segment_3;

   typedef  boost::variant<Point_3, Segment_3> variant_type;

   boost::optional<variant_type> operator() (Segment_3 s3, Segment_3 t3)
   {  Point_2 p2, q2;
      Point_3 p3, q3;

      // convert Segment_3 to Segment_2
      p3 = s3.source();
      q3 = s3.target();
      Segment_2 s2 (Point_2(p3.hx(),p3.hy(),p3.hw()),
                    Point_2(q3.hx(),q3.hy(),q3.hw()));

      p3 = t3.source();
      q3 = t3.target();
      Segment_2 t2 (Point_2(p3.hx(),p3.hy(),p3.hw()),
                    Point_2(q3.hx(),q3.hy(),q3.hw()));

      // convert intersection from Object_2 to Object_3
      // Note: there is not necessarily a spartial intersection,
      //       so all third components are faked!
      auto obj = intersection (s2,t2);
      if(! obj){
        return boost::none;
      }
      if (const Point_2* pi =  boost::get<Point_2>(&*obj))
      {
        return boost::make_optional(variant_type(Point_3(p2.hx(),p2.hy(),0,p2.hw())));
      }

      const Segment_2* si = boost::get<Segment_2>(&*obj);
      p2 = si->source();
      q2 = si->target();

      return boost::make_optional(variant_type(Segment_3(Point_3 (p2.hx(),p2.hy(),0,p2.hw()),
                                                         Point_3 (q2.hx(),q2.hy(),0,q2.hw())) ));
   }
};

template <class R>
struct Exact_triangulation_euclidean_traits_xy_3
 : public CGAL::Projection_traits_xy_3<R>
{
   typedef CGAL::Exact_triangulation_euclidean_traits_xy_3<R> Traits;
   typedef CGAL::Projection_traits_xy_3<R>       Tet;
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

   Exact_triangulation_euclidean_traits_xy_3(){}
   Exact_triangulation_euclidean_traits_xy_3(
      const Exact_triangulation_euclidean_traits_xy_3&){}
   Exact_triangulation_euclidean_traits_xy_3 &operator=(
      const Exact_triangulation_euclidean_traits_xy_3&) {return *this;}

   typedef Exact_intersect_xy_2<Rp,typename Rp::Kernel_tag> Intersect_2;

   Intersect_2 intersect_2_object () const
      { return Intersect_2(); }
};

} //namespace CGAL
#endif // CGAL_NEF3_EXACT_TRIANGULATION_EUCLIDEAN_TRAITS_XY_3_H
