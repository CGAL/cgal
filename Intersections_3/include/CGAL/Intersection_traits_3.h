// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Philipp Möller

#ifndef CGAL_INTERSECTION_TRAITS_3_H
#define CGAL_INTERSECTION_TRAITS_3_H

// to reuse the macros
#include <CGAL/Intersection_traits_2.h>
#include <CGAL/Bbox_3.h>

namespace CGAL  {

template<typename, typename, typename>
struct Intersection_traits_3;

template<typename K, typename A, typename B>
struct IT3 : public Intersection_traits_3<K, A, B> {};

CGAL_INTERSECTION_TRAITS_DIM_2(Line_3, Line_3, Point_3, Line_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Line_3, Plane_3, Point_3, Line_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Plane_3, Line_3, Point_3, Line_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Line_3, Ray_3, Point_3, Ray_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Ray_3, Line_3, Point_3, Ray_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Line_3, Segment_3, Point_3, Segment_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Segment_3, Line_3, Point_3, Segment_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Line_3, Triangle_3, Point_3, Segment_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Triangle_3, Line_3, Point_3, Segment_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Plane_3, Plane_3, Line_3, Plane_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Plane_3, Ray_3, Point_3, Ray_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Ray_3, Plane_3, Point_3, Ray_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Plane_3, Segment_3, Point_3, Segment_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Segment_3, Plane_3, Point_3, Segment_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Plane_3, Sphere_3, Point_3, Circle_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Sphere_3, Plane_3, Point_3, Circle_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_3(Plane_3, Triangle_3, Point_3, Segment_3, Triangle_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_3(Triangle_3, Plane_3, Point_3, Segment_3, Triangle_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_3(Ray_3, Ray_3, Point_3, Ray_3, Segment_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Ray_3, Segment_3, Point_3, Segment_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Segment_3, Ray_3, Point_3, Segment_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Ray_3, Triangle_3, Point_3, Segment_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Triangle_3, Ray_3, Point_3, Segment_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Segment_3, Segment_3, Point_3, Segment_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_2(Segment_3, Triangle_3, Point_3, Segment_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Triangle_3, Segment_3, Point_3, Segment_3, 3)

CGAL_INTERSECTION_TRAITS_DIM_3(Sphere_3, Sphere_3, Point_3, Circle_3, Sphere_3, 3)

template<typename K>
struct Intersection_traits_3<K, typename K::Triangle_3, typename K::Triangle_3>  {
  typedef typename 
  boost::variant< typename K::Point_3, typename K::Segment_3, typename K::Triangle_3,
                  typename std::vector< typename K::Point_3 > > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// !!! undocumented !!!

// Segment_3 Iso_cuboid_3
CGAL_INTERSECTION_TRAITS_DIM_2(Segment_3, Iso_cuboid_3, Point_3, Segment_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Iso_cuboid_3, Segment_3, Point_3, Segment_3, 3)

// Line_3 Iso_cuboid_3
CGAL_INTERSECTION_TRAITS_DIM_2(Line_3, Iso_cuboid_3, Point_3, Segment_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Iso_cuboid_3, Line_3, Point_3, Segment_3, 3)

// Ray_3 Iso_cuboid_3
CGAL_INTERSECTION_TRAITS_DIM_2(Ray_3, Iso_cuboid_3, Point_3, Segment_3, 3)
CGAL_INTERSECTION_TRAITS_DIM_2(Iso_cuboid_3, Ray_3, Point_3, Segment_3, 3)

// Iso_cuboid_3 Iso_cuboid_3, variant of one
template<typename K>
struct Intersection_traits_3<K, typename K::Iso_cuboid_3, typename K::Iso_cuboid_3>  {
  typedef typename 
  boost::variant< typename K::Iso_cuboid_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Intersections with BBox returns the same for Ray_3, Line_3 and Segment_3
// Bbox_3 Line_3
template<typename K>
struct Intersection_traits_3<K, CGAL::Bbox_3, typename K::Line_3>  {
  typedef typename 
  boost::variant< typename K::Segment_3, typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

template<typename K>
struct Intersection_traits_3<K, typename K::Line_3, CGAL::Bbox_3> : 
    public Intersection_traits_3<K, CGAL::Bbox_3, typename K::Line_3> {};

// Bbox_3 Segment_3
template<typename K>
struct Intersection_traits_3<K, typename K::Segment_3, CGAL::Bbox_3> : 
    public Intersection_traits_3<K, CGAL::Bbox_3, typename K::Line_3> {};

template<typename K>
struct Intersection_traits_3<K, CGAL::Bbox_3, typename K::Segment_3> : 
    public Intersection_traits_3<K, CGAL::Bbox_3, typename K::Line_3> {};

// Bbox_3 Ray_3
template<typename K>
struct Intersection_traits_3<K, typename K::Ray_3, CGAL::Bbox_3> : 
    public Intersection_traits_3<K, CGAL::Bbox_3, typename K::Line_3> {};

template<typename K>
struct Intersection_traits_3<K, CGAL::Bbox_3, typename K::Ray_3> : 
    public Intersection_traits_3<K, CGAL::Bbox_3, typename K::Line_3> {};


} // namespace

#endif /* CGAL_INTERSECTION_TRAITS_3_H */

