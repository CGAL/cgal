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

#include <CGAL/Intersection_traits.h>

#if !(CGAL_INTERSECTION_VERSION < 2)

namespace CGAL  {

CGAL_INTERSECTION_TRAITS_2(Line_3, Line_3, Point_3, Line_3)

CGAL_INTERSECTION_TRAITS_2(Line_3, Plane_3, Point_3, Line_3)
CGAL_INTERSECTION_TRAITS_2(Plane_3, Line_3, Point_3, Line_3)

CGAL_INTERSECTION_TRAITS_2(Line_3, Ray_3, Point_3, Ray_3)
CGAL_INTERSECTION_TRAITS_2(Ray_3, Line_3, Point_3, Ray_3)

CGAL_INTERSECTION_TRAITS_2(Line_3, Segment_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Segment_3, Line_3, Point_3, Segment_3)

CGAL_INTERSECTION_TRAITS_2(Line_3, Triangle_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Triangle_3, Line_3, Point_3, Segment_3)

CGAL_INTERSECTION_TRAITS_2(Plane_3, Plane_3, Line_3, Plane_3)

CGAL_INTERSECTION_TRAITS_2(Plane_3, Ray_3, Point_3, Ray_3)
CGAL_INTERSECTION_TRAITS_2(Ray_3, Plane_3, Point_3, Ray_3)

CGAL_INTERSECTION_TRAITS_2(Plane_3, Segment_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Segment_3, Plane_3, Point_3, Segment_3)

CGAL_INTERSECTION_TRAITS_2(Plane_3, Sphere_3, Point_3, Circle_3)
CGAL_INTERSECTION_TRAITS_2(Sphere_3, Plane_3, Point_3, Circle_3)

CGAL_INTERSECTION_TRAITS_3(Plane_3, Triangle_3, Point_3, Segment_3, Triangle_3)
CGAL_INTERSECTION_TRAITS_3(Triangle_3, Plane_3, Point_3, Segment_3, Triangle_3)

CGAL_INTERSECTION_TRAITS_3(Ray_3, Ray_3, Point_3, Ray_3, Segment_3)

CGAL_INTERSECTION_TRAITS_2(Ray_3, Segment_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Segment_3, Ray_3, Point_3, Segment_3)

CGAL_INTERSECTION_TRAITS_2(Ray_3, Triangle_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Triangle_3, Ray_3, Point_3, Segment_3)

CGAL_INTERSECTION_TRAITS_2(Segment_3, Segment_3, Point_3, Segment_3)

CGAL_INTERSECTION_TRAITS_2(Segment_3, Triangle_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Triangle_3, Segment_3, Point_3, Segment_3)

CGAL_INTERSECTION_TRAITS_3(Sphere_3, Sphere_3, Point_3, Circle_3, Sphere_3)

template<typename K>
struct Intersection_traits<K, typename K::Triangle_3, typename K::Triangle_3>  {
  typedef typename 
  boost::variant< typename K::Point_3, typename K::Segment_3, typename K::Triangle_3,
                  typename std::vector< typename K::Point_3 > > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};


// !!! undocumented !!!

// Segment_3 Iso_cuboid_3
CGAL_INTERSECTION_TRAITS_2(Segment_3, Iso_cuboid_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Iso_cuboid_3, Segment_3, Point_3, Segment_3)

// Line_3 Iso_cuboid_3
CGAL_INTERSECTION_TRAITS_2(Line_3, Iso_cuboid_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Iso_cuboid_3, Line_3, Point_3, Segment_3)

// Ray_3 Iso_cuboid_3
CGAL_INTERSECTION_TRAITS_2(Ray_3, Iso_cuboid_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Iso_cuboid_3, Ray_3, Point_3, Segment_3)

// Iso_cuboid_3 Iso_cuboid_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Iso_cuboid_3>  {
  typedef typename 
  boost::variant< typename K::Iso_cuboid_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};


// Intersections with BBox returns the same for Ray_3, Line_3 and Segment_3
// Bbox_3 Line_3
template<typename K>
struct Intersection_traits<K, CGAL::Bbox_3, typename K::Line_3>  {
  typedef typename 
  boost::variant< typename K::Segment_3, typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Line_3, CGAL::Bbox_3> : 
    public Intersection_traits<K, CGAL::Bbox_3, typename K::Line_3> {};

template<typename K>
inline
typename K::Intersect_3::template Result<CGAL::Bbox_3, typename K::Line_3>
intersection(const CGAL::Bbox_3& a, typename K::Line_3 const& b) {
  return K().intersect_3_object()(a, b);
}

template<typename K>
inline
typename K::Intersect_3::template Result<CGAL::Bbox_3, typename K::Line_3>
intersection(typename K::Line_3 const& b, const CGAL::Bbox_3& a) {
  return K().intersect_3_object()(b, a);
}

// Bbox_3 Segment_3
template<typename K>
struct Intersection_traits<K, typename K::Segment_3, CGAL::Bbox_3> : 
    public Intersection_traits<K, CGAL::Bbox_3, typename K::Line_3> {};

template<typename K>
struct Intersection_traits<K, CGAL::Bbox_3, typename K::Segment_3> : 
    public Intersection_traits<K, CGAL::Bbox_3, typename K::Line_3> {};

template<typename K>
inline
typename K::Intersect_3::template Result<CGAL::Bbox_3, typename K::Segment_3>
intersection(typename K::Segment_3 const& b, const CGAL::Bbox_3& a) {
  return K().intersect_3_object()(b, a);
}

template<typename K>
inline
typename K::Intersect_3::template Result<CGAL::Bbox_3, typename K::Segment_3>
intersection(typename CGAL::Bbox_3 const& b, typename K::Segment_3 const& a) {
  return K().intersect_3_object()(b, a);
}

// Bbox_3 Ray_3
template<typename K>
struct Intersection_traits<K, typename K::Ray_3, CGAL::Bbox_3> : 
    public Intersection_traits<K, CGAL::Bbox_3, typename K::Line_3> {};

template<typename K>
struct Intersection_traits<K, CGAL::Bbox_3, typename K::Ray_3> : 
    public Intersection_traits<K, CGAL::Bbox_3, typename K::Line_3> {};

template<typename K>
inline
typename K::Intersect_3::template Result<CGAL::Bbox_3, typename K::Ray_3>
intersection(typename K::Ray_3 const& b, const CGAL::Bbox_3& a) {
  return K().intersect_3_object()(b, a);
}

template<typename K>
inline
typename K::Intersect_3::template Result<CGAL::Bbox_3, typename K::Ray_3>
intersection(const CGAL::Bbox_3& a, typename K::Ray_3 const& b) {
  return K().intersect_3_object()(b, a);
}

// Intersection_traits for the circular kernel

// Sphere_3 Line_3
template<typename K>
struct Intersection_traits<K, typename K::Sphere_3, typename K::Line_3> {
  typedef boost::variant< std::pair< typename K::Circular_arc_point_3, unsigned int > > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Line_3, typename K::Sphere_3>
  : public Intersection_traits<K, typename K::Sphere_3, typename K::Line_3> {};

// Circle_3 Plane_3
template<typename K>
struct Intersection_traits<K, typename K::Circle_3, typename K::Plane_3> {
  typedef typename boost::variant< std::pair< typename K::Circular_arc_point_3,
                                              unsigned int >,
                                   typename K::Circle_3 > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Plane_3, typename K::Circle_3>
  : public Intersection_traits<K, typename K::Circle_3, typename K::Plane_3> {};

// Circle_3 Sphere_3
template<typename K>
struct Intersection_traits<K, typename K::Circle_3, typename K::Sphere_3> {
  typedef typename boost::variant< std::pair< typename K::Circular_arc_point_3,
                                              unsigned int >,
                                   typename K::Circle_3 > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Sphere_3, typename K::Circle_3>
  : public Intersection_traits<K, typename K::Circle_3, typename K::Sphere_3> {};

// Circle_3 Circle_3
template<typename K>
struct Intersection_traits<K, typename K::Circle_3, typename K::Circle_3> {
  typedef typename boost::variant< std::pair <typename K::Circular_arc_point_3, unsigned int >, 
                                   typename K::Circle_3 > result_type;
};

// Circle_3 Circle_3
template<typename K>
struct Intersection_traits<K, typename K::Circle_3, typename K::Line_3> {
  typedef typename boost::variant< 
    std::pair <typename K::Circular_arc_point_3, unsigned int > > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Line_3, typename K::Circle_3>
  : public Intersection_traits<K, typename K::Circle_3, typename K::Line_3> {};

// Circular_arc_3 Circular_arc_3 
template<typename K>
struct Intersection_traits<K, typename K::Circular_arc_3, typename K::Circular_arc_3> {
  typedef typename boost::variant< 
    typename K::Circle_3,
    std::pair <typename K::Circular_arc_point_3, unsigned int >,
    typename K::Circular_arc_3 > result_type;
};

// Circular_arc_3 Plane_3 
template<typename K>
struct Intersection_traits<K, typename K::Circular_arc_3, typename K::Plane_3> {
  typedef typename boost::variant< 
    std::pair <typename K::Circular_arc_point_3, unsigned int >,
    typename K::Circular_arc_3 > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Plane_3, typename K::Circular_arc_3>
  : public Intersection_traits<K, typename K::Circular_arc_3, typename K::Plane_3> {};


template<typename K>
struct Intersection_traits<K, typename K::Line_arc_3, typename K::Line_arc_3> {
  typedef typename boost::variant< 
    std::pair <typename K::Circular_arc_point_3, unsigned int >,
    typename K::Line_arc_3 > result_type;
};


} // namespace

#endif // !(CGAL_INTERSECTION_VERSION < 2)

#endif /* CGAL_INTERSECTION_TRAITS_3_H */

