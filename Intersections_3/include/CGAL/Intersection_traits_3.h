// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Philipp Möller

#ifndef CGAL_INTERSECTION_TRAITS_3_H
#define CGAL_INTERSECTION_TRAITS_3_H

#include <CGAL/Intersection_traits.h>
#include <vector>

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

CGAL_INTERSECTION_TRAITS_2(Line_3, Tetrahedron_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Tetrahedron_3, Line_3, Point_3, Segment_3)

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

CGAL_INTERSECTION_TRAITS_2(Ray_3, Tetrahedron_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Tetrahedron_3, Ray_3, Point_3, Segment_3)

CGAL_INTERSECTION_TRAITS_2(Segment_3, Segment_3, Point_3, Segment_3)

CGAL_INTERSECTION_TRAITS_2(Segment_3, Triangle_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Triangle_3, Segment_3, Point_3, Segment_3)

CGAL_INTERSECTION_TRAITS_2(Segment_3, Tetrahedron_3, Point_3, Segment_3)
CGAL_INTERSECTION_TRAITS_2(Tetrahedron_3, Segment_3, Point_3, Segment_3)

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

// Bbox_3 Segment_3
template<typename K>
struct Intersection_traits<K, typename K::Segment_3, CGAL::Bbox_3> :
    public Intersection_traits<K, CGAL::Bbox_3, typename K::Line_3> {};

template<typename K>
struct Intersection_traits<K, CGAL::Bbox_3, typename K::Segment_3> :
    public Intersection_traits<K, CGAL::Bbox_3, typename K::Line_3> {};

// Bbox_3 Ray_3
template<typename K>
struct Intersection_traits<K, typename K::Ray_3, CGAL::Bbox_3> :
    public Intersection_traits<K, CGAL::Bbox_3, typename K::Line_3> {};

template<typename K>
struct Intersection_traits<K, CGAL::Bbox_3, typename K::Ray_3> :
    public Intersection_traits<K, CGAL::Bbox_3, typename K::Line_3> {};

// Point_3 Iso_cuboid_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Point_3, typename K::Iso_cuboid_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Bbox_3 Point_3, variant of one
template<typename K>
struct Intersection_traits<K, Bbox_3, typename K::Point_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

  template<typename K>
struct Intersection_traits<K, typename K::Point_3, Bbox_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};


// Bbox_3 Iso_cuboid_3, variant of 1
template<typename K>
struct Intersection_traits<K, CGAL::Bbox_3, typename K::Iso_cuboid_3>  {
  typedef typename
  boost::variant<  typename K::Iso_cuboid_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Iso_cuboid_3, CGAL::Bbox_3>
    : public Intersection_traits<K, CGAL::Bbox_3, typename K::Iso_cuboid_3> {};

// Iso_cuboid_3 Point_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Point_3>  {
  typedef typename
    boost::variant< typename K::Point_3 > variant_type;
    typedef typename boost::optional< variant_type > result_type;
};

// Iso_cuboid_3 Plane_3, variant of 4
template<typename K>
struct Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Plane_3>  {
  typedef typename
  boost::variant< typename K::Point_3, typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;

  typedef typename boost::optional< variant_type > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Plane_3, typename K::Iso_cuboid_3>  {
  typedef typename
  boost::variant< typename K::Point_3, typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;

  typedef typename boost::optional< variant_type > result_type;
};

// Iso_cuboid_3 Triangle_3, variant of 4
template<typename K>
struct Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Triangle_3>  {
  typedef typename
  boost::variant< typename K::Point_3, typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;

  typedef typename boost::optional< variant_type > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Triangle_3, typename K::Iso_cuboid_3>  {
  typedef typename
  boost::variant< typename K::Point_3, typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;

  typedef typename boost::optional< variant_type > result_type;
};

// Bbox_3 Plane_3, variant of 4
template<typename K>
struct Intersection_traits<K, typename CGAL::Bbox_3, typename K::Plane_3>  {
  typedef typename
  boost::variant< typename K::Point_3, typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;

  typedef typename boost::optional< variant_type > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Plane_3, typename CGAL::Bbox_3>  {
  typedef typename
  boost::variant< typename K::Point_3, typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;

  typedef typename boost::optional< variant_type > result_type;
};

// Bbox_3 Triangle_3, variant of 4
template<typename K>
struct Intersection_traits<K, typename CGAL::Bbox_3, typename K::Triangle_3>  {
  typedef typename
  boost::variant< typename K::Point_3, typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;

  typedef typename boost::optional< variant_type > result_type;
};

template<typename K>
struct Intersection_traits<K, typename K::Triangle_3, typename CGAL::Bbox_3>  {
  typedef typename
  boost::variant< typename K::Point_3, typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;

  typedef typename boost::optional< variant_type > result_type;
};

// Point_3 Line_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Point_3, typename K::Line_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Line_3 Point_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Line_3, typename K::Point_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Point_3 Ray_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Point_3, typename K::Ray_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Ray_3 Point_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Ray_3, typename K::Point_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Point_3 Segment_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Point_3, typename K::Segment_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Segment_3 Point_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Segment_3, typename K::Point_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Point_3 Point_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Point_3, typename K::Point_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Point_3 Plane_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Point_3, typename K::Plane_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Plane_3 Point_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Plane_3, typename K::Point_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Point_3 Triangle_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Point_3, typename K::Triangle_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Triangle_3 Point_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Triangle_3, typename K::Point_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Point_3 Tetrahedron_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Point_3, typename K::Tetrahedron_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Tetrahedron_3 Point_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Tetrahedron_3, typename K::Point_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Point_3 Sphere_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Point_3, typename K::Sphere_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Sphere_3 Point_3, variant of one
template<typename K>
struct Intersection_traits<K, typename K::Sphere_3, typename K::Point_3>  {
  typedef typename
  boost::variant< typename K::Point_3 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

//Tetrahedron_3 Plane_3, variant of 4
template<class K>
struct Intersection_traits<K, typename K::Tetrahedron_3, typename K::Plane_3>
{
  typedef typename
  boost::variant< typename K::Point_3 , typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

//Plane_3 Tetrahedron_3, variant of 4
template<class K>
struct Intersection_traits<K, typename K::Plane_3,typename K::Tetrahedron_3>
{
  typedef typename
  boost::variant< typename K::Point_3 , typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

//Triangle_3 Tetrahedron_3, variant of 4
template<class K>
struct Intersection_traits<K, typename K::Triangle_3,typename K::Tetrahedron_3>
{
  typedef typename
  boost::variant< typename K::Point_3 , typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};

// Tetrahedron_3 Triangle_3, variant of 4
template<class K>
struct Intersection_traits<K, typename K::Tetrahedron_3,typename K::Triangle_3>
{
  typedef typename
  boost::variant< typename K::Point_3 , typename K::Segment_3,
                  typename K::Triangle_3, std::vector<typename K::Point_3> > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};


} // namespace

#endif /* CGAL_INTERSECTION_TRAITS_3_H */

