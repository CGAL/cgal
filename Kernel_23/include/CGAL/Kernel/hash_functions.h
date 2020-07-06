// Copyright (c) 2019
// GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KERNEL_HASH_FUNCTIONS_H
#define CGAL_KERNEL_HASH_FUNCTIONS_H

#include <boost/functional/hash.hpp>
#include <type_traits>

namespace CGAL
{

using boost::hash_value;

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Aff_transformation_2<K>& transform)
{
  std::size_t result = hash_value(transform.cartesian(0,0));
  for(int i=0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      if (!(i == 0 && j == 0))
        boost::hash_combine(result, hash_value(transform.cartesian(i,j)));
  return result;
}

inline std::size_t
hash_value (const Bbox_2& bbox)
{
  std::size_t result = hash_value(bbox.xmin());
  boost::hash_combine(result, hash_value(bbox.xmax()));
  boost::hash_combine(result, hash_value(bbox.ymin()));
  boost::hash_combine(result, hash_value(bbox.ymax()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Circle_2<K>& circle)
{
  std::size_t result = hash_value(circle.center());
  boost::hash_combine(result, hash_value(circle.squared_radius()));
  boost::hash_combine(result, hash_value(circle.orientation()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Iso_rectangle_2<K>& iso_rectangle)
{
  std::size_t result = hash_value((iso_rectangle.min)());
  boost::hash_combine(result, hash_value((iso_rectangle.max)()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Point_2<K>& point)
{
  std::size_t result = hash_value(point.x());
  boost::hash_combine(result, hash_value(point.y()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Segment_2<K>& segment)
{
  std::size_t result = hash_value(segment.source());
  boost::hash_combine(result, hash_value(segment.target()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Vector_2<K>& vector)
{
  std::size_t result = hash_value(vector.x());
  boost::hash_combine(result, hash_value(vector.y()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Weighted_point_2<K>& weighed_point)
{
  std::size_t result = hash_value(weighed_point.point());
  boost::hash_combine(result, hash_value(weighed_point.weight()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Aff_transformation_3<K>& transform)
{
  std::size_t result = hash_value(transform.cartesian(0,0));
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 4; ++j)
      if (!(i == 0 && j == 0))
        boost::hash_combine(result, hash_value(transform.cartesian(i,j)));
  return result;
}

inline std::size_t
hash_value (const Bbox_3& bbox)
{
  std::size_t result = hash_value(bbox.xmin());
  boost::hash_combine(result, hash_value(bbox.xmax()));
  boost::hash_combine(result, hash_value(bbox.ymin()));
  boost::hash_combine(result, hash_value(bbox.ymax()));
  boost::hash_combine(result, hash_value(bbox.zmin()));
  boost::hash_combine(result, hash_value(bbox.zmax()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Iso_cuboid_3<K>& iso_cuboid)
{
  std::size_t result = hash_value((iso_cuboid.min)());
  boost::hash_combine(result, hash_value((iso_cuboid.max)()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Point_3<K>& point)
{
  std::size_t result = hash_value(point.x());
  boost::hash_combine(result, hash_value(point.y()));
  boost::hash_combine(result, hash_value(point.z()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Segment_3<K>& segment)
{
  std::size_t result = hash_value(segment.source());
  boost::hash_combine(result, hash_value(segment.target()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Sphere_3<K>& sphere)
{
  std::size_t result = hash_value(sphere.center());
  boost::hash_combine(result, hash_value(sphere.squared_radius()));
  boost::hash_combine(result, hash_value(sphere.orientation()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Vector_3<K>& vector)
{
  std::size_t result = hash_value(vector.x());
  boost::hash_combine(result, hash_value(vector.y()));
  boost::hash_combine(result, hash_value(vector.z()));
  return result;
}

template <typename K>
inline std::enable_if_t<std::is_same<typename K::Rep_tag, Cartesian_tag>::value, std::size_t>
hash_value (const Weighted_point_3<K>& weighed_point)
{
  std::size_t result = hash_value(weighed_point.point());
  boost::hash_combine(result, hash_value(weighed_point.weight()));
  return result;
}

} //namespace CGAL

// overloads of std::hash used for using std::unordered_[set/map] on CGAL Kernel objects
namespace std
{

template <typename K> struct hash<CGAL::Aff_transformation_2<K> > {
  std::size_t operator() (const CGAL::Aff_transformation_2<K>& transform) const   {
    return CGAL::hash_value<K> (transform);
  }
};
template <> struct hash<CGAL::Bbox_2> {
  std::size_t operator() (const CGAL::Bbox_2& bbox) const   {
    return CGAL::hash_value (bbox);
  }
};
template <typename K> struct hash<CGAL::Circle_2<K> > {
  std::size_t operator() (const CGAL::Circle_2<K>& circle) const   {
    return CGAL::hash_value<K> (circle);
  }
};
template <typename K> struct hash<CGAL::Iso_rectangle_2<K> > {
  std::size_t operator() (const CGAL::Iso_rectangle_2<K>& iso_rectangle) const   {
    return CGAL::hash_value<K> (iso_rectangle);
  }
};
template <typename K> struct hash<CGAL::Point_2<K> > {
  std::size_t operator() (const CGAL::Point_2<K>& point) const   {
    return CGAL::hash_value<K> (point);
  }
};
template <typename K> struct hash<CGAL::Segment_2<K> > {
  std::size_t operator() (const CGAL::Segment_2<K>& segment) const   {
    return CGAL::hash_value<K> (segment);
  }
};
template <typename K> struct hash<CGAL::Vector_2<K> > {
  std::size_t operator() (const CGAL::Vector_2<K>& vector) const   {
    return CGAL::hash_value<K> (vector);
  }
};
template <typename K> struct hash<CGAL::Weighted_point_2<K> > {
  std::size_t operator() (const CGAL::Weighted_point_2<K>& weighted_point) const   {
    return CGAL::hash_value<K> (weighted_point);
  }
};
template <typename K> struct hash<CGAL::Aff_transformation_3<K> > {
  std::size_t operator() (const CGAL::Aff_transformation_3<K>& transform) const   {
    return CGAL::hash_value<K> (transform);
  }
};
template <> struct hash<CGAL::Bbox_3> {
  std::size_t operator() (const CGAL::Bbox_3& bbox) const   {
    return CGAL::hash_value (bbox);
  }
};
template <typename K> struct hash<CGAL::Iso_cuboid_3<K> > {
  std::size_t operator() (const CGAL::Iso_cuboid_3<K>& iso_cuboid) const   {
    return CGAL::hash_value<K> (iso_cuboid);
  }
};
template <typename K> struct hash<CGAL::Point_3<K> > {
  std::size_t operator() (const CGAL::Point_3<K>& point) const   {
    return CGAL::hash_value<K> (point);
  }
};
template <typename K> struct hash<CGAL::Segment_3<K> > {
  std::size_t operator() (const CGAL::Segment_3<K>& segment) const   {
    return CGAL::hash_value<K> (segment);
  }
};
template <typename K> struct hash<CGAL::Sphere_3<K> > {
  std::size_t operator() (const CGAL::Sphere_3<K>& sphere) const   {
    return CGAL::hash_value<K> (sphere);
  }
};
template <typename K> struct hash<CGAL::Vector_3<K> > {
  std::size_t operator() (const CGAL::Vector_3<K>& vector) const   {
    return CGAL::hash_value<K> (vector);
  }
};
template <typename K> struct hash<CGAL::Weighted_point_3<K> > {
  std::size_t operator() (const CGAL::Weighted_point_3<K>& weighted_point) const   {
    return CGAL::hash_value<K> (weighted_point);
  }
};

}

#endif  // CGAL_KERNEL_HASH_FUNCTIONS_H
