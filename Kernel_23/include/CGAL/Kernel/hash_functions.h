// Copyright (c) 2019,2025
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
//
// Test file: test/Kernel_23/test_hash_functions.cpp

#ifndef CGAL_KERNEL_HASH_FUNCTIONS_H
#define CGAL_KERNEL_HASH_FUNCTIONS_H

#include <boost/functional/hash.hpp>
#include <type_traits>

#include <CGAL/representation_tags.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Weighted_point_2.h>
#include <CGAL/Weighted_point_3.h>


namespace CGAL
{

using boost::hash_value;

template <typename K, typename = void>
inline constexpr bool has_rep_tag_v = false;

template <typename K>
inline constexpr bool has_rep_tag_v<K, std::void_t<typename K::Rep_tag>> = true;

template <typename K, typename = void>
struct Rep_tag {
  using type = void;
};

template <typename K>
struct Rep_tag<K, std::enable_if_t<has_rep_tag_v<K>>> {
  using type = typename K::Rep_tag;
};

template <typename K>
using Rep_tag_t = typename Rep_tag<K>::type;

template <typename K>
inline constexpr bool is_Cartesian_v = std::is_same<Rep_tag_t<K>, Cartesian_tag>::value;

template <typename K, typename = void>
struct Is_kernel_hashable : public std::false_type {};

template <typename K>
struct Is_kernel_hashable<K, std::void_t<decltype(hash_value(std::declval<typename K::FT>()))>> : public std::true_type {};

template <typename K>
inline constexpr bool is_kernel_hashable_v = Is_kernel_hashable<K>::value;

template <typename K, typename T>
using enable_if_Cartesian_and_hashable_t =
    std::enable_if_t<is_Cartesian_v<K> && is_kernel_hashable_v<K>, T>;

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Aff_transformation_2<K>& transform)
{
  std::size_t result = hash_value(transform.cartesian(0,0));
  for(int i=0; i < 3; ++i)
    for(int j=0; j < 3; ++j)
      // Skip (0,0) as it was already used to initialize the hash
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
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Circle_2<K>& circle)
{
  std::size_t result = hash_value(circle.center());
  boost::hash_combine(result, hash_value(circle.squared_radius()));
  boost::hash_combine(result, hash_value(circle.orientation()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Iso_rectangle_2<K>& iso_rectangle)
{
  std::size_t result = hash_value((iso_rectangle.min)());
  boost::hash_combine(result, hash_value((iso_rectangle.max)()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Point_2<K>& point)
{
  std::size_t result = hash_value(point.x());
  boost::hash_combine(result, hash_value(point.y()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Segment_2<K>& segment)
{
  std::size_t result = hash_value(segment.source());
  boost::hash_combine(result, hash_value(segment.target()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Vector_2<K>& vector)
{
  std::size_t result = hash_value(vector.x());
  boost::hash_combine(result, hash_value(vector.y()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Weighted_point_2<K>& weighed_point)
{
  std::size_t result = hash_value(weighed_point.point());
  boost::hash_combine(result, hash_value(weighed_point.weight()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Aff_transformation_3<K>& transform)
{
  std::size_t result = hash_value(transform.cartesian(0,0));
  for(int i = 0; i < 3; ++i)
    for(int j = (i == 0 ? 1 : 0); j < 4; ++j)
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
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Iso_cuboid_3<K>& iso_cuboid)
{
  std::size_t result = hash_value((iso_cuboid.min)());
  boost::hash_combine(result, hash_value((iso_cuboid.max)()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Point_3<K>& point)
{
  std::size_t result = hash_value(point.x());
  boost::hash_combine(result, hash_value(point.y()));
  boost::hash_combine(result, hash_value(point.z()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Segment_3<K>& segment)
{
  std::size_t result = hash_value(segment.source());
  boost::hash_combine(result, hash_value(segment.target()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Sphere_3<K>& sphere)
{
  std::size_t result = hash_value(sphere.center());
  boost::hash_combine(result, hash_value(sphere.squared_radius()));
  boost::hash_combine(result, hash_value(sphere.orientation()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Vector_3<K>& vector)
{
  std::size_t result = hash_value(vector.x());
  boost::hash_combine(result, hash_value(vector.y()));
  boost::hash_combine(result, hash_value(vector.z()));
  return result;
}

template <typename K>
inline enable_if_Cartesian_and_hashable_t<K, std::size_t>
hash_value (const Weighted_point_3<K>& weighed_point)
{
  std::size_t result = hash_value(weighed_point.point());
  boost::hash_combine(result, hash_value(weighed_point.weight()));
  return result;
}

struct Forward_to_hash_value {
  template <typename T>
  std::size_t operator()(T&& t) const {
    using boost::hash_value;
    return hash_value(std::forward<T>(t));
  }
};

template <typename K, typename = void>
struct Maybe_forward_to_hash_value {
  Maybe_forward_to_hash_value() = delete;
  Maybe_forward_to_hash_value(const Maybe_forward_to_hash_value&) = delete;
};

template <typename K>
struct Maybe_forward_to_hash_value<K, std::enable_if_t<is_kernel_hashable_v<K>>>
    : public Forward_to_hash_value {};


} //namespace CGAL

// overloads of std::hash used for using std::unordered_[set/map] on CGAL Kernel objects
namespace std {
template <typename K> struct hash<CGAL::Aff_transformation_2<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Circle_2<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Iso_rectangle_2<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Point_2<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Segment_2<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Vector_2<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Weighted_point_2<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <> struct hash<CGAL::Bbox_2> : CGAL::Forward_to_hash_value {};
template <typename K> struct hash<CGAL::Aff_transformation_3<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Iso_cuboid_3<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Point_3<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Segment_3<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Sphere_3<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Vector_3<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <typename K> struct hash<CGAL::Weighted_point_3<K>> : CGAL::Maybe_forward_to_hash_value<K> {};
template <> struct hash<CGAL::Bbox_3> : CGAL::Forward_to_hash_value {};
} // namespace std


#endif  // CGAL_KERNEL_HASH_FUNCTIONS_H
