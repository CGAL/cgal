// Copyright (c) 2022  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CARTESIAN_IS_TRIVIAL_CONSTRUCTION_H
#define CGAL_CARTESIAN_IS_TRIVIAL_CONSTRUCTION_H

#include <CGAL/type_traits/is_iterator.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/tags.h>
#include <CGAL/enum.h>
#include <CGAL/type_traits.h>

namespace CGAL {
namespace CartesianFunctors {

template <class Construction, typename ...Args>
struct Is_trivial_construction_base
{
  // If the return type of the construction, with the specified arguments, is a
  // reference or an iterator, them the construction is necessarily trivial.
  using return_type = decltype(std::declval<Construction>()(std::declval<Args>()...));
  enum {
    value = std::is_reference<return_type>::value || CGAL::is_iterator<return_type>::value
  };
};

template <class Construction, typename ...Args>
struct Is_trivial_construction : public Is_trivial_construction_base<Construction, Args...>
{};

template <typename K, typename... Args>
struct Is_trivial_construction<CommonKernelFunctors::Assign_2<K>, Args...>
    : public Tag_true
{};

template <typename K, typename... Args>
struct Is_trivial_construction<CommonKernelFunctors::Assign_3<K>, Args...>
    : public Tag_true
{};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_point_2<K>, Args...>
{
  typedef typename K::RT         RT;

  static Tag_true trivial(Return_base_tag, Origin);
  static Tag_true trivial(Return_base_tag, RT, RT);
  static Tag_true trivial(Origin);
  static Tag_true trivial(RT, RT);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_point_2<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_point_3<K>, Args...>
{
  typedef typename K::RT         RT;

  static Tag_true trivial(Return_base_tag, Origin);
  static Tag_true trivial(Return_base_tag, RT, RT, RT);
  static Tag_true trivial(Origin);
  static Tag_true trivial(RT, RT, RT);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_point_3<K>, Args...>::value
   };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_weighted_point_2<K>, Args...>
{
  typedef typename K::FT         FT;
  typedef typename K::Point_2    Point_2;

  static Tag_true trivial(Return_base_tag, Origin);
  static Tag_true trivial(Return_base_tag, Point_2, FT);
  static Tag_true trivial(Return_base_tag, FT, FT);
  static Tag_true trivial(Point_2, FT);
  static Tag_true trivial(Origin);
  static Tag_true trivial(FT, FT);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_weighted_point_2<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_weighted_point_3<K>, Args...>
{
  typedef typename K::FT         FT;
  typedef typename K::Point_3    Point_3;

  static Tag_true trivial(Return_base_tag, Origin);
  static Tag_true trivial(Return_base_tag, Point_3, FT);
  static Tag_true trivial(Return_base_tag, FT, FT, FT);
  static Tag_true trivial(Point_3, FT);
  static Tag_true trivial(Origin);
  static Tag_true trivial(FT, FT, FT);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_weighted_point_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_vector_2<K>, Args...>
{
  typedef typename K::RT         RT;
  typedef typename K::Point_2    Point_2;

  static Tag_true trivial(Return_base_tag, Null_vector);
  static Tag_true trivial(Return_base_tag, Origin, Point_2);
  static Tag_true trivial(Return_base_tag, RT, RT);
  static Tag_true trivial(Null_vector);
  static Tag_true trivial(Origin, Point_2);
  static Tag_true trivial(RT, RT);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_vector_2<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_vector_3<K>, Args...>
{
  typedef typename K::RT         RT;
  typedef typename K::Point_3    Point;
  typedef typename K::Line_3     Line;

  static Tag_true trivial(Return_base_tag, Null_vector);
  static Tag_true trivial(Return_base_tag, Origin, Point);
  static Tag_true trivial(Return_base_tag, RT, RT, RT);
  static Tag_true trivial(Return_base_tag, Line);
  static Tag_true trivial(Null_vector);
  static Tag_true trivial(Origin, Point);
  static Tag_true trivial(RT, RT, RT);
  static Tag_true trivial(Line);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_vector_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_direction_2<K>, Args...>
{
  typedef typename K::RT         RT;
  typedef typename K::Vector_2      Vector;

  static Tag_true trivial(Return_base_tag, RT, RT);
  static Tag_true trivial(Return_base_tag, Vector);
  static Tag_true trivial(RT, RT);
  static Tag_true trivial(Vector);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_direction_2<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_direction_3<K>, Args...>
{
  typedef typename K::RT         RT;
  typedef typename K::Vector_3   Vector;

  static Tag_true trivial(Return_base_tag, RT, RT, RT);
  static Tag_true trivial(Return_base_tag, Vector);
  static Tag_true trivial(RT, RT, RT);
  static Tag_true trivial(Vector);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_direction_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_line_2<K>, Args...>
{
  typedef typename K::RT         RT;

  static Tag_true trivial(Return_base_tag, RT, RT, RT);
  static Tag_true trivial(RT, RT, RT);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_line_2<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_line_3<K>, Args...>
{
  typedef typename K::Point_3     Point;
  typedef typename K::Vector_3    Vector;
  typedef typename K::Direction_3 Direction;

  static Tag_true trivial(Return_base_tag, Point, Vector);
  static Tag_true trivial(Return_base_tag, Point, Direction);
  static Tag_true trivial(Point, Vector);
  static Tag_true trivial(Point, Direction);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_line_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_segment_2<K>, Args...>
{
  typedef typename K::Point_2    Point;

  static Tag_true trivial(Return_base_tag, Point, Point);
  static Tag_true trivial(Point, Point);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CommonKernelFunctors::Construct_segment_2<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_segment_3<K>, Args...>
{
  typedef typename K::Point_3    Point;

  static Tag_true trivial(Return_base_tag, Point, Point);
  static Tag_true trivial(Point, Point);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CommonKernelFunctors::Construct_segment_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_circle_2<K>, Args...>
{
  typedef typename K::FT         FT;
  typedef typename K::Point_2    Point;

  static Tag_true trivial(Return_base_tag, Point, FT, Orientation);
  static Tag_true trivial(Return_base_tag, Point, Orientation);
  static Tag_true trivial(Point, FT, Orientation);
  static Tag_true trivial(Point, Orientation);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_circle_2<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Compute_squared_radius_2<K>, Args...>
{
  typedef typename K::Point_2     Point_2;
  typedef typename K::Circle_2    Circle_2;

  static Tag_true trivial(Point_2);
  static Tag_true trivial(Circle_2);
  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Compute_squared_radius_2<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Compute_squared_radius_3<K>, Args...>
{
  typedef typename K::Point_3     Point_3;
  typedef typename K::Circle_3    Circle_3;
  typedef typename K::Sphere_3    Sphere_3;

  static Tag_true trivial(Point_3);
  static Tag_true trivial(Circle_3);
  static Tag_true trivial(Sphere_3);
  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Compute_squared_radius_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_iso_rectangle_2<K>, Args...>
{
  typedef typename K::Point_2    Point;
  typedef typename K::RT         RT;

  static Tag_true trivial(Return_base_tag, Point, Point);
  static Tag_true trivial(Return_base_tag, Point, Point, int);
  static Tag_true trivial(Return_base_tag, Point, Point, Point, Point);
  static Tag_true trivial(Return_base_tag, RT, RT, RT, RT);
  static Tag_true trivial(Point, Point);
  static Tag_true trivial(Point, Point, int);
  static Tag_true trivial(Point, Point, Point, Point);
  static Tag_true trivial(RT, RT, RT, RT);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_iso_rectangle_2<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CartesianKernelFunctors::Construct_iso_cuboid_3<K>, Args...>
{
  typedef typename K::Point_3    Point;
  typedef typename K::RT         RT;

  static Tag_true trivial(Return_base_tag, Point, Point);
  static Tag_true trivial(Return_base_tag, Point, Point, int);
  static Tag_true trivial(Return_base_tag, Point, Point, Point, Point, Point, Point);
  static Tag_true trivial(Return_base_tag, RT, RT, RT, RT, RT, RT);
  static Tag_true trivial(Point, Point);
  static Tag_true trivial(Point, Point, int);
  static Tag_true trivial(Point, Point, Point, Point, Point, Point);
  static Tag_true trivial(RT, RT, RT, RT, RT, RT);
\
  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_iso_cuboid_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_ray_2<K>, Args...>
{
  typedef typename K::Point_2    Point;

  static Tag_true trivial(Return_base_tag, Point, Point);
  static Tag_true trivial(Point, Point);
\
  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CommonKernelFunctors::Construct_ray_2<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_ray_3<K>, Args...>
{
  typedef typename K::Point_3    Point;

  static Tag_true trivial(Return_base_tag, Point, Point);
  static Tag_true trivial(Point, Point);
\
  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CommonKernelFunctors::Construct_ray_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_triangle_2<K>, Args...> : public Tag_true
{};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_triangle_3<K>, Args...> : public Tag_true
{};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_plane_3<K>, Args...>
{
  typedef typename K::RT         RT;
  typedef typename K::Circle_3   Circle;

  static Tag_true trivial(Return_base_tag, RT, RT, RT, RT);
  static Tag_true trivial(Return_base_tag, Circle);
  static Tag_true trivial(RT, RT, RT, RT);
  static Tag_true trivial(Circle);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CommonKernelFunctors::Construct_plane_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_sphere_3<K>, Args...>
{
  typedef typename K::FT         FT;
  typedef typename K::Point_3    Point_3;
  typedef typename K::Circle_3   Circle_3;

  static Tag_true trivial(Return_base_tag, Point_3, FT, Orientation);
  static Tag_true trivial(Return_base_tag, Point_3, Orientation);
  static Tag_true trivial(Return_base_tag, Circle_3);
  static Tag_true trivial(Point_3, FT, Orientation);
  static Tag_true trivial(Point_3, Orientation);
  static Tag_true trivial(Circle_3);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_sphere_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_circle_3<K>, Args...>
{
  typedef typename K::Plane_3      Plane_3;
  typedef typename K::Sphere_3     Sphere_3;

  static Tag_true trivial(Return_base_tag, Plane_3, Sphere_3, int);
  static Tag_true trivial(Plane_3, Sphere_3, int);
  static Tag_true trivial(Sphere_3, Plane_3, int);

  static Tag_false trivial(...);

  enum { value = decltype(trivial(std::declval<CGAL::cpp20::remove_cvref_t<Args>>()...))::value ||
    Is_trivial_construction_base<CGAL::CartesianKernelFunctors::Construct_circle_3<K>, Args...>::value
  };
};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_second_point_3<K>, Args...> : public Tag_true
{};

template <typename K, typename... Args>
struct Is_trivial_construction<CGAL::CommonKernelFunctors::Construct_tetrahedron_3<K>, Args...> : public Tag_true
{};

} // end namespace CartesianFunctors
} // end namespace CGAL

#endif // CGAL_CARTESIAN_IS_TRIVIAL_CONSTRUCTION_H
