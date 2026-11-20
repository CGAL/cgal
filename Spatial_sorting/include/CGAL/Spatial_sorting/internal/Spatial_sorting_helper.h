
// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     :  Léo Valque

#ifndef CGAL_SPATIAL_SORTING_HELPER_3_H
#define CGAL_SPATIAL_SORTING_HELPER_3_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {
namespace Spatial_sorting {
namespace internal {

template<class Kernel, class = void>
struct Has_construct_bbox_2 : std::false_type {};

template<class Kernel>
struct Has_construct_bbox_2<Kernel,
  std::void_t<decltype(std::declval<const Kernel&>().construct_bbox_2_object())>>
  : std::true_type {};

template<class Kernel, class = void>
struct Has_construct_bbox_3 : std::false_type {};

template<class Kernel>
struct Has_construct_bbox_3<Kernel,
  std::void_t<decltype(std::declval<const Kernel&>().construct_bbox_3_object())>>
  : std::true_type {};

// Custom kernel containing only compute_x_2_object()
template<class Kernel, bool = Has_construct_bbox_2<Kernel>::value>
struct Kernel_wrapper_impl_2
{
  using type = Kernel_wrapper_impl_2<Kernel, false>;

  Kernel _k;
  Kernel_wrapper_impl_2(Kernel k) : _k(k) {}

  struct Construct_bbox_2
  {
    Kernel _k;
    explicit Construct_bbox_2(Kernel k) : _k(k){}

    template <typename T>
    Bbox_2 operator()(const T& obj) const
    {
    auto x = to_double(_k.compute_x_2_object()(obj));
    auto y = to_double(_k.compute_y_2_object()(obj));
    return Bbox_2(x, y, x, y);
    }
  };
  Construct_bbox_2 construct_bbox_2_object() const { return Construct_bbox_2(_k); }
};

// Custom kernel containing only compute_x_3_object()
template<class Kernel, bool = Has_construct_bbox_3<Kernel>::value>
struct Kernel_wrapper_impl_3
{
  using type = Kernel_wrapper_impl_3<Kernel, false>;

  Kernel _k;
  Kernel_wrapper_impl_3(Kernel k) : _k(k) {}

  struct Construct_bbox_3
  {
    Kernel _k;
    explicit Construct_bbox_3(Kernel k) : _k(k){}

    template <typename T>
    Bbox_3 operator()(const T& obj) const
    {
    auto x = to_double(_k.compute_x_3_object()(obj));
    auto y = to_double(_k.compute_y_3_object()(obj));
    auto z = to_double(_k.compute_z_3_object()(obj));
    return Bbox_3(x, y, z, x, y, z);
    }
  };
  Construct_bbox_3 construct_bbox_3_object() const { return Construct_bbox_3(_k); }
};

// If Kernel is complete, use it directly
template<class Kernel>
struct Kernel_wrapper_impl_2<Kernel, true>
{
  using type = Kernel;
};

template<class Kernel>
struct Kernel_wrapper_impl_3<Kernel, true>
{
  using type = Kernel;
};

template<class Kernel>
using Kernel_wrapper_2 = typename Kernel_wrapper_impl_2<Kernel>::type;

template<class Kernel>
using Kernel_wrapper_3 = typename Kernel_wrapper_impl_3<Kernel>::type;

} // namespace internal
} // namespace Spatial_sorting
} // namespace CGAL

#endif // CGAL_SPATIAL_SORTING_HELPER_3_H