// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_POINT_SET_ORACLE_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_POINT_SET_ORACLE_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_traits.h>
#include <CGAL/Alpha_wrap_3/internal/Oracle_base.h>

#include <CGAL/AABB_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>

#include <boost/range/value_type.hpp>

#include <algorithm>
#include <iostream>
#include <functional>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

// No longer used, but might find some purpose again in the future
template <class InputIterator, class PM>
struct Point_from_iterator_map
{
  using key_type = InputIterator;
  using value_type = typename boost::property_traits<PM>::value_type;
  using reference = typename boost::property_traits<PM>::reference;
  using category = boost::readable_property_map_tag;

  Point_from_iterator_map(const PM& pm = PM()) : pm(pm) { }

  inline friend reference get(const Point_from_iterator_map& map, const key_type it)
  {
    return get(map.pm, *it);
  }

  PM pm;
};

// Just some typedefs for readability
template <typename PointRange, typename GT_>
struct PS_oracle_traits
{
  using Point = typename boost::range_value<PointRange>::type;
  using Default_GT = typename Kernel_traits<Point>::Kernel;
  using Base_GT = typename Default::Get<GT_, Default_GT>::type; // = Kernel, usually
  using Geom_traits = Alpha_wrap_AABB_traits<Base_GT>; // Wrap the kernel to add Ball_3 + custom Do_intersect_3

  using PR_iterator = typename PointRange::const_iterator;
  using Primitive = AABB_primitive<PR_iterator,
                                   Input_iterator_property_map<PR_iterator> /*DPM*/,
                                   Input_iterator_property_map<PR_iterator> /*RPM*/,
                                   CGAL::Tag_false, // not external
                                   CGAL::Tag_false>; // no caching

  using AABB_traits = CGAL::AABB_traits<Geom_traits, Primitive>;
  using AABB_tree = CGAL::AABB_tree<AABB_traits>;
};

template <typename PointRange,
          typename GT_ = CGAL::Default,
          typename BaseOracle = int>
class Point_set_oracle
  : public AABB_tree_oracle<typename PS_oracle_traits<PointRange, GT_>::Geom_traits,
                            typename PS_oracle_traits<PointRange, GT_>::AABB_tree,
                            CGAL::Default, // Default_traversal_traits<AABB_traits>
                            BaseOracle>
{
  using PSOT = PS_oracle_traits<PointRange, GT_>;
  using Base_GT = typename PSOT::Base_GT;

public:
  using Geom_traits = typename PSOT::Geom_traits;

private:
  using AABB_tree = typename PSOT::AABB_tree;
  using Oracle_base = AABB_tree_oracle<Geom_traits, AABB_tree, CGAL::Default, BaseOracle>;

public:
  // Constructors
  Point_set_oracle()
    : Oracle_base(BaseOracle(), Base_GT())
  { }

  Point_set_oracle(const BaseOracle& base_oracle,
                   const Base_GT& gt = Base_GT())
    : Oracle_base(base_oracle, gt)
  { }

  Point_set_oracle(const Base_GT& gt,
                   const BaseOracle& base_oracle = BaseOracle())
    : Oracle_base(base_oracle, gt)
  { }

public:
  // adds a range of points to the oracle
  template <typename NamedParameters = parameters::Default_named_parameters>
  void add_point_set(const PointRange& points,
                     const NamedParameters& /*np*/ = CGAL::parameters::default_values())
  {
    if(points.empty())
    {
#ifdef CGAL_AW3_DEBUG
      std::cout << "Warning: Input is empty " << std::endl;
#endif
      return;
    }

#ifdef CGAL_AW3_DEBUG
    std::cout << "Insert into AABB tree (points)..." << std::endl;
#endif

    this->tree().insert(points.begin(), points.end());
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_POINT_SET_ORACLE_H
