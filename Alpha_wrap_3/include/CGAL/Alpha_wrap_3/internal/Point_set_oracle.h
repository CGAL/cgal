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

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_geom_traits.h>
#include <CGAL/Alpha_wrap_3/internal/Oracle_base.h>

#include <CGAL/AABB_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <functional>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

// Just some typedefs for readability
template <typename GT_>
struct PS_oracle_traits
{
  using Geom_traits = Alpha_wrap_AABB_geom_traits<GT_>; // Wrap the kernel to add Ball_3 + custom Do_intersect_3

  using Points = std::vector<typename GT_::Point_3>;
  using PR_iterator = typename Points::const_iterator;

  using Primitive = AABB_primitive<PR_iterator,
                                   Input_iterator_property_map<PR_iterator> /*DPM*/,
                                   Input_iterator_property_map<PR_iterator> /*RPM*/,
                                   CGAL::Tag_false, // not external
                                   CGAL::Tag_false>; // no caching

  using AABB_traits = CGAL::AABB_traits<Geom_traits, Primitive>;
  using AABB_tree = CGAL::AABB_tree<AABB_traits>;
};

template <typename GT_,
          typename BaseOracle = int>
class Point_set_oracle
  : public AABB_tree_oracle<typename PS_oracle_traits<GT_>::Geom_traits,
                            typename PS_oracle_traits<GT_>::AABB_tree,
                            CGAL::Default, // Default_traversal_traits<AABB_traits>
                            BaseOracle>
{
  using PSOT = PS_oracle_traits<GT_>;
  using Base_GT = GT_;

public:
  using Geom_traits = typename PSOT::Geom_traits;

private:
  using Points = typename PSOT::Points;
  using AABB_tree = typename PSOT::AABB_tree;
  using Oracle_base = AABB_tree_oracle<Geom_traits, AABB_tree, CGAL::Default, BaseOracle>;

private:
  Points m_points;

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
  template <typename PointRange,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_point_set(const PointRange& points,
                     const CGAL_NP_CLASS& /*np*/ = CGAL::parameters::default_values())
  {
    if(points.empty())
    {
#ifdef CGAL_AW3_DEBUG
      std::cout << "Warning: Input is empty " << std::endl;
#endif
      return;
    }

    const std::size_t old_size = m_points.size();
    m_points.insert(std::cend(m_points), std::cbegin(points), std::cend(points));

#ifdef CGAL_AW3_DEBUG
    std::cout << "Insert into AABB tree (points)..." << std::endl;
#endif

    this->tree().insert(std::next(std::cbegin(m_points), old_size), std::cend(m_points));

    CGAL_postcondition(this->tree().size() == m_points.size());
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_POINT_SET_ORACLE_H
