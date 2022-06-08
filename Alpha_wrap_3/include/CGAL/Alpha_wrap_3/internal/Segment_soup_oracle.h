// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : TBA
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_SEGMENT_SOUP_ORACLE_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_SEGMENT_SOUP_ORACLE_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_traits.h>
#include <CGAL/Alpha_wrap_3/internal/Oracle_base.h>

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_segment_primitive.h>
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

// Just some typedefs for readability
template <typename SegmentRange, typename GT_>
struct SS_oracle_traits
{
  using Segment = typename boost::range_value<SegmentRange>::type;
  using Default_GT = typename Kernel_traits<Segment>::Kernel;
  using Base_GT = typename Default::Get<GT_, Default_GT>::type; // = Kernel, usually
  using Geom_traits = Alpha_wrap_AABB_traits<Base_GT>; // Wrap the kernel to add Ball_3 + custom Do_intersect_3

  using SR_iterator = typename SegmentRange::const_iterator;
  using Primitive = AABB_primitive<SR_iterator,
                                   Input_iterator_property_map<SR_iterator>, // DPM
                                   CGAL::internal::Source_of_segment_3_iterator_property_map<Geom_traits, SR_iterator>, // RPM
                                   CGAL::Tag_false, // not external
                                   CGAL::Tag_false>; // no caching

  using AABB_traits = CGAL::AABB_traits<Geom_traits, Primitive>;
  using AABB_tree = CGAL::AABB_tree<AABB_traits>;
};

template <typename SegmentRange,
          typename GT_ = CGAL::Default,
          typename BaseOracle = int>
class Segment_soup_oracle
  : public AABB_tree_oracle<typename SS_oracle_traits<SegmentRange, GT_>::Geom_traits,
                            typename SS_oracle_traits<SegmentRange, GT_>::AABB_tree,
                            CGAL::Default, // Default_traversal_traits<AABB_traits>
                            BaseOracle>
{
  using SSOT = SS_oracle_traits<SegmentRange, GT_>;
  using Base_GT = typename SSOT::Base_GT;

public:
  using Geom_traits = typename SSOT::Geom_traits;

private:
  using AABB_tree = typename SSOT::AABB_tree;
  using Oracle_base = AABB_tree_oracle<Geom_traits, AABB_tree, CGAL::Default, BaseOracle>;

public:
  // Constructors
  Segment_soup_oracle()
    : Oracle_base(BaseOracle(), Base_GT())
  { }

  Segment_soup_oracle(const BaseOracle& base_oracle,
                      const Base_GT& gt = Base_GT())
    : Oracle_base(base_oracle, gt)
  { }

  Segment_soup_oracle(const Base_GT& gt,
                      const BaseOracle& base_oracle = BaseOracle())
    : Oracle_base(base_oracle, gt)
  { }

public:
  template <typename NamedParameters = parameters::Default_named_parameters>
  void add_segment_soup(const SegmentRange& segments,
                        const NamedParameters& /*np*/ = CGAL::parameters::default_values())
  {
    if(segments.empty())
    {
#ifdef CGAL_AW3_DEBUG
      std::cout << "Warning: Input is empty " << std::endl;
#endif
      return;
    }

#ifdef CGAL_AW3_DEBUG
    std::cout << "Insert into AABB tree (segments)..." << std::endl;
#endif
    this->tree().insert(segments.begin(), segments.end());

    CGAL_postcondition(this->tree().size() == segments.size());
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_SEGMENT_SOUP_ORACLE_H
