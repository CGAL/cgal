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
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_SEGMENT_SOUP_ORACLE_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_SEGMENT_SOUP_ORACLE_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_geom_traits.h>
#include <CGAL/Alpha_wrap_3/internal/Oracle_base.h>

#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_segment_primitive_3.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <functional>
#include <memory>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

// Just some typedefs for readability
template <typename GT_>
struct SS_oracle_traits
{
  using Geom_traits = Alpha_wrap_AABB_geom_traits<GT_>; // Wrap the kernel to add Ball_3 + custom Do_intersect_3

  using Segment = typename GT_::Segment_3;
  using Segments = std::vector<Segment>;
  using Segments_ptr = std::shared_ptr<Segments>;
  using SR_iterator = typename Segments::const_iterator;

  using Primitive = AABB_primitive<SR_iterator,
                                   Input_iterator_property_map<SR_iterator>, // DPM
                                   CGAL::internal::Source_of_segment_3_iterator_property_map<Geom_traits, SR_iterator>, // RPM
                                   CGAL::Tag_false, // not external
                                   CGAL::Tag_false>; // no caching

  using AABB_traits = CGAL::AABB_traits_3<Geom_traits, Primitive>;
  using AABB_tree = CGAL::AABB_tree<AABB_traits>;
};

template <typename GT_,
          typename BaseOracle = int>
class Segment_soup_oracle
  : public AABB_tree_oracle<typename SS_oracle_traits<GT_>::Geom_traits,
                            typename SS_oracle_traits<GT_>::AABB_tree,
                            CGAL::Default, // Default_traversal_traits<AABB_traits>
                            BaseOracle>
{
  using SSOT = SS_oracle_traits<GT_>;
  using Base_GT = GT_;

public:
  using Geom_traits = typename SSOT::Geom_traits;

private:
  using Segment = typename SSOT::Segment;
  using Segments = typename SSOT::Segments;
  using Segments_ptr = typename SSOT::Segments_ptr;
  using AABB_tree = typename SSOT::AABB_tree;
  using Oracle_base = AABB_tree_oracle<Geom_traits, AABB_tree, CGAL::Default, BaseOracle>;

private:
  Segments_ptr m_segments_ptr;

public:
  // Constructors
  Segment_soup_oracle(const BaseOracle& base_oracle,
                      const Base_GT& gt = Base_GT())
    : Oracle_base(base_oracle, gt)
  {
    m_segments_ptr = std::make_shared<Segments>();
  }

  Segment_soup_oracle(const Base_GT& gt,
                      const BaseOracle& base_oracle = BaseOracle())
    : Segment_soup_oracle(base_oracle, gt)
  { }

  Segment_soup_oracle()
    : Segment_soup_oracle(BaseOracle(), Base_GT())
  { }

public:
  template <typename SegmentRange,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_segment_soup(const SegmentRange& segments,
                        const CGAL_NP_CLASS& /*np*/ = CGAL::parameters::default_values())
  {
    if(segments.empty())
    {
#ifdef CGAL_AW3_DEBUG
      std::cout << "Warning: Input is empty (SS)" << std::endl;
#endif
      return;
    }

    typename Geom_traits::Is_degenerate_3 is_degenerate = this->geom_traits().is_degenerate_3_object();

    const std::size_t old_size = m_segments_ptr->size();

    for(const Segment& s : segments)
    {
      if(is_degenerate(s))
      {
#ifdef CGAL_AW3_DEBUG
        std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
        continue;
      }

      m_segments_ptr->push_back(s);
    }

#ifdef CGAL_AW3_DEBUG
    std::cout << "Insert into AABB tree (segments)..." << std::endl;
#endif
    this->tree().insert(std::next(std::cbegin(*m_segments_ptr), old_size), std::cend(*m_segments_ptr));

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of a facet that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

    CGAL_postcondition(this->tree().size() == m_segments_ptr->size());
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_SEGMENT_SOUP_ORACLE_H
