// Copyright (c) 2019-2022 Google LLC (USA).
// Copyright (c) 2025 GeometryFactory (France).
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
#ifndef CGAL_ALPHA_WRAP_2_INTERNAL_SEGMENT_SOUP_ORACLE_H
#define CGAL_ALPHA_WRAP_2_INTERNAL_SEGMENT_SOUP_ORACLE_H

#include <CGAL/license/Alpha_wrap_2.h>

#include <CGAL/Alpha_wrap_2/internal/Alpha_wrap_AABB_geom_traits.h>
#include <CGAL/Alpha_wrap_2/internal/Oracle_base.h>

#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/property_map.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <functional>
#include <memory>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

// Just some typedefs for readability
template <typename GT_>
struct SS_oracle_traits
{
  using Geom_traits = Alpha_wrap_AABB_geom_traits<GT_>; // Wrap the kernel to add Disk_2 + custom Do_intersect_2

  using Point = typename GT_::Point_2;
  using Segment = typename GT_::Segment_2;
  using Segments = std::vector<Segment>;
  using Segments_ptr = std::shared_ptr<Segments>;
  using SR_iterator = typename Segments::const_iterator;

  using Primitive = AABB_primitive<SR_iterator,
                                   Input_iterator_property_map<SR_iterator>, // DPM
                                   CGAL::internal::Source_of_segment_2_iterator_property_map<Geom_traits, SR_iterator>, // RPM
                                   CGAL::Tag_false, // not external
                                   CGAL::Tag_false>; // no caching

  using AABB_traits = CGAL::AABB_traits_2<Geom_traits, Primitive>;
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
  using Point = typename SSOT::Point;
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
  void clear()
  {
    m_segments_ptr->clear();
    Oracle_base::clear();
  }

  template <typename SegmentRange,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_segments(const SegmentRange& segments,
                    const CGAL_NP_CLASS& /*np*/ = CGAL::parameters::default_values())
  {
#ifdef CGAL_AW2_DEBUG
    std::cout << "Insert into AABB tree (" << segments.size() << " segments)..." << std::endl;
#endif

    if(segments.empty())
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Warning: Input is empty (SS)" << std::endl;
#endif
      return;
    }

    const std::size_t old_size = m_segments_ptr->size();

    typename Geom_traits::Is_degenerate_2 is_degenerate = this->geom_traits().is_degenerate_2_object();

    for(const Segment& s : segments)
    {
      if(is_degenerate(s))
      {
#ifdef CGAL_AW2_DEBUG
        std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
        continue;
      }

      m_segments_ptr->push_back(s);
    }

    this->tree().rebuild(std::cbegin(*m_segments_ptr), std::cend(*m_segments_ptr));

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of an edge that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW2_DEBUG
    std::cout << "SS Tree: " << this->tree().size() << " primitives" << std::endl;
#endif
  }

  template <typename TriangleRange,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_triangles(const TriangleRange& triangles,
                     const CGAL_NP_CLASS& /*np*/ = CGAL::parameters::default_values())
  {
#ifdef CGAL_AW2_DEBUG
    std::cout << "Insert into AABB Tree (" << triangles.size()  << " triangles)..." << std::endl;
#endif

    if(triangles.empty())
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Warning: Input is empty (TS)" << std::endl;
#endif
      return;
    }

    const std::size_t old_size = m_segments_ptr->size();

    typename Geom_traits::Construct_segment_2 segment = this->geom_traits().construct_segment_2_object();
    typename Geom_traits::Is_degenerate_2 is_degenerate = this->geom_traits().is_degenerate_2_object();

    for(const auto& tr : triangles)
    {
      for(int i=0; i<3; ++i)
      {
        Segment s = segment(tr[i], tr[(i+1)%3]);
        if(is_degenerate(s))
        {
#ifdef CGAL_AW2_DEBUG
          std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
          continue;
        }

        m_segments_ptr->push_back(s);
      }
    }

    this->tree().rebuild(std::cbegin(*m_segments_ptr), std::cend(*m_segments_ptr));

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of an edge that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW2_DEBUG
    std::cout << "SS Tree: " << this->tree().size() << " primitives" << std::endl;
#endif
  }

  template <typename PointRange, typename FaceRange,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_polygon_soup(const PointRange& points,
                        const FaceRange& faces,
                        const CGAL_NP_CLASS& np = CGAL::parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    using PPM = typename GetPointMap<PointRange, CGAL_NP_CLASS>::const_type;
    using Point_ref = typename boost::property_traits<PPM>::reference;

    using Face = typename boost::range_value<FaceRange>::type;

#ifdef CGAL_AW2_DEBUG
    std::cout << "Insert into AABB tree (" << faces.size() << " polygons)..." << std::endl;
#endif

    if(points.empty() || faces.empty())
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Warning: Input is empty (PS)" << std::endl;
#endif
      return;
    }

    const std::size_t old_size = m_segments_ptr->size();

    PPM pm = choose_parameter<PPM>(get_parameter(np, internal_np::point_map));
    static_assert(std::is_same<typename boost::property_traits<PPM>::value_type, Point>::value);

    typename Geom_traits::Construct_segment_2 segment = this->geom_traits().construct_segment_2_object();
    typename Geom_traits::Is_degenerate_2 is_degenerate = this->geom_traits().is_degenerate_2_object();

    for(const Face& f : faces)
    {
      if(f.size() < 2)
        continue;

      for(std::size_t i=0,n=f.size(); i<n; ++i)
      {
        Segment s = segment(get(pm, points[f[i]]), get(pm, points[f[(i+1)%n]]));
        if(is_degenerate(s))
        {
#ifdef CGAL_AW2_DEBUG
          std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
          continue;
        }

        m_segments_ptr->push_back(s);
      }
    }

    this->tree().rebuild(std::cbegin(*m_segments_ptr), std::cend(*m_segments_ptr));

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of an edge that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW2_DEBUG
    std::cout << "SS Tree: " << this->tree().size() << " primitives" << std::endl;
#endif
  }

  template <typename MultiLineString,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_multilinestring(const MultiLineString& mls,
                           const CGAL_NP_CLASS& /*np*/ = CGAL::parameters::default_values())
  {
    using LineString = typename boost::range_value<MultiLineString>::type;

    const std::size_t old_size = m_segments_ptr->size();

    typename Geom_traits::Construct_segment_2 segment = this->geom_traits().construct_segment_2_object();
    typename Geom_traits::Is_degenerate_2 is_degenerate = this->geom_traits().is_degenerate_2_object();

#ifdef CGAL_AW2_DEBUG
    std::cout << "Insert into AABB tree (multi-linestring)..." << std::endl;
#endif

    if(mls.empty())
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Warning: Input is empty (multi-linestring)" << std::endl;
#endif
      return;
    }

    for(const LineString& ls : mls)
    {
      for(std::size_t i=0; i<ls.size()-1; ++i)
      {
        const Segment s = segment(ls[i], ls[i+1]);
        if(is_degenerate(s))
        {
#ifdef CGAL_AW2_DEBUG
          std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
          continue;
        }
        m_segments_ptr->push_back(s);
      }
    }

    this->tree().rebuild(std::cbegin(*m_segments_ptr), std::cend(*m_segments_ptr));

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of an edge that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW2_DEBUG
    std::cout << "SS Tree: " << this->tree().size() << " primitives" << std::endl;
#endif
  }

  template <typename PolygonalChains,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_polygonal_chains(const PolygonalChains& pcs,
                            const CGAL_NP_CLASS& np = CGAL::parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    using Polygonal_chain = typename boost::range_value<PolygonalChains>::type;
    using value_type = boost::range_value<Polygonal_chain>::type;

    using IPM = Identity_property_map<value_type>;
    using Point_map = typename internal_np::Lookup_named_param_def<
                        internal_np::point_t,
                        CGAL_NP_CLASS,
                        IPM>::type;

    Point_map point_map = parameters::choose_parameter<Point_map>(
                            parameters::get_parameter(np, internal_np::point_map));

    const bool close_pc = choose_parameter(get_parameter(np, internal_np::close_chains), false);

    typename Geom_traits::Construct_segment_2 segment = this->geom_traits().construct_segment_2_object();
    typename Geom_traits::Is_degenerate_2 is_degenerate = this->geom_traits().is_degenerate_2_object();

#ifdef CGAL_AW2_DEBUG
    std::cout << "Insert into AABB tree (polygonal chain)..." << std::endl;
#endif

    if(pcs.empty())
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Warning: Input is empty (polygonal chain)" << std::endl;
#endif
      return;
    }

    const std::size_t old_size = m_segments_ptr->size();

    for(const Polygonal_chain& pc : pcs)
    {
      for(std::size_t i=0; i<pc.size()-1; ++i)
      {
        const Segment s = segment(get(point_map, pc[i]),
                                  get(point_map, pc[i+1]));
        if(is_degenerate(s))
        {
#ifdef CGAL_AW2_DEBUG
          std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
          continue;
        }
        m_segments_ptr->push_back(s);
      }

      if(close_pc && pc.size() > 1)
      {
        const Segment s = segment(get(point_map, pc[pc.size() - 1]),
                                  get(point_map, pc[0]));
        if(is_degenerate(s))
        {
#ifdef CGAL_AW2_DEBUG
          std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
        }
        else
        {
          m_segments_ptr->push_back(s);
        }
      }
    }

    this->tree().rebuild(std::cbegin(*m_segments_ptr), std::cend(*m_segments_ptr));

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of an edge that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW2_DEBUG
    std::cout << "SS Tree: " << this->tree().size() << " primitives" << std::endl;
#endif
  }

  template <typename MultipolygonWithHoles,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_multipolygon(const MultipolygonWithHoles& mp,
                        const CGAL_NP_CLASS& /*np*/ = CGAL::parameters::default_values())
  {
    using Polygon_with_holes_2 = typename MultipolygonWithHoles::Polygon_with_holes_2;

    typename Geom_traits::Is_degenerate_2 is_degenerate = this->geom_traits().is_degenerate_2_object();

#ifdef CGAL_AW2_DEBUG
    std::cout << "Insert into AABB tree (multi-polygon)..." << std::endl;
#endif

    if(mp.polygons_with_holes().empty())
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Warning: Input is empty (multi-polygon)" << std::endl;
#endif
      return;
    }

    const std::size_t old_size = m_segments_ptr->size();

    for(const Polygon_with_holes_2& polygon : mp.polygons_with_holes()) {
      for(const Segment& s : polygon.outer_boundary().edges()) {
        if(is_degenerate(s))
        {
#ifdef CGAL_AW2_DEBUG
          std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
          continue;
        }
        m_segments_ptr->push_back(s);
      }
      for(const auto& hole : polygon.holes()) {
        for(const Segment& s : hole.edges()) {
          if(is_degenerate(s))
          {
#ifdef CGAL_AW2_DEBUG
            std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
            continue;
          }
          m_segments_ptr->push_back(s);
        }
      }
    }

    this->tree().rebuild(std::cbegin(*m_segments_ptr), std::cend(*m_segments_ptr));

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of an edge that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW2_DEBUG
    std::cout << "SS Tree: " << this->tree().size() << " primitives" << std::endl;
#endif
  }

  template <typename PolygonWithHoles,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_polygon_with_holes(const PolygonWithHoles& pwh,
                              const CGAL_NP_CLASS& /*np*/ = CGAL::parameters::default_values())
  {
    typename Geom_traits::Is_degenerate_2 is_degenerate = this->geom_traits().is_degenerate_2_object();

#ifdef CGAL_AW2_DEBUG
    std::cout << "Insert into AABB tree (polygon with holes)..." << std::endl;
#endif

    if(pwh.outer_boundary().is_empty() && pwh.holes().empty())
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Warning: Input is empty (polygon with holes)" << std::endl;
#endif
      return;
    }

    const std::size_t old_size = m_segments_ptr->size();

    for(const Segment& s : pwh.outer_boundary().edges()) {
      if(is_degenerate(s))
      {
#ifdef CGAL_AW2_DEBUG
        std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
        continue;
      }
      m_segments_ptr->push_back(s);
    }

    for(const auto& hole : pwh.holes()) {
      for(const Segment& s : hole.edges()) {
        if(is_degenerate(s))
        {
#ifdef CGAL_AW2_DEBUG
          std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
          continue;
        }
        m_segments_ptr->push_back(s);
      }
    }

    this->tree().rebuild(std::cbegin(*m_segments_ptr), std::cend(*m_segments_ptr));

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of an edge that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW2_DEBUG
    std::cout << "SS Tree: " << this->tree().size() << " primitives" << std::endl;
#endif
  }

  template <typename Traits_, typename Container_,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_polygon(const CGAL::Polygon_2<Traits_, Container_>& p,
                   const CGAL_NP_CLASS& /*np*/ = CGAL::parameters::default_values())
  {
    typename Geom_traits::Is_degenerate_2 is_degenerate = this->geom_traits().is_degenerate_2_object();

#ifdef CGAL_AW2_DEBUG
    std::cout << "Insert into AABB tree (polyon)..." << std::endl;
#endif

    if(p.is_empty())
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Warning: Input is empty (polyon)" << std::endl;
#endif
      return;
    }

    const std::size_t old_size = m_segments_ptr->size();

    for(const Segment& s : p.edges()) {
      if(is_degenerate(s))
      {
#ifdef CGAL_AW2_DEBUG
        std::cerr << "Warning: ignoring degenerate segment " << s << std::endl;
#endif
        continue;
      }
      m_segments_ptr->push_back(s);
    }

    this->tree().rebuild(std::cbegin(*m_segments_ptr), std::cend(*m_segments_ptr));

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of an edge that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW2_DEBUG
    std::cout << "SS Tree: " << this->tree().size() << " primitives" << std::endl;
#endif
  }
};

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_INTERNAL_SEGMENT_SOUP_ORACLE_H
