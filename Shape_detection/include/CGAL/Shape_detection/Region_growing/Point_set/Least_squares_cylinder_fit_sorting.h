// Copyright (c) 2021 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_CYLINDER_FIT_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_CYLINDER_FIT_SORTING_H

#include <CGAL/license/Shape_detection.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

  /*!
    \ingroup PkgShapeDetectionRGOnPoints

    \brief Sorting of 3D points with respect to the local cylinder fit quality.

    Indices of 3D input points are sorted with respect to the quality of the
    least squares cylinder fit applied to the neighboring points of each point.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam Item_
    a descriptor representing a given point. Must be a model of `Hashable`.

    \tparam NeighborQuery
    a model of `NeighborQuery`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type `Item` and value type is `Kernel::Point_3`

    \tparam NormalMap
    a model of `ReadablePropertyMap` whose key type is `Item` and value type is `Kernel::Vector_3`
  */
  template<
  typename GeomTraits,
  typename Item_,
  typename NeighborQuery,
  typename PointMap,
  typename NormalMap>
  class Least_squares_cylinder_fit_sorting {

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Neighbor_query = NeighborQuery;
    using Point_map = PointMap;
    using Normal_map = NormalMap;
    /// \endcond

    /// Item type.
    using Item = Item_;

    /// Seed range.
    using Seed_range = std::vector<Item>;

    /// @}

  private:
    using FT = typename Traits::FT;
    using Compare_scores = internal::Compare_scores<FT>;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam InputRange
      a model of `ConstRange` whose iterator type is `RandomAccessIterator`

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      an instance of `InputRange` with 3D points

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access point's neighbors

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{item_map}
          \cgalParamDescription{an instance of a model of `ReadablePropertyMap` with `InputRange::const_iterator`
                                as key type and `Item` as value type.}
          \cgalParamDefault{A default is provided when `Item` is `InputRange::const_iterator` or its value type.}
        \cgalParamNEnd
        \cgalParamNBegin{point_map}
          \cgalParamDescription{an instance of `PointMap` that maps an item to `Kernel::Point_3`}
          \cgalParamDefault{`PointMap()`}
        \cgalParamNEnd
        \cgalParamNBegin{normal_map}
          \cgalParamDescription{ an instance of `NormalMap` that maps an item to `Kernel::Vector_3`}
          \cgalParamDefault{`NormalMap()`}
        \cgalParamNEnd
        \cgalParamNBegin{geom_traits}
          \cgalParamDescription{an instance of `GeomTraits`}
          \cgalParamDefault{`GeomTraits()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
    */
    template<typename InputRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
    Least_squares_cylinder_fit_sorting(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      const CGAL_NP_CLASS& np = parameters::default_values()) :
      m_neighbor_query(neighbor_query),
      m_point_map(parameters::choose_parameter<PointMap>(parameters::get_parameter(np, internal_np::point_map))),
      m_normal_map(parameters::choose_parameter<NormalMap>(parameters::get_parameter(np, internal_np::normal_map))),
      m_traits(parameters::choose_parameter<GeomTraits>(parameters::get_parameter(np, internal_np::geom_traits))) {

      CGAL_precondition(input_range.size() > 0);

      using NP_helper = internal::Default_property_map_helper<CGAL_NP_CLASS, Item, typename InputRange::const_iterator, internal_np::item_map_t>;
      using Item_map = typename NP_helper::type;
      Item_map item_map = NP_helper::get(np);

      m_ordered.resize(input_range.size());

      std::size_t index = 0;
      for (auto it = input_range.begin(); it != input_range.end(); it++)
        m_ordered[index++] = get(item_map, it);

      m_scores.resize(input_range.size());
    }

    /// @}

    /// \name Sorting
    /// @{

    /*!
      \brief sorts indices of input points.
    */
    void sort() {
      std::size_t seed_cutoff = compute_scores();
      CGAL_postcondition(m_scores.size() > 0);
      Compare_scores cmp(m_scores);

      std::vector<std::size_t> order(m_ordered.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), cmp);

      order.resize(m_ordered.size() - seed_cutoff);

      std::vector<Item> tmp(order.size());
      for (std::size_t i = 0; i < order.size(); i++)
        tmp[i] = m_ordered[order[i]];

      m_ordered.swap(tmp);
    }
    /// @}

    /// \name Access
    /// @{

    /*!
      \brief returns an instance of `Seed_range` to access the ordered `Items`
      of input points.
    */
    const Seed_range ordered() {
      return m_ordered;
    }
    /// @}

  private:
    Neighbor_query& m_neighbor_query;
    const Point_map m_point_map;
    const Normal_map m_normal_map;
    const Traits m_traits;
    Seed_range m_ordered;
    std::vector<FT> m_scores;

    std::size_t compute_scores() {

      std::vector<Item> neighbors;
      std::size_t idx = 0;
      std::size_t seed_cutoff = 0;
      for (const Item& item : m_ordered) {
        neighbors.clear();
        m_neighbor_query(item, neighbors);
        neighbors.push_back(item);

        m_scores[idx] = -internal::create_cylinder(
          neighbors, m_point_map, m_normal_map, m_traits).second;

        if (m_scores[idx] == -(std::numeric_limits<double>::max)())
          seed_cutoff++;

        idx++;
      }

      return seed_cutoff;
    }
  };

  /*!
      \ingroup PkgShapeDetectionRGOnPointSet3
      shortcut to ease the definition of the class when using `CGAL::Point_set_3`.
      To be used together with `make_least_squares_cylinder_fit_sorting()`.
      \relates Least_squares_cylinder_fit_sorting
   */
  template <class PointSet3, class NeighborQuery>
  using Least_squares_cylinder_fit_sorting_for_point_set =
    Least_squares_cylinder_fit_sorting<typename Kernel_traits<typename PointSet3::Point_3>::Kernel,
                                       typename PointSet3::Index,
                                       NeighborQuery,
                                       typename PointSet3::Point_map,
                                       typename PointSet3::Vector_map>;

  /*!
      \ingroup PkgShapeDetectionRGOnPointSet3
      returns an instance of the sorting class to be used with `CGAL::Point_set_3`, with point and normal maps added to `np`.
      \relates Least_squares_cylinder_fit_sorting
   */
  template <class PointSet3, class NeighborQuery, typename CGAL_NP_TEMPLATE_PARAMETERS>
  Least_squares_cylinder_fit_sorting_for_point_set<PointSet3,NeighborQuery>
  make_least_squares_cylinder_fit_sorting(const PointSet3& ps,
                                          NeighborQuery& neighbor_query,
                                          const CGAL_NP_CLASS np = parameters::default_values())
  {
    return Least_squares_cylinder_fit_sorting_for_point_set<PointSet3,NeighborQuery>
      (ps, neighbor_query, np.point_map(ps.point_map()).normal_map(ps.normal_map()));
  }

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_LEAST_SQUARES_CYLINDER_FIT_SORTING_H
