// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_K_NEIGHBOR_QUERY_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_K_NEIGHBOR_QUERY_H

#include <CGAL/license/Shape_detection.h>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>

namespace CGAL {
namespace Shape_detection {
namespace Point_set {

  /*!
    \ingroup PkgShapeDetectionRGOnPoints

    \brief K nearest neighbors search in a set of 2D or 3D points.

    This class returns the K nearest neighbors of a query point in a point set.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam Item_
    a descriptor representing a given point. Must be a model of `Hashable`.

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_2` or `Kernel::Point_3`

    \cgalModels{NeighborQuery}
  */
  template<
  typename GeomTraits,
  typename Item_,
  typename PointMap>
  class K_neighbor_query {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Point_map = PointMap;
    using Point_type = typename boost::property_traits<Point_map>::value_type;
    /// \endcond

    /// Item type. Needs to be compatible to Region_type::Item.
    using Item = Item_;
    using Region = std::vector<Item>;

  private:
    using Search_base = typename std::conditional<
      std::is_same<typename Traits::Point_2, Point_type>::value,
      CGAL::Search_traits_2<Traits>,
      CGAL::Search_traits_3<Traits> >::type;

    using Search_traits =
      CGAL::Search_traits_adapter<Item, PointMap, Search_base>;

    using Distance =
      CGAL::Distance_adapter<
      Item,
      PointMap,
      CGAL::Euclidean_distance<Search_base> >;

    using Splitter =
      CGAL::Sliding_midpoint<Search_traits>;

    using Search_tree =
      CGAL::Kd_tree<Search_traits, Splitter, CGAL::Tag_true, CGAL::Tag_true>;

    using Neighbor_search =
      CGAL::Orthogonal_k_neighbor_search<
      Search_traits,
      Distance,
      Splitter,
      Search_tree>;

    using Tree =
      typename Neighbor_search::Tree;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes a Kd-tree with input points.

      \tparam InputRange
      a model of `ConstRange` whose iterator type is `RandomAccessIterator`

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param input_range
      an instance of `InputRange` with 2D or 3D points

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{item_map}
          \cgalParamDescription{an instance of a model of `ReadablePropertyMap` with `InputRange::const_iterator`
                                as key type and `Item` as value type.}
          \cgalParamDefault{A default is provided when `Item` is `InputRange::const_iterator` or its value type.}
        \cgalParamNEnd
        \cgalParamNBegin{k_neighbors}
          \cgalParamDescription{the number of returned neighbors per each query point}
          \cgalParamType{`std::size_t`}
          \cgalParamDefault{12}
        \cgalParamNEnd
        \cgalParamNBegin{point_map}
          \cgalParamDescription{an instance of `PointMap` that maps an item from `input_range`
          to `Kernel::Point_2` or to `Kernel::Point_3`}
          \cgalParamDefault{`PointMap()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
      \pre `K > 0`
    */
    template<typename InputRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
    K_neighbor_query(
      const InputRange& input_range,
      const CGAL_NP_CLASS& np = parameters::default_values()) :
    m_point_map(parameters::choose_parameter<PointMap>(parameters::get_parameter(np, internal_np::point_map))),
    m_distance(m_point_map)
    {
      CGAL_precondition(input_range.size() > 0);

      using NP_helper = internal::Default_property_map_helper<CGAL_NP_CLASS, Item, typename InputRange::const_iterator, internal_np::item_map_t>;
      using Item_map = typename NP_helper::type;
      Item_map item_map = NP_helper::get(np);

      m_tree_ptr.reset( new Tree(make_transform_iterator_from_property_map(make_prevent_deref(input_range.begin()), item_map),
                                 make_transform_iterator_from_property_map(make_prevent_deref(input_range.end()), item_map),
                                 Splitter(),
                                 Search_traits(m_point_map)) );

      const std::size_t K = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::k_neighbors), 12);
      CGAL_precondition(K > 0);
      m_number_of_neighbors = K;
      m_tree_ptr->build();
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `NeighborQuery::operator()()`.

      This operator finds the `K` closest points to the point around the query point `query`
      using a Kd-tree. These references to the points are returned in `neighbors`.

      \param query
      `Item` of the query point

      \param neighbors
      `Items` of points, which are neighbors of the query point

      \pre `query` is a valid const_iterator of InputRange
    */
    void operator()(
      const Item query,
      std::vector<Item>& neighbors) const {

      neighbors.clear();
      Neighbor_search neighbor_search(
        *m_tree_ptr,
        get(m_point_map, query),
        static_cast<unsigned int>(m_number_of_neighbors),
        0, true, m_distance);
      for (auto it = neighbor_search.begin(); it != neighbor_search.end(); ++it)
        neighbors.push_back(it->first);
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    void set_k(const std::size_t k) {
      m_number_of_neighbors = k;
    }
    /// \endcond

  private:
    Point_map m_point_map;
    std::size_t m_number_of_neighbors;
    Distance m_distance;
    std::shared_ptr<Tree> m_tree_ptr; // TODO: move constructor?
  };

  /*!
      \ingroup PkgShapeDetectionRGOnPointSet3
      shortcut to ease the definition of the class when using `CGAL::Point_set_3`.
      To be used together with `make_k_neighbor_query()`.
      \relates K_neighbor_query
   */
  template <class PointSet3>
  using K_neighbor_query_for_point_set =
    K_neighbor_query<typename Kernel_traits<typename PointSet3::Point_3>::Kernel,
                     typename PointSet3::Index,
                     typename PointSet3::Point_map>;

  /*!
      \ingroup PkgShapeDetectionRGOnPointSet3
      returns an instance of the sorting class to be used with `CGAL::Point_set_3`, with point and normal maps added to `np`.
      \relates K_neighbor_query
   */
  template <class PointSet3, typename CGAL_NP_TEMPLATE_PARAMETERS>
  K_neighbor_query_for_point_set<PointSet3>
  make_k_neighbor_query(const PointSet3& ps, CGAL_NP_CLASS np = parameters::default_values())
  {
    return K_neighbor_query_for_point_set<PointSet3>(
      ps, np.point_map(ps.point_map()));
  }

} // namespace Point_set
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_K_NEIGHBOR_QUERY_H
