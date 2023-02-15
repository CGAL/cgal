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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <queue>
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Shape_detection {

  /*!
    \ingroup PkgShapeDetectionRG

    \brief Main class/entry point for running the region growing algorithm.

    This version of the region growing algorithm enables to detect regions in a set
    of user-defined items
    - given a way to access neighbors of each item via the `NeighborQuery` parameter class and
    - control if items form a valid region type via the `RegionType` parameter class,
    - the `SeedMap` property map enables to define the seeding order of items and skip unnecessary items.

    \tparam InputRange
    must be a model of `ConstRange`.

    \tparam NeighborQuery
    must be a model of `NeighborQuery`.

    \tparam RegionType
    must be a model of `RegionType`.

    \tparam SeedMap
    must be an `LvaluePropertyMap` whose key and value types are `std::size_t`.
    %Default is `CGAL::Identity_property_map`.
  */
  template<
  typename InputRange,
  typename NeighborQuery,
  typename RegionType,
  typename SeedMap = CGAL::Identity_property_map<std::size_t> >
  class Region_growing {

  public:

    /// \cond SKIP_IN_MANUAL
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Region_type = RegionType;
    using Seed_map = SeedMap;

    using Visited_items = std::vector<bool>;
    using Running_queue = std::queue<std::size_t>;
    using Indices       = std::vector<std::size_t>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes the region growing algorithm.

      \param input_range
      a range of input items for region growing

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access item's neighbors

      \param region_type
      an instance of `RegionType` that is used internally to
      control if items form a valid region type

      \param seed_map
      an instance of `SeedMap` property map that is used internally to
      set the order of items in the region growing processing queue. If it maps
      to `std::size_t(-1)`, the corresponding item is skipped.

      \pre `input_range.size() > 0`
    */
    Region_growing(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      RegionType& region_type,
      const SeedMap seed_map = SeedMap()) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_region_type(region_type),
    m_seed_map(seed_map) {

      CGAL_precondition(input_range.size() > 0);
      clear();
    }

    /// @}

    /// \name Detection
    /// @{

    /*!
      \brief runs the region growing algorithm and fills an output iterator
      with the found regions.

      \tparam OutputIterator
      must be an output iterator whose value type is `std::vector<std::size_t>`.

      \param regions
      an output iterator that stores regions, where each region is returned
      as a vector of indices of the items, which belong to this region

      \return past-the-end position in the output sequence
    */
    template<typename OutputIterator>
    OutputIterator detect(OutputIterator regions) {

      clear();
      Indices region;

      // Grow regions.
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const std::size_t seed_index = get(m_seed_map, i);

        // Skip items that user does not want to use.
        if (seed_index == std::size_t(-1))
          continue;

        CGAL_precondition(
          seed_index < m_input_range.size());

        // Try to grow a new region from the index of the seed item.
        if (!m_visited[seed_index]) {
          propagate(seed_index, region);

          // Check global conditions.
          if (!m_region_type.is_valid_region(region))
            revert(region);
          else
            *(regions++) = region;
        }
      }
      return regions;
    }

    /// @}

    /// \name Unassigned Items
    /// @{

    /*!
      \brief fills an output iterator with indices of all unassigned items.

      \tparam OutputIterator
      must be an output iterator whose value type is `std::size_t`.

      \param output
      an output iterator that stores indices of all items, which are not assigned
      to any region

      \return past-the-end position in the output sequence
    */
    template<typename OutputIterator>
    OutputIterator unassigned_items(OutputIterator output) const {

      // Return indices of all unassigned items.
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const std::size_t seed_index = get(m_seed_map, i);

        // Skip items that user does not want to use.
        if (seed_index == std::size_t(-1))
          continue;

        CGAL_precondition(
          seed_index < m_input_range.size());

        if (!m_visited[seed_index])
          *(output++) = seed_index;
      }
      return output;
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    void clear() {

      m_visited.clear();
      m_visited.resize(m_input_range.size(), false);
    }

    void release_memory() {

      m_visited.clear();
      m_visited.shrink_to_fit();
    }
    /// \endcond

  private:

    void propagate(const std::size_t seed_index, Indices& region) {
      region.clear();

      // Use two queues, while running on this queue, push to the other queue;
      // When the queue is done, update the shape of the current region and swap to the other queue;
      // depth_index is the index of the queue we are using.
      Running_queue running_queue[2];
      bool depth_index = 0;

      // Once the index of an item is pushed to the queue, it is pushed to the region too.
      m_visited[seed_index] = true;
      running_queue[depth_index].push(seed_index);
      region.push_back(seed_index);

      // Update internal properties of the region.
      m_region_type.update(region);

      Indices neighbors;
      while (
        !running_queue[depth_index].empty() ||
        !running_queue[!depth_index].empty()) {

        // Call the next item index of the queue and remove it from the queue.
        const std::size_t item_index = running_queue[depth_index].front();
        running_queue[depth_index].pop();

        // Get neighbors of the current item.
        neighbors.clear();
        m_neighbor_query(item_index, neighbors);

        // Visit all found neighbors.
        for (const std::size_t neighbor_index : neighbors) {

          // Skip items that user does not want to use.
          if (neighbor_index == std::size_t(-1))
            continue;

          CGAL_precondition(
            neighbor_index < m_input_range.size());

          if (!m_visited[neighbor_index] &&
            m_region_type.is_part_of_region(item_index, neighbor_index, region)) {

            // Add this neighbor to the other queue so that we can visit it later.
            m_visited[neighbor_index] = true;
            running_queue[!depth_index].push(neighbor_index);
            region.push_back(neighbor_index);
          }
        }

        // Update internal properties of the region.
        if (running_queue[depth_index].empty()) {

          m_region_type.update(region);
          depth_index = !depth_index;
        }
      }
    }

    void revert(const Indices& region) {
      for (const std::size_t item_index : region)
        m_visited[item_index] = false;
    }

    // Fields.
    const Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    Region_type& m_region_type;
    const Seed_map m_seed_map;

    Visited_items m_visited;
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_H
