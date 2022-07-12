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

// Boost includes.
#include <boost/unordered_map.hpp>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/is_iterator.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>

namespace CGAL {
namespace Shape_detection {

  /*!
    \ingroup PkgShapeDetectionRG

    \brief Main class/entry point for running the region growing algorithm.

    This version of the region growing algorithm enables to detect regions in a set
    of user-defined items
    - given a way to access neighbors of each item via the `NeighborQuery` parameter class and
    - control if items form a valid region type via the `RegionType` parameter class,
    - optional `SeedRange` defining the seeding order of items and skipping unnecessary items.

    \tparam NeighborQuery
    a model of `NeighborQuery`

    \tparam RegionType
    a model of `RegionType`
  */
  template<
    typename NeighborQuery,
    typename RegionType >
  class Region_growing {

  public:
    /// \name Types
    /// \cond SKIP_IN_MANUAL
    using Neighbor_query = NeighborQuery;
    using Region_type = RegionType;
    /// \endcond

    /// Item type.
    using Item = typename RegionType::Item;
    using Region = std::vector<Item>;

    /// Primitive and region type
    using Primitive_and_region = std::pair<typename Region_type::Primitive, Region>;

    /// Item to region property map.
    using Region_map = typename Region_type::Region_index_map;

  private:
    using Running_queue = std::queue<Item>;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes the region growing algorithm.

      \tparam InputRange
        a model of `ConstRange`

      \param input_range
      a range of input items for region growing

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access item's neighbors

      \param region_type
      an instance of `RegionType` that is used internally to
      control if items form a valid region type

      \pre `input_range.size() > 0`
    */
    template<typename InputRange>
    Region_growing(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      RegionType& region_type) :
      m_neighbor_query(neighbor_query),
      m_region_type(region_type),
      m_region_map(region_type.region_index_map()),
      m_visited(m_visited_map) {

      CGAL_precondition(input_range.size() > 0);
      m_seed_range.resize(input_range.size());

      std::size_t idx = 0;
      for (auto it = input_range.begin(); it != input_range.end(); it++)
        m_seed_range[idx++] = internal::conditional_deref<typename InputRange::const_iterator, Item>()(it);

      clear(input_range);
    }

    /*!
      \brief initializes the region growing algorithm.

      \tparam InputRange
        a model of `ConstRange

      \param input_range
      a range of input items for region growing

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access item's neighbors

      \param region_type
      an instance of `RegionType` that is used internally to
      control if items form a valid region type

      \param seed_range
      a vector of `Item` that is used as seeds for the region growing.
      Defaults to the full input_range.

      \pre `input_range.size() > 0`
    */
    template<typename InputRange, typename SeedRange>
    Region_growing(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      RegionType& region_type,
      SeedRange& seed_range) :
      m_neighbor_query(neighbor_query),
      m_region_type(region_type),
      m_region_map(region_type.region_index_map()),
      m_visited(m_visited_map) {

      CGAL_precondition(input_range.size() > 0);
      CGAL_precondition(seed_range.size() > 0);

      m_seed_range.resize(seed_range.size());

      std::size_t idx = 0;
      for (auto it = seed_range.begin(); it != seed_range.end(); it++)
        m_seed_range[idx++] = internal::conditional_deref<typename SeedRange::const_iterator, Item>()(it);

      clear(input_range);
    }

    /// @}

    /// \name Detection
    /// @{

    /*!
      \brief runs the region growing algorithm and fills an output iterator
      with the fitted primitive and their region.

      \tparam PrimitiveAndRegionOutputIterator
      a model of `OutputIterator` whose value type is `Primitive_and_region`

      \param regions
      an iterator of type `PrimitiveAndRegionOutputIterator`.

      \return past-the-end position in the output sequence
    */
    template<typename PrimitiveAndRegionOutputIterator>
    PrimitiveAndRegionOutputIterator detect(PrimitiveAndRegionOutputIterator regions) {
//      clear(); TODO: this is not valid to comment this clear()
      m_visited_map.clear(); // tmp replacement for the line above

      Region region;
      std::size_t idx = 0;

      // Grow regions.
      for (auto it = m_seed_range.begin(); it != m_seed_range.end(); it++) {
        const Item seed = *it;

        // Try to grow a new region from the index of the seed item.
        if (!get(m_visited, seed)) {
          typename Region_type::Primitive primitive = m_region_type.primitive();
          const bool is_success = propagate(seed, region, primitive);

          // Check global conditions.
          if (!is_success || !m_region_type.is_valid_region(region)) {
            revert(region);
          } else {
            *(regions++) = std::pair<typename RegionType::Primitive, Region>(primitive, region);
            fill_region_map(idx++, region);
          }
        }
      }

      return regions;
    }

    /*!
      \brief provides a property map that provides the region index (or std::size_t(-1)) for each input element.

      \return Property map that maps each iterator of the input range to a region index.
    */

    const Region_map &region_map() {
      return m_region_map;
    }

    /// @}

    /// \name Unassigned Items
    /// @{

    /*!
      \brief fills an output iterator with all unassigned items.

      \tparam ItemOutputIterator
      a model of `OutputIterator` whose value type is `Item`

      \param output
      an iterator of type `PrimitiveAndRegionOutputIterator`.

      \return past-the-end position in the output sequence
    */
    template<typename InputRange, typename ItemOutputIterator>
    ItemOutputIterator unassigned_items(const InputRange& input_range, ItemOutputIterator output) const {
      for (auto it = input_range.begin(); it != input_range.end(); it++) {
        Item i = internal::conditional_deref<typename InputRange::const_iterator, Item>()(it);
        if (!get(m_visited, i))
          *(output++) = i;
      }
      return output;
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    template <class InputRange>
    void clear(const InputRange& input_range) {
      for (auto it = input_range.begin(); it != input_range.end(); it++) {
        put(m_region_map, internal::conditional_deref<typename InputRange::const_iterator, typename Region_map::key_type>()(it), std::size_t(-1));
      }
      m_visited_map.clear();
    }
    /// \endcond

  private:
    Neighbor_query& m_neighbor_query;
    Region_type& m_region_type;
    Region_map m_region_map;
    std::vector<Item> m_seed_range;

    using VisitedMap = std::unordered_set<typename Region_type::Item, internal::hash_item<typename Region_type::Item> >;
    VisitedMap m_visited_map;
    Boolean_property_map<VisitedMap> m_visited;

    void fill_region_map(std::size_t idx, const Region& region) {
      for (auto item : region) {
        put(m_region_map, internal::conditional_deref<Item, typename Region_map::key_type>()(item), idx);
      }
    }

    bool propagate(const Item &seed, Region& region, typename RegionType::Primitive &primitive) {
      region.clear();

      // Use two queues, while running on this queue, push to the other queue;
      // When the queue is done, update the shape of the current region and swap to the other queue;
      // depth_index is the index of the queue we are using.
      Running_queue running_queue[2];
      bool depth_index = 0;

      // Once the index of an item is pushed to the queue, it is pushed to the region too.
      put(m_visited, seed, true);
      running_queue[depth_index].push(seed);
      region.push_back(seed);

      // Update internal properties of the region.
      const bool is_well_created = m_region_type.update(region);
      primitive = m_region_type.primitive();
      if (!is_well_created) return false;

      bool grown = true;
      std::vector<std::pair<const Item, const Item> > rejected, rejected_swap;
      while (grown) {
        grown = false;

        Region neighbors;
        while (
          !running_queue[depth_index].empty() ||
          !running_queue[!depth_index].empty()) {

          while (!running_queue[depth_index].empty()) {

            // Call the next item index of the queue and remove it from the queue.
            const Item item = running_queue[depth_index].front();
            running_queue[depth_index].pop();

            // Get neighbors of the current item.
            neighbors.clear();
            m_neighbor_query(item, neighbors);

            // Visit all found neighbors.
            for (const Item neighbor : neighbors) {

              if (!get(m_visited, neighbor)) {
                if (m_region_type.is_part_of_region(item, neighbor, region)) {

                  // Add this neighbor to the other queue so that we can visit it later.
                  put(m_visited, neighbor, true);
                  running_queue[!depth_index].push(neighbor);
                  region.push_back(neighbor);
                  grown = true;
                }
                else {
                  // Add this neighbor to the rejected queue so I won't be checked again before refitting the primitive.
                  put(m_visited, neighbor, true);
                  rejected.push_back(std::pair<const Item, const Item>(item, neighbor));
                }
              }
            }
          }
          depth_index = !depth_index;
        }

        // Update internal properties of the region.
        // The region expanded with the current primitive to its largest extent.
        // After refitting the growing may continue, but it is only continued if the refitted primitive still fits all elements of the region.
        if (grown) {
          m_region_type.update(region);

          // Verify that associated elements are still within the tolerance.
          bool fits = true;
          Item former = region.front();
          for (const Item item : region) {
            if (!m_region_type.is_part_of_region(former, item, region)) {
              fits = false;
              break;
            }
            former = item;
          }

          // The refitted primitive does not fit all elements of the region, so the growing stops here.
          if (!fits) {
            // Reset visited flags for items that were rejected
            for (const std::pair<const Item, const Item>& p : rejected)
              put(m_visited, p.second, false);
            return true;
          }

          // If it fits, update the primitive.
          primitive = m_region_type.primitive();

          // Try to continue growing the region by considering formerly rejected elements.
          for (const std::pair<const Item, const Item>& p : rejected) {
            if (m_region_type.is_part_of_region(p.first, p.second, region)) {

              // Add this neighbor to the other queue so that we can visit it later.
              put(m_visited, p.second, true);
              running_queue[depth_index].push(p.second);
              region.push_back(p.second);
            }
            else rejected_swap.push_back(p);
          }
          rejected.clear();
          rejected.swap(rejected_swap);
        }
      }

      // Reset visited flags for items that were rejected
      for (const std::pair<const Item, const Item>& p : rejected)
        put(m_visited, p.second, false);

      return true;
    }

    void revert(const Region& region) {
      for (const Item item : region)
        put(m_visited, item, false);
    }
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_H
