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
#include <unordered_set>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/type_traits/is_iterator.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {

namespace internal {
  template <typename RegionType, typename RegionMap,
            bool b=std::is_same<RegionMap, typename RegionType::Region_index_map>::value>
  struct RM_creator{
    static RegionMap create(RegionType) { return RegionMap(); }
  };

  template <typename RegionType, typename RegionMap>
  struct RM_creator<RegionType, RegionMap, true>{
    static RegionMap create(RegionType& r ) { return r.region_index_map(); }
  };
}

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

    \tparam RegionMap a model of `ReadWritePropertyMap` whose key type is `Item`
    and value type is `std::size_t`.
  */
  template<
    typename NeighborQuery,
    typename RegionType,
    typename RegionMap = typename RegionType::Region_index_map>
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
    using Region_map = RegionMap;

  private:
    using Running_queue = std::queue<Item>;

  public:
    /// \name Initialization (RegionMap is the default type)
    /// @{

    /*!
      \brief initializes the region growing algorithm.

      \tparam InputRange
        a model of `ConstRange`

      \tparam ItemMap
        a model of `ReadablePropertyMap` with `InputRange::const_iterator` as key type and `Item` as value type.
        A default can be deduced using the value type of `InputRange` and `Item` to be
        either `CGAL::Dereference_property_map` or `CGAL::Identity_property_map`.

      \param input_range
      a range of input items for region growing.

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access item's neighbors

      \param region_type
      an instance of `RegionType` that is used internally to
      control if items form a valid region type

      \param item_map
        an instance of the property map to retrieve items from input values

      \pre `input_range.size() > 0`
    */
    template<typename InputRange, typename ItemMap = Default>
    Region_growing(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      RegionType& region_type,
      ItemMap item_map = ItemMap()
#ifndef DOXYGEN_RUNNING
      , std::enable_if_t<!std::is_same<ItemMap, RegionMap>::value>* = 0
#endif
      ) :
      m_neighbor_query(neighbor_query),
      m_region_type(region_type),
      m_region_map(internal::RM_creator<Region_type, RegionMap>::create(region_type)),
      m_visited(m_visited_map)
    {
      CGAL_precondition(input_range.size() > 0);
      m_seed_range.resize(input_range.size());

      using Item_helper = internal::Item_map_helper<ItemMap, Item, typename InputRange::const_iterator>;
      using Item_map = typename Item_helper::type;
      Item_map item_map_ = Item_helper::get(item_map);

      std::size_t idx = 0;
      for (auto it = input_range.begin(); it != input_range.end(); it++)
        m_seed_range[idx++] = get(item_map_, it);

      clear(input_range, item_map_);
    }

    /*!
      \brief initializes the region growing algorithm.

      \tparam InputRange
        a model of `ConstRange`

      \tparam SeedRange
        a model of `ConstRange` with `Item` as value type

      \tparam ItemMap
        a model of `ReadablePropertyMap` with `InputRange::const_iterator` as key type and `Item` as value type.
        A default can be deduced using the value type of `InputRange` and `Item` to be
        either `CGAL::Dereference_property_map` or `CGAL::Identity_property_map`.

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

      \param item_map
        an instance of the property map to retrieve items from input values

      \pre `input_range.size() > 0`
    */
    template<typename InputRange, typename SeedRange, typename ItemMap = Default>
    Region_growing(
      const InputRange& input_range,
      SeedRange& seed_range,
      NeighborQuery& neighbor_query,
      RegionType& region_type,
      ItemMap item_map = ItemMap()
#ifndef DOXYGEN_RUNNING
      , std::enable_if_t<!std::is_same<ItemMap, RegionMap>::value>* = 0
#endif
      ) :
      m_neighbor_query(neighbor_query),
      m_region_type(region_type),
      m_region_map(internal::RM_creator<Region_type, RegionMap>::create(region_type)),
      m_visited(m_visited_map) {

      CGAL_precondition(input_range.size() > 0);
      CGAL_precondition(seed_range.size() > 0);
      m_seed_range.resize(seed_range.size());

      using Item_helper = internal::Item_map_helper<ItemMap, Item, typename InputRange::const_iterator>;
      using Item_map = typename Item_helper::type;
      Item_map item_map_ = Item_helper::get(item_map);

      std::size_t idx = 0;
      for (auto it = seed_range.begin(); it != seed_range.end(); it++)
        m_seed_range[idx++] = *it;

      clear(input_range, item_map_);
    }
    /// @}

    /// \name Initialization (RegionMap is a user provided type)
    /// @{

    /*!
      \brief initializes the region growing algorithm.

      \tparam InputRange
        a model of `ConstRange`

      \tparam ItemMap
        a model of `ReadablePropertyMap` with `InputRange::const_iterator` as key type and `Item` as value type.
        A default can be deduced using the value type of `InputRange` and `Item` to be
        either `CGAL::Dereference_property_map` or `CGAL::Identity_property_map`.

      \param input_range
      a range of input items for region growing.

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access item's neighbors

      \param region_type
      an instance of `RegionType` that is used internally to
      control if items form a valid region type

      \param item_map
        an instance of the property map to retrieve items from input values

      \param rm external property map that will be filled when calling `detect()`

      \pre `input_range.size() > 0`
    */
    template<typename InputRange, typename ItemMap = Default>
    Region_growing(
      const InputRange& input_range,
      NeighborQuery& neighbor_query,
      RegionType& region_type,
      Region_map rm,
      ItemMap item_map = ItemMap()) :
      m_neighbor_query(neighbor_query),
      m_region_type(region_type),
      m_region_map(rm),
      m_visited(m_visited_map)
    {
      CGAL_precondition(input_range.size() > 0);
      m_seed_range.resize(input_range.size());

      using Item_helper = internal::Item_map_helper<ItemMap, Item, typename InputRange::const_iterator>;
      using Item_map = typename Item_helper::type;
      Item_map item_map_ = Item_helper::get(item_map);

      std::size_t idx = 0;
      for (auto it = input_range.begin(); it != input_range.end(); it++)
        m_seed_range[idx++] = get(item_map_, it);

      clear(input_range, item_map_);
    }

    /*!
      \brief initializes the region growing algorithm.

      \tparam InputRange
        a model of `ConstRange`

      \tparam SeedRange
        a model of `ConstRange` with `Item` as value type

      \tparam ItemMap
        a model of `ReadablePropertyMap` with `InputRange::const_iterator` as key type and `Item` as value type.
        A default can be deduced using the value type of `InputRange` and `Item` to be
        either `CGAL::Dereference_property_map` or `CGAL::Identity_property_map`.

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

      \param item_map
        an instance of the property map to retrieve items from input values

      \param rm external property map that will be filled when calling `detect()`

      \pre `input_range.size() > 0`
    */
    template<typename InputRange, typename SeedRange, typename ItemMap = Default>
    Region_growing(
      const InputRange& input_range,
      SeedRange& seed_range,
      NeighborQuery& neighbor_query,
      RegionType& region_type,
      Region_map rm,
      ItemMap item_map = ItemMap()) :
      m_neighbor_query(neighbor_query),
      m_region_type(region_type),
      m_region_map(rm),
      m_visited(m_visited_map) {

      CGAL_precondition(input_range.size() > 0);
      CGAL_precondition(seed_range.size() > 0);
      m_seed_range.resize(seed_range.size());

      using Item_helper = internal::Item_map_helper<ItemMap, Item, typename InputRange::const_iterator>;
      using Item_map = typename Item_helper::type;
      Item_map item_map_ = Item_helper::get(item_map);

      std::size_t idx = 0;
      for (auto it = seed_range.begin(); it != seed_range.end(); it++)
        m_seed_range[idx++] = *it;

      clear(input_range, item_map_);
    }

    /// @}

    /// \name Detection
    /// @{

    /*!
      \brief runs the region growing algorithm and fills an output iterator
      with the fitted primitive and their region.

      \tparam PrimitiveAndRegionOutputIterator
      a model of `OutputIterator` whose value type is `Primitive_and_region`

      \param region_out
      an output iterator of type `PrimitiveAndRegionOutputIterator`.

      \return past-the-end position in the output sequence
    */
    template<typename PrimitiveAndRegionOutputIterator = Emptyset_iterator>
    PrimitiveAndRegionOutputIterator detect(PrimitiveAndRegionOutputIterator region_out = PrimitiveAndRegionOutputIterator()) {
      //      clear(); TODO: this is not valid to comment this clear()
      m_visited_map.clear(); // tmp replacement for the line above

      Region region;
      m_nb_regions = 0;

      // Grow regions.
      for (auto it = m_seed_range.begin(); it != m_seed_range.end(); it++) {
        const Item seed = *it;

        // Try to grow a new region from the index of the seed item.
        if (!get(m_visited, seed)) {
          const bool is_success = propagate(seed, region);

          // Check global conditions.
          if (!is_success || !m_region_type.is_valid_region(region)) {
            revert(region);
          }
          else {
            fill_region_map(m_nb_regions++, region);
            if (!std::is_same<PrimitiveAndRegionOutputIterator, Emptyset_iterator>::value)
              *region_out++ = std::make_pair(m_region_type.primitive(), std::move(region));
          }
        }
      }

      return region_out;
    }

    /*!
      \brief provides a property map that provides the region index (or std::size_t(-1)) for each input element.

      \return Property map that maps each iterator of the input range to a region index.
    */
    const Region_map& region_map() {
      return m_region_map;
    }

    /// @}

    /// \name Unassigned Items
    /// @{

    /*!
      \brief fills an output iterator with all unassigned items.

      \tparam ItemOutputIterator
      a model of `OutputIterator` whose value type is `Item`

      \tparam InputRange
        a model of `ConstRange`

      \tparam ItemMap
        a model of `ReadablePropertyMap` with `InputRange::const_iterator` as key type and `Item` as value type.
        A default can be deduced using the value type of `InputRange` and `Item` to be
        either `CGAL::Dereference_property_map` or `CGAL::Identity_property_map`.

      \param input_range
      a range of input items for region growing

      \param output
      an iterator of type `PrimitiveAndRegionOutputIterator`.

      \param item_map
        an instance of the property map to retrieve items from input values

      \return past-the-end position in the output sequence
    */
    template<typename InputRange, typename ItemOutputIterator, typename ItemMap = Default>
    ItemOutputIterator unassigned_items(const InputRange& input_range, ItemOutputIterator output, ItemMap item_map = ItemMap()) const
    {
      using Item_helper = internal::Item_map_helper<ItemMap, Item, typename InputRange::const_iterator>;
      using Item_map = typename Item_helper::type;
      Item_map item_map_ = Item_helper::get(item_map);

      for (auto it = input_range.begin(); it != input_range.end(); it++) {
        Item i = get(item_map_,it);
        if (!get(m_visited, i))
          *(output++) = i;
      }
      return output;
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    template <class InputRange, class ItemMap = Default>
    void clear(const InputRange& input_range, ItemMap item_map = ItemMap()) {
      using Item_helper = internal::Item_map_helper<ItemMap, Item, typename InputRange::const_iterator>;
      using Item_map = typename Item_helper::type;
      Item_map item_map_ = Item_helper::get(item_map);

      m_nb_regions = 0;
      typename boost::property_traits<Region_map>::value_type init_value(-1);
      for (auto it = input_range.begin(); it != input_range.end(); it++) {
        Item item = get(item_map_, it);
        put(m_region_map, item, init_value);
      }
      // TODO if we want to allow subranges while NeighborQuery operates on the full range
      // (like for faces in a PolygonMesh) we should fill a non-visited map rather than a visited map
      m_visited_map.clear();
    }

    std::size_t number_of_regions_detected() const
    {
      return m_nb_regions;
    }
    /// \endcond

  private:
    Neighbor_query& m_neighbor_query;
    Region_type& m_region_type;
    Region_map m_region_map;

    std::vector<Item> m_seed_range;
    std::size_t m_nb_regions = 0;

    using VisitedMap = std::unordered_set<typename Region_type::Item, internal::hash_item<typename Region_type::Item> >;
    VisitedMap m_visited_map;
    Boolean_property_map<VisitedMap> m_visited;

    void fill_region_map(std::size_t idx, const Region& region) {
      typedef typename boost::property_traits<Region_map>::value_type Id;
      for (auto item : region) {
        put(m_region_map, item, static_cast<Id>(idx));
      }
    }

    bool propagate(const Item &seed, Region& region) {
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
            for (Item neighbor : neighbors) {

              if (!get(m_visited, neighbor)) {
                if (m_region_type.is_part_of_region(neighbor, region)) {

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
          for (Item item : region) {
            if (!m_region_type.is_part_of_region(item, region)) {
              fits = false;
              break;
            }
          }

          // The refitted primitive does not fit all elements of the region, so the growing stops here.
          if (!fits) {
            // Reset visited flags for items that were rejected
            for (const std::pair<const Item, const Item>& p : rejected)
              put(m_visited, p.second, false);
            return true;
          }

          // Try to continue growing the region by considering formerly rejected elements.
          for (const std::pair<const Item, const Item>& p : rejected) {
            if (m_region_type.is_part_of_region(p.second, region)) {

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
      for (Item item : region)
        put(m_visited, item, false);
    }
  };

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_H
