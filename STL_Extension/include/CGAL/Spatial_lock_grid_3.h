// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Clement Jamin

#ifndef CGAL_STL_EXTENSION_SPATIAL_LOCK_GRID_3_H
#define CGAL_STL_EXTENSION_SPATIAL_LOCK_GRID_3_H

#ifdef CGAL_LINKED_WITH_TBB

#include <CGAL/Bbox_3.h>
#include <CGAL/number_utils.h>

#include <atomic>
#include <thread>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/cache_aligned_allocator.h>

#include <algorithm>
#include <array>
#include <limits>
#include <vector>

namespace CGAL {

struct Tag_no_lock {};
struct Tag_non_blocking {};
struct Tag_priority_blocking {};

//*****************************************************************************
// class Spatial_lock_grid_base_3
// (Uses Curiously recurring template pattern)
//*****************************************************************************

template <typename Derived>
class Spatial_lock_grid_base_3
{
private:
  static bool *init_TLS_grid(int num_cells_per_axis)
  {
    int num_cells = num_cells_per_axis*
      num_cells_per_axis*num_cells_per_axis;
    bool *local_grid = new bool[num_cells];
    for (int i = 0 ; i < num_cells ; ++i)
      local_grid[i] = false;
    return local_grid;
  }

public:
  bool *get_thread_local_grid()
  {
    return m_tls_grids.local();
  }

  void set_bbox(const Bbox_3 &bbox)
  {
    // Compute resolutions
    m_bbox = bbox;
    double n = static_cast<double>(m_num_grid_cells_per_axis);
    m_resolution_x = n / (bbox.xmax() - bbox.xmin());
    m_resolution_y = n / (bbox.ymax() - bbox.ymin());
    m_resolution_z = n / (bbox.zmax() - bbox.zmin());

#ifdef CGAL_CONCURRENT_MESH_3_VERBOSE
    std::cerr << "Locking data structure Bounding Box = "
      << "[" << bbox.xmin() << ", " << bbox.xmax() << "], "
      << "[" << bbox.ymin() << ", " << bbox.ymax() << "], "
      << "[" << bbox.zmin() << ", " << bbox.zmax() << "]"
      << std::endl;
#endif
  }

  const Bbox_3 &get_bbox() const
  {
    return m_bbox;
  }

  bool is_locked_by_this_thread(int cell_index)
  {
    return get_thread_local_grid()[cell_index];
  }

  template <typename P3>
  bool is_locked(const P3 &point)
  {
    return is_cell_locked(grid_index(point));
  }

  template <typename P3>
  bool is_locked_by_this_thread(const P3 &point)
  {
    return get_thread_local_grid()[grid_index(point)];
  }

  template <bool no_spin = false>
  bool try_lock(int cell_index)
  {
    return get_thread_local_grid()[cell_index]
        || try_lock_cell<no_spin>(cell_index);
  }

  template <bool no_spin = false>
  bool try_lock(int index_x, int index_y, int index_z, int lock_radius)
  {
    if (lock_radius == 0)
    {
      int index_to_lock = grid_index(index_x, index_y, index_z);
      return try_lock<no_spin>(index_to_lock);
    }
    else
    {
      // We have to lock the square
      std::vector<int> locked_cells_tmp;

      // For each cell inside the square
      for (int i = (std::max)(0, index_x-lock_radius) ;
           i <= (std::min)(m_num_grid_cells_per_axis - 1, index_x+lock_radius) ;
           ++i)
      {
        for (int j = (std::max)(0, index_y-lock_radius) ;
             j <= (std::min)(m_num_grid_cells_per_axis - 1, index_y+lock_radius) ;
             ++j)
        {
          for (int k = (std::max)(0, index_z-lock_radius) ;
               k <= (std::min)(m_num_grid_cells_per_axis - 1, index_z+lock_radius) ;
               ++k)
          {
            int index_to_lock = grid_index(i, j, k);
            // Try to lock it
            if (try_lock<no_spin>(index_to_lock))
            {
              locked_cells_tmp.push_back(index_to_lock);
            }
            else
            {
              // failed => we unlock already locked cells and return false
              for(int cell_index : locked_cells_tmp)
              {
                unlock(cell_index);
              }
              return false;
            }
          }
        }
      }

      return true;
    }
  }


  template <bool no_spin = false>
  bool try_lock(int cell_index, int lock_radius)
  {
    if (lock_radius == 0)
    {
      return try_lock<no_spin>(cell_index);
    }
    else
    {
      int index_z = cell_index/(m_num_grid_cells_per_axis*m_num_grid_cells_per_axis);
      cell_index -= index_z*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis;
      int index_y = cell_index/m_num_grid_cells_per_axis;
      cell_index -= index_y*m_num_grid_cells_per_axis;
      int index_x = cell_index;

      return try_lock<no_spin>(index_x, index_y, index_z, lock_radius);
    }
  }

  // P3 must provide .x(), .y(), .z()
  template <bool no_spin = false, typename P3>
  bool try_lock(const P3 &point, int lock_radius = 0)
  {
    auto [index_x, index_y, index_z] = get_grid_indices(point);

    if (lock_radius == 0)
    {
      return try_lock<no_spin>(grid_index(index_x, index_y, index_z));
    }
    else
    {
      return try_lock<no_spin>(index_x, index_y, index_z, lock_radius);
    }
  }

  void unlock(int cell_index)
  {
    // Unlock lock and shared grid
    unlock_cell(cell_index);
    get_thread_local_grid()[cell_index] = false;
  }

  void unlock_all_points_locked_by_this_thread()
  {
    std::vector<int> &tls_locked_cells = m_tls_locked_cells.local();
    for(int cell_index : tls_locked_cells)
    {
      // If we still own the lock
      if (get_thread_local_grid()[cell_index] == true)
        unlock(cell_index);
    }
    tls_locked_cells.clear();
  }

  void unlock_all_tls_locked_cells_but_one(int cell_index_to_keep_locked)
  {
    std::vector<int> &tls_locked_cells = m_tls_locked_cells.local();
    bool cell_to_keep_found = false;
    for(int cell_index : tls_locked_cells)
    {
      // If we still own the lock
      if (get_thread_local_grid()[cell_index] == true)
      {
        if (cell_index == cell_index_to_keep_locked)
          cell_to_keep_found = true;
        else
          unlock(cell_index);
      }
    }
    tls_locked_cells.clear();
    if (cell_to_keep_found)
      tls_locked_cells.push_back(cell_index_to_keep_locked);
  }

  template <typename P3>
  void unlock_all_tls_locked_locations_but_one_point(const P3 &point)
  {
    unlock_all_tls_locked_cells_but_one(grid_index(point));
  }

  bool check_if_all_cells_are_unlocked()
  {
    int num_cells = m_num_grid_cells_per_axis*
      m_num_grid_cells_per_axis*m_num_grid_cells_per_axis;
    bool unlocked = true;
    for (int i = 0 ; unlocked && i < num_cells ; ++i)
      unlocked = !is_cell_locked(i);
    return unlocked;
  }

  bool check_if_all_tls_cells_are_unlocked()
  {
    int num_cells = m_num_grid_cells_per_axis*
      m_num_grid_cells_per_axis*m_num_grid_cells_per_axis;
    bool unlocked = true;
    for (int i = 0 ; unlocked && i < num_cells ; ++i)
      unlocked = (get_thread_local_grid()[i] == false);
    return unlocked;
  }

protected:

  // Constructor
  Spatial_lock_grid_base_3(const Bbox_3 &bbox,
                                          int num_grid_cells_per_axis)
    : m_num_grid_cells_per_axis(num_grid_cells_per_axis),
      m_tls_grids([num_grid_cells_per_axis](){ return init_TLS_grid(num_grid_cells_per_axis); })
  {
    set_bbox(bbox);
  }

  /// Destructor
  ~Spatial_lock_grid_base_3()
  {
    for( auto* grid : m_tls_grids )
    {
      delete [] grid;
    }
  }

  template <typename P3>
  std::array<int, 3> get_grid_indices(const P3& point) const
  {
    int index_x = static_cast<int>( (CGAL::to_double(point.x()) - m_bbox.xmin()) * m_resolution_x);
    index_x = std::clamp(index_x, 0, m_num_grid_cells_per_axis - 1);
    int index_y = static_cast<int>( (CGAL::to_double(point.y()) - m_bbox.ymin()) * m_resolution_y);
    index_y = std::clamp(index_y, 0, m_num_grid_cells_per_axis - 1);
    int index_z = static_cast<int>( (CGAL::to_double(point.z()) - m_bbox.zmin()) * m_resolution_z);
    index_z = std::clamp(index_z, 0, m_num_grid_cells_per_axis - 1);

    return {index_x, index_y, index_z};
  }

  int grid_index(int index_x, int index_y, int index_z) const
  {
    return
      index_z*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis
      + index_y*m_num_grid_cells_per_axis
      + index_x;
  }

  template <typename P3>
  int grid_index(const P3& point) const
  {
    auto [index_x, index_y, index_z] = get_grid_indices(point);

    return grid_index(index_x, index_y, index_z);
  }

  auto* derived() { return static_cast<Derived*>(this); }

  bool is_cell_locked(int cell_index)
  {
    return derived()->is_cell_locked_impl(cell_index);
  }

  template <bool no_spin = false>
  bool try_lock_cell(int cell_index)
  {
    return derived()
      ->template try_lock_cell_impl<no_spin>(cell_index);
  }

  void unlock_cell(int cell_index)
  {
    derived()->unlock_cell_impl(cell_index);
  }

  int                                             m_num_grid_cells_per_axis;
  Bbox_3                                          m_bbox;
  double                                          m_resolution_x;
  double                                          m_resolution_y;
  double                                          m_resolution_z;

  // TLS
  typedef tbb::enumerable_thread_specific<
    bool*,
    tbb::cache_aligned_allocator<bool*>,
    tbb::ets_key_per_instance>                               TLS_grid;
  typedef tbb::enumerable_thread_specific<std::vector<int> > TLS_locked_cells;

  TLS_grid                                        m_tls_grids;
  TLS_locked_cells                                m_tls_locked_cells;
};



//*****************************************************************************
// class Spatial_lock_grid_3
//*****************************************************************************
template <typename Grid_lock_tag = Tag_priority_blocking>
class Spatial_lock_grid_3;


//*****************************************************************************
// class Spatial_lock_grid_3<Tag_non_blocking>
//*****************************************************************************
template <>
class Spatial_lock_grid_3<Tag_non_blocking>
  : public Spatial_lock_grid_base_3<
      Spatial_lock_grid_3<Tag_non_blocking> >
{
  typedef Spatial_lock_grid_base_3<
    Spatial_lock_grid_3<Tag_non_blocking> > Base;

public:
  // Constructors
  Spatial_lock_grid_3(const Bbox_3 &bbox, int num_grid_cells_per_axis)
  : Base(bbox, num_grid_cells_per_axis),
    m_grid(num_grid_cells_per_axis*num_grid_cells_per_axis*num_grid_cells_per_axis)
  {
    int num_cells =
      num_grid_cells_per_axis*num_grid_cells_per_axis*num_grid_cells_per_axis;

    // Initialize grid (useless?)
    for (int i = 0 ; i < num_cells ; ++i)
      m_grid[i] = false;
  }

  ~Spatial_lock_grid_3()
  {
  }

  bool is_cell_locked_impl(int cell_index)
  {
    return (m_grid[cell_index] == true);
  }

  template <bool no_spin>
  bool try_lock_cell_impl(int cell_index)
  {
    bool v1 = true, v2 = false;
    if(m_grid[cell_index].compare_exchange_strong(v2,v1))
    {
      get_thread_local_grid()[cell_index] = true;
      m_tls_locked_cells.local().push_back(cell_index);
      return true;
    }
    return false;
  }

  void unlock_cell_impl(int cell_index)
  {
    m_grid[cell_index] = false;
  }

protected:

  std::vector<std::atomic<bool> > m_grid;
};


//*****************************************************************************
// class Spatial_lock_grid_3<Tag_priority_blocking>
//*****************************************************************************

template <>
class Spatial_lock_grid_3<Tag_priority_blocking>
  : public Spatial_lock_grid_base_3<Spatial_lock_grid_3<Tag_priority_blocking> >
{
  using Self = Spatial_lock_grid_3<Tag_priority_blocking>;
  using Base = Spatial_lock_grid_base_3<Self>;

  using priority_t = unsigned int;

public:
  // Constructors

  Spatial_lock_grid_3(const Bbox_3 &bbox, int num_grid_cells_per_axis)
  : Base(bbox, num_grid_cells_per_axis),
    m_grid(num_grid_cells_per_axis*num_grid_cells_per_axis*num_grid_cells_per_axis),
    m_tls_thread_priorities(init_TLS_thread_priorities)
  {
    // Explicitly initialize the atomics
    for ( auto& cell : m_grid )
      cell = 0;
  }

  /// Destructor
  ~Spatial_lock_grid_3()
  {
  }

  bool is_cell_locked_impl(int cell_index)
  {
    return (m_grid[cell_index] != 0);
  }

  template <bool no_spin>
  bool try_lock_cell_impl(int cell_index)
  {
    const priority_t this_thread_priority = m_tls_thread_priorities.local();

    // NO SPIN
    if (no_spin)
    {
      priority_t old_value = 0;
      if(m_grid[cell_index].compare_exchange_strong(old_value, this_thread_priority))
      {
        get_thread_local_grid()[cell_index] = true;
        m_tls_locked_cells.local().push_back(cell_index);
        return true;
      }
    }
    // SPIN
    else
    {
      for(;;)
      {
        priority_t old_value = 0;
        if(m_grid[cell_index].compare_exchange_weak(old_value, this_thread_priority))
        {
          get_thread_local_grid()[cell_index] = true;
          m_tls_locked_cells.local().push_back(cell_index);
          return true;
        }
        else if (old_value > this_thread_priority)
        {
          // Another "more priority" thread owns the lock, we back off
          return false;
        }
        else
        {
          std::this_thread::yield();
        }
      }
    }

    return false;
  }

  void unlock_cell_impl(int cell_index)
  {
    m_grid[cell_index] = 0;
  }

private:
  static priority_t init_TLS_thread_priorities()
  {
    static std::atomic<priority_t> last_id;
    priority_t id = ++last_id;
    // Ensure it is > 0
    return (1 + id%((std::numeric_limits<priority_t>::max)()));
  }

protected:

  std::vector<std::atomic<priority_t>>               m_grid;

  using TLS_thread_uint_ids = tbb::enumerable_thread_specific<priority_t>;
  TLS_thread_uint_ids                                m_tls_thread_priorities;
};

} //namespace CGAL

#else // !CGAL_LINKED_WITH_TBB

namespace CGAL {

template <typename Grid_lock_tag = void>
class Spatial_lock_grid_3
{
};

}

#endif // CGAL_LINKED_WITH_TBB

#endif // CGAL_STL_EXTENSION_SPATIAL_LOCK_GRID_3_H
