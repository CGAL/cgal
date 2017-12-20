// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Clement Jamin

#ifndef CGAL_STL_EXTENSION_SPATIAL_LOCK_GRID_3_H
#define CGAL_STL_EXTENSION_SPATIAL_LOCK_GRID_3_H

#ifdef CGAL_LINKED_WITH_TBB

#include <CGAL/Bbox_3.h>

#include <boost/bind.hpp>

#include <tbb/atomic.h>
#if TBB_IMPLEMENT_CPP0X
# include <tbb/compat/thread>
#else
# include <thread>
#endif
#include <tbb/enumerable_thread_specific.h>

#include <algorithm>
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
    return is_cell_locked(get_grid_index(point));
  }

  template <typename P3>
  bool is_locked_by_this_thread(const P3 &point)
  {
    return get_thread_local_grid()[get_grid_index(point)];
  }

  bool try_lock(int cell_index)
  {
    return try_lock<false>(cell_index);
  }

  template <bool no_spin>
  bool try_lock(int cell_index)
  {
    return get_thread_local_grid()[cell_index]
        || try_lock_cell<no_spin>(cell_index);
  }


  bool try_lock(int index_x, int index_y, int index_z, int lock_radius)
  {
    return try_lock<false>(index_x, index_y, index_z, lock_radius);
  }

  template <bool no_spin>
  bool try_lock(int index_x, int index_y, int index_z, int lock_radius)
  {
    if (lock_radius == 0)
    {
      int index_to_lock =
        index_z*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis
        + index_y*m_num_grid_cells_per_axis
        + index_x;
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
            int index_to_lock =
              k*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis
              + j*m_num_grid_cells_per_axis
              + i;
            // Try to lock it
            if (try_lock<no_spin>(index_to_lock))
            {
              locked_cells_tmp.push_back(index_to_lock);
            }
            else
            {
              // failed => we unlock already locked cells and return false
              std::vector<int>::const_iterator it = locked_cells_tmp.begin();
              std::vector<int>::const_iterator it_end = locked_cells_tmp.end();
              for( ; it != it_end ; ++it)
              {
                unlock(*it);
              }
              return false;
            }
          }
        }
      }

      return true;
    }
  }


  bool try_lock(int cell_index, int lock_radius)
  {
    return try_lock<false>(cell_index, lock_radius);
  }

  template <bool no_spin>
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
  template <typename P3>
  bool try_lock(const P3 &point, int lock_radius = 0)
  {
    return try_lock<false, P3>(point, lock_radius);
  }

  // P3 must provide .x(), .y(), .z()
  template <bool no_spin, typename P3>
  bool try_lock(const P3 &point, int lock_radius = 0)
  {
    // Compute index on grid
    int index_x = static_cast<int>( (CGAL::to_double(point.x()) - m_bbox.xmin()) * m_resolution_x);
    //index_x = std::max( 0, std::min(index_x, m_num_grid_cells_per_axis - 1) );
    index_x =
      (index_x < 0 ?
        0
        : (index_x >= m_num_grid_cells_per_axis ?
            m_num_grid_cells_per_axis - 1
            : index_x
          )
      );
    int index_y = static_cast<int>( (CGAL::to_double(point.y()) - m_bbox.ymin()) * m_resolution_y);
    //index_y = std::max( 0, std::min(index_y, m_num_grid_cells_per_axis - 1) );
    index_y =
      (index_y < 0 ?
        0
        : (index_y >= m_num_grid_cells_per_axis ?
            m_num_grid_cells_per_axis - 1
            : index_y
          )
      );
    int index_z = static_cast<int>( (CGAL::to_double(point.z()) - m_bbox.zmin()) * m_resolution_z);
    //index_z = std::max( 0, std::min(index_z, m_num_grid_cells_per_axis - 1) );
    index_z =
      (index_z < 0 ?
        0
        : (index_z >= m_num_grid_cells_per_axis ?
            m_num_grid_cells_per_axis - 1
            : index_z
          )
      );

    if (lock_radius == 0)
    {
      int index =
        index_z*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis
        + index_y*m_num_grid_cells_per_axis
        + index_x;
      return try_lock<no_spin>(index);
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
    std::vector<int>::const_iterator it = tls_locked_cells.begin();
    std::vector<int>::const_iterator it_end = tls_locked_cells.end();
    for( ; it != it_end ; ++it)
    {
      // If we still own the lock
      int cell_index = *it;
      if (get_thread_local_grid()[cell_index] == true)
        unlock(cell_index);
    }
    tls_locked_cells.clear();
  }

  void unlock_all_tls_locked_cells_but_one(int cell_index_to_keep_locked)
  {
    std::vector<int> &tls_locked_cells = m_tls_locked_cells.local();
    std::vector<int>::const_iterator it = tls_locked_cells.begin();
    std::vector<int>::const_iterator it_end = tls_locked_cells.end();
    bool cell_to_keep_found = false;
    for( ; it != it_end ; ++it)
    {
      // If we still own the lock
      int cell_index = *it;
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
    unlock_all_tls_locked_cells_but_one(get_grid_index(point));
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
      m_tls_grids(boost::bind(init_TLS_grid, num_grid_cells_per_axis))
  {
    set_bbox(bbox);
  }

  /// Destructor
  ~Spatial_lock_grid_base_3()
  {
    for( TLS_grid::iterator it_grid = m_tls_grids.begin() ;
             it_grid != m_tls_grids.end() ;
             ++it_grid )
    {
      delete [] *it_grid;
    }
  }

  template <typename P3>
  int get_grid_index(const P3& point) const
  {
    // Compute indices on grid
    int index_x = static_cast<int>( (CGAL::to_double(point.x()) - m_bbox.xmin()) * m_resolution_x);
    //index_x = std::max( 0, std::min(index_x, m_num_grid_cells_per_axis - 1) );
    index_x =
      (index_x < 0 ?
        0
        : (index_x >= m_num_grid_cells_per_axis ?
            m_num_grid_cells_per_axis - 1
            : index_x
          )
      );
    int index_y = static_cast<int>( (CGAL::to_double(point.y()) - m_bbox.ymin()) * m_resolution_y);
    //index_y = std::max( 0, std::min(index_y, m_num_grid_cells_per_axis - 1) );
    index_y =
      (index_y < 0 ?
        0
        : (index_y >= m_num_grid_cells_per_axis ?
            m_num_grid_cells_per_axis - 1
            : index_y
          )
      );
    int index_z = static_cast<int>( (CGAL::to_double(point.z()) - m_bbox.zmin()) * m_resolution_z);
    //index_z = std::max( 0, std::min(index_z, m_num_grid_cells_per_axis - 1) );
    index_z =
      (index_z < 0 ?
        0
        : (index_z >= m_num_grid_cells_per_axis ?
            m_num_grid_cells_per_axis - 1
            : index_z
          )
      );

    return
      index_z*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis
      + index_y*m_num_grid_cells_per_axis
      + index_x;
  }

  bool is_cell_locked(int cell_index)
  {
    return static_cast<Derived*>(this)->is_cell_locked_impl(cell_index);
  }

  bool try_lock_cell(int cell_index)
  {
    return try_lock_cell<false>(cell_index);
  }

  template <bool no_spin>
  bool try_lock_cell(int cell_index)
  {
    return static_cast<Derived*>(this)
      ->template try_lock_cell_impl<no_spin>(cell_index);
  }
  void unlock_cell(int cell_index)
  {
    static_cast<Derived*>(this)->unlock_cell_impl(cell_index);
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
  : Base(bbox, num_grid_cells_per_axis)
  {
    int num_cells =
      num_grid_cells_per_axis*num_grid_cells_per_axis*num_grid_cells_per_axis;

    m_grid.resize(num_cells);
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
    bool old_value = m_grid[cell_index].compare_and_swap(true, false);
    if (old_value == false)
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

  std::vector<tbb::atomic<bool> > m_grid;
};


//*****************************************************************************
// class Spatial_lock_grid_3<Tag_priority_blocking>
//*****************************************************************************

template <>
class Spatial_lock_grid_3<Tag_priority_blocking>
  : public Spatial_lock_grid_base_3<Spatial_lock_grid_3<Tag_priority_blocking> >
{
  typedef Spatial_lock_grid_base_3<
    Spatial_lock_grid_3<Tag_priority_blocking> > Base;

public:
  // Constructors

  Spatial_lock_grid_3(const Bbox_3 &bbox, int num_grid_cells_per_axis)
  : Base(bbox, num_grid_cells_per_axis),
    m_tls_thread_priorities(init_TLS_thread_priorities)
  {
    int num_cells =
      num_grid_cells_per_axis*num_grid_cells_per_axis*num_grid_cells_per_axis;
    m_grid.resize(num_cells);
    // Explicitly initialize the atomics
    std::vector<tbb::atomic<unsigned int> >::iterator it     = m_grid.begin();
    std::vector<tbb::atomic<unsigned int> >::iterator it_end = m_grid.end();
    for ( ; it != it_end ; ++it)
      *it = 0;
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
    unsigned int this_thread_priority = m_tls_thread_priorities.local();

    // NO SPIN
    if (no_spin)
    {
      unsigned int old_value =
        m_grid[cell_index].compare_and_swap(this_thread_priority, 0);
      if (old_value == 0)
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
        unsigned int old_value =
          m_grid[cell_index].compare_and_swap(this_thread_priority, 0);
        if (old_value == 0)
        {
          get_thread_local_grid()[cell_index] = true;
          m_tls_locked_cells.local().push_back(cell_index);
          return true;
        }
        else if (old_value > this_thread_priority)
        {
          // Another "more prioritary" thread owns the lock, we back off
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
  static unsigned int init_TLS_thread_priorities()
  {
    static tbb::atomic<unsigned int> last_id;
    unsigned int id = ++last_id;
    // Ensure it is > 0
    return (1 + id%((std::numeric_limits<unsigned int>::max)()));
  }

protected:

  std::vector<tbb::atomic<unsigned int> >               m_grid;

  typedef tbb::enumerable_thread_specific<unsigned int> TLS_thread_uint_ids;
  TLS_thread_uint_ids                                   m_tls_thread_priorities;
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
