// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Clement Jamin

#ifdef CONCURRENT_MESH_3

#ifndef CGAL_MESH_3_LOCKING_DATA_STRUCTURES_H
#define CGAL_MESH_3_LOCKING_DATA_STRUCTURES_H

#include <CGAL/Mesh_3/Concurrent_mesher_config.h>

#include <CGAL/Bbox_3.h>

#include <boost/bind.hpp>

#include <tbb/tbb.h>
#include <tbb/compat/thread>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/recursive_mutex.h>

#include <algorithm>

namespace CGAL {
namespace Mesh_3 {

//******************************************************************************
// class Grid_locking_ds_base
// (Uses Curiously recurring template pattern)
//******************************************************************************

bool *init_TLS_grid(int num_cells_per_axis)
{
  int num_cells = num_cells_per_axis*
    num_cells_per_axis*num_cells_per_axis;
  bool *local_grid = new bool[num_cells];
  for (int i = 0 ; i < num_cells ; ++i)
    local_grid[i] = false;
  return local_grid;
}

template <typename Derived>
class Grid_locking_ds_base
{
private:
  
public:
  bool try_lock(int cell_index)
  {
    bool ret = false;
    // Already locked by this thread?
    if (m_tls_grids.local()[cell_index])
    {
      ret = true;
    }
    // Otherwise, try to lock it
    else
    {
      if (try_lock_cell(cell_index))
      {
        ret = true;
        m_tls_grids.local()[cell_index] = true;
        m_tls_locked_cells.local().push_back(cell_index);
      }
    }
    return ret;
  }
  
  bool try_lock(int index_x, int index_y, int index_z, int lock_radius)
  {
    if (lock_radius == 0)
    {
      int index_to_lock = 
        index_z*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis
        + index_y*m_num_grid_cells_per_axis 
        + index_x;
      return try_lock(index_to_lock);
    }
    else
    {
      // We have to lock the square
      std::vector<int> locked_cells_tmp;

      // For each cell inside the square
      for (int i = std::max(0, index_x-lock_radius) ; 
           i <= std::min(m_num_grid_cells_per_axis - 1, index_x+lock_radius) ; 
           ++i)
      {
        for (int j = std::max(0, index_y-lock_radius) ; 
             j <= std::min(m_num_grid_cells_per_axis - 1, index_y+lock_radius) ; 
             ++j)
        {
          for (int k = std::max(0, index_z-lock_radius) ; 
               k <= std::min(m_num_grid_cells_per_axis - 1, index_z+lock_radius) ;
               ++k)
          {
            int index_to_lock = 
              k*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis
              + j*m_num_grid_cells_per_axis 
              + i;
            // Try to lock it
            if (try_lock(index_to_lock))
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
    if (lock_radius == 0)
    {
      return try_lock(cell_index);
    }
    else
    {
      int index_z = cell_index/(m_num_grid_cells_per_axis*m_num_grid_cells_per_axis);
      cell_index -= index_z*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis;
      int index_y = cell_index/m_num_grid_cells_per_axis;
      cell_index -= index_y*m_num_grid_cells_per_axis;
      int index_x = cell_index;

      return try_lock(index_x, index_y, index_z, lock_radius);
    }
  }

  /// P3 must provide .x(), .y(), .z()
  /// Returns a pair "success or not + index of the grid cell"
  template <typename P3>
  std::pair<bool, int> try_lock(const P3 &point, int lock_radius = 0)
  {
    // Compute indices on grid
    int index_x = static_cast<int>( (to_double(point.x()) - m_xmin) * m_resolution_x);
    index_x = std::max( 0, std::min(index_x, m_num_grid_cells_per_axis - 1) );
    int index_y = static_cast<int>( (to_double(point.y()) - m_ymin) * m_resolution_y);
    index_y = std::max( 0, std::min(index_y, m_num_grid_cells_per_axis - 1) );
    int index_z = static_cast<int>( (to_double(point.z()) - m_zmin) * m_resolution_z);
    index_z = std::max( 0, std::min(index_z, m_num_grid_cells_per_axis - 1) );
    
    int index = 
      index_z*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis
      + index_y*m_num_grid_cells_per_axis 
      + index_x;

    if (lock_radius == 0)
    {
      return std::make_pair(try_lock(index), index);
    }
    else
    {
      return std::make_pair(
        try_lock(index_x, index_y, index_z, lock_radius), 
        index);
    }
  }

  template <typename P3>
  void unlock(const P3 &point)
  {
    // Compute indices on grid
    int index_x = static_cast<int>( (point.x() - m_xmin) * m_resolution_x);
    index_x = std::max( 0, std::min(index_x, m_num_grid_cells_per_axis - 1) );
    int index_y = static_cast<int>( (point.y() - m_ymin) * m_resolution_y);
    index_y = std::max( 0, std::min(index_y, m_num_grid_cells_per_axis - 1) );
    int index_z = static_cast<int>( (point.z() - m_zmin) * m_resolution_z);
    index_z = std::max( 0, std::min(index_z, m_num_grid_cells_per_axis - 1) );
    
    int index = 
      index_z*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis
      + index_y*m_num_grid_cells_per_axis 
      + index_x;

    unlock(index);
  }

  void unlock(int cell_index)
  {
    // Unlock lock and shared grid
    unlock_cell(cell_index);
    m_tls_grids.local()[cell_index] = false;
  }
  
  void unlock_all_tls_locked_cells()
  {
    std::vector<int> &tls_locked_cells = m_tls_locked_cells.local();
    std::vector<int>::const_iterator it = tls_locked_cells.begin();
    std::vector<int>::const_iterator it_end = tls_locked_cells.end();
    for( ; it != it_end ; ++it)
    {
      // If we still own the lock
      int cell_index = *it;
      if (m_tls_grids.local()[cell_index] == true)
        unlock(cell_index);
    }
    tls_locked_cells.clear();
  }
  
  bool check_if_all_tls_cells_are_unlocked()
  {
    int num_cells = m_num_grid_cells_per_axis*
      m_num_grid_cells_per_axis*m_num_grid_cells_per_axis;
    bool unlocked = true;
    for (int i = 0 ; unlocked && i < num_cells ; ++i)
      unlocked = (m_tls_grids.local()[i] == false);
    return unlocked; 
  }

protected:
  
  // Constructor
  Grid_locking_ds_base(const Bbox_3 &bbox, 
                       int num_grid_cells_per_axis)
    : m_num_grid_cells_per_axis(num_grid_cells_per_axis),
      m_tls_grids(boost::bind(init_TLS_grid, num_grid_cells_per_axis))
  {
    // Keep mins and resolutions
    m_xmin = bbox.xmin();
    m_ymin = bbox.ymin();
    m_zmin = bbox.zmin();
    double n = static_cast<double>(num_grid_cells_per_axis);
    m_resolution_x = n / (bbox.xmax() - m_xmin);
    m_resolution_y = n / (bbox.ymax() - m_ymin);
    m_resolution_z = n / (bbox.zmax() - m_zmin);
  }

  /// Destructor
  ~Grid_locking_ds_base()
  {
    for( TLS_grid::iterator it_grid = m_tls_grids.begin() ; 
             it_grid != m_tls_grids.end() ; 
             ++it_grid )
    {
      delete [] *it_grid;
    }
  }

  bool try_lock_cell(int cell_index) 
  { 
    return static_cast<Derived*>(this)->try_lock_cell_impl(cell_index);
  }
  void unlock_cell(int cell_index) 
  { 
    static_cast<Derived*>(this)->unlock_cell_impl(cell_index);
  }

  int                                             m_num_grid_cells_per_axis;
  double                                          m_xmin;
  double                                          m_ymin;
  double                                          m_zmin;
  double                                          m_resolution_x;
  double                                          m_resolution_y;
  double                                          m_resolution_z;

  // TLS
  typedef tbb::enumerable_thread_specific<bool*>              TLS_grid;
  typedef tbb::enumerable_thread_specific<std::vector<int> >  TLS_locked_cells;

  TLS_grid                                        m_tls_grids;
  TLS_locked_cells                                m_tls_locked_cells;
};



//******************************************************************************
// class Simple_grid_locking_ds
//******************************************************************************

class Simple_grid_locking_ds
  : public Grid_locking_ds_base<Simple_grid_locking_ds>
{
public:
  
  // Constructors
  Simple_grid_locking_ds(const Bbox_3 &bbox, 
                         int num_grid_cells_per_axis)
  : Grid_locking_ds_base(bbox, num_grid_cells_per_axis)
  {
    int num_cells = 
      num_grid_cells_per_axis*num_grid_cells_per_axis*num_grid_cells_per_axis;
    
    m_grid = new tbb::atomic<bool>[num_cells];
    // Initialize grid
    for (int i = 0 ; i < num_cells ; ++i)
      m_grid[i] = false;
  }

  ~Simple_grid_locking_ds()
  {
    delete [] m_grid;
  }
  
  bool try_lock_cell_impl(int cell_index)
  {
    bool old_value = m_grid[cell_index].compare_and_swap(true, false);
    return (old_value == false);
  }
  
  void unlock_cell_impl(int cell_index)
  {
    m_grid[cell_index] = false;
  }

protected:

  tbb::atomic<bool> *                             m_grid;
};


//******************************************************************************
// class Simple_grid_locking_ds_with_thread_ids
//******************************************************************************

class Simple_grid_locking_ds_with_thread_ids
  : public Grid_locking_ds_base<Simple_grid_locking_ds_with_thread_ids>
{
public:
  // Constructors
  
  Simple_grid_locking_ds_with_thread_ids(const Bbox_3 &bbox, 
                                         int num_grid_cells_per_axis)
  : Grid_locking_ds_base(bbox, num_grid_cells_per_axis),
    m_tls_thread_ids(
      [=]() -> unsigned int // CJTODO: lambdas OK?
      {
        static unsigned int last_id = 0;
        return ++last_id;
      }
    )
  {
    int num_cells = 
      num_grid_cells_per_axis*num_grid_cells_per_axis*num_grid_cells_per_axis;
    m_grid = new tbb::atomic<unsigned int>[num_cells];
    // Initialize grid
    for (int i = 0 ; i < num_cells ; ++i)
      m_grid[i] = 0;
  }

  /// Destructor
  ~Simple_grid_locking_ds_with_thread_ids()
  {
    delete [] m_grid;
  }
  
  bool try_lock_cell_impl(int cell_index)
  {
    bool ret = false;

    unsigned int this_thread_id = m_tls_thread_ids.local();
    unsigned int old_value;
    do
    {
      old_value = m_grid[cell_index].compare_and_swap(
        this_thread_id, 0);
      if (old_value == 0)
      {
        ret = true;
        m_tls_grids.local()[cell_index] = true;
        m_tls_locked_cells.local().push_back(cell_index);
      }
      else
      {
        std::this_thread::yield();
      }
    } while (ret == false && old_value < this_thread_id);

    return ret;
  }
  
  void unlock_cell_impl(int cell_index)
  {
    m_grid[cell_index] = 0;
  }
  
protected:
  tbb::atomic<unsigned int> *                           m_grid;

  typedef tbb::enumerable_thread_specific<unsigned int> TLS_thread_uint_ids;
  TLS_thread_uint_ids                                   m_tls_thread_ids;
};

//******************************************************************************
// class Simple_grid_locking_ds_with_mutex
//******************************************************************************

class Simple_grid_locking_ds_with_mutex
  : public Grid_locking_ds_base<Simple_grid_locking_ds_with_mutex>
{
public:
  
  // Constructors
  Simple_grid_locking_ds_with_mutex(const Bbox_3 &bbox,
                                    int num_grid_cells_per_axis)
  : Grid_locking_ds_base(bbox, num_grid_cells_per_axis)
  {
    int num_cells = 
      num_grid_cells_per_axis*num_grid_cells_per_axis*num_grid_cells_per_axis;
    m_grid = new tbb::recursive_mutex[num_cells];
  }

  /// Destructor
  ~Simple_grid_locking_ds_with_mutex()
  {
    delete [] m_grid;
  }
  
  bool try_lock_cell_impl(int cell_index)
  {
    return m_grid[cell_index].try_lock();
  }
  
  void unlock_cell_impl(int cell_index)
  {
    m_grid[cell_index].unlock();
  }

protected:
  
  tbb::recursive_mutex *                          m_grid;
};


//typedef Simple_grid_locking_ds LockDataStructureType;
//typedef Simple_grid_locking_ds_with_mutex LockDataStructureType;
typedef Simple_grid_locking_ds_with_thread_ids LockDataStructureType;


} //namespace Mesh_3
} //namespace CGAL

#endif // CGAL_MESH_3_LOCKING_DATA_STRUCTURES_H
#endif // CONCURRENT_MESH_3
