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
// $URL: $
// $Id: $
//
// Author(s)     : Clement Jamin

#ifdef CONCURRENT_MESH_3

#ifndef CGAL_MESH_3_WORKSHARING_DATA_STRUCTURES_H
#define CGAL_MESH_3_WORKSHARING_DATA_STRUCTURES_H

#include <CGAL/Bbox_3.h>

#include <tbb/concurrent_queue.h>
#include <tbb/task.h>

// CJTODO TEMP: not thread-safe => move it to Mesher_3
extern CGAL::Bbox_3 g_bbox;

namespace CGAL {
namespace Mesh_3 {

// Forward declarations
class Dynamic_load_based_worksharing_ds;
// Typedef
typedef Dynamic_load_based_worksharing_ds Worksharing_ds_type;



class Work_statistics
{
public:
  // Constructors
  
  Work_statistics(const Bbox_3 &bbox, 
                  int num_grid_cells_per_axis)
    : m_num_grid_cells_per_axis(num_grid_cells_per_axis)
  {
    m_laziest_cell_index = 0;
    m_laziest_cell_occupation = 1000;

    int num_cells =
      num_grid_cells_per_axis*num_grid_cells_per_axis*num_grid_cells_per_axis;
    m_occupation_grid = new tbb::atomic<int>[num_cells];
    // Initialize grid
    for (int i = 0 ; i < num_cells ; ++i)
      m_occupation_grid[i] = 0;

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
  ~Work_statistics()
  {
    delete [] m_occupation_grid;
  }

  void add_occupation(int cell_index, int to_add, int num_items_in_work_queue)
  {
    int new_occupation = 
      (m_occupation_grid[cell_index].fetch_and_add(to_add)) 
      + to_add;

    // If this cell is the current most lazy, update the value
    if (cell_index == m_laziest_cell_index)
    {
      if (num_items_in_work_queue == 0)
        // So that it won't stay long the laziest
        m_laziest_cell_occupation = 999999;
      else
        m_laziest_cell_occupation = new_occupation;
    }
    else if (num_items_in_work_queue > 0 
      && new_occupation <= m_laziest_cell_occupation)
    {
      m_laziest_cell_index = cell_index;
      m_laziest_cell_occupation = new_occupation;
    }
  }
  
  void add_occupation(int index_x, int index_y, int index_z, 
                      int to_add, int num_items_in_work_queue)
  {
    int index = 
      index_z*m_num_grid_cells_per_axis*m_num_grid_cells_per_axis
      + index_y*m_num_grid_cells_per_axis 
      + index_x;
    return add_occupation(index, to_add, num_items_in_work_queue);
  }
  
  /// P3 must provide .x(), .y(), .z()
  template <typename P3>
  int compute_index(const P3 &point)
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

    return index;
  }

  /// P3 must provide .x(), .y(), .z()
  // Returns index in grid
  template <typename P3>
  int add_occupation(const P3 &point, int to_add, int num_items_in_work_queue)
  {
    int index = compute_index(point);
    add_occupation(index, to_add, num_items_in_work_queue);
    return index;
  }

  int get_laziest_cell_index()
  {
    return m_laziest_cell_index;
  }
  
protected:
  int                                             m_num_grid_cells_per_axis;
  double                                          m_xmin;
  double                                          m_ymin;
  double                                          m_zmin;
  double                                          m_resolution_x;
  double                                          m_resolution_y;
  double                                          m_resolution_z;
  tbb::atomic<int> *                              m_occupation_grid;

  tbb::atomic<int>                                m_laziest_cell_index;
  tbb::atomic<int>                                m_laziest_cell_occupation;
};


/* 
 * ==============
 * class WorkItem
 * Abstract base class for a piece of work.
 * ==============
 */
class WorkItem 
{
public:
  WorkItem() {}
  // Derived class defines the actual work.
  virtual void run() = 0;
  virtual void set_index(int) = 0;
  virtual int get_index() const = 0;
};

template<typename Func>
class ConcreteWorkItem
  : public WorkItem
{
public:
  ConcreteWorkItem(const Func& func)
    : m_func(func), m_index(-1)
  {}
  
  void run() 
  {
    m_func();
    delete this;
  }
  
  void set_index(int index)
  {
    m_index = index;
  }

  int get_index() const
  {
    return m_index;
  }

private:
  Func  m_func;
  int   m_index;
};



/* 
 * =================
 * class RunWorkItem
 * =================
 */
class RunWorkItem
  : public tbb::task 
{
public:
  RunWorkItem() {}

private:
  /*override*/inline tbb::task* execute();
};



/* 
 * =======================================
 * class Dynamic_load_based_worksharing_ds
 * =======================================
 */
class Dynamic_load_based_worksharing_ds
{
public:
  // Constructors
  Dynamic_load_based_worksharing_ds()
    : m_stats(g_bbox, MESH_3_WORK_STATS_GRID_NUM_CELLS_PER_AXIS)
  {
    for (int i = 0 ; i < MESH_3_WORK_STATS_GRID_NUM_CELLS ; ++i)
      m_num_items[i] = 0;
  }

  /// Destructor
  ~Dynamic_load_based_worksharing_ds()
  {
  }

  template <typename P3>
  void add(WorkItem * p_item, const P3 &point, tbb::task &parent_task)
  {
    int index = m_stats.compute_index(point);
    p_item->set_index(index);
    m_work_items[index].push(p_item);
    ++m_num_items[index];
    // CJTODO: try "spawn" instead of enqueue
    tbb::task::enqueue(*new(parent_task.allocate_child()) RunWorkItem);
  }

  void run_next_work_item()
  {
    WorkItem *p_item = 0;
    int index = m_stats.get_laziest_cell_index();
    bool popped = m_work_items[index].try_pop(p_item);
    // If queue is empty
    if (!popped)
    {
      // Look for an non-empty queue
      for (index = 0 ; !popped ; ++index)
      {
        CGAL_assertion(index < MESH_3_WORK_STATS_GRID_NUM_CELLS);
        popped = m_work_items[index].try_pop(p_item);
      }

      --index;
    }
    --m_num_items[index];
    CGAL_assertion(p_item != 0);
    m_stats.add_occupation(index, 1, m_num_items[index]);
    p_item->run();
    m_stats.add_occupation(index, -1, m_num_items[index]);
  }

protected:
  Work_statistics                   m_stats; 
  tbb::concurrent_queue<WorkItem*>  m_work_items[MESH_3_WORK_STATS_GRID_NUM_CELLS];
  tbb::atomic<int>                  m_num_items [MESH_3_WORK_STATS_GRID_NUM_CELLS];
};


} //namespace Mesh_3
} //namespace CGAL

extern CGAL::Mesh_3::Worksharing_ds_type g_worksharing_ds;

namespace CGAL
{
namespace Mesh_3
{

inline tbb::task* RunWorkItem::execute()
{
  g_worksharing_ds.run_next_work_item();
  return NULL;
}

/* 
 * =====================
 * function enqueue_work
 * =====================
 */
template<typename Func, typename P3>
void enqueue_work(Func f, tbb::task &parent_task, const P3 &point)
{
  g_worksharing_ds.add(new ConcreteWorkItem<Func>(f), 
                       point,
                       parent_task);
}

} //namespace Mesh_3
} //namespace CGAL

#endif // CGAL_MESH_3_WORKSHARING_DATA_STRUCTURES_H
#endif // CONCURRENT_MESH_3
