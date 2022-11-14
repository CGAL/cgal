// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Clement Jamin

#ifndef CGAL_MESH_3_WORKSHARING_DATA_STRUCTURES_H
#define CGAL_MESH_3_WORKSHARING_DATA_STRUCTURES_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#ifdef CGAL_LINKED_WITH_TBB

#include <CGAL/Mesh_3/Concurrent_mesher_config.h>

#include <CGAL/Bbox_3.h>

#include <tbb/concurrent_queue.h>
#include <tbb/task_group.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/concurrent_vector.h>
#include <tbb/scalable_allocator.h>

#include <atomic>

#include <vector>

namespace CGAL { namespace Mesh_3 {

// Forward declarations
class Load_based_worksharing_ds;
class Auto_worksharing_ds;

// Typedef

// Load-based
#ifdef CGAL_MESH_3_LOAD_BASED_WORKSHARING
  typedef Load_based_worksharing_ds WorksharingDataStructureType;
// Task-scheduler with TLS work buffers
// => 1 work-buffer / thread
#else
  typedef Auto_worksharing_ds WorksharingDataStructureType;
#endif



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

    m_num_cells =
      num_grid_cells_per_axis*num_grid_cells_per_axis*num_grid_cells_per_axis;
    m_occupation_grid = new std::atomic<int>[m_num_cells];
    m_num_batches_grid = new std::atomic<int>[m_num_cells];
    // Initialize grid
    for (int i = 0 ; i < m_num_cells ; ++i)
    {
      m_occupation_grid[i] = 0;
      m_num_batches_grid[i] = 0;
    }

    set_bbox(bbox);
  }

  /// Destructor
  virtual ~Work_statistics()
  {
    delete [] m_occupation_grid;
    delete [] m_num_batches_grid;
  }

  void set_bbox(const Bbox_3 &bbox)
  {
    // Keep mins and resolutions
    m_xmin = bbox.xmin();
    m_ymin = bbox.ymin();
    m_zmin = bbox.zmin();
    double n = static_cast<double>(m_num_grid_cells_per_axis);
    m_resolution_x = n / (bbox.xmax() - m_xmin);
    m_resolution_y = n / (bbox.ymax() - m_ymin);
    m_resolution_z = n / (bbox.zmax() - m_zmin);

#ifdef CGAL_CONCURRENT_MESH_3_VERBOSE
    std::cerr << "Worksharing data structure Bounding Box = ["
      << bbox.xmin() << ", " << bbox.xmax() << "], "
      << bbox.ymin() << ", " << bbox.ymax() << "], "
      << bbox.zmin() << ", " << bbox.zmax() << "]"
      << std::endl;
#endif
  }

  void add_batch(int cell_index, int to_add)
  {
    m_num_batches_grid[cell_index].fetch_add(to_add);
  }

  void add_occupation(int cell_index, int to_add, int)
  {
    m_occupation_grid[cell_index].fetch_add(to_add);

    /*int new_occupation =
      (m_occupation_grid[cell_index].fetch_add(to_add))
      + to_add;
    //m_num_batches_grid[cell_index] = num_items_in_work_queue;

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
    }*/
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
  int compute_index(const P3 &point) const
  {
    // Compute indices on grid
    int index_x = static_cast<int>( (to_double(point.x()) - m_xmin) * m_resolution_x);
    index_x = (std::max)( 0, (std::min)(index_x, m_num_grid_cells_per_axis - 1) );
    int index_y = static_cast<int>( (to_double(point.y()) - m_ymin) * m_resolution_y);
    index_y = (std::max)( 0, (std::min)(index_y, m_num_grid_cells_per_axis - 1) );
    int index_z = static_cast<int>( (to_double(point.z()) - m_zmin) * m_resolution_z);
    index_z = (std::max)( 0, (std::min)(index_z, m_num_grid_cells_per_axis - 1) );

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

  int get_laziest_cell_index() const
  {
    //return m_laziest_cell_index;


    /*
    // Look for best occupation/work ratio
    int laziest_index = 0;
    float laziest_ratio = 200000.f;
    for (int i = 0 ; i < m_num_cells ; ++i)
    {
      if (m_num_batches_grid[i] > 0)
      {
        float ratio =
          static_cast<float>(m_occupation_grid[i])
          / m_num_batches_grid[i];
        if (ratio < laziest_ratio)
        {
          laziest_index = i;
          laziest_ratio = ratio;
        }
      }
    }
    return laziest_index;*/

    // Look for the least occupied
    /*int laziest_index = 0;
    int smallest_occupation = 99999;
    for (int i = 0 ; i < m_num_cells ; ++i)
    {
      if (m_num_batches_grid[i] > 1)
      {
        if (m_occupation_grid[i] < smallest_occupation)
        {
          laziest_index = i;
          smallest_occupation = m_occupation_grid[i];
        }
      }
    }
    //std::cerr << "Occ=" << m_occupation_grid[laziest_index]
    //  << " / Bat=" << m_num_batches_grid[laziest_index]
    //  << std::endl;
    return laziest_index;*/


    // Rotate
    static std::atomic<int> last_cell_index;
    //std::cerr << "last=" << last_cell_index << std::endl;
    int i = (last_cell_index + 1) % m_num_cells;
    for ( ; i != last_cell_index ; i = (i + 1) % m_num_cells)
    {
      //std::cerr << "#" << i << "=" << m_num_batches_grid[i] << std::endl;
      if (m_num_batches_grid[i] > 0)
      {
        break;
      }
    }
    last_cell_index = i;
    return i;
  }

protected:
  double                                          m_xmin;
  double                                          m_ymin;
  double                                          m_zmin;
  double                                          m_resolution_x;
  double                                          m_resolution_y;
  double                                          m_resolution_z;

  int                                             m_num_grid_cells_per_axis;
  int                                             m_num_cells;
  std::atomic<int> *                              m_occupation_grid;
  std::atomic<int> *                              m_num_batches_grid;

  std::atomic<int>                                m_laziest_cell_index;
  std::atomic<int>                                m_laziest_cell_occupation;
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
  virtual ~WorkItem() { }

  // Derived class defines the actual work.
  virtual void operator()() const = 0;
  virtual bool less_than(const WorkItem &) const = 0;
};

struct CompareTwoWorkItems
{
  bool operator()(const WorkItem *p1, const WorkItem *p2) const
  {
    return p1->less_than(*p2);
  }
};

/*
 * ==============
 * class MeshRefinementWorkItem
 * Concrete class for a piece of work in the refinement process.
 * ==============
 */
template<typename Func, typename Quality>
class MeshRefinementWorkItem
  : public WorkItem
{
public:
  MeshRefinementWorkItem(const Func& func, const Quality &quality)
    : m_func(func), m_quality(quality)
  {}

  ~MeshRefinementWorkItem() override
  {}

  void operator()() const override
  {
    m_func();
    tbb::scalable_allocator<MeshRefinementWorkItem>().deallocate(
        const_cast<MeshRefinementWorkItem *>(this), 1);
  }

  bool less_than (const WorkItem &other) const override
  {
    /*try
    {
      const MeshRefinementWorkItem& other_cwi = dynamic_cast<const MeshRefinementWorkItem<Func,Quality>&>(other);
      return m_quality < other_cwi.m_quality;
    }
    catch (const std::bad_cast&)
    {
      return false;
    }*/
    const MeshRefinementWorkItem& other_cwi = static_cast<const MeshRefinementWorkItem<Func,Quality>&>(other);
    return m_quality < other_cwi.m_quality;
  }

private:
  Func      m_func;
  Quality   m_quality;
};




/*
 * ==============
 * class SimpleFunctorWorkItem
 * Concrete class for a work item embedding a simple functor
 * ==============
 */
template<typename Func>
class SimpleFunctorWorkItem
  : public WorkItem
{
public:
  SimpleFunctorWorkItem(const Func& func)
    : m_func(func)
  {}

  ~SimpleFunctorWorkItem() override = default;

  void operator()() const override
  {
    m_func();
    tbb::scalable_allocator<SimpleFunctorWorkItem>().deallocate(
        const_cast<SimpleFunctorWorkItem *>(this), 1);
  }

  // Irrelevant here
  bool less_than (const WorkItem &other) const override
  {
    // Just compare addresses
    return this < &other;
  }

private:
  Func      m_func;
};


/*
 * ===============
 * class WorkBatch
 * ===============
 */
class WorkBatch
{
public:
  typedef std::vector<WorkItem *> Batch;
  typedef Batch::const_iterator   BatchConstIterator;
  typedef Batch::iterator         BatchIterator;

  WorkBatch() {}

  void add_work_item(WorkItem *p_item)
  {
    m_batch.push_back(p_item);
  }

  void operator()() const
  {
    std::sort(m_batch.begin(), m_batch.end(), CompareTwoWorkItems());
    BatchIterator it = m_batch.begin();
    BatchIterator it_end = m_batch.end();
    for ( ; it != it_end ; ++it)
      (*it)->operator()();
  }

  size_t size() const
  {
    return m_batch.size();
  }

  void clear()
  {
    m_batch.clear();
  }

protected:
  mutable Batch m_batch;
};




/*
 * ===================
 * class WorkItemTask
 * ===================
 */
class WorkItemTask
{
public:
  WorkItemTask(WorkItem *pwi)
    : m_pwi(pwi)
  {
  }

private:
  inline void operator()() const;

  WorkItem *m_pwi;
};


/*
 * =======================================
 * class Simple_worksharing_ds
 * =======================================
 */
class Simple_worksharing_ds
{
  template <typename Func>
  void enqueue_work(Func f, tbb::task_group &task_group) const
  {
    task_group.run(f);
  }
};

/*
 * ==================
 * class TokenTask
 * ==================
 */
class TokenTask
{
public:
  TokenTask(Load_based_worksharing_ds *p_wsds)
    : m_worksharing_ds(p_wsds) {}

  inline void operator()() const;

private:
  Load_based_worksharing_ds *m_worksharing_ds;
};

/*
 * =======================================
 * class Load_based_worksharing_ds
 * =======================================
 */
class Load_based_worksharing_ds
{
public:
  // Constructors
  Load_based_worksharing_ds(const Bbox_3 &bbox)
    : NUM_WORK_ITEMS_PER_BATCH(
        Concurrent_mesher_config::get().num_work_items_per_batch),
      m_num_cells_per_axis(
        Concurrent_mesher_config::get().work_stats_grid_num_cells_per_axis),
      m_num_cells(m_num_cells_per_axis*m_num_cells_per_axis*m_num_cells_per_axis),
      m_stats(bbox, m_num_cells_per_axis)
  {
    m_tls_work_buffers = new TLS_WorkBuffer[m_num_cells];
    m_work_batches = new tbb::concurrent_queue<WorkBatch>[m_num_cells];
    m_num_batches = new std::atomic<int>[m_num_cells];

    for (int i = 0 ; i < m_num_cells ; ++i)
      m_num_batches[i] = 0;

    set_bbox(bbox);
  }

  /// Destructor
  virtual ~Load_based_worksharing_ds()
  {
    delete [] m_tls_work_buffers;
    delete [] m_work_batches;
    delete [] m_num_batches;
  }

  void set_bbox(const Bbox_3 &bbox)
  {
    m_stats.set_bbox(bbox);
  }

  template <typename P3, typename Func, typename Quality>
  void enqueue_work(Func f, const Quality &quality, tbb::task_group &task_group, const P3 &point)
  {
    WorkItem *p_item = new MeshRefinementWorkItem<Func, Quality>(f, quality);
    int index = m_stats.compute_index(point);
    WorkBatch &wb = m_tls_work_buffers[index].local();
    wb.add_work_item(p_item);
    if (wb.size() >= NUM_WORK_ITEMS_PER_BATCH)
    {
      add_batch_and_enqueue_task(wb, index, task_group);
      wb.clear();
    }
  }

  // Returns true if some items were flushed
  bool flush_work_buffers(tbb::task_group &task_group)
  {
    int num_flushed_items = 0;

    for (int i = 0 ; i < m_num_cells ; ++i)
    {
      for (TLS_WorkBuffer::iterator it_buffer = m_tls_work_buffers[i].begin() ;
           it_buffer != m_tls_work_buffers[i].end() ;
           ++it_buffer )
      {
        if (it_buffer->size() > 0)
        {
          add_batch(*it_buffer, i);
          it_buffer->clear();
          ++num_flushed_items;
        }
      }
    }

    for (int i = 0 ; i < num_flushed_items ; ++i)
      enqueue_task(task_group);

    return (num_flushed_items > 0);
  }

  void run_next_work_item()
  {
    WorkBuffer wb;
    int index = m_stats.get_laziest_cell_index();
    bool popped = m_work_batches[index].try_pop(wb);

    if (!popped)
    {
      // Look for an non-empty queue
      for (index = 0 ; !popped ; ++index)
      {
        CGAL_assertion(index < m_num_cells);
        popped = m_work_batches[index].try_pop(wb);
      }

      --index;
    }

    CGAL_assertion(index < m_num_cells);

    --m_num_batches[index];
    m_stats.add_batch(index, -1);
    add_occupation(index, 1);

#ifdef CGAL_CONCURRENT_MESH_3_VERY_VERBOSE
    std::cerr << "Running a batch of " << wb.size() <<
      " elements on cell #" << index << std::endl;
#endif
    wb();
    add_occupation(index, -1);
  }

protected:

  // TLS
  typedef WorkBatch                                   WorkBuffer;
  typedef tbb::enumerable_thread_specific<WorkBuffer> TLS_WorkBuffer;

  void add_batch(const WorkBuffer &wb, int index)
  {
    m_work_batches[index].push(wb);
    ++m_num_batches[index];
    m_stats.add_batch(index, 1);
  }

  void enqueue_task(tbb::task_group &task_group)
  {
    // Warning: when using "enqueue", the system will use up to two threads
    // even if you told task_scheduler_init to use only one
    // (see http://software.intel.com/en-us/forums/showthread.php?t=101669)
    task_group.run(TokenTask(this));
  }

  void add_batch_and_enqueue_task(const WorkBuffer &wb, int index,
                                  tbb::task_group &task_group)
  {
    add_batch(wb, index);
    enqueue_task(task_group);
  }

  void add_occupation(int cell_index, int to_add, int occupation_radius = 1)
  {
    int index_z = cell_index/(m_num_cells_per_axis*
                              m_num_cells_per_axis);
    cell_index -= index_z*
                  m_num_cells_per_axis*
                  m_num_cells_per_axis;
    int index_y = cell_index/m_num_cells_per_axis;
    cell_index -= index_y*m_num_cells_per_axis;
    int index_x = cell_index;

    // For each cell inside the square
    for (int i = (std::max)(0, index_x-occupation_radius) ;
          i <= (std::min)(m_num_cells_per_axis - 1, index_x+occupation_radius) ;
          ++i)
    {
      for (int j = (std::max)(0, index_y-occupation_radius) ;
            j <= (std::min)(m_num_cells_per_axis - 1, index_y+occupation_radius) ;
            ++j)
      {
        for (int k = (std::max)(0, index_z-occupation_radius) ;
              k <= (std::min)(m_num_cells_per_axis - 1, index_z+occupation_radius) ;
              ++k)
        {
          int index =
            k*m_num_cells_per_axis*m_num_cells_per_axis
            + j*m_num_cells_per_axis
            + i;

          int weight =
             (occupation_radius + 1 - std::abs(i - index_x))
            *(occupation_radius + 1 - std::abs(j - index_y))
            *(occupation_radius + 1 - std::abs(k - index_z));

          m_stats.add_occupation(index, to_add*weight, m_num_batches[index]);
        }
      }
    }

  }

  const size_t                      NUM_WORK_ITEMS_PER_BATCH;

  int                               m_num_cells_per_axis;
  int                               m_num_cells;
  Work_statistics                   m_stats;
  TLS_WorkBuffer                   *m_tls_work_buffers;
  tbb::concurrent_queue<WorkBatch> *m_work_batches;
  std::atomic<int>                 *m_num_batches;
};








/*
 * ===================
 * class WorkBatchTask
 * ===================
 */
class WorkBatchTask
{
public:
  WorkBatchTask(const WorkBatch &wb)
    : m_wb(wb)
  {
  }

  inline void operator()() const;

private:
  WorkBatch m_wb;
};

/*
 * =======================================
 * class Auto_worksharing_ds
 * =======================================
 */
class Auto_worksharing_ds
{
public:
  // Constructors
  Auto_worksharing_ds(const Bbox_3 &bbox)
    : NUM_WORK_ITEMS_PER_BATCH(
        Concurrent_mesher_config::get().num_work_items_per_batch)
  {
    set_bbox(bbox);
  }

  /// Destructor
  virtual ~Auto_worksharing_ds()
  {
  }

  void set_bbox(const Bbox_3 &/*bbox*/)
  {
    // We don't need it.
  }

  template <typename Func>
  void enqueue_work(Func f, tbb::task_group &task_group)
  {
    //WorkItem *p_item = new SimpleFunctorWorkItem<Func>(f);
    WorkItem *p_item =
      tbb::scalable_allocator<SimpleFunctorWorkItem<Func> >().allocate(1);
    new (p_item) SimpleFunctorWorkItem<Func>(f);
    WorkBatch &workbuffer = m_tls_work_buffers.local();
    workbuffer.add_work_item(p_item);
    if (workbuffer.size() >= NUM_WORK_ITEMS_PER_BATCH)
    {
      add_batch_and_enqueue_task(workbuffer, task_group);
      workbuffer.clear();
    }
  }

  template <typename Func, typename Quality>
  void enqueue_work(Func f, const Quality &quality, tbb::task_group &task_group)
  {
    WorkItem *p_item =
      tbb::scalable_allocator<MeshRefinementWorkItem<Func, Quality> >()
      .allocate(1);
    new (p_item) MeshRefinementWorkItem<Func, Quality>(f, quality);
    WorkBatch &workbuffer = m_tls_work_buffers.local();
    workbuffer.add_work_item(p_item);
    if (workbuffer.size() >= NUM_WORK_ITEMS_PER_BATCH)
    {
      add_batch_and_enqueue_task(workbuffer, task_group);
      workbuffer.clear();
    }
  }

  // Returns true if some items were flushed
  bool flush_work_buffers(tbb::task_group &task_group)
  {
    int num_flushed_items = 0;

    std::vector<WorkBatchTask> tasks;
    tasks.reserve(m_tls_work_buffers.size());

    for (TLS_WorkBuffer::iterator it_buffer = m_tls_work_buffers.begin() ;
          it_buffer != m_tls_work_buffers.end() ;
          ++it_buffer )
    {
      if (it_buffer->size() > 0)
      {
        tasks.push_back(create_task(*it_buffer));
        it_buffer->clear();
        ++num_flushed_items;
      }
    }

    for (auto it = tasks.begin() ;
      it != tasks.end() ; ++it)
    {
      enqueue_task(*it, task_group);
    }

    return (num_flushed_items > 0);
  }

  [[deprecated]] int approximate_number_of_enqueued_element() const {
    return 0;
  }

protected:

  // TLS
  typedef WorkBatch                                        WorkBuffer;
  typedef tbb::enumerable_thread_specific<WorkBuffer>      TLS_WorkBuffer;

  WorkBatchTask create_task(const WorkBuffer &wb) const
  {
    return { wb };
  }

  void enqueue_task(const WorkBatchTask& task,
                    tbb::task_group &task_group) const
  {
    task_group.run(task);
  }

  void add_batch_and_enqueue_task(const WorkBuffer &wb,
                                  tbb::task_group &task_group) const
  {
    enqueue_task(create_task(wb), task_group);
  }

  const size_t                      NUM_WORK_ITEMS_PER_BATCH;
  TLS_WorkBuffer                    m_tls_work_buffers;
};




inline void TokenTask::operator()() const
{
  m_worksharing_ds->run_next_work_item();
}

inline void WorkItemTask::operator()() const
{
  m_pwi->operator()();
}

inline void WorkBatchTask::operator()() const
{
  m_wb.operator()();
}

} } //namespace CGAL::Mesh_3

#else // !CGAL_LINKED_WITH_TBB

namespace CGAL { namespace Mesh_3 {
  typedef void WorksharingDataStructureType;
} } //namespace CGAL::Mesh_3

#endif // CGAL_LINKED_WITH_TBB

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_WORKSHARING_DATA_STRUCTURES_H
