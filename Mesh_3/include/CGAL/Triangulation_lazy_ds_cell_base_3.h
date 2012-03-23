// Copyright (c) 1999-2005  INRIA Sophia-Antipolis (France).
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

// cell of a triangulation data structure of any dimension <=3 for lazy compact container

#ifndef CGAL_TRIANGULATION_LAZY_DS_CELL_BASE_3_H
#define CGAL_TRIANGULATION_LAZY_DS_CELL_BASE_3_H

#include <CGAL/Triangulation_ds_cell_base_3.h>

#ifdef CONCURRENT_MESH_3
  #include <tbb/tbb.h>

  #include <CGAL/Mesh_3/Locking_data_structures.h> // CJODO TEMP?
  #include <CGAL/BBox_3.h>
  // CJTODO TEMP: not thread-safe => move it to Mesher_3
  extern CGAL::Bbox_3 g_bbox;
# ifdef CGAL_MESH_3_LOCKING_STRATEGY_SIMPLE_GRID_LOCKING
  extern CGAL::Mesh_3::Refinement_grid_type g_lock_grid;
#elif defined(CGAL_MESH_3_LOCKING_STRATEGY_CELL_LOCK)
# include <utility>
# include <vector>
# include <tbb/enumerable_thread_specific.h>
  extern tbb::enumerable_thread_specific<std::vector<std::pair<void*, unsigned int> > > g_tls_locked_cells;
# endif

#endif

namespace CGAL {
  
template < typename TDS = void >
class Triangulation_lazy_ds_cell_base_3
  : public Triangulation_ds_cell_base_3<TDS>
{
public:
  typedef Triangulation_lazy_ds_cell_base_3<TDS> Self;

  typedef TDS                          Triangulation_data_structure;
  typedef typename TDS::Vertex_handle  Vertex_handle;
  typedef typename TDS::Cell_handle    Cell_handle;
  typedef typename TDS::Vertex         Vertex;
  typedef typename TDS::Cell           Cell;
  typedef typename TDS::Cell_data      TDS_data;

#ifdef CONCURRENT_MESH_3
  bool try_lock()
  {
    bool success = true;
    
# ifdef CGAL_MESH_3_LOCKING_STRATEGY_SIMPLE_GRID_LOCKING
    // Lock the element area on the grid
    for (int iVertex = 0 ; success && iVertex < 4 ; ++iVertex)
    {
      Vertex_handle vh = vertex(iVertex);
      //if (vh != infinite_vertex()) // CJTODO: à tester?
      std::pair<bool, int> r = g_lock_grid.try_lock(vh->point());
      success = r.first;
    }
# elif defined(CGAL_MESH_3_LOCKING_STRATEGY_CELL_LOCK)
    success = m_mutex.try_lock();
    if (success) 
      g_tls_locked_cells.local().push_back(std::make_pair(this, m_erase_counter));
# endif

    return success;
  }

# ifdef CGAL_MESH_3_LOCKING_STRATEGY_CELL_LOCK
  void lock()
  {
    m_mutex.lock();
    g_tls_locked_cells.local().push_back(std::make_pair(this, m_erase_counter));
  }

  void unlock()
  {
    m_mutex.unlock();
  }
#elif defined(CGAL_MESH_3_LOCKING_STRATEGY_SIMPLE_GRID_LOCKING)
  void lock()
  {
    // Active wait
    while (!try_lock())
      tbb::this_tbb_thread::yield(); 
  }

  void unlock()
  {
    for (int iVertex = 0 ; iVertex < 4 ; ++iVertex)
    {
      Vertex_handle vh = vertex(iVertex);
      g_lock_grid.unlock(vh->point());
    }
  }
# endif

  typedef tbb::atomic<unsigned int> Erase_counter_type;
#else
  typedef unsigned int Erase_counter_type;
#endif
  
  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_lazy_ds_cell_base_3<TDS2> Other; };

  // Constructors
  // We DO NOT INITIALIZE m_erase_counter since it is managed by the Compact_container

  Triangulation_lazy_ds_cell_base_3() {}

  Triangulation_lazy_ds_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2, Vertex_handle v3)
    : Triangulation_ds_cell_base_3(v0, v1, v2, v3) {}

  Triangulation_lazy_ds_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2, Vertex_handle v3,
                            Cell_handle   n0, Cell_handle   n1,
                            Cell_handle   n2, Cell_handle   n3)
    : Triangulation_ds_cell_base_3(v0, v1, v2, v3, n0, n1, n2, n3) {}
  
  // Erase counter (cf. Compact_container)

  unsigned int get_erase_counter() const
  {
    return m_erase_counter; 
  }
  void set_erase_counter(unsigned int c) 
  {
    m_erase_counter = c;
  }
  void increment_erase_counter()
  {
    ++m_erase_counter;
  }

protected:
  Erase_counter_type                m_erase_counter;
#ifdef CGAL_MESH_3_LOCKING_STRATEGY_CELL_LOCK
  mutable Cell_mutex_type           m_mutex;
#endif
};

// Specialization for void.
template <>
class Triangulation_lazy_ds_cell_base_3<void>
{
public:
  typedef internal::Dummy_tds_3                         Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Cell_handle     Cell_handle;
  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_lazy_ds_cell_base_3<TDS2> Other; };
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_LAZY_DS_CELL_BASE_3_H
