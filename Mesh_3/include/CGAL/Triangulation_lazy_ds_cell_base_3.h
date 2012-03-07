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
  Erase_counter_type m_erase_counter;
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
