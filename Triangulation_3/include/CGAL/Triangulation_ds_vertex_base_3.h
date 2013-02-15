// Copyright (c) 1999,2000,2001,2002,2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

#ifndef CGAL_TRIANGULATION_DS_VERTEX_BASE_3_H
#define CGAL_TRIANGULATION_DS_VERTEX_BASE_3_H

#include <CGAL/basic.h>
#include <CGAL/internal/Dummy_tds_3.h>

namespace CGAL {


// Without erase counter
template <bool Use_erase_counter>
class Triangulation_ds_vertex_base_3_base
{
public:
  // Dummy
  unsigned int get_erase_counter() const { return 0; }
  void set_erase_counter(unsigned int) {}
  void increment_erase_counter() {}
};

#ifdef CGAL_LINKED_WITH_TBB
// Specialized version (with erase counter)
template <>
class Triangulation_ds_vertex_base_3_base<true>
{
public:
  
  // Erase counter (cf. Compact_container)
  unsigned int get_erase_counter() const
  {
    return this->m_erase_counter;
  }
  void set_erase_counter(unsigned int c)
  {
	  this->m_erase_counter = c;
  }
  void increment_erase_counter()
  {
    ++this->m_erase_counter;
  }
  
protected:
  typedef tbb::atomic<unsigned int> Erase_counter_type;
  Erase_counter_type                m_erase_counter;

};
#endif // CGAL_LINKED_WITH_TBB



template < typename TDS = void >
class Triangulation_ds_vertex_base_3
: public Triangulation_ds_vertex_base_3_base<
    TDS::Vertex_container_strategy::Uses_erase_counter>
{
public:
  typedef TDS                          Triangulation_data_structure;
  typedef typename TDS::Vertex_handle  Vertex_handle;
  typedef typename TDS::Cell_handle    Cell_handle;

  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_ds_vertex_base_3<TDS2> Other; };

  
  Triangulation_ds_vertex_base_3()
    : _c()
  {
  }

  Triangulation_ds_vertex_base_3(Cell_handle c)
    : _c(c) 
  {
  }

  Cell_handle cell() const 
  { return _c; }  

  void set_cell(Cell_handle c)
  {
    _c = c;
  }

  // the following trivial is_valid allows
  // the user of derived cell base classes
  // to add their own purpose checking
  bool is_valid(bool = false, int = 0) const
  {
    return cell() != Cell_handle();
  }

  // For use by the Compact_container.
  void *   for_compact_container() const
  { return _c.for_compact_container(); }
  void * & for_compact_container()
  { return _c.for_compact_container(); }

private:
  Cell_handle _c;
};

template < class TDS >
inline
std::istream&
operator>>(std::istream &is, Triangulation_ds_vertex_base_3<TDS> &)
  // no combinatorial information.
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os, const Triangulation_ds_vertex_base_3<TDS> &)
  // no combinatorial information.
{
  return os;
}

// Specialization for void.
template <>
class Triangulation_ds_vertex_base_3<void>
{
public:
  typedef internal::Dummy_tds_3                         Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Cell_handle     Cell_handle;
  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_ds_vertex_base_3<TDS2> Other; };
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_DS_VERTEX_BASE_3_H
