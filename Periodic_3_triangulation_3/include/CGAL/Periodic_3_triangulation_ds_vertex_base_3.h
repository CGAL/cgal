// Copyright (c) 1999-2003,2007-2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_TRIANGULATION_DS_VERTEX_BASE_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_DS_VERTEX_BASE_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/basic.h>
#include <CGAL/internal/Dummy_tds_3.h>
#include <CGAL/Periodic_3_offset_3.h>

namespace CGAL {

template < typename TDS = void >
class Periodic_3_triangulation_ds_vertex_base_3
{
public:
  typedef TDS                          Triangulation_data_structure;
  typedef typename TDS::Vertex_handle  Vertex_handle;
  typedef typename TDS::Cell_handle    Cell_handle;
  typedef CGAL::Periodic_3_offset_3    Offset;

  template <typename TDS2>
  struct Rebind_TDS {
    typedef Periodic_3_triangulation_ds_vertex_base_3<TDS2> Other;
  };

  Periodic_3_triangulation_ds_vertex_base_3()
    : _c(), _off(), offset_flag(false)
#ifdef CGAL_PERIODIC_TRIANGULATION_USE_VISITED_VERTEX_BOOLEAN
      , visited_for_vertex_extractor(false)
#endif
  {}

  Periodic_3_triangulation_ds_vertex_base_3(const Cell_handle& c)
    : _c(c), _off(), offset_flag(false)
#ifdef CGAL_PERIODIC_TRIANGULATION_USE_VISITED_VERTEX_BOOLEAN
      , visited_for_vertex_extractor(false)
#endif
  {}

  const Cell_handle& cell() const
  { return _c; }

  void set_cell(const Cell_handle& c)
  { _c = c; }

  const Offset& offset() const
  { return _off; }

  void set_offset(const Offset& off)
  { _off = off; offset_flag=true; }

  void clear_offset() {
    offset_flag=false;
    _off = Offset();
  }

  bool get_offset_flag() const { return offset_flag; }

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
  void for_compact_container(void *p)
  { _c.for_compact_container(p); }

private:
  Cell_handle _c;
  Offset _off;
  bool offset_flag;

public:
  // Using 'visited_for_vertex_extractor' below allows to avoid using sets while
  // gathering incident/adjacent elements, instead simply marking vertices when
  // they are visited.
  // IMPORTANT: this should only be used when sure that the triangulation is
  // _always_ a 1-cover periodic triangulation. Otherwise, the same vertex might
  // appear multiple times with different offsets but will be ignored because
  // it will have been marked as already visited, and bugs appear...
#ifdef CGAL_PERIODIC_TRIANGULATION_USE_VISITED_VERTEX_BOOLEAN
  // The typedef and the bool are used by Triangulation_data_structure::Vertex_extractor
  // The names are chosen complicated so that we do not have to document them
  // (privacy by obfuscation)
  typedef bool Has_visited_for_vertex_extractor;
  bool visited_for_vertex_extractor;
#endif
};

template < class TDS >
inline
std::istream&
operator>>(std::istream &is, Periodic_3_triangulation_ds_vertex_base_3<TDS> &)
  // no combinatorial information.
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os,
    const Periodic_3_triangulation_ds_vertex_base_3<TDS> &)
  // no combinatorial information.
{
  return os;
}

// Specialization for void.
template <>
class Periodic_3_triangulation_ds_vertex_base_3<void>
{
public:
  typedef internal::Dummy_tds_3             Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Cell_handle     Cell_handle;
  template <typename TDS2>
  struct Rebind_TDS {
    typedef Periodic_3_triangulation_ds_vertex_base_3<TDS2> Other;
  };
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_DS_VERTEX_BASE_3_H
