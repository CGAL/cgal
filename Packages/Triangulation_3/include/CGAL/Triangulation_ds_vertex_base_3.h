// ============================================================================
//
// Copyright (c) 1999,2000,2001,2002,2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_ds_vertex_base_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_VERTEX_BASE_3_H
#define CGAL_TRIANGULATION_DS_VERTEX_BASE_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Dummy_tds_3.h>

CGAL_BEGIN_NAMESPACE

template < typename TDS = void >
class Triangulation_ds_vertex_base_3
{
public:
  typedef TDS                          Triangulation_data_structure;
  typedef typename TDS::Vertex_handle  Vertex_handle;
  typedef typename TDS::Cell_handle    Cell_handle;

  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_ds_vertex_base_3<TDS2> Other; };

  Triangulation_ds_vertex_base_3()
    : _c(NULL) {}
  
  Triangulation_ds_vertex_base_3(Cell_handle c)
    : _c(c) {}

  Cell_handle cell() const
  { return _c; }

  void set_cell(Cell_handle c)
  { _c = c; }

  // the following trivial is_valid allows
  // the user of derived cell base classes 
  // to add their own purpose checking
  bool is_valid(bool, int ) const
  { 
    return cell() != NULL;
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
  typedef Dummy_tds_3             Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Cell_handle     Cell_handle;
  template <typename TDS2>
  struct Rebind_TDS { typedef Triangulation_ds_vertex_base_3<TDS2> Other; };
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_DS_VERTEX_BASE_3_H
