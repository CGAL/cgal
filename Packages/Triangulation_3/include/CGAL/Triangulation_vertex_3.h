// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_vertex_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_VERTEX_3_H
#define CGAL_TRIANGULATION_VERTEX_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template <class Tr>
class Triangulation_vertex_3
  : public Tr::Triangulation_data_structure::Vertex
{
  typedef typename Tr::Triangulation_data_structure::Vertex Vtds;
  typedef typename Tr::Triangulation_data_structure::Cell Ctds;

  typedef typename Tr::Vertex_handle    Vertex_handle;
  typedef typename Tr::Cell_handle      Cell_handle;
  typedef typename Tr::Cell             Cell;

public:
 
  Triangulation_vertex_3()
     : Vtds() {}

  void set_cell(Cell_handle c)
  {
    Vtds::set_cell(&(*c));
  }

  Cell_handle cell() const
  {
    return (Cell *) Vtds::cell();
  }

  Vertex_handle handle() const
  {
    return Vertex_handle(this);
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_VERTEX_3_H
