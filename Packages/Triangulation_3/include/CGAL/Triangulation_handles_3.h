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
// file          : include/CGAL/Triangulation_handles_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_HANDLES_3_H
#define CGAL_TRIANGULATION_HANDLES_3_H

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Pointer.h>

CGAL_BEGIN_NAMESPACE

template <class Tr>
class Triangulation_vertex_handle_3
  : public Pointer<CGAL_TYPENAME_MSVC_NULL Tr::Vertex>
{
  typedef typename Tr::Vertex                        Vertex;
  typedef typename Tr::Vertex_iterator               Vertex_iterator;
  typedef Pointer<Vertex>                            Ptr;
  typedef Triangulation_vertex_handle_3<Tr>          Vertex_handle;
public:

  Triangulation_vertex_handle_3()
    : Ptr(NULL)
  {}

  Triangulation_vertex_handle_3(Vertex * v)
    : Ptr(v)
  {}

  Triangulation_vertex_handle_3(const Vertex_iterator & vit)
    : Ptr(&(*vit))
  {}

  Vertex_handle & operator=(Vertex * v)
  {
    ptr() = v;
    return *this;
  }

  Vertex_handle & operator=(const Vertex_handle & v)
  {
    ptr() = v.ptr();
    return *this;
  }
};

template <class Tr>
class Triangulation_cell_handle_3
  : public Pointer<CGAL_TYPENAME_MSVC_NULL Tr::Cell>
{
  typedef typename Tr::Cell                         Cell;
  typedef typename Tr::Cell_iterator                Cell_iterator;
  typedef typename Tr::Cell_circulator              Cell_circulator;
  typedef Pointer<Cell>                             Ptr;
  typedef Triangulation_cell_handle_3<Tr>           Cell_handle;
public:

  Triangulation_cell_handle_3()
    : Ptr(NULL)
  {}

  Triangulation_cell_handle_3(Cell * c)
    : Ptr(c)
  {}

  Triangulation_cell_handle_3(const Cell_iterator & cit)
    : Ptr(&(*cit))
  {}

  Triangulation_cell_handle_3(const Cell_circulator & ccir)
    : Ptr(&(*ccir))
  {}

  Cell_handle & operator=(Cell * c)
  {
    ptr() = c;
    return *this;
  }

  Cell_handle & operator=(const Cell_handle & c)
  {
    ptr() = c.ptr();
    return *this;
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_HANDLES_3_H
