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
#include <CGAL/Triangulation_vertex_3.h>
#include <CGAL/Triangulation_cell_3.h>
#include <CGAL/Triangulation_iterators_3.h>
#include <CGAL/Triangulation_circulators_3.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds > class Triangulation_cell_3;
template < class Gt, class Tds > class Triangulation_vertex_3;
template < class Gt, class Tds > class Triangulation_cell_iterator_3;
template < class Gt, class Tds > class Triangulation_vertex_iterator_3;
template < class Gt, class Tds > class Triangulation_cell_circulator_3;

template < class Gt, class Tds>
class Triangulation_vertex_handle_3
  : public Pointer<Triangulation_vertex_3<Gt,Tds> >
{
public:
  typedef Pointer<Triangulation_vertex_3<Gt,Tds> >   Ptr;
  typedef Triangulation_vertex_3<Gt,Tds>             Vertex;
  typedef Triangulation_vertex_handle_3<Gt,Tds>      Vertex_handle;

  typedef Triangulation_vertex_iterator_3<Gt,Tds>    Vertex_iterator;

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

template <class Gt, class Tds>
Triangulation_vertex_3<Gt,Tds> * 
handle2pointer(const Triangulation_vertex_handle_3<Gt,Tds> v)
{
  return v.ptr();
}

template <class Gt, class Tds>
class Triangulation_cell_handle_3
  : public Pointer<Triangulation_cell_3<Gt,Tds> >
{
public:
  typedef Pointer<Triangulation_cell_3<Gt,Tds> >    Ptr;
  typedef Triangulation_cell_3<Gt,Tds>              Cell;
  typedef Triangulation_cell_handle_3<Gt,Tds>       Cell_handle;

  typedef Triangulation_cell_iterator_3<Gt,Tds>     Cell_iterator;
  typedef Triangulation_cell_circulator_3<Gt,Tds>   Cell_circulator;

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

template <class Gt, class Tds>
Triangulation_cell_3<Gt,Tds> *
handle2pointer(const Triangulation_cell_handle_3<Gt,Tds> c)
{
  return c.ptr();
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_HANDLES_3_H
