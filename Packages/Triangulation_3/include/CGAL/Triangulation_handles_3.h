// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_HANDLES_3_H
#define CGAL_TRIANGULATION_HANDLES_3_H

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
class Triangulation_cell_3;

template <  class Gt, class Tds >
class Triangulation_vertex_3;

template <  class Gt, class Tds>
class Triangulation_cell_iterator_3;

template <  class Gt, class Tds>
class Triangulation_vertex_iterator_3;

template <  class Gt, class Tds>
class Triangulation_cell_circulator_3;

//template <  class Gt, class Tds>
//class Triangulation_vertex_circulator_3;


template < class Gt, class Tds>
class Triangulation_vertex_handle_3
  :public Pointer<Triangulation_vertex_3<Gt,Tds> > 
{
public:
  typedef Pointer<Triangulation_vertex_3<Gt,Tds> > Pointer;
  typedef Triangulation_vertex_3<Gt,Tds> Vertex;
  typedef Triangulation_vertex_handle_3<Gt,Tds> Vertex_handle;
  
  typedef Triangulation_vertex_iterator_3<Gt,Tds>      Vertex_iterator;
  
  inline 
  Triangulation_vertex_handle_3()
    : Pointer(NULL)
  {}

  inline  
  Triangulation_vertex_handle_3(const Vertex* v)
    : Pointer((Vertex*)v)
    {}

  inline  
  Triangulation_vertex_handle_3(const Vertex_iterator & vit)
    : Pointer(&(*vit))
    {}
  
  inline Vertex_handle & operator=(const Vertex* & v)
  {
    ptr() = v ;
    return *this;
  }
    
  inline Vertex_handle & operator=(const Vertex_handle & v)
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
  :public Pointer<Triangulation_cell_3<Gt,Tds> > 
{
public:
  typedef Pointer<Triangulation_cell_3<Gt,Tds> > Pointer;
  typedef Triangulation_cell_3<Gt,Tds> Cell;
  typedef Triangulation_cell_handle_3<Gt,Tds> Cell_handle;
  
  typedef Triangulation_cell_iterator_3<Gt,Tds> Cell_iterator;
  typedef Triangulation_cell_circulator_3<Gt,Tds> Cell_circulator;
  
  inline 
  Triangulation_cell_handle_3()
    : Pointer(NULL)
  {}

  inline  
  Triangulation_cell_handle_3(const Cell* c)
    : Pointer((Cell*)c)
  {}

  inline  
  Triangulation_cell_handle_3(const Cell_iterator & cit)
    : Pointer(&(*cit))
  {}
  
  inline  
  Triangulation_cell_handle_3(const Cell_circulator & ccir)
    : Pointer(&(*ccir))
  {}

  inline Cell_handle & operator=(const Cell* & c)
  {
    ptr() = c ;
    return *this;
  }
    
  inline Cell_handle & operator=(const Cell_handle & c)
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
