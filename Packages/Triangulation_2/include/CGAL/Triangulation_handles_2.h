// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : Triangulation/include/CGAL/Triangulation_handles_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  < Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_HANDLES_2_H
#define CGAL_TRIANGULATION_HANDLES_2_H

#include <CGAL/Triangulation_short_names_2.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
class Triangulation_face_2;

template <  class Gt, class Tds >
class Triangulation_vertex_2;

template <  class Gt, class Tds>
class Triangulation_face_iterator_2;

template <  class Gt, class Tds>
class Triangulation_vertex_iterator_2;

template <  class Gt, class Tds>
class Triangulation_face_circulator_2;

template <  class Gt, class Tds>
class Triangulation_vertex_circulator_2;


template <  class Gt, class Tds>
class Triangulation_face_handle_2
  :public Pointer<Triangulation_face_2<Gt,Tds> > 
{
public:
  typedef Pointer<Triangulation_face_2<Gt,Tds> > Pointer;
  typedef Triangulation_face_2<Gt,Tds> Face;
  typedef Triangulation_face_handle_2<Gt,Tds> Face_handle;
  
  typedef Triangulation_face_iterator_2<Gt,Tds>      Face_iterator;
  typedef Triangulation_face_circulator_2<Gt,Tds>    Face_circulator;
  
  inline 
  Triangulation_face_handle_2()
    : Pointer(NULL)
  {}

  inline  
  Triangulation_face_handle_2(const Face* p)
    : Pointer((Face*)p)
  {}

  inline Face_handle& operator=(const Face*& p)
    {
        ptr() = p ;
        return *this;
    }
    
    inline Face_handle& operator=(const Face_handle& p)
    {
        ptr() = p.ptr();
        return *this;
    }
  
   inline  
    Triangulation_face_handle_2(const Face_iterator& fit)
        : Pointer(&(*fit))
    {}
  

  inline  
   Triangulation_face_handle_2(const Face_circulator& fc)
        : Pointer(&(*fc))
    {}
};

template < class Gt, class Tds>
class Triangulation_vertex_handle_2
  :public Pointer<Triangulation_vertex_2<Gt,Tds> > 
{
public:
  typedef Pointer<Triangulation_vertex_2<Gt,Tds> > Pointer;
  typedef Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  
  typedef Triangulation_vertex_iterator_2<Gt,Tds>      Vertex_iterator;
  typedef Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;
  
  inline 
  Triangulation_vertex_handle_2()
    : Pointer(NULL)
  {}

  inline  
  Triangulation_vertex_handle_2(const Vertex* p)
        : Pointer((Vertex*)p)
    {}

  inline Vertex_handle& operator=(const Vertex*& p)
    {
        ptr() = p ;
        return *this;
    }
    
    inline Vertex_handle& operator=(const Vertex_handle& p)
    {
        ptr() = p.ptr();
        return *this;
    }
  
   inline  
   Triangulation_vertex_handle_2(const Vertex_iterator& vit)
        : Pointer(&(*vit))
    {}

  
  inline  
   Triangulation_vertex_handle_2(const Vertex_circulator& vc)
        : Pointer(&(*vc))
    {}
  
};

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_HANDLES_2_H
