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
// source        : $Source$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  < Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_HANDLES_2_H
#define CGAL_TRIANGULATION_HANDLES_2_H

#include <CGAL/Triangulation_short_names_2.h>

template < class Gt, class Tds >
class CGAL_Triangulation_face_2;

template <  class Gt, class Tds >
class CGAL_Triangulation_vertex_2;

template <  class Gt, class Tds>
class CGAL_Triangulation_face_iterator_2;

template <  class Gt, class Tds>
class CGAL_Triangulation_vertex_iterator_2;

template <  class Gt, class Tds>
class CGAL_Triangulation_face_circulator_2;

template <  class Gt, class Tds>
class CGAL_Triangulation_vertex_circulator_2;


template <  class Gt, class Tds>
class CGAL_Triangulation_face_handle_2
  :public CGAL_Pointer<CGAL_Triangulation_face_2<Gt,Tds> > 
{
public:
  typedef CGAL_Pointer<CGAL_Triangulation_face_2<Gt,Tds> > Pointer;
  typedef CGAL_Triangulation_face_2<Gt,Tds> Face;
  typedef CGAL_Triangulation_face_handle_2<Gt,Tds> Face_handle;
  
  typedef CGAL_Triangulation_face_iterator_2<Gt,Tds>      Face_iterator;
  typedef CGAL_Triangulation_face_circulator_2<Gt,Tds>    Face_circulator;
  
  inline 
  CGAL_Triangulation_face_handle_2()
    : Pointer(NULL)
  {}

  inline  
  CGAL_Triangulation_face_handle_2(const Face* p)
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
    CGAL_Triangulation_face_handle_2(const Face_iterator& fit)
        : Pointer(&(*fit))
    {}
  

  inline  
   CGAL_Triangulation_face_handle_2(const Face_circulator& fc)
        : Pointer(&(*fc))
    {}
};

template < class Gt, class Tds>
class CGAL_Triangulation_vertex_handle_2
  :public CGAL_Pointer<CGAL_Triangulation_vertex_2<Gt,Tds> > 
{
public:
  typedef CGAL_Pointer<CGAL_Triangulation_vertex_2<Gt,Tds> > Pointer;
  typedef CGAL_Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef CGAL_Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  
  typedef CGAL_Triangulation_vertex_iterator_2<Gt,Tds>      Vertex_iterator;
  typedef CGAL_Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;
  
  inline 
  CGAL_Triangulation_vertex_handle_2()
    : Pointer(NULL)
  {}

  inline  
  CGAL_Triangulation_vertex_handle_2(const Vertex* p)
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
   CGAL_Triangulation_vertex_handle_2(const Vertex_iterator& vit)
        : Pointer(&(*vit))
    {}

  
  inline  
   CGAL_Triangulation_vertex_handle_2(const Vertex_circulator& vc)
        : Pointer(&(*vc))
    {}
  
};

#endif CGAL_TRIANGULATION_HANDLES_2_H
