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
class Triangulation_all_faces_iterator_2;

template <  class Gt, class Tds>
class Triangulation_all_vertices_iterator_2;

template <  class Gt, class Tds>
class Triangulation_finite_faces_iterator_2;

template <  class Gt, class Tds>
class Triangulation_finite_vertices_iterator_2;

template <  class Gt, class Tds>
class Triangulation_face_circulator_2;

template <  class Gt, class Tds>
class Triangulation_vertex_circulator_2;


template <  class Gt, class Tds>
class Triangulation_face_handle_2
  :public Pointer<Triangulation_face_2<Gt,Tds> > 
{
public:
  typedef Pointer<Triangulation_face_2<Gt,Tds> >        Pointer;
  typedef Triangulation_face_2<Gt,Tds>                  Face;
  typedef Triangulation_face_handle_2<Gt,Tds>           Face_handle;
  
  typedef Triangulation_all_faces_iterator_2<Gt,Tds>    All_faces_iterator;
  typedef Triangulation_finite_faces_iterator_2<Gt,Tds> Finite_faces_iterator;
  typedef Triangulation_face_circulator_2<Gt,Tds>       Face_circulator;
  
  Triangulation_face_handle_2()
    : Pointer(NULL)
    {}

  Triangulation_face_handle_2(const Face* p)
    : Pointer(p)
    {}

  Triangulation_face_handle_2(const Face_handle& p)
    : Pointer(p.ptr())
    {}

  Triangulation_face_handle_2(const Finite_faces_iterator& fit)
    : Pointer(&(*fit))
    {}
  
  Triangulation_face_handle_2(const All_faces_iterator& fit)
    : Pointer(&(*fit))
    {}

  Triangulation_face_handle_2(const Face_circulator& fc)
    : Pointer(&(*fc))
    {}

  Face_handle& operator=(Face* p)
    {
      ptr() = p ;
      return *this;
    }
    
  Face_handle& operator=(const Face_handle& p)
    {
      ptr() = p.ptr();
      return *this;
    }
};

template < class Gt, class Tds>
class Triangulation_vertex_handle_2
  :public Pointer<Triangulation_vertex_2<Gt,Tds> > 
{
public:
  typedef Pointer<Triangulation_vertex_2<Gt,Tds> >  Pointer;
  typedef Triangulation_vertex_2<Gt,Tds>            Vertex;
  typedef Triangulation_vertex_handle_2<Gt,Tds>     Vertex_handle;
  
  typedef Triangulation_all_vertices_iterator_2<Gt,Tds> 
                                                    All_vertices_iterator;
  typedef Triangulation_finite_vertices_iterator_2<Gt,Tds> 
                                                    Finite_vertices_iterator;
  typedef Triangulation_vertex_circulator_2<Gt,Tds> Vertex_circulator;
  
  Triangulation_vertex_handle_2()
    : Pointer(NULL)
    {}

  Triangulation_vertex_handle_2(const Vertex* p)
    : Pointer(p)
    {}

  Triangulation_vertex_handle_2(const Vertex_handle& p)
    : Pointer(p.ptr())
    {}

  Triangulation_vertex_handle_2(const Finite_vertices_iterator& vit)
    : Pointer(&(*vit))
    {}
  
  Triangulation_vertex_handle_2(const All_vertices_iterator& vit)
    : Pointer(&(*vit))
    {}

  Triangulation_vertex_handle_2(const Vertex_circulator& vc)
    : Pointer(&(*vc))
    {}
  
  Vertex_handle& operator=(Vertex* p)
    {
        ptr() = p ;
        return *this;
    }
    
  Vertex_handle& operator=(const Vertex_handle& p)
    {
        ptr() = p.ptr();
        return *this;
    }
  
 };

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_HANDLES_2_H
