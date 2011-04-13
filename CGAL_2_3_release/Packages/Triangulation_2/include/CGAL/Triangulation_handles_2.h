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

#include <CGAL/Pointer.h>

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
class Triangulation_line_face_circulator_2;


template <  class Gt, class Tds>
class Triangulation_face_handle_2
  :public Pointer<Triangulation_face_2<Gt,Tds> > 
{
public:
  typedef Pointer<Triangulation_face_2<Gt,Tds> >        Pointer_;
  typedef Triangulation_face_2<Gt,Tds>                  Face;
  typedef Triangulation_face_handle_2<Gt,Tds>           Face_handle;
  
  typedef Triangulation_all_faces_iterator_2<Gt,Tds>    All_faces_iterator;
  typedef Triangulation_finite_faces_iterator_2<Gt,Tds> Finite_faces_iterator;
  typedef Triangulation_face_circulator_2<Gt,Tds>       Face_circulator;
  typedef Triangulation_line_face_circulator_2<Gt,Tds> Line_face_circulator;

  Triangulation_face_handle_2()
    : Pointer_(NULL)
    {}

  Triangulation_face_handle_2( Face* p)
    : Pointer_(p)
    {}

  Triangulation_face_handle_2( const Face_handle& p)
    : Pointer_(p.ptr())
    {}

  Triangulation_face_handle_2(const Finite_faces_iterator& fit)
    : Pointer_(&(*fit))
    {}
  
  Triangulation_face_handle_2(const All_faces_iterator& fit)
    : Pointer_(&(*fit))
    {}

  Triangulation_face_handle_2(const Face_circulator& fc)
    : Pointer_(&(*fc))
    {}

  Triangulation_face_handle_2(const Line_face_circulator& fc)
    : Pointer_(&(*fc))
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
  typedef Pointer<Triangulation_vertex_2<Gt,Tds> >  Pointer_;
  typedef Triangulation_vertex_2<Gt,Tds>            Vertex;
  typedef Triangulation_vertex_handle_2<Gt,Tds>     Vertex_handle;
  
  typedef Triangulation_all_vertices_iterator_2<Gt,Tds> 
                                                    All_vertices_iterator;
  typedef Triangulation_finite_vertices_iterator_2<Gt,Tds> 
                                                    Finite_vertices_iterator;
  typedef Triangulation_vertex_circulator_2<Gt,Tds> Vertex_circulator;
  
  Triangulation_vertex_handle_2()
    : Pointer_(NULL)
    {}

  Triangulation_vertex_handle_2( Vertex* p)
    : Pointer_(p)
    {}

  Triangulation_vertex_handle_2(const Vertex_handle& p)
    : Pointer_(p.ptr())
    {}

  Triangulation_vertex_handle_2(const Finite_vertices_iterator& vit)
    : Pointer_(&(*vit))
    {}
  
  Triangulation_vertex_handle_2(const All_vertices_iterator& vit)
    : Pointer_(&(*vit))
    {}

  Triangulation_vertex_handle_2(const Vertex_circulator& vc)
    : Pointer_(&(*vc))
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
