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
// file          : Triangulation/include/CGAL/Triangulation_ds_handles_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  < Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_HANDLES_2_H
#define CGAL_TRIANGULATION_DS_HANDLES_2_H

#include <CGAL/Pointer.h>

CGAL_BEGIN_NAMESPACE

// template <class Tds >
// class Triangulation_face_iterator_2;

// template <class Tds >
// class Triangulation_face_circulator_2;

// template <class Tds >
// class Triangulation_vertex_iterator_2;

template <class Gt, class Tds>
class Triangulation_line_face_circulator_2;

template <class Tds>
class Triangulation_ds_face_handle_2
  : public Pointer<typename Tds::Face> 
{
public:
  typedef Triangulation_ds_face_handle_2<Tds>          Face_handle;
  typedef Pointer<typename Tds::Face>                  Pointer_;
  typedef typename Tds::Face                           Face;
  typedef typename Tds::Face_iterator                  Face_iterator;
  typedef typename Tds::Face_circulator                Face_circulator;
  //typedef Triangulation_line_face_circulator_2<Gt,Tds> Line_face_circulator;

  Triangulation_ds_face_handle_2()
    : Pointer_(NULL)
    {}

  Triangulation_ds_face_handle_2( Face* p)
    : Pointer_(p)
    {}

  Triangulation_ds_face_handle_2( const Face_handle& p)
    : Pointer_(p.ptr())
    {}

  Triangulation_ds_face_handle_2(const Face_iterator& fit)
    : Pointer_(&(*fit))
    {}
  
  Triangulation_ds_face_handle_2(const Face_circulator& fc)
    : Pointer_(&(*fc))
    {}

  template<class Gt>
  Triangulation_ds_face_handle_2(
	      const Triangulation_line_face_circulator_2<Gt,Tds>& lfc)
    : Pointer_(lfc->handle())
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

template <class Tds>
class Triangulation_ds_vertex_handle_2
  : public Pointer<typename Tds::Vertex > 
{
public:
  typedef Triangulation_ds_vertex_handle_2<Tds>     Vertex_handle;
  typedef Pointer< typename Tds::Vertex>            Pointer_;
  typedef typename Tds::Vertex                      Vertex;
  typedef typename Tds::Vertex_iterator             Vertex_iterator;
  typedef typename Tds::Vertex_circulator           Vertex_circulator ;
  
  Triangulation_ds_vertex_handle_2()
    : Pointer_(NULL)
    {}

  Triangulation_ds_vertex_handle_2( Vertex* p)
    : Pointer_(p)
    {}

  Triangulation_ds_vertex_handle_2(const Vertex_handle& p)
    : Pointer_(p.ptr())
    {}

  Triangulation_ds_vertex_handle_2(const Vertex_iterator& vit)
    : Pointer_(&(*vit))
    {}
  
  Triangulation_ds_vertex_handle_2(const Vertex_circulator& vc)
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

#endif //CGAL_TRIANGULATION_DS_HANDLES_2_H
