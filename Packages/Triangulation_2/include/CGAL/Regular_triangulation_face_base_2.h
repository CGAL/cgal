// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Regular_triangulation_face_base_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Frederic Fichel, Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_FACE_BASE_2_H
#define CGAL_REGULAR_TRIANGULATION_FACE_BASE_2_H

#include <list>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>

CGAL_BEGIN_NAMESPACE


template <class Gt>
class Regular_triangulation_face_base_2
  :  public Triangulation_face_base_2<Gt>
{
public:
  typedef Gt Geom_traits;
  typedef Triangulation_face_base_2<Gt> Fbase;
  typedef Regular_triangulation_face_base_2<Gt> Regular_face_base;
  typedef typename Gt::Weighted_point   Weighted_point;
  typedef std::list<void*>              Vertex_list;
protected:
  Vertex_list vlist;

public:
 Regular_triangulation_face_base_2()
   : Fbase(),  vlist()
  {}

  Regular_triangulation_face_base_2(void* v0, void* v1, void* v2)
    : Fbase(v0,v1,v2), vlist()
  { }

  Regular_triangulation_face_base_2(void* v0, void* v1, void* v2,
				    void* n0, void* n1, void* n2)
    : Fbase(v0,v1,v2,n0,n1,n2),  vlist()
  { }

  ~Regular_triangulation_face_base_2()
  { 
    vlist.clear();
  }


  Vertex_list& vertex_list()
  {
    return vlist;
  }


};

CGAL_END_NAMESPACE 

#endif // CGAL_REGULAR_TRIANGULATION_FACE_BASE_2_H

