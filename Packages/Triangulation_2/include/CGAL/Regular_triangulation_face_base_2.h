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
// package       : Triangulation_2
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


template <class Fb = Triangulation_ds_face_base_2 <> >
class Regular_triangulation_face_base_2
  :  public Fb
{
  typedef Fb                                            Fbase;
  typedef typename Fbase::Triangulation_data_structure  TDS;
public:
  typedef TDS                                  Triangulation_data_structure;
  typedef typename TDS::Vertex_handle          Vertex_handle;
  typedef typename TDS::Face_handle            Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other   Fb2;
    typedef Regular_triangulation_face_base_2<Fb2>             Other;
  };

  typedef std::list<Vertex_handle>             Vertex_list;

protected:
  Vertex_list vlist;

public:
 Regular_triangulation_face_base_2()
   : Fbase(),  vlist()
  {}

  Regular_triangulation_face_base_2(Vertex_handle v0, 
				    Vertex_handle v1, 
				    Vertex_handle v2)
    : Fbase(v0,v1,v2), vlist()
  { }

  Regular_triangulation_face_base_2(Vertex_handle v0, 
				    Vertex_handle v1, 
				    Vertex_handle v2,
				    Face_handle n0, 
				    Face_handle n1, 
				    Face_handle n2)
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

