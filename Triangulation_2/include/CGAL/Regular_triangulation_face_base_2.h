// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Frederic Fichel, Mariette Yvinec

#ifndef CGAL_REGULAR_TRIANGULATION_FACE_BASE_2_H
#define CGAL_REGULAR_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/license/Triangulation_2.h>


#include <list>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_face_base_2.h>

namespace CGAL {


template <class Gt, class Fb = Triangulation_face_base_2<Gt> >
class Regular_triangulation_face_base_2
  :  public Fb
{
  typedef Fb                                            Fbase;
  typedef typename Fbase::Triangulation_data_structure  TDS;
public:
  typedef Gt                                   Geom_traits;
  typedef TDS                                  Triangulation_data_structure;
  typedef typename TDS::Vertex_handle          Vertex_handle;
  typedef typename TDS::Face_handle            Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other   Fb2;
    typedef Regular_triangulation_face_base_2<Gt,Fb2>             Other;
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

} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_FACE_BASE_2_H
