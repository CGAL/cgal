// Copyright (c) 2019  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_TRIANGULATION_FACE_BASE_WITH_ID_2_H
#define CGAL_TRIANGULATION_FACE_BASE_WITH_ID_2_H

#include <CGAL/Triangulation_face_base_2.h>

namespace CGAL {

template < typename GT,
           typename Fb = Triangulation_face_base_2<GT> >
class Triangulation_face_base_with_id_2
  : public Fb
{
public:
  typedef typename Fb::Vertex_handle                 Vertex_handle;
  typedef typename Fb::Face_handle                   Face_handle;
  typedef typename Fb::Point                         Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Triangulation_face_base_with_id_2<GT, Fb2>   Other;
  };

  Triangulation_face_base_with_id_2() : Fb() { }

  Triangulation_face_base_with_id_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2)
    : Fb(v0, v1, v2)
  { }

  Triangulation_face_base_with_id_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
                                    Face_handle n0, Face_handle n1, Face_handle n2)
    : Fb(v0, v1, v2, n0, n1, n2)
  { }

  int& id() { return _id; }
  int id() const { return _id; }

private:
  int _id;
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_FACE_BASE_WITH_ID_2_H
