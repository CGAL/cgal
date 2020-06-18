// Copyright (c) 2019  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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

  int& id() { return face_id; }
  int id() const { return face_id; }

  int& edge_id(const std::size_t i) { return edge_ids[i]; }
  int edge_id(const std::size_t i) const { return edge_ids[i]; }

private:
  int face_id;
  std::array<int, 3> edge_ids;
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_FACE_BASE_WITH_ID_2_H
