// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_POINT_WITH_VERTEX_H
#define CGAL_POINT_WITH_VERTEX_H

#include <CGAL/license/Boolean_set_operations_2.h>


namespace CGAL {

template <class Arrangement_>
class Point_with_vertex
{
protected:
  typedef typename Arrangement_::Vertex_handle        Vertex_handle;
  typedef typename Arrangement_::Vertex_const_handle  Vertex_const_handle;

public:
  Vertex_handle  m_v;

  Point_with_vertex()
  {};

  Point_with_vertex(Vertex_handle v) : m_v(v)
  {}

  Vertex_handle vertex() const
  {
    return (m_v);
  }

  Vertex_handle vertex()
  {
    return (m_v);
  }

  void set_halfedge(Vertex_handle v)
  {
    m_v = v;
  }
};

} //namespace CGAL
#endif
