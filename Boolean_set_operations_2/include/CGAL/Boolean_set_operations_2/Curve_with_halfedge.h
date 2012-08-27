// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_CURVE_WITH_HALFEDGE_H
#define CGAL_CURVE_WITH_HALFEDGE_H

namespace CGAL {

template <class Arrangement_>
class Curve_with_halfedge
{
protected:
  typedef typename Arrangement_::Halfedge_handle        Halfedge_handle;
  typedef typename Arrangement_::Halfedge_const_handle  Halfedge_const_handle;
  
public:
  Halfedge_handle  m_he;

  Curve_with_halfedge()
  {};

  Curve_with_halfedge(Halfedge_handle he) : m_he(he)
  {}

  Halfedge_handle halfedge() const
  {
    return (m_he);
  }

  Halfedge_handle halfedge()
  {
    return (m_he);
  }

  void set_halfedge(Halfedge_handle he)
  {
    m_he = he;
  }
};

} //namespace CGAL
#endif
