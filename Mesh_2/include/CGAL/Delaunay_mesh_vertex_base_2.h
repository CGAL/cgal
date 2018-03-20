// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Jane Tournois

#ifndef CGAL_DELAUNAY_VERTEX_BASE_2_H
#define CGAL_DELAUNAY_VERTEX_BASE_2_H

#include <CGAL/license/Mesh_2.h>


#include <CGAL/Triangulation_vertex_base_2.h>

namespace CGAL {

template <class Gt,
          class Vb = Triangulation_vertex_base_2<Gt> >
class Delaunay_mesh_vertex_base_2 : public Vb
{
public:
  typedef typename Gt::FT          FT;
  typedef typename Vb::Point       Point;
  typedef typename Vb::Face_handle Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
    typedef Delaunay_mesh_vertex_base_2<Gt,Vb2> Other;
  };

protected:
  FT sizing_info_;

public:
  Delaunay_mesh_vertex_base_2()
    : Vb()
    , sizing_info_(0.)
  {}

  Delaunay_mesh_vertex_base_2(Point p)
    : Vb(p)
    , sizing_info_(0.)
  {}

  Delaunay_mesh_vertex_base_2(Point p,
                              Face_handle f)
    : Vb(p, f)
    , sizing_info_(0.)
  {}

  void set_sizing_info(const FT& s)
  { 
    sizing_info_ = s;
  }
  const FT& sizing_info() const { return sizing_info_; }

};

} // namespace CGAL

#endif
