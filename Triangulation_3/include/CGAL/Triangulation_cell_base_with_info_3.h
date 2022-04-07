// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

// cell of a triangulation of any dimension <=3

#ifndef CGAL_TRIANGULATION_CELL_BASE_WITH_INFO_3_H
#define CGAL_TRIANGULATION_CELL_BASE_WITH_INFO_3_H

#include <CGAL/license/Triangulation_3.h>


#include <CGAL/Triangulation_cell_base_3.h>

namespace CGAL {

template < typename Info_, typename GT,
           typename Cb = Triangulation_cell_base_3<GT> >
class Triangulation_cell_base_with_info_3
  : public Cb
{
  Info_ _info;
public:
  typedef typename Cb::Vertex_handle                   Vertex_handle;
  typedef typename Cb::Cell_handle                     Cell_handle;
  typedef Info_                                        Info;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other       Cb2;
    typedef Triangulation_cell_base_with_info_3<Info, GT, Cb2>  Other;
  };

  Triangulation_cell_base_with_info_3()
    : Cb() {}

  Triangulation_cell_base_with_info_3(Vertex_handle v0, Vertex_handle v1,
                                      Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {}

  Triangulation_cell_base_with_info_3(Vertex_handle v0, Vertex_handle v1,
                                      Vertex_handle v2, Vertex_handle v3,
                                      Cell_handle   n0, Cell_handle   n1,
                                      Cell_handle   n2, Cell_handle   n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3) {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_CELL_BASE_WITH_INFO_3_H
