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

#ifndef CGAL_TRIANGULATION_VERTEX_BASE_WITH_INFO_3_H
#define CGAL_TRIANGULATION_VERTEX_BASE_WITH_INFO_3_H

#include <CGAL/license/Triangulation_3.h>


#include <CGAL/Triangulation_vertex_base_3.h>

namespace CGAL {

template < typename Info_, typename GT,
           typename Vb = Triangulation_vertex_base_3<GT> >
class Triangulation_vertex_base_with_info_3
  : public Vb
{
  Info_ _info;
public:

  typedef typename Vb::Cell_handle                   Cell_handle;
  typedef typename Vb::Point                         Point;
  typedef Info_                                      Info;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other          Vb2;
    typedef Triangulation_vertex_base_with_info_3<Info, GT, Vb2>   Other;
  };

  Triangulation_vertex_base_with_info_3()
    : Vb() {}

  Triangulation_vertex_base_with_info_3(const Point & p)
    : Vb(p) {}

  Triangulation_vertex_base_with_info_3(const Point & p, Cell_handle c)
    : Vb(p, c) {}

  Triangulation_vertex_base_with_info_3(Cell_handle c)
    : Vb(c) {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_VERTEX_BASE_WITH_INFO_3_H
