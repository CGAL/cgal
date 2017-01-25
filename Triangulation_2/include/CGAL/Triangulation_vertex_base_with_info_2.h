// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec, Sylvain Pion

#ifndef CGAL_TRIANGULATION_VERTEX_BASE_WITH_INFO_2_H
#define CGAL_TRIANGULATION_VERTEX_BASE_WITH_INFO_2_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/Triangulation_vertex_base_2.h>

namespace CGAL {

template < typename Info_, typename GT,
           typename Vb = Triangulation_vertex_base_2<GT> >
class Triangulation_vertex_base_with_info_2
  : public Vb
{
  Info_ _info;

public:
  typedef typename Vb::Face_handle                   Face_handle;
  typedef typename Vb::Point                         Point;
  typedef Info_                                      Info;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other          Vb2;
    typedef Triangulation_vertex_base_with_info_2<Info, GT, Vb2>   Other;
  };

  Triangulation_vertex_base_with_info_2()
    : Vb() {}

  Triangulation_vertex_base_with_info_2(const Point & p)
    : Vb(p) {}

  Triangulation_vertex_base_with_info_2(const Point & p, Face_handle c)
    : Vb(p, c) {}

  Triangulation_vertex_base_with_info_2(Face_handle c)
    : Vb(c) {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_VERTEX_BASE_WITH_INFO_2_H
