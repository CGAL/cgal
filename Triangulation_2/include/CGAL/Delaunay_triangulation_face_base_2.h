// Copyright (c) 2025 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_DELAUNAY_TRIANGULATION_FACE_BASE_2_H
#define CGAL_DELAUNAY_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/config.h>
#include <CGAL/assertions.h>
#include <CGAL/Triangulation_face_base_2.h>

namespace CGAL {

template < typename Gt, typename Fb = Triangulation_face_base_2<Gt> >
class Delaunay_triangulation_face_base_2
  : public Fb
{
public:
  typedef Gt                                           Geom_traits;
  typedef typename Fb::Vertex_handle                   Vertex_handle;
  typedef typename Fb::Face_handle                     Face_handle;

  typedef typename Geom_traits::Point_2                Point_2;
  typedef typename Geom_traits::Point_2                Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Delaunay_triangulation_face_base_2<Gt, Fb2>    Other;
  };

public:
Delaunay_triangulation_face_base_2()
       : Fb() {}

       Delaunay_triangulation_face_base_2(Vertex_handle v0,
                                          Vertex_handle v1,
                                          Vertex_handle v2)
    : Fb(v0,v1,v2) {}

    Delaunay_triangulation_face_base_2(Vertex_handle v0,
                                       Vertex_handle v1,
                                       Vertex_handle v2,
                                       Face_handle n0,
                                       Face_handle n1,
                                       Face_handle n2)
    : Fb(v0,v1,v2,n0,n1,n2) {}

  static int ccw(int i) {return Triangulation_cw_ccw_2::ccw(i);}
  static int  cw(int i) {return Triangulation_cw_ccw_2::cw(i);}

  template <typename GT_>
  Point_2 circumcenter(const GT_& gt) const
  {
    return gt.construct_circumcenter_2_object()(this->vertex(0)->point(),
                                                this->vertex(1)->point(),
                                                this->vertex(2)->point());
  }

  Point_2 circumcenter() const
  {
    return circumcenter(Geom_traits());
  }
};

} // namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_FACE_BASE_2_H
