// Copyright (c) 1999-2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion

// cell of a triangulation of any dimension <=3,
// storing its circumcenter lazily.

#ifndef CGAL_TRIANGULATION_CELL_BASE_WITH_CIRCUMCENTER_3_H
#define CGAL_TRIANGULATION_CELL_BASE_WITH_CIRCUMCENTER_3_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_ds_cell_base_3.h>

namespace CGAL {

template < typename GT, typename Cb = Triangulation_ds_cell_base_3<> >
class Triangulation_cell_base_with_circumcenter_3
  : public Cb
{
  typedef typename GT::Point_3                         Point_3;

  mutable Point_3 * circumcenter_;

  void invalidate_circumcenter()
  {
      if (circumcenter_) {
          delete circumcenter_;
          circumcenter_ = NULL;
      }
  }

public:
  typedef typename Cb::Vertex_handle                   Vertex_handle;
  typedef typename Cb::Cell_handle                     Cell_handle;

  typedef GT                                           Geom_traits;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other         Cb2;
    typedef Triangulation_cell_base_with_circumcenter_3<GT, Cb2>  Other;
  };

  Triangulation_cell_base_with_circumcenter_3()
    : Cb(), circumcenter_(NULL) {}

  Triangulation_cell_base_with_circumcenter_3
        (const Triangulation_cell_base_with_circumcenter_3 &c)
    : Cb(c), circumcenter_(c.circumcenter_ != NULL ? new Point_3(*(c.circumcenter_)) : NULL)
  {}

  Triangulation_cell_base_with_circumcenter_3&
  operator=(const Triangulation_cell_base_with_circumcenter_3 &c)
  {
      Triangulation_cell_base_with_circumcenter_3 tmp=c;
      std::swap(tmp, *this);
      return *this;
  }

  Triangulation_cell_base_with_circumcenter_3(
	                    Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3), circumcenter_(NULL) {}

  Triangulation_cell_base_with_circumcenter_3(
	                    Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2, Vertex_handle v3,
                            Cell_handle   n0, Cell_handle   n1,
                            Cell_handle   n2, Cell_handle   n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3), circumcenter_(NULL) {}

  ~Triangulation_cell_base_with_circumcenter_3()
  {
      delete circumcenter_;
  }

  // We must override the functions that modify the vertices.
  // And if the point inside a vertex is modified, we fail,
  // but there's not much we can do for this now.
  void set_vertex(int i, Vertex_handle v)
  {
      invalidate_circumcenter();
      Cb::set_vertex(i, v);
  }

  void set_vertices()
  {
      invalidate_circumcenter();
      Cb::set_vertices();
  }

  void set_vertices(Vertex_handle v0, Vertex_handle v1,
                    Vertex_handle v2, Vertex_handle v3)
  {
      invalidate_circumcenter();
      Cb::set_vertices(v0, v1, v2, v3);
  }

  const Point_3 &
  circumcenter(const Geom_traits& gt = Geom_traits()) const
  {
      if (circumcenter_ == NULL) {
	  circumcenter_ = new Point_3(
	      gt.construct_circumcenter_3_object()(this->vertex(0)->point(),
	                                           this->vertex(1)->point(),
						   this->vertex(2)->point(),
						   this->vertex(3)->point()));
      }
      return *circumcenter_;
  }
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_CELL_BASE_WITH_CIRCUMCENTER_3_H
