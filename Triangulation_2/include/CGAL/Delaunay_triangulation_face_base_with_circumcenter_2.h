// Copyright (c) 1999-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion

// face of a triangulation of any dimension <=2,
// storing its circumcenter lazily.

#ifndef CGAL_DELAUNAY_TRIANGULATION_FACE_BASE_WITH_CIRCUMCENTER_2_H
#define CGAL_DELAUNAY_TRIANGULATION_FACE_BASE_WITH_CIRCUMCENTER_2_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/Delaunay_triangulation_face_base_2.h>

namespace CGAL {

template < typename GT, typename Fb = Delaunay_triangulation_face_base_2<GT> >
class Delaunay_triangulation_face_base_with_circumcenter_2
  : public Fb
{
  typedef typename GT::Point_2                         Point;

  mutable Point * circumcenter_;

public:
  void invalidate_circumcenter()
  {
      if (circumcenter_) {
          delete circumcenter_;
          circumcenter_ = nullptr;
      }
  }

public:
  typedef typename Fb::Vertex_handle                   Vertex_handle;
  typedef typename Fb::Face_handle                     Face_handle;

  typedef GT                                           Geom_traits;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other         Fb2;
    typedef Delaunay_triangulation_face_base_with_circumcenter_2<GT, Fb2>  Other;
  };

  Delaunay_triangulation_face_base_with_circumcenter_2()
    : Fb(), circumcenter_(nullptr) {}

  Delaunay_triangulation_face_base_with_circumcenter_2
        (const Delaunay_triangulation_face_base_with_circumcenter_2 &c)
    : Fb(c), circumcenter_(c.circumcenter_ != nullptr ? new Point(*(c.circumcenter_)) : nullptr)
  {}

  Delaunay_triangulation_face_base_with_circumcenter_2
        (Delaunay_triangulation_face_base_with_circumcenter_2 &&c)
    : Fb(std::move(c)), circumcenter_(std::exchange(c.circumcenter_, nullptr))
  {
  }

  Delaunay_triangulation_face_base_with_circumcenter_2&
  operator=(const Delaunay_triangulation_face_base_with_circumcenter_2 &c)
  {
      Delaunay_triangulation_face_base_with_circumcenter_2 tmp=c;
      std::swap(tmp, *this);
      return *this;
  }

  Delaunay_triangulation_face_base_with_circumcenter_2&
  operator=(Delaunay_triangulation_face_base_with_circumcenter_2 &&c)
  {
      Fb::operator=(std::move(c));
      circumcenter_ = std::exchange(c.circumcenter_, nullptr);
      return *this;
  }

  Delaunay_triangulation_face_base_with_circumcenter_2(
                            Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2)
    : Fb(v0, v1, v2), circumcenter_(nullptr) {}

  Delaunay_triangulation_face_base_with_circumcenter_2(
                            Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2,
                            Face_handle   n0, Face_handle   n1,
                            Face_handle   n2)
    : Fb(v0, v1, v2, n0, n1, n2), circumcenter_(nullptr) {}

  ~Delaunay_triangulation_face_base_with_circumcenter_2()
  {
      delete circumcenter_;
  }

  // We must override the functions that modify the vertices.
  // And if the point inside a vertex is modified, we fail,
  // but there's not much we can do for this now.
  void set_vertex(int i, Vertex_handle v)
  {
      invalidate_circumcenter();
      Fb::set_vertex(i, v);
  }

  void set_vertices()
  {
      invalidate_circumcenter();
      Fb::set_vertices();
  }

  void set_vertices(Vertex_handle v0, Vertex_handle v1,
                    Vertex_handle v2)
  {
      invalidate_circumcenter();
      Fb::set_vertices(v0, v1, v2);
  }

  void set_circumcenter(const Point& p) const
  {
      if (circumcenter_ == nullptr) {
        circumcenter_ = new Point(p);
      } else {
        *circumcenter_ = p;
      }
  }

  const Point& circumcenter(const Geom_traits& gt = Geom_traits()) const
  {
      if (circumcenter_ == nullptr) {
        circumcenter_ = new Point(this->Fb::circumcenter(gt));
      } else {
        CGAL_expensive_assertion(
          this->Fb::circumcenter(gt) == *circumcenter_);
      }

      return *circumcenter_;
  }
};

} //namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_FACE_BASE_WITH_CIRCUMCENTER_2_H
