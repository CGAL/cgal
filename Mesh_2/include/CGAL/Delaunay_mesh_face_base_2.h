// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DELAUNAY_FACE_BASE_2_H
#define CGAL_DELAUNAY_FACE_BASE_2_H

#include <CGAL/license/Mesh_2.h>


#include <CGAL/Constrained_Delaunay_triangulation_face_base_2.h>
#include <CGAL/Has_timestamp.h>

namespace CGAL {

template <class Gt,
          class Fb = Constrained_Delaunay_triangulation_face_base_2<Gt> >
class Delaunay_mesh_face_base_2 : public Fb
{
public:
  typedef Gt Geom_traits;
  typedef typename Fb::Vertex_handle Vertex_handle;
  typedef typename Fb::Face_handle Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other Fb2;
    typedef Delaunay_mesh_face_base_2<Gt,Fb2> Other;
  };

protected:
  bool in_domain;

public:
  Delaunay_mesh_face_base_2(): Fb(), in_domain(false) {}

  Delaunay_mesh_face_base_2(Vertex_handle v0,
                            Vertex_handle v1,
                            Vertex_handle v2)
    : Fb(v0,v1,v2), in_domain(false) {}

  Delaunay_mesh_face_base_2(Vertex_handle v0,
                            Vertex_handle v1,
                            Vertex_handle v2,
                            Face_handle n0,
                            Face_handle n1,
                            Face_handle n2)
    : Fb(v0,v1,v2,n0,n1,n2), in_domain(false) {}

  inline
  bool is_in_domain() const { return in_domain; }

  inline
  void set_in_domain(const bool b) { in_domain=b; }

  /** compatibility with CGAL-3.2 */
  inline
  bool is_marked() const { return in_domain; }

  /** compatibility with CGAL-3.2 */
  inline
  void set_marked(const bool b) { in_domain=b; }

  typedef Tag_true Has_timestamp;

  std::size_t time_stamp() const { return time_stamp_; }

  void set_time_stamp(const std::size_t& ts) { time_stamp_ = ts; }

  std::size_t time_stamp_;
};

} // namespace CGAL

#endif
