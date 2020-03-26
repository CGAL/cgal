// Copyright (c) 2013 INRIA Sophia-Antipolis (France),
//               2014-2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Jane Tournois, Raul Gallegos
//

#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_FACE_BASE_2_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/license/Mesh_2.h>


#include <CGAL/triangulation_assertions.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

namespace CGAL {

template <typename Kernel,
          class Fb = CGAL::Constrained_triangulation_face_base_2<Kernel> >
class Constrained_Delaunay_triangulation_face_base_2
  : public Fb
{
public:
  typedef Fb Base;
  typedef typename Base::Vertex_handle  Vertex_handle;
  typedef typename Base::Face_handle    Face_handle;
  typedef std::pair<Face_handle, int>   Edge_cdt;
  typedef Constrained_Delaunay_triangulation_face_base_2 CDT_face_base;

  typedef typename Kernel::Point_2      Point;
  typedef typename Kernel::Segment_2    Segment;

  struct Edge
  {
    Face_handle face;
    int index;
    Edge()
      : face(), index(-1)
    {}
    Edge(Face_handle face_, int index_)
      : face(face_), index(index_)
    {}
  };

private:
  bool m_blind;
  Edge m_blinding_constraint;

public:
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Base::template Rebind_TDS<TDS2>::Other Fb2;
    typedef Constrained_Delaunay_triangulation_face_base_2<Kernel,Fb2>
      Other;
  };

public:
  Constrained_Delaunay_triangulation_face_base_2()
    : Base(),
      m_blind(false)
  {
  }
  Constrained_Delaunay_triangulation_face_base_2( Vertex_handle v1,
                                                  Vertex_handle v2,
                                                  Vertex_handle v3)
    : Base(v1,v2,v3),
      m_blind(false)
  {
  }
  Constrained_Delaunay_triangulation_face_base_2( Vertex_handle v1,
                                                  Vertex_handle v2,
                                                  Vertex_handle v3,
                                                  Face_handle f1,
                                                  Face_handle f2,
                                                  Face_handle f3)
    : Base(v1,v2,v3,f1,f2,f3),
      m_blind(false)
  {
  }
  Constrained_Delaunay_triangulation_face_base_2(Face_handle f)
    : Base(f),
      m_blind(false)
  {
  }

  // sees its circumcenter or not?
  bool is_blind() const { return m_blind; }
  void set_blind(const bool b){ m_blind = b; }

  // if blind, the constrained edge that prevents the face
  // to see its circumcenter
  Edge_cdt blinding_constraint() const
  {
    CGAL_precondition(this->is_blind());
    return std::make_pair(m_blinding_constraint.face,
                          m_blinding_constraint.index);
  }
  void set_blinding_constraint(const Edge_cdt& e)
  {
    CGAL_precondition(this->is_blind());
    CGAL_precondition(e.first->is_constrained(e.second));
    m_blinding_constraint.face = e.first;
    m_blinding_constraint.index = e.second;
  }
};

} //namespace CGAL

#endif //CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_FACE_BASE_2_2_H
