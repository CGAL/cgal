// Copyright (c) 1997, 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>

#ifndef CGAL_ALPHA_SHAPE_FACE_BASE_2_H
#define CGAL_ALPHA_SHAPE_FACE_BASE_2_H

#include <CGAL/license/Alpha_shapes_2.h>


#include <CGAL/utility.h>
#include <CGAL/internal/Lazy_alpha_nt_2.h>
#include <CGAL/Default.h>
#include <CGAL/Triangulation_face_base_2.h>

namespace CGAL {


template <class Gt,
          class Fb_ = Default,
          class ExactAlphaComparisonTag = Tag_false,
          class Weighted_tag = Tag_false>
class Alpha_shape_face_base_2
  : public Default::Get<Fb_, Triangulation_face_base_2<Gt> >::type
{
  typedef typename Default::Get<Fb_, Triangulation_face_base_2<Gt> >::type Fb;

public:
  typedef typename Fb::Vertex_handle                          Vertex_handle;
  typedef typename Fb::Face_handle                            Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other     Fb2;
    typedef Alpha_shape_face_base_2<
      Gt, Fb2, ExactAlphaComparisonTag, Weighted_tag>         Other;
  };

  typedef typename internal::Alpha_nt_selector_2<
    Gt, ExactAlphaComparisonTag, Weighted_tag>::Type_of_alpha Type_of_alpha;
  typedef Type_of_alpha                                       NT;

  typedef Type_of_alpha FT;
  typedef Triple<Type_of_alpha, Type_of_alpha, Type_of_alpha> Interval_3;

private:
  Interval_3 vec_edge[3];
  Type_of_alpha A;

public:
  Alpha_shape_face_base_2()  : Fb()     {}

  Alpha_shape_face_base_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2)
    : Fb(v0, v1, v2)     {}

  Alpha_shape_face_base_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
                          Face_handle n0, Face_handle n1, Face_handle n2)
    : Fb(v0, v1, v2, n0, n1, n2)
    {}

  const Type_of_alpha & get_alpha() const
    {
      return A;
    }

  void set_alpha(const Type_of_alpha & AA)
    {
      A = AA;
    }

  const Interval_3 & get_ranges(int i) const
    {
      return vec_edge[i];
    }

  void set_ranges(int i, const Interval_3& Inter)
    {
      vec_edge[i]=Inter;
    }


};

} //namespace CGAL

#endif //ALPHA_SHAPE_FACE_BASE_2_H
