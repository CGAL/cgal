// Copyright (c) 2017 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_FUNCTOR_WITH_OFFSET_POINTS_ADAPTOR_2_H
#define CGAL_FUNCTOR_WITH_OFFSET_POINTS_ADAPTOR_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

namespace CGAL {

template < class K_, class Functor_ >
class Functor_with_offset_points_adaptor_2
  : public Functor_
{
  typedef K_                                    Kernel;
  typedef Functor_                              Functor;

  typedef typename Kernel::Point_2              Point;
  typedef typename Kernel::Offset               Offset;

  typedef typename Kernel::Construct_point_2    Construct_point_2;

public:
  typedef typename Functor::result_type         result_type;

  Functor_with_offset_points_adaptor_2(const Functor& functor,
                                       const Construct_point_2& cp)
    : Functor_(functor), cp(cp)
  { }

  // gives access to function calls without offset
  using Functor::operator();

  result_type operator()(const Point& p0, const Point& p1,
                         const Offset& o0, const Offset& o1) const {
    return operator()(cp(p0,o0), cp(p1,o1));
  }
  result_type operator()(const Point& p0, const Point& p1, const Point& p2,
                         const Offset& o0, const Offset& o1, const Offset& o2) const {
    return operator()(cp(p0,o0), cp(p1,o1), cp(p2,o2));
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3,
                         const Offset& o0, const Offset& o1,
                         const Offset& o2, const Offset& o3) const {
    return operator()(cp(p0,o0), cp(p1,o1), cp(p2,o2), cp(p3,o3));
  }
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3, const Point& p4,
                         const Offset& o0, const Offset& o1, const Offset& o2,
                         const Offset& o3, const Offset& o4) const {
    return operator()(cp(p0,o0), cp(p1,o1), cp(p2,o2), cp(p3,o3), cp(p4,o4));
  }

  const Construct_point_2 cp;
};

}  // namespace CGAL

#endif /* CGAL_FUNCTOR_WITH_OFFSET_POINTS_ADAPTOR_2_H */
