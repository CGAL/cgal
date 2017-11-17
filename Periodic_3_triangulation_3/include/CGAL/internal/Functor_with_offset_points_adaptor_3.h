// Copyright (c) 1999-2004,2006-2009,2014-2015   INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Aymeric Pellé <Aymeric.Pelle@sophia.inria.fr>
//                 Mael Rouxel-Labbé
#ifndef CGAL_FUNCTOR_WITH_OFFSET_POINTS_ADAPTOR_3_H
#define CGAL_FUNCTOR_WITH_OFFSET_POINTS_ADAPTOR_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

namespace CGAL {

template < class K_, class Functor_ >
class Functor_with_offset_points_adaptor_3
  : public Functor_
{
  typedef K_                                    Kernel;
  typedef Functor_                              Functor;

  typedef typename Kernel::Point_3              Point;
  typedef typename Kernel::Offset               Offset;

  typedef typename Kernel::Construct_point_3    Construct_point_3;

public:
  typedef typename Functor::result_type         result_type;

  Functor_with_offset_points_adaptor_3(const Functor& functor,
                                       const Construct_point_3& cp)
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

  const Construct_point_3 cp;
};

}  // namespace CGAL

#endif /* CGAL_FUNCTOR_WITH_OFFSET_POINTS_ADAPTOR_3_H */
