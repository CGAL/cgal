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
#ifndef CGAL_PERIODIC_3_CONSTRUCT_POINT_3_H
#define CGAL_PERIODIC_3_CONSTRUCT_POINT_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

namespace CGAL {

template < typename K, typename Construct_point_3_base>
class Periodic_3_construct_point_3
  : public Construct_point_3_base
{
  typedef Construct_point_3_base             Base;

  typedef K                                  Kernel;

  typedef typename Kernel::Point_3           Point;
  typedef typename Kernel::Weighted_point_3  Weighted_point;
  typedef typename Kernel::Offset            Offset;
  typedef typename Kernel::Iso_cuboid_3      Iso_cuboid_3;

public:
  Periodic_3_construct_point_3(const Iso_cuboid_3* dom,
                               const Construct_point_3_base& cp)
    : Base(cp), _dom(dom)
  { }

  using Base::operator();

  Point operator() ( const Point& p, const Offset& o ) const {
    return operator()(p.x() + (_dom->xmax() - _dom->xmin()) * o.x(),
                      p.y() + (_dom->ymax() - _dom->ymin()) * o.y(),
                      p.z() + (_dom->zmax() - _dom->zmin()) * o.z());
  }

  Point operator() ( const Weighted_point& p, const Offset& o ) const {
    return operator()(p.x() + (_dom->xmax() - _dom->xmin()) * o.x(),
                      p.y() + (_dom->ymax() - _dom->ymin()) * o.y(),
                      p.z() + (_dom->zmax() - _dom->zmin()) * o.z());
  }

private:
  const Iso_cuboid_3* _dom;
};

} // namespace CGAL

#endif // CGAL_PERIODIC_3_CONSTRUCT_POINT_3_H
