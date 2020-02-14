// Copyright (c) 2017   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_3_CONSTRUCT_WEIGHTED_POINT_3_H
#define CGAL_PERIODIC_3_CONSTRUCT_WEIGHTED_POINT_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

namespace CGAL {

template <typename K_, typename Construct_weighted_point_3_base_>
class Periodic_3_construct_weighted_point_3
  : public Construct_weighted_point_3_base_
{
  typedef Construct_weighted_point_3_base_      Base;

  typedef K_                                    Kernel;

  typedef typename Kernel::Point_3              Point_3;
  typedef typename Kernel::Weighted_point_3     Weighted_point_3;
  typedef typename Kernel::Offset               Offset;
  typedef typename Kernel::Iso_cuboid_3         Iso_cuboid_3;

public:
  Periodic_3_construct_weighted_point_3(const Iso_cuboid_3* dom,
                                        const Base& wp)
    : Base(wp), _dom(dom)
  { }

  using Base::operator();

  Weighted_point_3 operator() (const Weighted_point_3& p, const Offset& o) const
  {
    return Weighted_point_3(Point_3(p.x() + (_dom->xmax() - _dom->xmin()) * o.x(),
                                    p.y() + (_dom->ymax() - _dom->ymin()) * o.y(),
                                    p.z() + (_dom->zmax() - _dom->zmin()) * o.z()),
                            p.weight());
  }

private:
  const Iso_cuboid_3* _dom;
};

} // namespace CGAL

#endif // CGAL_PERIODIC_3_CONSTRUCT_WEIGHTED_POINT_3_H
