// Copyright (c) 2024 Max-Planck-Institute Saarbruecken (Germany), GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//                 Andreas Fabri
// =============================================================================

#ifndef CGAL_FRECHET_DISTANCE_INTERNAL_FRECHET_DISTANCE_TRAITS_H
#define CGAL_FRECHET_DISTANCE_INTERNAL_FRECHET_DISTANCE_TRAITS_H

#include <CGAL/license/Frechet_distance.h>

#include <CGAL/Bbox_d.h>
#include <array>

namespace CGAL {
namespace Frechet_distance {
namespace internal {

template <class NT, int dimension>
class Frechet_distance_traits
{
public:
  using Dimension = Dimension_tag<dimension>;

  using FT = NT;
  using Point_d = std::array<FT, dimension>;
  using Cartesian_const_iterator_d = typename Point_d::const_iterator;

  struct Construct_bbox_d
  {
    Bbox_d<Dimension> operator()(const Point_d& p) const
    {
      Bbox_d<Dimension> bb;
      for (int i=0;i<dimension; ++i)
      {
        (bb.min)(i)=to_interval(p[i]).first;
        (bb.max)(i)=to_interval(p[i]).second;
      }

      return bb;
    }
  };

  struct Construct_cartesian_const_iterator_d
  {
    Cartesian_const_iterator_d operator()(const Point_d& p) const
    {
      return p.begin();
    }
    Cartesian_const_iterator_d operator()(const Point_d& p, int) const
    {
      return p.end();
    }
  };

  Construct_bbox_d construct_bbox_d_object() const { return Construct_bbox_d(); }

  Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const
  {
     return Construct_cartesian_const_iterator_d();
  }
};

}
}
}  // end of namespace CGAL

#endif  // CGAL_FRECHET_DISTANCE_INTERNAL_FRECHET_DISTANCE_TRAITS_H
