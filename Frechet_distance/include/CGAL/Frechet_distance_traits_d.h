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

#ifndef CGAL_Frechet_distance_TRAITS_D_H
#define CGAL_Frechet_distance_TRAITS_D_H

#include <CGAL/license/Frechet_distance.h>

#include <CGAL/Dimension.h>

namespace CGAL
{
/*!
 * \ingroup PkgFrechetDistanceTraits
 *
* \cgalModels{FrechetDistanceTraits}
* \tparam K geometric traits class
*/
template <class K>
class Frechet_distance_traits_d
{
public:
  using Dimension = Dimension_tag<K::Dimension::value>;

  using Kernel = K;
  using FT = typename Kernel::FT;
  using Point_d = typename Kernel::Point_d;
  using Construct_bbox_d = typename Kernel::Construct_bbox_d;
  using Cartesian_const_iterator_d = typename Kernel::Cartesian_const_iterator_d;
  using Construct_cartesian_const_iterator_d = typename Kernel::Construct_cartesian_const_iterator_d;

  Construct_bbox_d construct_bbox_d_object() const {
     return Construct_bbox_d();
  }

  Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const {
     return Construct_cartesian_const_iterator_d();
  }
};

}  // end of namespace CGAL

#endif  // CGAL_Frechet_distance_TRAITS_D_H
