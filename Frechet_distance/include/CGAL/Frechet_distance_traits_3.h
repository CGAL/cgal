// Copyright (c) 2024 Max-Planck-Institute Saarbruecken (Germany), CNRS (France), GeometryFactory (France)
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

#ifndef CGAL_Frechet_distance_TRAITS_3_H
#define CGAL_Frechet_distance_TRAITS_3_H

#include <CGAL/license/Frechet_distance.h>

#include <CGAL/Dimension.h>

namespace CGAL
{
/*!
 * \ingroup PkgFrechetDistanceTraits
 *
* \cgalModels{FrechetDistanceTraits}
* \tparam K  geometric traits class
*/
template <class K>
class Frechet_distance_traits_3
{
public:
  using Dimension = Dimension_tag<3>;

  using Kernel = K;
  using FT = typename Kernel::FT;
  using Point_d = typename Kernel::Point_3;
  using Construct_bbox_d = typename Kernel::Construct_bbox_3;
  using Cartesian_const_iterator_d = typename Kernel::Cartesian_const_iterator_3;
  using Construct_cartesian_const_iterator_d = typename Kernel::Construct_cartesian_const_iterator_3;
  using Compare_squared_distance_d = typename Kernel::Compare_squared_distance_3;

  Construct_bbox_d construct_bbox_d_object() const {
     return Construct_bbox_d();
  }

  Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const {
     return Construct_cartesian_const_iterator_d();
  }
};

}  // end of namespace CGAL

#endif  // CGAL_Frechet_distance_TRAITS_3_H
