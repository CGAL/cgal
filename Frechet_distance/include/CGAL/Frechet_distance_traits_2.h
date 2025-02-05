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

#ifndef CGAL_Frechet_distance_TRAITS_2_H
#define CGAL_Frechet_distance_TRAITS_2_H

#include <CGAL/license/Frechet_distance.h>

#include <CGAL/Dimension.h>

namespace CGAL
{
/*!
 * \ingroup PkgFrechetDistanceTraits
 *
* \cgalModels{FrechetDistanceTraits}
* \tparam K  geometric traits class. Must be a model of `Kernel`
*/
template <class K>
class Frechet_distance_traits_2
{
public:
  using Dimension = Dimension_tag<2>;

  using Kernel = K;
  using FT = typename Kernel::FT;
  using Point_d = typename Kernel::Point_2;
  using Construct_bbox_d = typename Kernel::Construct_bbox_2;
  using Cartesian_const_iterator_d = typename Kernel::Cartesian_const_iterator_2;
  using Construct_cartesian_const_iterator_d = typename Kernel::Construct_cartesian_const_iterator_2;
  using Compare_squared_distance_d = typename Kernel::Compare_squared_distance_2;

  Construct_bbox_d construct_bbox_d_object() const {
     return Construct_bbox_d();
  }

  Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const {
     return Construct_cartesian_const_iterator_d();
  }
  
  Compare_squared_distance_d construct_compare_squared_distance_d_object() const {
   return Compare_squared_distance_d();
  }
};

}  // end of namespace CGAL

#endif  // CGAL_Frechet_distance_TRAITS_2_H
