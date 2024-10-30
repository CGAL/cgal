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
#include <CGAL/Frechet_distance/internal/id.h>
#include <CGAL/Frechet_distance/internal/Get_exact_kernel.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Bbox.h>

// TODO: is it too restrictive to use vectors by default?
#include <vector>

namespace CGAL
{
/*!
 * \ingroup PkgFrechetDistanceRef
 *
* \cgalModels{FrechetDistanceTraits}
* \tparam GT  geometric traits class
*/
template <class GT>
class Frechet_distance_traits_d
{
public:
  using Dimension = Dimension_tag<GT::Dimension::value>;

  using Kernel = GT;
  using FT = typename Kernel::FT;
  using Point = typename Kernel::Point_d;

  static constexpr bool is_filtered = CGAL::Frechet_distance_::internal::Get_exact_kernel<Kernel>::is_filtered;
  static constexpr bool  is_floating_point = CGAL::Frechet_distance_::internal::Get_exact_kernel<Kernel>::is_floating_point;
/*
  using distance_t = Interval_nt<false>;
  using Approximate_kernel = Kernel_d_interface<Cartesian_base_d<distance_t,Dimension_tag<dimension>>>;
  using Approximate_point = typename Approximate_kernel::Point_d;
*/
  using Squared_distance = typename Kernel::Squared_distance_d;
/*
  using Exact_kernel = typename CGAL::Frechet_distance_::internal::Get_exact_kernel<Kernel>::type;
  using Exact_point = typename Exact_kernel::Point_d;

  using D2D = NT_converter<distance_t,double>;
  using A2E = KernelD_converter<Approximate_kernel, Exact_kernel, Default, D2D>;

  using FT2I = NT_converter<typename Kernel::FT,distance_t>;
  using K2A = KernelD_converter<Kernel, Approximate_kernel, Default, FT2I>;
*/
};

}  // end of namespace CGAL

#endif  // CGAL_Frechet_distance_TRAITS_D_H
