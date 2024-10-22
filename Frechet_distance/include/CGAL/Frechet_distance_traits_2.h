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
#include <CGAL/Frechet_distance/internal/id.h>
#include <CGAL/Frechet_distance/internal/Get_exact_kernel.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian_converter.h>

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
class Frechet_distance_traits_2
{
public:
    static const int dimension = 2;
    using Kernel = GT;
    using FT = typename Kernel::FT;
    using Point = typename Kernel::Point_2;
    using Bbox = Bbox_2;

    static constexpr bool is_filtered = CGAL::Frechet_distance_::internal::Get_exact_kernel<Kernel>::is_filtered;
    static constexpr bool  is_floating_point = CGAL::Frechet_distance_::internal::Get_exact_kernel<Kernel>::is_floating_point;

    using distance_t = CGAL::Interval_nt<false>;
    using Approximate_kernel = CGAL::Simple_cartesian<distance_t>;
    using Approximate_point = typename Approximate_kernel::Point_2;
    using Construct_bbox = typename Approximate_kernel::Construct_bbox_2;
    using Squared_distance = typename Approximate_kernel::Compute_squared_distance_2;
    using Difference_of_points = typename Approximate_kernel::Construct_vector_2;
    using Scaled_vector = typename Approximate_kernel::Construct_scaled_vector_2;
    using Translated_point = typename Approximate_kernel::Construct_translated_point_2;

    using Exact_kernel = typename CGAL::Frechet_distance_::internal::Get_exact_kernel<Kernel>::type;
    using Exact_point = typename Exact_kernel::Point_2;

    using D2D = NT_converter<distance_t,double>;
    using F2E = Cartesian_converter<Approximate_kernel, Exact_kernel, D2D>;

    using FT2I = NT_converter<typename Kernel::FT,distance_t>;
    using K2F = Cartesian_converter<Kernel, Approximate_kernel, FT2I>;

    /// @todo remove?
    using Polyline = std::vector<Point>;
    using Polylines = std::vector<Polyline>;
    using PolylineID = std::size_t;
    using PolylineIDs = std::vector<PolylineID>;
};

}  // end of namespace CGAL

#endif  // CGAL_Frechet_distance_TRAITS_2_H
