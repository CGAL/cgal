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
#include <CGAL/Interval_nt.h>
#include <CGAL/Simple_cartesian.h>
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
    static constexpr int dimension = GT::Dimension::value;

    using Kernel = GT;
    using FT = typename Kernel::FT;
    using Point = typename Kernel::Point_d;

    /// @todo  remove?
    using Polyline = std::vector<Point>;
    using Polylines = std::vector<Polyline>;
    using PolylineID = std::size_t;
    using PolylineIDs = std::vector<PolylineID>;

};

}  // end of namespace CGAL

#endif  // CGAL_Frechet_distance_TRAITS_D_H
