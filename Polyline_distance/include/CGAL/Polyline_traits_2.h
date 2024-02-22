// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//
// =============================================================================

#ifndef CGAL_POLYLINE_TRAITS_2_H
#define CGAL_POLYLINE_TRAITS_2_H

// TODO: add doxygen code

/// \cgalModels{PolylineDistanceTraits}


// TODO: is it too restrictive to use vectors by default?
#include <vector>

namespace CGAL
{

// TODO: This is just a starter to enable using the same types all over the package.
template <class T, class CNT>
class Polyline_traits_2
{
public:
	using BaseTraits = T;
	// TODO: remove FT?
	using FT = typename BaseTraits::FT;
	using Point = typename BaseTraits::Point_2;
	using Compare_squared_distance = typename BaseTraits::Compare_squared_distance_2;
	using Construct_midpoint = typename BaseTraits::Construct_midpoint_2;
	using Construction_number_type = CNT;
	// TODO: remove?
	using Polyline = std::vector<Point>;
	using Polylines = std::vector<Polyline>;
	using PolylineID = std::size_t;
	using PolylineIDs = std::vector<PolylineID>;
};

} // end of namespace CGAL

#endif // CGAL_POLYLINE_TRAITS_2_H
