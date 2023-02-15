// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_WEIGHTS_H
#define CGAL_WEIGHTS_H

/**
* \ingroup PkgWeightsRef
* \file CGAL/Weights.h
* A convenience header that includes all weights.
*/

/// \cond SKIP_IN_MANUAL
#include <CGAL/Weights/utils.h>
/// \endcond

#include <CGAL/Weights/uniform_weights.h>

#include <CGAL/Weights/shepard_weights.h>
#include <CGAL/Weights/inverse_distance_weights.h>

#include <CGAL/Weights/three_point_family_weights.h>
#include <CGAL/Weights/wachspress_weights.h>
#include <CGAL/Weights/authalic_weights.h>
#include <CGAL/Weights/mean_value_weights.h>
#include <CGAL/Weights/tangent_weights.h>
#include <CGAL/Weights/discrete_harmonic_weights.h>
#include <CGAL/Weights/cotangent_weights.h>

#include <CGAL/Weights/uniform_region_weights.h>
#include <CGAL/Weights/triangular_region_weights.h>
#include <CGAL/Weights/barycentric_region_weights.h>
#include <CGAL/Weights/voronoi_region_weights.h>
#include <CGAL/Weights/mixed_voronoi_region_weights.h>

#endif // CGAL_WEIGHTS_H
