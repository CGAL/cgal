// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_H

/// \cond SKIP_IN_MANUAL
#include <CGAL/license/Shape_detection.h>
/// \endcond

/**
* \ingroup PkgShapeDetectionRef
* \file CGAL/Shape_detection/Region_growing/Point_set.h
* A convenience header that includes all classes related to the region growing algorithm on a point set.
*/

#include <CGAL/Shape_detection/Region_growing/Point_set/K_neighbor_query.h>
#include <CGAL/Shape_detection/Region_growing/Point_set/Sphere_neighbor_query.h>

#include <CGAL/Shape_detection/Region_growing/Point_set/Least_squares_line_fit_region.h>
#include <CGAL/Shape_detection/Region_growing/Point_set/Least_squares_circle_fit_region.h>
#include <CGAL/Shape_detection/Region_growing/Point_set/Least_squares_plane_fit_region.h>
#include <CGAL/Shape_detection/Region_growing/Point_set/Least_squares_sphere_fit_region.h>
#include <CGAL/Shape_detection/Region_growing/Point_set/Least_squares_cylinder_fit_region.h>

#include <CGAL/Shape_detection/Region_growing/Point_set/Least_squares_line_fit_sorting.h>
#include <CGAL/Shape_detection/Region_growing/Point_set/Least_squares_circle_fit_sorting.h>
#include <CGAL/Shape_detection/Region_growing/Point_set/Least_squares_plane_fit_sorting.h>
#include <CGAL/Shape_detection/Region_growing/Point_set/Least_squares_sphere_fit_sorting.h>
#include <CGAL/Shape_detection/Region_growing/Point_set/Least_squares_cylinder_fit_sorting.h>

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_H
