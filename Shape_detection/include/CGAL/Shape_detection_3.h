// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez, Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_3_DEPRECATED_H
#define CGAL_SHAPE_DETECTION_3_DEPRECATED_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>

// Deprecated -->
#include <CGAL/Shape_detection/deprecated/Region_growing.h>
#include <CGAL/Shape_detection/deprecated/Shape_detection_traits.h>

#define CGAL_DEPRECATED_MESSAGE_DETAILS \
  "CGAL::Shape_detection_3 namespace has been replaced by the namespace "\
  "CGAL::Shape_detection."

#include <CGAL/internal/deprecation_warning.h>

namespace CGAL {
  namespace Shape_detection_3 = Shape_detection;
}
// --<

#endif // CGAL_SHAPE_DETECTION_3_DEPRECATED_H
