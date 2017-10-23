// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_H
#define CGAL_CLASSIFICATION_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/classify.h>
#include <CGAL/Classification/Sum_of_weighted_features_classifier.h>
#include <CGAL/Classification/ETHZ_random_forest_classifier.h>

#ifdef CGAL_LINKED_WITH_OPENCV
#include <CGAL/Classification/OpenCV_random_forest_classifier.h>
#endif

#include <CGAL/Classification/Color.h>
#include <CGAL/Classification/Evaluation.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label.h>
#include <CGAL/Classification/Label_set.h>
#include <CGAL/Classification/Local_eigen_analysis.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Point_set_feature_generator.h>
#include <CGAL/Classification/Point_set_neighborhood.h>

#include <CGAL/Classification/Feature/Distance_to_plane.h>
#include <CGAL/Classification/Feature/Echo_scatter.h>
#include <CGAL/Classification/Feature/Eigen.h>
#include <CGAL/Classification/Feature/Elevation.h>
#include <CGAL/Classification/Feature/Hsv.h>
#include <CGAL/Classification/Feature/Simple_feature.h>
#include <CGAL/Classification/Feature/Vertical_dispersion.h>
#include <CGAL/Classification/Feature/Verticality.h>

#endif // CGAL_CLASSIFICATION_H
