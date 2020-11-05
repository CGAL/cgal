// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_H
#define CGAL_CLASSIFICATION_H

#include <CGAL/license/Classification.h>

#include <CGAL/Classification/classify.h>
#include <CGAL/Classification/Sum_of_weighted_features_classifier.h>
#include <CGAL/Classification/ETHZ/Random_forest_classifier.h>

#ifdef CGAL_LINKED_WITH_OPENCV
#include <CGAL/Classification/OpenCV/Random_forest_classifier.h>
#endif

#ifdef CGAL_LINKED_WITH_TENSORFLOW
#include <CGAL/Classification/TensorFlow/Neural_network_classifier.h>
#endif

#include <CGAL/Classification/Cluster.h>
#include <CGAL/Classification/Evaluation.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label.h>
#include <CGAL/Classification/Label_set.h>
#include <CGAL/Classification/Local_eigen_analysis.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Point_set_feature_generator.h>
#include <CGAL/Classification/Point_set_neighborhood.h>
#include <CGAL/Classification/Mesh_feature_generator.h>
#include <CGAL/Classification/Mesh_neighborhood.h>
#include <CGAL/Classification/property_maps.h>

#include <CGAL/Classification/Feature/Cluster_mean_of_feature.h>
#include <CGAL/Classification/Feature/Cluster_size.h>
#include <CGAL/Classification/Feature/Cluster_variance_of_feature.h>
#include <CGAL/Classification/Feature/Cluster_vertical_extent.h>
#include <CGAL/Classification/Feature/Color_channel.h>
#include <CGAL/Classification/Feature/Distance_to_plane.h>
#include <CGAL/Classification/Feature/Echo_scatter.h>
#include <CGAL/Classification/Feature/Eigenvalue.h>
#include <CGAL/Classification/Feature/Elevation.h>
#include <CGAL/Classification/Feature/Simple_feature.h>
#include <CGAL/Classification/Feature/Vertical_dispersion.h>
#include <CGAL/Classification/Feature/Verticality.h>

#endif // CGAL_CLASSIFICATION_H
