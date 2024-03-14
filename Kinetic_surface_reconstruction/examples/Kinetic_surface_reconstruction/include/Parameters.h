// Copyright (c) 2023 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Florent Lafarge, Dmitry Anisimov, Simon Giraudot

#ifndef CGAL_KSR_ALL_PARAMETERS_EXAMPLES_H
#define CGAL_KSR_ALL_PARAMETERS_EXAMPLES_H

// STL includes.
#include <string>

namespace CGAL {
namespace KSR {

  template<typename FT>
  struct All_parameters {

    // Path to the input data file.
    std::string data; // required!

    // Boolean tags.
    bool with_normals; // do we use normals
    bool verbose;// verbose basic info
    bool debug; // verbose more info

    // Shape detection / Shape regularization.
    // See the corresponding CGAL packages.
    std::size_t k_neighbors;
    FT maximum_distance;
    FT maximum_angle;
    std::size_t min_region_size;

    bool regparallel, regcoplanar, regorthogonal, regsymmetric;
    FT angle_tolerance;
    FT maximum_offset;

    // Partitioning.
    // See KSR/parameters.h
    unsigned int k_intersections;
    FT enlarge_bbox_ratio;
    bool reorient;

    std::size_t max_octree_depth;
    std::size_t max_octree_node_size;

    // Reconstruction.
    FT graphcut_beta; // magic parameter between 0 and 1
    bool use_ground;
    //std::map<typename KSP::Face_support, bool> external_nodes;

    // Constructor.
    All_parameters() :
      data(""),
      // boolean tags
      with_normals(true),
      verbose(false),
      debug(false),
      // shape detection / shape regularization
      k_neighbors(12),
      min_region_size(0),
      max_octree_node_size(40),
      max_octree_depth(3),
      maximum_distance(0),
      maximum_angle(10),
      regparallel(false),
      regcoplanar(true),
      regorthogonal(false),
      regsymmetric(false),
      angle_tolerance(5),
      maximum_offset(0),
      // partition
      k_intersections(1),
      enlarge_bbox_ratio(FT(11) / FT(10)),
      reorient(false),
      use_ground(false),
      // reconstruction
      graphcut_beta(FT(1) / FT(2))
    { }

  };

} // KSR
} // CGAL

#endif // CGAL_KSR_ALL_PARAMETERS_EXAMPLES_H
