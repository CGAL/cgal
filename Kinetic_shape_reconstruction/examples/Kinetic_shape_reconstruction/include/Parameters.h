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

    // Label indices defined in the ply header:
    // ground (gi),
    // building boundary (bi),
    // building interior (ii),
    // vegetation (vi).
    std::string gi, bi, ii, vi;

    // Main parameters.
    FT scale; // meters
    FT noise; // meters

    // Boolean tags.
    const bool with_normals; // do we use normals
    bool verbose;// verbose basic info
    bool debug; // verbose more info

    // Shape detection / Shape regularization.
    // See the corresponding CGAL packages.
    std::size_t k_neighbors;
    FT distance_threshold;
    FT angle_threshold;
    std::size_t min_region_size;
    bool regularize;

    // Partitioning.
    // See KSR/parameters.h
    unsigned int k_intersections;
    const unsigned int n_subdivisions;
    const FT enlarge_bbox_ratio;
    const bool reorient;

    // Reconstruction.
    FT graphcut_beta; // magic parameter between 0 and 1

    // Constructor.
    All_parameters() :
    data(""),
    gi("0"), bi("1"), ii("2"), vi("3"), // semantic labels mapping
    // main parameters
    scale(FT(4)),
    noise(FT(2)),
    // boolean tags
    with_normals(true),
    verbose(false),
    debug(false),
    // shape detection / shape regularization
    k_neighbors(12),
    distance_threshold(noise / FT(2)),
    angle_threshold(FT(25)),
    min_region_size(100),
    regularize(false),
    // partition
    k_intersections(1),
    n_subdivisions(0),
    enlarge_bbox_ratio(FT(11) / FT(10)),
    reorient(false),
    // reconstruction
    graphcut_beta(FT(1) / FT(2))
    { }

    // Update all parameters, which depend on scale and noise.
    void update_dependent() {
      distance_threshold = noise / FT(2);
    }
  };

} // KSR
} // CGAL

#endif // CGAL_KSR_ALL_PARAMETERS_EXAMPLES_H
