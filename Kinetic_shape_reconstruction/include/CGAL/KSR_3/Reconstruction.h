// Copyright (c) 2019 GeometryFactory SARL (France).
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

#ifndef CGAL_KSR_3_RECONSTRUCTION_H
#define CGAL_KSR_3_RECONSTRUCTION_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR_3/Data_structure.h>

namespace CGAL {
namespace KSR_3 {

template<typename GeomTraits>
class Reconstruction {

public:
  using Kernel = GeomTraits;

private:
  using FT        = typename Kernel::FT;
  using Point_3   = typename Kernel::Point_3;
  using Segment_3 = typename Kernel::Segment_3;

  using Data_structure = KSR_3::Data_structure<Kernel>;

public:

  template<
  typename InputRange,
  typename PointMap,
  typename VectorMap,
  typename LabelMap>
  Reconstruction(
    const InputRange& input_range,
    const PointMap point_map,
    const VectorMap normal_map,
    const LabelMap label_map,
    Data_structure& data) :
  m_data(data) {
    CGAL_assertion_msg(false, "TODO: RECONSTRUCTION, INITIALIZE!");
  }

  template<typename NamedParameters>
  void detect_planar_shapes(
    const NamedParameters& np) {
    CGAL_assertion_msg(false, "TODO: RECONSTRUCTION, DETECT PLANAR SHAPES!");
  }

  template<typename NamedParameters>
  void partition(
    const NamedParameters& np) {
    CGAL_assertion_msg(false, "TODO: RECONSTRUCTION, PARTITION!");
  }

  template<typename NamedParameters>
  void compute_model(
    const NamedParameters& np) {
    CGAL_assertion_msg(false, "TODO: RECONSTRUCTION, COMPUTE MODEL!");
  }

private:
  Data_structure& m_data;
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_RECONSTRUCTION_H
