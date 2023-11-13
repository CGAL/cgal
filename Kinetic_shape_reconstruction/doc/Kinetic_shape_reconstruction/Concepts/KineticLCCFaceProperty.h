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

/*!
\ingroup PkgKineticShapePartitionConcepts
\cgalConcept

The concept `KineticLCCFaceProperty` is part of the item properties of the `Linear_Cell_Complex` used to store the resulting partition of `CGAL::Kinetic_shape_partition_3`.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Kinetic_shape_partition_3::LCC_Base_Properties::Face_property`}
\cgalHasModelsEnd

\sa `CGAL::Kinetic_shape_partition_3`
\sa `LinearCellComplexItems`
*/

struct KineticLCCFaceProperty {
  /// Stores the index of the input polygon the provided that support plane of this face. Indices -1 till -6 correspond to bbox faces, -7 to faces from octree
  int input_polygon_index;
  /// Support plane of the face derived from the corresponding input polygon or from octree nodes used for subdivision.
  typename Intersection_kernel::Plane_3 plane;
  /// Does this face overlap with the corresponding input polygon.
  bool part_of_initial_polygon;
};
