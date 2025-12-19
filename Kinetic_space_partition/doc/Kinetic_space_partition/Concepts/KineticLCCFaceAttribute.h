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
\ingroup PkgKineticSpacePartitionConcepts
\cgalConcept

The concept `KineticLCCFaceAttribute` refines `CellAttribute` to store additional information for each face related to its associated input polygon.

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Cell_attribute `Cell_attribute<LCC, CGAL::Kinetic_space_partition_3::Linear_cell_complex_min_items::Face_attribute>`\endlink}
\cgalHasModelsEnd

\sa `CGAL::Kinetic_space_partition_3`
\sa `LinearCellComplexItems`
*/

struct KineticLCCFaceAttribute {
/// \name Access Members
/// @{
  /// 3D plane type compatible with `Kinetic_space_partition_3::Intersection_kernel`
  typedef unspecified_type Plane_3;
  /// Stores the index of the input polygon the provided that support plane of this face. Negative numbers correspond to the values defined in the enum `Kinetic_space_partition_3::Face_support`.
  int input_polygon_index;
  /// Support plane of the face derived from the corresponding input polygon or from octree nodes used for subdivision.
  Plane_3 plane;
  /// Does this face overlap with the corresponding input polygon.
  bool part_of_initial_polygon;
/// @}
};
