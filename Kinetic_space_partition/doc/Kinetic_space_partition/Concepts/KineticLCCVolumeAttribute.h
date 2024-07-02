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

The concept `KineticLCCVolumeAttribute` refines `CellAttribute` to store the barycenter and an id.

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Cell_attribute `Cell_attribute<LCC, CGAL::Kinetic_space_partition_3::Linear_cell_complex_min_items::Volume_attribute>`\endlink}
\cgalHasModelsEnd

\sa `CGAL::Kinetic_space_partition_3`
\sa `LinearCellComplexItems`
*/

struct KineticLCCVolumeAttribute {
/// \name Access Members
/// @{
  /// 3D point type compatible with `Kinetic_space_partition_3::Intersection_kernel`
  typedef unspecified_type Point_3;
  /// Contains the barycenter of the volume.
  Point_3 barycenter;
  /// 0-based volume id.
  std::size_t volume_id;
/// @}
};
