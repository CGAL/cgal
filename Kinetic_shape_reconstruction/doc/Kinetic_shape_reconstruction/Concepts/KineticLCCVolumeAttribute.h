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

The concept `KineticLCCVolumeAttribute` refines `CellAttribute` to store the barycenter and an id.

\cgalHasModelsBegin
\cgalHasModelsBare{`CellAttribute<LCC, CGAL::Kinetic_shape_partition_3::LCC_Base_Attribute::Volume_attribute>`}
\cgalHasModelsEnd

\sa `CGAL::Kinetic_shape_partition_3`
\sa `LinearCellComplexItems`
*/

struct KineticLCCVolumeAttribute {
/// \name Access Members
/// @{
  /// Contains the bary_cernter of the volume.
  typename Intersection_kernel::Point_3 barycenter;
  /// 0-based volume id.
  std::size_t volume_id;
/// @}
};
