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

The concept `KineticLCCItems` refines the concept of `LinearCellComplexItems` by adding 2-attributes and 3-attributes to store information from the kinetic partition.

The third type in Attributes tuple must be a model of the `KineticLCCFaceAttribute` concept.
The fourth type in Attributes tuple must be a model of the `KineticLCCVolumeAttribute` concept.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Kinetic_shape_partition_3::Lcc_min_items`}
\cgalHasModelsEnd

\sa `CGAL::Kinetic_shape_partition_3`
\sa `KineticLCCFaceAttribute`
\sa `KineticLCCVolumeAttribute`
\sa `LinearCellComplexItems`
*/

struct KineticLCCItems {
/// \name Types
/// @{
  /// Using the index-based version of the `LinearCellComplex` concept.
  typedef CGAL::Tag_true Use_index;
  typedef std::uint32_t Index_type;
/// @}
};
