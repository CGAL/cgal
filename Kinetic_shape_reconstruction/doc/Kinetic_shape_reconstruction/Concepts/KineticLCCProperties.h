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

The concept `KineticLCCProperties` refines the concept of `LinearCellComplexItems`. It describes the item properties of the Linear_Cell_Complex used to store the resulting partition of `CGAL::Kinetic_shape_partition_3`.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Kinetic_shape_partition_3::LCC_Base_Properties`}
\cgalHasModelsEnd

\sa `CGAL::Kinetic_shape_partition_3`
\sa `LinearCellComplexItems`
*/

struct KineticLCCProperties {
  /// Using the index-based version of the `Linear_Cell_Complex`.
  typedef CGAL::Tag_true Use_index;
  typedef std::uint32_t Index_type;

  /// Model of `KineticLCCFaceProperty` concept.
  struct Face_property {} Face_property;
  /// Model of `KineticLCCVolumeProperty` concept.
  struct Volume_property {} Volume_property;
  /// Model of `KineticLCCDartWrapper` concept
  template<class LCC>
  struct Dart_wrapper {} Dart_wrapper; .
};
