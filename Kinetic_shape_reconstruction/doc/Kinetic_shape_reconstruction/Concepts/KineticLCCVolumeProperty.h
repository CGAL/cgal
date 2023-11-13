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

The concept `KineticLCCVolumeProperty` is part of the item properties of the `Linear_Cell_Complex` used to store the resulting partition of `CGAL::Kinetic_shape_partition_3`.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Kinetic_shape_partition_3::LCC_Base_Properties::Volume_property`}
\cgalHasModelsEnd

\sa `CGAL::Kinetic_shape_partition_3`
\sa `LinearCellComplexItems`
*/

struct KineticLCCVolumeProperty {
    typename Intersection_kernel::Point_3 barycenter; // Contains the bary_cernter of the volume.
    std::size_t volume_id; // 0-based volume id.
};
