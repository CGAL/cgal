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

The concept `KineticLCCDartWrapper` is part of the item properties of the `Linear_Cell_Complex` used to store the resulting partition of `CGAL::Kinetic_shape_partition_3`.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Kinetic_shape_partition_3::LCC_Base_Properties::Dart_wrapper`}
\cgalHasModelsEnd

\sa `CGAL::Kinetic_shape_partition_3`
\sa `LinearCellComplexItems`
*/

struct KineticLCCDartWrapper {
    /// A Cell_attribute_with_point is required as the vertex property, but an additional type can be specified as the second template parameter.
    typedef CGAL::Cell_attribute_with_point< LCC, void> Vertex_attribute;
    /// Using a model of the concept `KineticLCCFaceProperty`
    typedef CGAL::Cell_attribute< LCC, Face_property > Face_attribute;
    /// Using a model of the concept `KineticLCCVolumeProperty`
    typedef CGAL::Cell_attribute< LCC, Volume_property > Volume_attribute;

    /// The Attributes tuple contains the Properties for the 0,1,2 and 3-cells, i.e., vertex, dart, face and volume properties.
    typedef std::tuple<Vertex_attribute, void, Face_attribute, Volume_attribute> Attributes;
};
