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
  // Using the index-based version of the `Linear_Cell_Complex`.
  typedef CGAL::Tag_true Use_index;
  typedef std::uint32_t Index_type;

  struct Face_property {
    int input_polygon_index; // Stores the index of the input polygon the provided that support plane of this face. Indices -1 till -6 correspond to bbox faces, -7 to faces from octree
    typename Intersection_kernel::Plane_3 plane; // Support plane of the face derived from the corresponding input polygon or from octree nodes used for subdivision.
    bool part_of_initial_polygon; // Does this face overlap with the corresponding input polygon.
  } Face_property;

  struct Volume_property {
    typename Intersection_kernel::Point_3 barycenter; // Contains the bary_cernter of the volume.
    std::size_t volume_id; // 0-based volume id.
  } Volume_property;

  template<class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute_with_point< LCC, void> Vertex_attribute;
    typedef CGAL::Cell_attribute< LCC, Face_property > Face_attribute;
    typedef CGAL::Cell_attribute< LCC, Volume_property > Volume_attribute;

    typedef std::tuple<Vertex_attribute, void, Face_attribute, Volume_attribute> Attributes;
  } Dart_wrapper;
};