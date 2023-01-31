// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot, Dmitry Anisimov, Sven Oesau

#ifndef CGAL_KINETIC_SHAPE_PARTITIONING_TRAITS_3_H
#define CGAL_KINETIC_SHAPE_PARTITIONING_TRAITS_3_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/Cartesian_converter.h>

namespace CGAL {

/*!
  \ingroup PkgKineticPartition
  \brief %Default traits class for the `CGAL::Kinetic_shape_partition_3`.

  \cgalModels `KineticShapePartitioningTraits_3`

  \tparam GeomTraits must be a model of the concept `Kernel`. This Kernel is used for non critical calculations and assumed for the input data.

  \tparam IntersectionTraits must be a model of the concept `Kernel`. A Kernel with exact constructions is advised.

  \tparam InputRange must be a model of `Range` with random access iterators, providing input points through the following property map.

  \tparam PointMap must be a model of `ReadablePropertyMap` with `std::iterator_traits<Input_range::iterator>::%value_type` as key type and `Geom_traits::Point_3` as value type.
*/
template<typename GeomTraits, typename IntersectionTraits, typename InputRange, typename PointMap>
class Kinetic_shape_partitioning_traits_3 {

public:
  using Kernel = GeomTraits;
  using Intersection_Kernel = IntersectionTraits;
  using Input_range = InputRange; /// must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`and value type is `Point_3`
  using Point_map = PointMap;

  using FT = typename Kernel::FT;

  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;

  using Vector_2 = typename Kernel::Vector_2;
  using Vector_3 = typename Kernel::Vector_3;

  using Direction_2 = typename Kernel::Direction_2;

  using Segment_2 = typename Kernel::Segment_2;
  using Segment_3 = typename Kernel::Segment_3;

  using Line_2 = typename Kernel::Line_2;
  using Line_3 = typename Kernel::Line_3;

  using Plane_3 = typename Kernel::Plane_3;

  using IK_Point_2 = typename Intersection_Kernel::Point_2;
  using IK_Point_3 = typename Intersection_Kernel::Point_3;
  using IK_Segment_3 = typename Intersection_Kernel::Segment_3;
  using IK_Line_3 = typename Intersection_Kernel::Line_3;

  using Transform_3 = typename Kernel::Aff_transformation_3;

  using Triangle_2 = typename Kernel::Triangle_2;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Tetrahedron_3 = typename Kernel::Tetrahedron_3;

  using From_exact = CGAL::Cartesian_converter<Intersection_Kernel, Kernel>;
  using To_exact = CGAL::Cartesian_converter<Kernel, Intersection_Kernel>;

};

} // namespace CGAL

#endif // CGAL_KINETIC_SHAPE_PARTITIONING_TRAITS_3_H
