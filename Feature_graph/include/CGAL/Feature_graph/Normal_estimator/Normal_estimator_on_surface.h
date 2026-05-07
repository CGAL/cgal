// Copyright (c) 2026 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ange Clement

#ifndef CGAL_FEATURE_GRAPH_NORMAL_ESTIMATOR_ON_SURFACE_H
#define CGAL_FEATURE_GRAPH_NORMAL_ESTIMATOR_ON_SURFACE_H

namespace CGAL {

namespace Feature_graph {

namespace Normal_estimator
{

/*!
* \ingroup PkgFeatureGraphNormalEstimator
*
* \brief Functor that assign a normal on elements of a surface.
*
* \tparam Vector_3 the type of the normal vector model of `Kernel::Vector_3`.
*
* \cgalModels{NormalEstimator}
*/
template <typename Vector_3>
struct Normal_estimator_on_surface
{
public:
  /// \name Types
  /// @{

  /*!
  * The type of the normal vector.
  */
  typedef Vector_3 Normal_type;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  * Constructor that pre-computes the normals on the surface.
  *
  * \tparam Surface a model of `FaceListGraph` that represents a surface mesh.
  *
  * \param surface the surface where the normals are evaluated.
  */
  template <typename Surface>
  Normal_estimator_on_surface(const Surface& surface);

  /// @}

  /// \name Functor
  /// @{

  /*!
  * returns the normal vector of the surface element described by a type and an index.
  *
  * \tparam DimensionTag a tag that represent the element type.
  *         Can be `CGAL::Dimension_tag<0>`, `CGAL::Dimension_tag<1>` or `CGAL::Dimension_tag<2>`
  * \tparam Index the type of index of the element to evaluate.
  *
  * \param element_index the index of the element to evaluate.
  */
  template <typename DimensionTag, typename Index>
  Normal_type operator()(const Index& element_index) const;

  /// @}
};

} /* namespace Normal_estimator */

} /* namespace Feature_graph */

} /* namespace CGAL */

#endif // CGAL_FEATURE_GRAPH_NORMAL_ESTIMATOR_ON_SURFACE_H