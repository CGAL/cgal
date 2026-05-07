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

#ifndef CGAL_FEATURE_GRAPH_AMBROSIOTORTORELLI_ON_IMAGE_H
#define CGAL_FEATURE_GRAPH_AMBROSIOTORTORELLI_ON_IMAGE_H

namespace CGAL {

namespace Feature_graph {

/*!
* \ingroup PkgFeatureGraphRef
*
* \brief Class that evaluates the Ambriosio-Tortorelli normals and sharpness measure from an image.
* Two functors can then be retrieved to access the estimations.
*
* \tparam Vector_3 the type of the normal vector model of `Kernel::Vector_3`.
*
* \cgalModels{NormalEstimator}
*/
template <typename Vector_3>
struct AmbrosioTortorelli_on_image
{
public:
  /// \name Types
  /// @{

  /*!
  * The type of the normal vector.
  */
  typedef Vector_3 Normal_type;

  /*!
  * The type of the functor that allows to retrieve the sharpness values.
  * \cgalModels{SharpnessEstimator}
  */
  typedef unspecified_type Sharpness_functor;
  /*!
  * The type of the functor that allows to retrieve the normals.
  * \cgalModels{NormalEstimator}
  */
  typedef unspecified_type Normal_functor;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  * evaluates the normal and sharpness values using the Ambrosio-Tortorelli energy optimization.
  *
  * \tparam Image the image type model of `FeatureImage_3`
  * \tparam FT a model of `RealEmbeddable`
  *
  * \param image the image.
  * \param selection_threshold a threshold on the sharpness value.
  *     Elements with a sharpness value lower than this threshold are considered flat
  *     and will be given a negative value.
  */
  template <typename Image, typename FT = Sharpness_functor::Sharpness_value_type>
  AmbrosioTortorelli_on_image(const Image& image, const FT& selection_threshold = FT(0.25));

  /// @}

  /// \name Functor Accessors
  /// @{

  /*!
  * returns the functor that allows to retrieve the sharpness values.
  */
  Sharpness_functor get_sharpness_functor() const;
  /*!
  * returns the functor that allows to retrieve the normals.
  */
  Normal_functor get_normal_functor() const;

  /// @}
};

} /* namespace Feature_graph */

} /* namespace CGAL */

#endif // CGAL_FEATURE_GRAPH_AMBROSIOTORTORELLI_ON_IMAGE_H