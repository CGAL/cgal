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

#ifndef CGAL_FEATURE_GRAPH_IMAGE_AMBROSIOTORTORELLI_H
#define CGAL_FEATURE_GRAPH_IMAGE_AMBROSIOTORTORELLI_H

#include <CGAL/license/Feature_graph.h>
#include <CGAL/Feature_graph/internal/DGtal_interface.h>

namespace CGAL {

namespace Feature_graph {

namespace internal
{

template <typename SharpnessContainer, typename SharpnessType>
struct Image_container_sharpness_estimator
{
  typedef SharpnessType Sharpness_number_type;

  const SharpnessContainer& sharpness_container;
  Sharpness_number_type selection_threshold;

  Image_container_sharpness_estimator(const SharpnessContainer& sharpness_container, const Sharpness_number_type& selection_threshold)
    : sharpness_container(sharpness_container), selection_threshold(selection_threshold)
  { }

  template <typename DimensionTag, typename Descriptor, typename Domain>
  Sharpness_number_type operator()(const Descriptor& element_descriptor, const Domain& domain) const
  {
    return sharpness_container.at(element_descriptor) - selection_threshold;
  }
};

template <typename NormalContainer, typename NormalType>
struct Image_container_normal_estimator
{
  typedef NormalType Normal_type;

  const NormalContainer& normal_container;

  Image_container_normal_estimator(const NormalContainer& normal_container)
    : normal_container(normal_container)
  { }

  // requires element_descriptor to describe a surfel
  template <typename DimensionTag, typename Descriptor, typename Domain>
  Normal_type operator()(const Descriptor& element_descriptor, const Domain& domain) const
  {
    return normal_container.at(element_descriptor);
  }
};

} /* namespace internal */


template <typename K>
struct Image_AmbrosioTortorelli
{
public:
  typedef typename K::Vector_3 Normal_type;
  typedef double Sharpness_type;
  typedef internal::Image_container_sharpness_estimator<
      std::unordered_map<std::size_t, Sharpness_type>,
    Sharpness_type> Sharpness_estimator;
  typedef internal::Image_container_normal_estimator<
      std::unordered_map<std::size_t, Normal_type>,
    Normal_type> Normal_estimator;

  Image_AmbrosioTortorelli() { }

  template <typename Image>
  Image_AmbrosioTortorelli(const Image& image,
    const Sharpness_type& selection_threshold = Sharpness_type(0.25),
    const float& value_outside = 0,
    const double& alpha = 0.05, const double& lambda = 0.02,
    const double& width_start = 2.0, const double& width_end = 0.5, const double& width_divide = 2.0)
    : selection_threshold(selection_threshold)
  {
    typedef boost::associative_property_map<std::unordered_map<std::size_t, Sharpness_type>> SharpnessPropertyMap;
    typedef boost::associative_property_map<std::unordered_map<std::size_t, Normal_type>> NormalPropertyMap;
    compute_AmbrosioTortorelli_on_labelled_image<K, Image, SharpnessPropertyMap, Sharpness_type, NormalPropertyMap, Normal_type>(
        image, SharpnessPropertyMap(sharpness_values), NormalPropertyMap(normals),
        alpha, lambda, width_start, width_end, width_divide,
        value_outside);
  }

  Sharpness_estimator sharpness_estimator() const
  {
    return Sharpness_estimator(sharpness_values, selection_threshold);
  }
  Normal_estimator normal_estimator() const
  {
    return Normal_estimator(normals);
  }

private:
  std::unordered_map<std::size_t, Sharpness_type> sharpness_values; // on pointels, linels and surfels
  std::unordered_map<std::size_t, Normal_type> normals; // only on surfels
  Sharpness_type selection_threshold;
};

} /* namespace Feature_graph */

} /* namespace CGAL */

#endif // CGAL_FEATURE_GRAPH_IMAGE_AMBROSIOTORTORELLI_H