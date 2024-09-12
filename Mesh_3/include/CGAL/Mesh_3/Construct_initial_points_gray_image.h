// Copyright (c) 2015,2016 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Jane Tournois and Ange Clement

#ifndef CGAL_MESH_3_CONSTRUCT_INITIAL_POINTS_GRAY_IMAGE_H
#define CGAL_MESH_3_CONSTRUCT_INITIAL_POINTS_GRAY_IMAGE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/Construct_initial_points_labeled_image.h>

#include <CGAL/Image_3.h>

namespace CGAL
{

/*!
 * \ingroup PkgMesh3Initializers
 *
 * Functor for generating initial points in gray images.
 * This functor is a model of the `InitialPointsGenerator` concept,
 * and can be passed as a parameter to `CGAL::make_mesh_3` using the
 * `CGAL::parameters::initial_points_generator()` parameter function.
 *
 * On images that contain multiple non-connected components,
 * this functor will scan the full image and
 * output points on every component.
 *
 * \sa `CGAL::parameters::initial_points_generator()`
 * \sa `CGAL::make_mesh_3()`
 * \sa `CGAL::Construct_initial_points_labeled_image`
 */
template <typename Functor = CGAL::Null_functor>
struct Construct_initial_points_gray_image
{
  const CGAL::Image_3 & image_;
  const double iso_value_;
  Functor image_values_to_subdomain_indices_;

  template <typename FT>
  Construct_initial_points_gray_image(const CGAL::Image_3 & image,
         const FT& iso_value,
         const Functor image_values_to_subdomain_indices = CGAL::Null_functor())
      : image_(image)
      , iso_value_(iso_value)
      , image_values_to_subdomain_indices_(image_values_to_subdomain_indices)
  { }

  /*!
  * \brief Constructs points by collecting them on the surface of all objects
  * in the image,
  * even if they are non-connected components.
  * Using this functor guarantees to initialize each connected component.
  *
  * \tparam OutputIterator model of `OutputIterator`, collecting points of type
  * `MeshDomain::Intersection`
  * \tparam MeshDomain model of `MeshDomain_3`
  * \tparam C3t3 model of `MeshComplex_3InTriangulation_3`
  *
  * \param n a lower bound on the number of points to be constructed
  */
  template <typename OutputIterator, typename MeshDomain, typename C3t3>
  OutputIterator operator()(OutputIterator pts, const MeshDomain& domain, const C3t3& c3t3, int n = 20) const
  {
    using CGAL::Mesh_3::internal::Create_gray_image_values_to_subdomain_indices;
    typedef Create_gray_image_values_to_subdomain_indices<Functor> C_i_v_t_s_i;
    typedef typename C_i_v_t_s_i::type Image_values_to_subdomain_indices;
    Image_values_to_subdomain_indices transform_fct =
      C_i_v_t_s_i()(image_values_to_subdomain_indices_, iso_value_);
    Construct_initial_points_labeled_image(image_).operator()(pts, domain, transform_fct, c3t3, n);
    return pts;
  }
};

} // end namespace CGAL

#endif // CGAL_MESH_3_CONSTRUCT_INITIAL_POINTS_GRAY_IMAGE_H
