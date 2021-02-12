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
// Author(s)     : Laurent Rineau, Jane Tournois

#ifndef CGAL_MESH_3_INITIALIZE_TRIANGULATION_FROM_GRAY_IMAGE_H
#define CGAL_MESH_3_INITIALIZE_TRIANGULATION_FROM_GRAY_IMAGE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/initialize_triangulation_from_labeled_image.h>

#include <CGAL/tags.h>

template<class C3T3, class MeshDomain, class MeshCriteria,
                typename FT, typename Image_word_type>
void initialize_triangulation_from_gray_image(C3T3& c3t3,
    const MeshDomain& domain,
    const CGAL::Image_3& image,
    const MeshCriteria& criteria,
    const FT& iso_value,
    Image_word_type,
    bool protect_features = false)
{
  typedef CGAL::Null_functor Functor;
  Functor image_values_to_subdomain_indices;

  using CGAL::Mesh_3::internal::Create_gray_image_values_to_subdomain_indices;
  typedef Create_gray_image_values_to_subdomain_indices<Functor> C_i_v_t_s_i;
  typedef typename C_i_v_t_s_i::type Image_values_to_subdomain_indices;
  Image_values_to_subdomain_indices transform_fct =
    C_i_v_t_s_i()(image_values_to_subdomain_indices, iso_value);

  initialize_triangulation_from_labeled_image(c3t3, domain, image, criteria,
                                                                  Image_word_type(),
                                                                  protect_features,
                                                                  transform_fct);
}

#endif // CGAL_MESH_3_INITIALIZE_TRIANGULATION_FROM_GRAY_IMAGE_H
