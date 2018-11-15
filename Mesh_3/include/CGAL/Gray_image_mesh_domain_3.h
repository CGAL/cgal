// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2012  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Stephane Tayeb, Laurent Rineau
//

#ifndef CGAL_GRAY_IMAGE_MESH_DOMAIN_3_H
#define CGAL_GRAY_IMAGE_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Random.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/Image_to_labeled_function_wrapper.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {

/**
 * @class Gray_image_mesh_domain_3
 *
 *
 */
template<class Image,
         class BGT,
         typename Image_word_type_ = float,
         typename Transform = Mesh_3::internal::Greater_than<double>,
         typename Subdomain_index = int>
class
CGAL_DEPRECATED_MSG
( "The class template `CGAL::Gray_image_mesh_domain_3` is now deprecated. "
  "Use the static member function template "
  "`Labeled_mesh_domain_3<K>::create_gray_image_mesh_domain` instead.")
Gray_image_mesh_domain_3
  : public Labeled_mesh_domain_3<BGT, Subdomain_index>
{
public:
  typedef Image_word_type_ Image_word_type;
  typedef Mesh_3::Image_to_labeled_function_wrapper<Image_word_type,
                                                    double,
                                                    Subdomain_index,
                                                    false>           Wrapper;

  typedef Labeled_mesh_domain_3<BGT, Subdomain_index>                Base;

  typedef typename Base::Sphere_3 Sphere_3;
  typedef typename Base::FT FT;
  typedef BGT Geom_traits;
  typedef CGAL::Bbox_3 Bbox_3;

  /// Constructor
  Gray_image_mesh_domain_3(const Image& image,
                           const Image_word_type iso_value,
                           const Image_word_type value_outside = 0.,
                           const FT& error_bound = FT(1e-3),
                           CGAL::Random* p_rng = NULL)
    : Base(Wrapper(image, 
                   Transform(iso_value),
                   Transform(iso_value)(value_outside)),
           Mesh_3::internal::compute_bounding_box(image),
           error_bound,
           p_rng)
  {
    CGAL_assertion(Transform(iso_value)(value_outside) == 0);
  }

  Gray_image_mesh_domain_3(const Image& image,
                           const Transform& transform,
                           const Image_word_type value_outside = 0.,
                           const FT& error_bound = FT(1e-3),
                           CGAL::Random* p_rng = NULL)
    : Base(Wrapper(image, transform, transform(value_outside)),
           Mesh_3::internal::compute_bounding_box(image),
           error_bound,
           p_rng)
  {
    CGAL_assertion(transform(value_outside) == 0);
  }

  /// Destructor
  virtual ~Gray_image_mesh_domain_3() {}
};  // end class Gray_image_mesh_domain_3

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>


#endif // CGAL_GRAY_IMAGE_MESH_DOMAIN_3_H
