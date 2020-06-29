// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************

#ifndef CGAL_LABELED_IMAGE_MESH_DOMAIN_3_H
#define CGAL_LABELED_IMAGE_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Random.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/Image_to_labeled_function_wrapper.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Default.h>

namespace CGAL {

/**
 * @class Labeled_image_mesh_domain_3
 *
 *
 */
template<class Image,
         class BGT,
         typename Image_word_type_ = unsigned char,
         typename Subdomain_index = int,
         class Null_subdomain_index = Default,
         class Wrapper_ = Default >
class
CGAL_DEPRECATED_MSG
( "The class template `CGAL::Labeled_image_mesh_domain_3` is now deprecated. "
  "Use the static member function template "
  "`Labeled_mesh_domain_3<K>::create_labeled_image_mesh_domain` instead.")
Labeled_image_mesh_domain_3
  : public Labeled_mesh_domain_3<BGT, Subdomain_index>
{
public:
  typedef Image_word_type_ Image_word_type;
  typedef typename Default::Get
    <Wrapper_,
     Mesh_3::Image_to_labeled_function_wrapper<Image_word_type,
                                               int,
                                               Subdomain_index>
     >::type Wrapper;
  typedef typename Default::Get<Null_subdomain_index,
                                CGAL::Null_subdomain_index>::type Null;

  typedef Labeled_mesh_domain_3<BGT, Subdomain_index> Base;

  typedef typename Base::Sphere_3 Sphere_3;
  typedef typename Base::FT FT;
  typedef BGT Geom_traits;
  typedef CGAL::Bbox_3 Bbox_3;
  typedef CGAL::Identity<Subdomain_index> Identity;

  /// Constructor
  Labeled_image_mesh_domain_3(const Image& image,
                              const FT& error_bound = FT(1e-3),
                              Subdomain_index value_outside = 0,
                              Null null = Null(),
                              CGAL::Random* p_rng = nullptr)
    : Base(Wrapper(image, Identity(), value_outside),
           compute_bounding_box(image),
           error_bound,
           parameters::null_subdomain_index = null,
           parameters::p_rng = p_rng)
  {}

  Labeled_image_mesh_domain_3(const Image& image,
                              const FT error_bound,
                              CGAL::Random* p_rng)
    : Base(Wrapper(image),
           compute_bounding_box(image),
           error_bound,
           p_rng)
  {}

  /// Destructor
  virtual ~Labeled_image_mesh_domain_3() {}

  using Base::bbox;

private:
  /// Returns a box enclosing image \c im
  Bbox_3 compute_bounding_box(const Image& im) const
  {
    return Bbox_3(-im.vx()+im.tx(),
                  -im.vy()+im.ty(),
                  -im.vz()+im.tz(),
                  double(im.xdim()+1)*im.vx()+im.tx(),
                  double(im.ydim()+1)*im.vy()+im.ty(),
                  double(im.zdim()+1)*im.vz()+im.tz());
  }
};  // end class Labeled_image_mesh_domain_3



}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_LABELED_IMAGE_MESH_DOMAIN_3_H
