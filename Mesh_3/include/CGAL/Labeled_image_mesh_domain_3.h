// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
class Labeled_image_mesh_domain_3
: public Labeled_mesh_domain_3
<typename Default::Get
   <Wrapper_,
    Mesh_3::Image_to_labeled_function_wrapper<Image, BGT,
                                              Image_word_type_,
                                              Subdomain_index>
    >::type,
 BGT,
 Null_subdomain_index
 >
{
public:
  typedef Image_word_type_ Image_word_type;
  typedef typename Default::Get
    <Wrapper_,
     Mesh_3::Image_to_labeled_function_wrapper<Image, BGT,
                                               Image_word_type,
                                               Subdomain_index>
     >::type Wrapper;
  typedef typename Default::Get<Null_subdomain_index,
                                CGAL::Null_subdomain_index>::type Null;

  typedef Labeled_mesh_domain_3<Wrapper, BGT, Null_subdomain_index> Base;

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
                              CGAL::Random* p_rng = NULL)
    : Base(Wrapper(image, Identity(), value_outside),
           compute_bounding_box(image),
           error_bound,
           null,
           p_rng)
  {}

  Labeled_image_mesh_domain_3(const Image& image,
                              const CGAL::Bbox_3& bbox,
                              const FT& error_bound = FT(1e-3),
                              Subdomain_index value_outside = 0,
                              Null null = Null(),
                              CGAL::Random* p_rng = NULL)
    : Base(Wrapper(image, Identity(), value_outside),
           bbox,
           error_bound,
           null,
           p_rng)
  {}

  Labeled_image_mesh_domain_3(const Image& image,
                              const FT error_bound,
                              CGAL::Random* p_rng)
    : Base(Wrapper(image),
           compute_bounding_box(image),
           error_bound,
           Null(),
           p_rng)
  {}

  Labeled_image_mesh_domain_3(const Image& image,
                              const CGAL::Bbox_3& bbox,
                              const FT error_bound,
                              CGAL::Random* p_rng)
    : Base(Wrapper(image),
           bbox,
           error_bound,
           Null(),
           p_rng)
  {}

  /// Destructor
  virtual ~Labeled_image_mesh_domain_3() {}

  using Base::bbox;

private:
  /// Returns a box enclosing image \c im
  Bbox_3 compute_bounding_box(const Image& im) const
  {
    return Bbox_3(-im.vx(),
                  -im.vy(),
                  -im.vz(),
                  double(im.xdim()+1)*im.vx(),
                  double(im.ydim()+1)*im.vy(),
                  double(im.zdim()+1)*im.vz());
  }

private:
  // Disabled copy constructor & assignment operator
  typedef Labeled_image_mesh_domain_3<Image, BGT> Self;
  Labeled_image_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Labeled_image_mesh_domain_3



}  // end namespace CGAL



#endif // CGAL_LABELED_IMAGE_MESH_DOMAIN_3_H
