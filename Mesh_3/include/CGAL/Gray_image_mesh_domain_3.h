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
//
//
// Author(s)     : Stephane Tayeb, Laurent Rineau
//

#ifndef CGAL_GRAY_IMAGE_MESH_DOMAIN_3_H
#define CGAL_GRAY_IMAGE_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Random.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/Image_to_labeled_function_wrapper.h>

#include <functional>

namespace CGAL {

namespace internal {

  template<typename T>
  struct Greater_than {
    typedef T argument_type;
    Greater_than(const T& second) : second(second) {}
    bool operator()(const T& first) const {
      return std::greater<T>()(first, second);
    }
    T second;
  };

}

/**
 * @class Gray_image_mesh_domain_3
 *
 *
 */
template<class Image,
         class BGT,
         typename Image_word_type = float,
         typename Transform = internal::Greater_than<double>,
         typename Subdomain_index = int>
class Gray_image_mesh_domain_3
  : public Labeled_mesh_domain_3<
  Mesh_3::Image_to_labeled_function_wrapper<Image, BGT,
                                            Image_word_type,
                                            Subdomain_index,
                                            Transform,
                                            false> ,
  BGT>
{
public:
  typedef Mesh_3::Image_to_labeled_function_wrapper<Image, BGT,
                                                    Image_word_type,
                                                    Subdomain_index,
                                                    Transform,
                                                    false>           Wrapper;

  typedef Labeled_mesh_domain_3<Wrapper, BGT>                        Base;

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
           compute_bounding_box(image),
           error_bound,
           Null_subdomain_index(),
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
           compute_bounding_box(image),
           error_bound,
           Null_subdomain_index(),
           p_rng)
  {
    CGAL_assertion(transform(value_outside) == 0);
  }

  /// Destructor
  virtual ~Gray_image_mesh_domain_3() {}


private:
  /// Returns a box enclosing image \c im
  Bbox_3 compute_bounding_box(const Image& im) const
  {
    return Bbox_3(-1,-1,-1,
                  double(im.xdim())*im.vx()+1,
                  double(im.ydim())*im.vy()+1,
                  double(im.zdim())*im.vz()+1);
  }

private:
  // Disabled copy constructor & assignment operator
  typedef Gray_image_mesh_domain_3<Image, BGT> Self;
  Gray_image_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Gray_image_mesh_domain_3



}  // end namespace CGAL



#endif // CGAL_GRAY_IMAGE_MESH_DOMAIN_3_H
