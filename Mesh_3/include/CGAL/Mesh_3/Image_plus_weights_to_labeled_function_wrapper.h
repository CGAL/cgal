// Copyright (c) 2016  GeometryFactory Sarl (France).
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
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_MESH_3_IMAGE_PLUS_WEIGHTS_TO_LABELED_FUNCTION_WRAPPER_H
#define CGAL_MESH_3_IMAGE_PLUS_WEIGHTS_TO_LABELED_FUNCTION_WRAPPER_H

#include <CGAL/function_objects.h>


namespace CGAL {

namespace ImageIO {

template <typename Image_word_type,
          typename Weights_type>
class Weighted_indicator_factory
{
  const Image_3& image;
  const Image_3& weights;

public:

  Weighted_indicator_factory(const Image_3& image,
                             const Image_3& weights)
    : image(image)
    , weights(weights)
  {}

  class Indicator : public std::unary_function<Image_word_type, double>
  {
    const Weighted_indicator_factory& f;
    const Image_word_type label;
  public:
    Indicator(const Weighted_indicator_factory& f,
              const Image_word_type label) : f(f), label(label) {}

    double operator()(const Image_word_type& x) const
    {
      const std::ptrdiff_t offset = &x - (Image_word_type *)f.image.data();
      const Weights_type w = (std::max)(
          Weights_type(128), ((Weights_type *)f.weights.data())[offset]);
      return (x == label) ? w : (255 - w);
    }
  }; // end nested class Indicator

  Indicator indicator(const Image_word_type label) const
  {
    return Indicator(*this, label);
  }
}; // end template class ImageIO::Weighted_indicator_factory

} // end namespace CGAL::ImageIO

namespace Mesh_3 {


/**
 * @class Image_plus_weights_to_labeled_function_wrapper
 *
 * Wraps a pair of images into a labeled function which takes his values into
 * N. Uses weighted trilinear interpolation.
 */
template<class Image_,
         class BGT,
         typename Image_word_type = unsigned char,
         typename Weights_type = unsigned char,
         typename Return_type = int,
         typename Transform = Identity<Return_type>
         >
class Image_plus_weights_to_labeled_function_wrapper
{
public:
  // Types
  typedef Return_type return_type;
  typedef Image_word_type word_type;
  typedef typename BGT::Point_3   Point_3;

  /// Constructor
  Image_plus_weights_to_labeled_function_wrapper
  (const Image_& image,
   const Image_& weights_image,
   const Transform& transform = Transform(),
   const Return_type value_outside = 0)
    : r_im_(image)
    , r_weights_im_(weights_image)
    , transform(transform)
    , value_outside(value_outside)
    , indicator_factory(image, weights_image)
  {
  }

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Image_plus_weights_to_labeled_function_wrapper() {}

  /**
   * Returns an int corresponding to the label at point \c p
   * @param p the input point
   * @return the label at point \c p
   */
  return_type operator()(const Point_3& p, const bool = true) const
  {

    return transform(
      r_im_.template labellized_trilinear_interpolation<Image_word_type>(
          CGAL::to_double(p.x()),
          CGAL::to_double(p.y()),
          CGAL::to_double(p.z()),
          value_outside,
          indicator_factory));
  }

private:
  /// Labeled image to wrap
  const Image_& r_im_;
  const Image_& r_weights_im_;
  const Transform transform;
  const Return_type value_outside;
  CGAL::ImageIO::Weighted_indicator_factory<Image_word_type,
                                            Weights_type> indicator_factory;
};  // end class Image_plus_weights_to_labeled_function_wrapper


}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_IMAGE_PLUS_WEIGHTS_TO_LABELED_FUNCTION_WRAPPER_H
